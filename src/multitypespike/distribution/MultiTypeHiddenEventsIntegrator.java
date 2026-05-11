package multitypespike.distribution;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

public class MultiTypeHiddenEventsIntegrator implements Loggable {

    private final ContinuousOutputModel[] p0geComArray;
    private final double[][] b;
    private final double[][][] M, b_ij;

    private final double[][] storedResults;

    private final int nTypes;
    private final double[] intervalEndTimes;
    private final Tree tree;

    protected double integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance;

    private final boolean isParallelizedCalculation;
    private final Executor pool;
    private final double minimalProportionForParallelization;
    private double parallelizationThreshold;

    private final double[] weightOfNodeSubTree;

    // Optional storage of π continuous output trajectories for testing purposes
    private final boolean storePiTrajectories;
    private final ContinuousOutputModel[] piIntegrationResults;

    public MultiTypeHiddenEventsIntegrator(Parameterization parameterization, Tree tree, ContinuousOutputModel[] p0geComArray,
                                           double absoluteTolerance, double relativeTolerance, boolean storePiTrajectories,
                                           boolean isParallelizedCalculation, Executor pool, double[] weightOfNodeSubTree,
                                           double minimalProportionForParallelization) {

        this.tree = tree;
        this.p0geComArray = p0geComArray;
        this.b = parameterization.getBirthRates();
        this.M = parameterization.getMigRates();
        this.b_ij = parameterization.getCrossBirthRates();
        this.nTypes = parameterization.getNTypes();
        int nodeCount = tree.getNodeCount();
        this.intervalEndTimes = parameterization.getIntervalEndTimes();

        this.storedResults = new double[nodeCount][2 * nTypes];

        integrationMinStep = parameterization.getTotalProcessLength() * 1e-100;
        integrationMaxStep = parameterization.getTotalProcessLength() / 10;
        this.absoluteTolerance = absoluteTolerance;
        this.relativeTolerance = relativeTolerance;

        this.storePiTrajectories = storePiTrajectories;
        this.piIntegrationResults = storePiTrajectories ? new ContinuousOutputModel[nodeCount] : null;

        this.isParallelizedCalculation = isParallelizedCalculation;
        this.pool = pool;
        this.weightOfNodeSubTree = weightOfNodeSubTree;
        this.minimalProportionForParallelization = minimalProportionForParallelization;

    }


    private final class BranchIntegrator implements FirstOrderDifferentialEquations {

        private int interval;
        private final int nodeNr;
        private final double[] p0geBuffer;
        private final double[] inv_ge;
        private static final double EPS = 1e-3; // Small value to avoid division by zero

        BranchIntegrator(int interval, int nodeNr) {
            this.interval   = interval;
            this.nodeNr     = nodeNr;
            this.p0geBuffer = new double[2 * nTypes];
            this.inv_ge     = new double[nTypes];
        }

        void setInterval(int interval) { this.interval = interval; }

        @Override public int getDimension() { return 2 * nTypes; }

        @Override
        public void computeDerivatives(double t, double[] y, double[] yDot) {

            ContinuousOutputModel com = p0geComArray[nodeNr];
            com.setInterpolatedTime(t);
            double[] state = com.getInterpolatedState();
            System.arraycopy(state, 0, p0geBuffer, 0, 2 * nTypes);

            final double[] lambda = b[interval];
            final double[][] lambda_ij = b_ij[interval];
            final double[][] migRate = M[interval];

            // Precompute inverses
            for (int k = 0; k < nTypes; k++) {
                inv_ge[k] = 1.0 / Math.max(p0geBuffer[nTypes + k], EPS);
            }

            Arrays.fill(yDot, 0.0);

            for (int i = 0; i < nTypes; i++) {
                final double ge_i  = p0geBuffer[nTypes + i];
                final double p0_i  = p0geBuffer[i];
                final double inv_i = inv_ge[i];

                double dyi = 0.0;
                for (int j = 0; j < nTypes; j++) {
                    if (j == i) continue;
                    final double ge_j = p0geBuffer[nTypes + j];
                    final double p0_j = p0geBuffer[j];
                    dyi += (lambda_ij[j][i] * p0_j + migRate[j][i]) * (ge_i * inv_ge[j]) * y[j];
                    dyi -= (lambda_ij[i][j] * p0_i + migRate[i][j]) * (ge_j * inv_i)     * y[i];
                }

                /*  π equations: (0 .. dim-1)  */
                yDot[i] = dyi;

                /*  hidden events equations: (dim .. 2*dim-1)  */
                yDot[nTypes + i] = 2.0 * y[i] * lambda[i] * p0_i;
            }
        }

        void integrate(double tStart, double tEnd, double[] state,
                       ContinuousOutputModel segment) {

            FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                    integrationMinStep, integrationMaxStep,
                    absoluteTolerance, relativeTolerance);

            if (segment != null) integrator.addStepHandler(segment);

            integrator.integrate(this, tStart, state, tEnd, state);
        }
    }

    public void integrateAlongEdge(Node node, double tStart, Parameterization parameterization,
                                   double finalSampleOffset, double[] initialState) {

        final int nodeNr  = node.getNr();
        final double tEnd = parameterization.getNodeTime(node, finalSampleOffset);

        double thisTime     = tStart;
        int    thisInterval = parameterization.getIntervalIndex(thisTime);
        final int endInterval = parameterization.getNodeIntervalIndex(node, finalSampleOffset);

        double[] state = initialState.clone();

        BranchIntegrator system = new BranchIntegrator(thisInterval, nodeNr);

        ContinuousOutputModel fullModel = storePiTrajectories ? new ContinuousOutputModel() : null;

        while (thisInterval < endInterval) {
            final double nextTime = intervalEndTimes[thisInterval];

            if (Utils.lessThanWithPrecision(thisTime, nextTime)) {
                system.setInterval(thisInterval);

                if (storePiTrajectories) {
                    ContinuousOutputModel segment = new ContinuousOutputModel();
                    system.integrate(thisTime, nextTime, state, segment);
                    fullModel.append(segment);
                } else {
                    system.integrate(thisTime, nextTime, state, null);
                }
            }

            thisTime = nextTime;
            thisInterval++;
        }

        if (Utils.lessThanWithPrecision(thisTime, tEnd)) {
            system.setInterval(thisInterval);

            if (storePiTrajectories) {
                ContinuousOutputModel segment = new ContinuousOutputModel();
                system.integrate(thisTime, tEnd, state, segment);
                fullModel.append(segment);
            } else {
                system.integrate(thisTime, tEnd, state, null);
            }
        }

        storeResultsAtNode(state, nodeNr);

        if (storePiTrajectories)
            piIntegrationResults[nodeNr] = fullModel;
    }


    public void integrateAtNode(Node node,
                                double parentTime,
                                Parameterization parameterization,
                                double finalSampleOffset) {

        final double nodeTime = parameterization.getNodeTime(node, finalSampleOffset);

        if (!node.isRoot() && !node.isDirectAncestor()) {
            final int parentNr = node.getParent().getNr();
            double[] state = getInitialConditionsAtNode(parentNr);
            integrateAlongEdge(node, parentTime, parameterization, finalSampleOffset, state);
        }

        if (node.isLeaf()) return;

        Node child1 = node.getChild(0);
        Node child2 = node.getChild(1);

        Node firstChild  = child1;
        Node secondChild = child2;

        // Heavier subtree first
        if (weightOfNodeSubTree[child2.getNr()] > weightOfNodeSubTree[child1.getNr()]) {
            firstChild  = child2;
            secondChild = child1;
        }

        // Decide if we parallelise the second child only
        boolean parallel = isParallelizedCalculation &&
                        weightOfNodeSubTree[firstChild.getNr()]  > parallelizationThreshold &&
                        weightOfNodeSubTree[secondChild.getNr()] > parallelizationThreshold;

        if (parallel && pool != null) {
            final Node second = secondChild;

            // Submit asynchronously with CompletableFuture
            CompletableFuture<Void> secondFuture = CompletableFuture.runAsync(
                    () -> integrateAtNode(second, nodeTime, parameterization, finalSampleOffset), pool);

            // Process first child on current thread
            integrateAtNode(firstChild, nodeTime, parameterization, finalSampleOffset);

            // Wait for second child (non-blocking join)
            secondFuture.join();
        } else {
            integrateAtNode(firstChild,  nodeTime, parameterization, finalSampleOffset);
            integrateAtNode(secondChild, nodeTime, parameterization, finalSampleOffset);
        }
    }

    /**
     * Performs recursive integration of π and expected hidden events along the tree.
     */
    public void integrateHiddenEvents(double[] startTypePriorProbs,
                                      Parameterization parameterization,
                                      double finalSampleOffset) {

        Node root   = this.tree.getRoot();
        int  rootNr = root.getNr();

        double rootTime = parameterization.getNodeTime(root, finalSampleOffset);

        // Set initial conditions at root
        double[] state = getInitialConditionsAtRoot(startTypePriorProbs, rootTime, rootNr);
        storeResultsAtNode(state, rootNr);

        updateParallelizationThreshold();

        // Start pre-order traversal integration from root
        integrateAtNode(root, rootTime, parameterization, finalSampleOffset);
    }

    // Public accessors

    /** Returns the interpolated p0 and ge values at a given time for an edge. */
    public double[] getP0Ge(int nodeNr, double time) {
        //  p0:  (0 .. dim-1)
        //  ge: (dim .. 2*dim-1)
        ContinuousOutputModel p0geCom = p0geComArray[nodeNr];
        p0geCom.setInterpolatedTime(time);
        return p0geCom.getInterpolatedState();
    }

    public ContinuousOutputModel getPiIntegrationResultsForNode(int nodeNr) {
        if (!storePiTrajectories)
            throw new IllegalStateException("π trajectories not stored; set storePiTrajectories to true.");
        return piIntegrationResults[nodeNr];
    }

    /**
     *  Returns the expected number of hidden events per type for an edge.
     */
    public double[] getExpNrHiddenEventsForNode(int nodeNr) {
        double[] expHiddenEvents = new double[nTypes];
        System.arraycopy(storedResults[nodeNr], nTypes, expHiddenEvents, 0, nTypes);
        return expHiddenEvents;
    }

    /**
     *  Returns the array of π values at a node.
     */
    public double[] getPiAtNode(int nodeNr) {
        double[] pi = new double[nTypes];
        System.arraycopy(storedResults[nodeNr], 0, pi, 0, nTypes);
        return pi;
    }

    // Initial conditions helpers

    private double[] getInitialConditionsAtRoot(double[] startTypePriorProbs, double rootTime, int rootNr) {
        double[] state    = new double[2 * nTypes];
        double   total    = 0.0;
        double[] p0geInit = getP0Ge(rootNr, rootTime);

        for (int i = 0; i < nTypes; i++) {
            state[i] = p0geInit[i + nTypes] * startTypePriorProbs[i];
            total += state[i];
        }
        for (int i = 0; i < nTypes; i++) {
            state[i] /= total;
            state[nTypes + i] = 0.0;
        }
        return state;
    }

    public double[] getInitialConditionsAtNode(int nodeNr) {
        double[] state = new double[2 * nTypes];

        // Copy π from parent
        System.arraycopy(storedResults[nodeNr], 0, state, 0, nTypes);

        // Hidden events set to zero
        Arrays.fill(state, nTypes, 2 * nTypes, 0.0);

        return state;
    }

    private void storeResultsAtNode(double[] state, int nodeNr) {
        System.arraycopy(state, 0, storedResults[nodeNr], 0, 2 * nTypes);
    }


    /**
     * Integrates π and hidden events along a single-lineage tree (for testing only).
     */
    public void integrateSingleLineage(double[] startTypePriorProbs,
                                       Parameterization parameterization,
                                       double startTime,
                                       double endTime) {

        Node root = this.tree.getRoot();

        if (!root.isLeaf()) {
            throw new RuntimeException("Tree has more than one lineage! Only single-lineage trees supported.");
        }

        final int rootNr = root.getNr();

        // Build initial state vector (2 * nTypes)
        double[] state = new double[2 * nTypes];
        double[] p0geInit = getP0Ge(rootNr, startTime);
        double total = 0.0;

        for (int i = 0; i < nTypes; i++) {
            state[i] = p0geInit[i + nTypes] * startTypePriorProbs[i];
            total += state[i];
        }
        for (int i = 0; i < nTypes; i++) {
            state[i] /= total;
            state[nTypes + i] = 0.0; // hidden events start at zero
        }

        BranchIntegrator system = new BranchIntegrator(parameterization.getIntervalIndex(startTime), rootNr);

        // Optional: store continuous π trajectories
        ContinuousOutputModel com = storePiTrajectories ? new ContinuousOutputModel() : null;

        FirstOrderIntegrator integrator =
                new DormandPrince54Integrator(
                        integrationMinStep,
                        integrationMaxStep,
                        absoluteTolerance,
                        relativeTolerance);

        if (storePiTrajectories) integrator.addStepHandler(com);

        // Integrate along the single branch
        integrator.integrate(system, startTime, state, endTime, state);

        storeResultsAtNode(state, rootNr);

        if (storePiTrajectories)
            piIntegrationResults[rootNr] = com;
    }

    // Parallelization helpers

    /**
     * Perform an initial traversal of the tree to get the 'weights' (sum of all its edges lengths) of all subtrees
     * Useful for performing parallelized calculations on the tree.
     * The weights of the subtrees tell us the depth at which parallelization should stop, so as to not parallelize on subtrees that are too small.
     * Results are stored in 'weightOfNodeSubTree' array
     *
     * @param tree tree whose subtree to compute weights of
     */
    private void getAllSubTreesWeights(TreeInterface tree) {
        Node   root   = tree.getRoot();
        double weight = 0;
        for (final Node child : root.getChildren())
            weight += getSubTreeWeight(child);
        weightOfNodeSubTree[root.getNr()] = weight;
    }

    /**
     * Perform an initial traversal of the subtree to get its 'weight': sum of all its edges.
     *
     * @param node root of subtree
     * @return weight
     */
    private double getSubTreeWeight(Node node) {

        // If leaf, stop recursion, get length of branch above and return
        if (node.isLeaf()) {
            weightOfNodeSubTree[node.getNr()] = node.getLength();
            return node.getLength();
        }

        // Else, iterate over the children of the node
        double weight = 0;
        for (final Node child : node.getChildren())
            weight += getSubTreeWeight(child);

        // Add length of parental branch
        weight += node.getLength();

        // Store the value
        weightOfNodeSubTree[node.getNr()] = weight;

        return weight;
    }

    private void updateParallelizationThreshold() {
        if (isParallelizedCalculation) {
            getAllSubTreesWeights(tree);
            // set 'parallelizationThreshold' to a fraction of the whole tree weight.
            // The size of this fraction is determined by a tuning parameter.
            // This parameter should be adjusted (increased) if more computation cores are available
            parallelizationThreshold =
                    weightOfNodeSubTree[tree.getRoot().getNr()] * minimalProportionForParallelization;
        }
    }


    @Override public void init(PrintStream out) { }
    @Override public void log(long sample, PrintStream out) { }
    @Override public void close(PrintStream out) { }
}