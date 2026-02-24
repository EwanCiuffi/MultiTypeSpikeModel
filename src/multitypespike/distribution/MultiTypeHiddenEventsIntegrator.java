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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;


public class MultiTypeHiddenEventsIntegrator implements Loggable {

    private final ContinuousOutputModel[] p0geComArray;
    private final double[][] b;
    private final double[][][] M, b_ij;

    public final double[][] storedResults;

    private final int nTypes;
    public double[] intervalEndTimes;
    private final Tree tree;

    protected double integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance;
    private static final double EPS = 1e-6; // Small value to avoid division by zero

    // Optional storage of π continuous output trajectories for testing purposes
    private final boolean storePiTrajectories;
    private final ContinuousOutputModel[] piIntegrationResults;


    private final boolean isParallelizedCalculation;
    private final ThreadPoolExecutor pool;
    private final double minimalProportionForParallelization;
    private double parallelizationThreshold;


    private final double[] weightOfNodeSubTree;

    public MultiTypeHiddenEventsIntegrator(Parameterization parameterization, Tree tree, ContinuousOutputModel[] p0geComArray,
                                           double absoluteTolerance, double relativeTolerance, boolean storePiTrajectories,
                                           boolean isParallelizedCalculation, ThreadPoolExecutor pool, double[] weightOfNodeSubTree, double minimalProportionForParallelization) {

        this.tree = tree;
        this.p0geComArray = p0geComArray;
        this.b = parameterization.getBirthRates();
        this.M = parameterization.getMigRates();
        this.b_ij = parameterization.getCrossBirthRates();
        this.nTypes = parameterization.getNTypes();
        int nodeCount = tree.getNodeCount();
        this.intervalEndTimes = parameterization.getIntervalEndTimes();

        this.storedResults = new double[nodeCount][2 *nTypes];

        integrationMinStep = parameterization.getTotalProcessLength() * 1e-12;
        integrationMaxStep= parameterization.getTotalProcessLength() / 10;
        this.absoluteTolerance = absoluteTolerance;
        this.relativeTolerance = relativeTolerance;

        this.storePiTrajectories = storePiTrajectories;
        this.piIntegrationResults = storePiTrajectories ? new ContinuousOutputModel[nodeCount] : null;

        this.isParallelizedCalculation = isParallelizedCalculation;
        this.pool = pool;
        this.weightOfNodeSubTree = weightOfNodeSubTree;
        this.minimalProportionForParallelization = minimalProportionForParallelization;
    }


    private ContinuousOutputModel getP0GeIntegrationResults(int nodeNr) {
    return p0geComArray[nodeNr];
}

    /**
     * Returns the interpolated p0 and ge values at a given time for an edge.
     */
    public double[] getP0Ge(int nodeNr, double time) {
        //  p0:  (0 .. dim-1)
        //  ge: (dim .. 2*dim-1)
        ContinuousOutputModel p0geCom = getP0GeIntegrationResults(nodeNr);
        p0geCom.setInterpolatedTime(time);
        double[] p0ge = p0geCom.getInterpolatedState();

        // Trim away small negative values due to numerical integration errors
        for (int i=0; i < p0ge.length; i++) {
            if (p0ge[i] < 0)
                p0ge[i] = 0.0;
        }
        return p0ge;
    }

    private final class BranchIntegrator implements FirstOrderDifferentialEquations {

        private final int nodeNr;
        private final int interval;

        BranchIntegrator(int nodeNr, int interval) {
            this.nodeNr = nodeNr;
            this.interval = interval;
        }

        @Override
        public int getDimension() {
            return 2 * nTypes;
        }

        @Override
        public void computeDerivatives(double t, double[] y, double[] yDot) {

            double[] p0ge = getP0Ge(nodeNr, t);

            for (int i = 0; i < nTypes; i++) {

                final double ge_i = p0ge[nTypes + i];
                final double p0_i = p0ge[i];

                yDot[i] = 0;

                for (int j = 0; j < nTypes; j++) {

                    if (j == i) continue;

                    final double ge_j = p0ge[nTypes + j];
                    final double p0_j = p0ge[j];

                    yDot[i] += ((b_ij[interval][j][i] *  p0_j + M[interval][j][i]) * (ge_i
                            / Math.max(ge_j, EPS))) * y[j];

                    yDot[i] -= ((b_ij[interval][i][j] * p0_i + M[interval][i][j]) * (ge_j
                            / Math.max(ge_i, EPS))) * y[i];
                }

                yDot[nTypes + i] = 2.0 * y[i] * b[interval][i] * p0_i;
            }
        }

        public void integrate(double tStart, double tEnd, double[] state, ContinuousOutputModel segment) {

            FirstOrderIntegrator integrator = threadLocalIntegrator.get();

            integrator.clearStepHandlers();

            if (segment != null) integrator.addStepHandler(segment);

            integrator.integrate(this, tStart, state, tEnd, state);
        }

    }

    private final ThreadLocal<FirstOrderIntegrator> threadLocalIntegrator =
        ThreadLocal.withInitial(() ->
                new DormandPrince54Integrator(integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance)
    );


    public void integrateAlongEdge(Node node, double tStart, Parameterization parameterization,
                                   double finalSampleOffset, double[] initialState) {

        final int nodeNr = node.getNr();
        double thisTime = tStart;
        final double tEnd = parameterization.getNodeTime(node, finalSampleOffset);
        int thisInterval = parameterization.getIntervalIndex(thisTime);
        final int endInterval = parameterization.getNodeIntervalIndex(node, finalSampleOffset);
        double[] state = initialState.clone();

//        if (storedResults[node.getNr()][0] != 0)
//            throw new RuntimeException("Branch integrated twice for node " + node.getNr());

        ContinuousOutputModel fullModel = storePiTrajectories ? new ContinuousOutputModel() : null;

        while (thisInterval < endInterval) {

            final double nextTime = intervalEndTimes[thisInterval];

            if (Utils.lessThanWithPrecision(thisTime, nextTime)) {

                BranchIntegrator system = new BranchIntegrator(nodeNr, thisInterval);

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

            BranchIntegrator system = new BranchIntegrator(nodeNr, thisInterval);

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

        final double nodeTime =
                parameterization.getNodeTime(node, finalSampleOffset);

        if (!node.isRoot() && !node.isDirectAncestor()) {
            final int parentNr = node.getParent().getNr();
            double[] state = getInitialConditionsAtNode(parentNr);
            integrateAlongEdge(node, parentTime, parameterization, finalSampleOffset, state);
        }

        if (node.isLeaf())
            return;

        Node child1 = node.getChild(0);
        Node child2 = node.getChild(1);

        Node firstChild = child1;
        Node secondChild = child2;

        // Heavier subtree first
        if (weightOfNodeSubTree[child2.getNr()] >
                weightOfNodeSubTree[child1.getNr()]) {
            firstChild = child2;
            secondChild = child1;
        }

        Future<?> secondFuture = null;

        // Decide if we parallelise the second child only
        boolean parallel =
                isParallelizedCalculation &&
                        weightOfNodeSubTree[firstChild.getNr()] > parallelizationThreshold &&
                        weightOfNodeSubTree[secondChild.getNr()] > parallelizationThreshold;

        if (parallel) {
            final Node second = secondChild;
            secondFuture = pool.submit(() -> {
                integrateAtNode(second, nodeTime, parameterization, finalSampleOffset);
                return null;
            });
        }

        // Process first child locally
        integrateAtNode(firstChild, nodeTime, parameterization, finalSampleOffset);

        // If we didn’t parallelise, do second child now
        if (!parallel) {
            integrateAtNode(secondChild, nodeTime, parameterization, finalSampleOffset);
        } else {
            try {
                secondFuture.get();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
    }

    /**
     * Performs recursive integration of π and expected hidden events along the tree.
     */
    public void integrateHiddenEvents(double[] startTypePriorProbs, Parameterization parameterization, double finalSampleOffset) {

        Node root = this.tree.getRoot();
        int rootNr = root.getNr();

        double rootTime = parameterization.getNodeTime(root, finalSampleOffset);

        // Set initial conditions at root
        double[] state = getInitialConditionsAtRoot(startTypePriorProbs, rootTime, rootNr);

        // Store initial conditions at root
        storeResultsAtNode(state, rootNr);

        updateParallelizationThreshold();

        // Start pre-order traversal integration from root
        integrateAtNode(root, rootTime, parameterization, finalSampleOffset);
    }


    public ContinuousOutputModel getPiIntegrationResultsForNode(int nodeNr) {
        if (!storePiTrajectories) {
            throw new IllegalStateException("π trajectories not stored; set storePiTrajectories to true.");
        }
        return piIntegrationResults[nodeNr];
    }

    /**
     *  Retrieves the expected number of hidden events per type for a branch
     */
    public double[] getExpNrHiddenEventsForNode(int nodeNr) {
        double[] expHiddenEvents = new double[nTypes];
        System.arraycopy(storedResults[nodeNr], nTypes, expHiddenEvents, 0, nTypes);
        return expHiddenEvents;
    }

    /**
     *  Retrieves array of π values at a node
     */
    public double[] getPiAtNode(int nodeNr) {
        double[] pi = new double[nTypes];
        System.arraycopy(storedResults[nodeNr], 0, pi, 0, nTypes);
        return pi;
    }


    private double[] getInitialConditionsAtRoot(double[] startTypePriorProbs, double rootTime, int rootNr) {

        double[] state = new double[2 * nTypes];
        double total = 0.0;
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
            state[nTypes + i] = 0.0;   // hidden events start at zero
        }

        // Optional: store continuous π trajectories
        ContinuousOutputModel com = storePiTrajectories ? new ContinuousOutputModel() : null;

        if (storePiTrajectories) {
            FirstOrderIntegrator integrator = threadLocalIntegrator.get();
            integrator.clearStepHandlers();
            integrator.addStepHandler(com);
        }

        // Integrate along the single branch
        BranchIntegrator system = new BranchIntegrator(rootNr,
                parameterization.getIntervalIndex(startTime));

        FirstOrderIntegrator integrator = threadLocalIntegrator.get();
        integrator.integrate(system, startTime, state, endTime, state);

        // Store results
        storeResultsAtNode(state, rootNr);

        if (storePiTrajectories) {
            piIntegrationResults[rootNr] = com;
        }
    }

    /**
     * Perform an initial traversal of the tree to get the 'weights' (sum of all its edges lengths) of all subtrees
     * Useful for performing parallelized calculations on the tree.
     * The weights of the subtrees tell us the depth at which parallelization should stop, so as to not parallelize on subtrees that are too small.
     * Results are stored in 'weightOfNodeSubTree' array
     *
     * @param tree tree whose subtree to compute weights of
     */
    private void getAllSubTreesWeights(TreeInterface tree) {
        Node root = tree.getRoot();
        double weight = 0;
        for (final Node child : root.getChildren()) {
            weight += getSubTreeWeight(child);
        }
        weightOfNodeSubTree[root.getNr()] = weight;
    }

    /**
     * Perform an initial traversal of the subtree to get its 'weight': sum of all its edges.
     *
     * @param node root of subtree
     * @return weight
     */
    private double getSubTreeWeight(Node node) {

        // if leaf, stop recursion, get length of branch above and return
        if (node.isLeaf()) {
            weightOfNodeSubTree[node.getNr()] = node.getLength();
            return node.getLength();
        }

        // else, iterate over the children of the node
        double weight = 0;
        for (final Node child : node.getChildren()) {
            weight += getSubTreeWeight(child);
        }
        // add length of parental branch
        weight += node.getLength();
        // store the value
        weightOfNodeSubTree[node.getNr()] = weight;

        return weight;
    }

    private void updateParallelizationThreshold() {
        if (isParallelizedCalculation) {
            getAllSubTreesWeights(tree);
            // set 'parallelizationThreshold' to a fraction of the whole tree weight.
            // The size of this fraction is determined by a tuning parameter. This parameter should be adjusted (increased) if more computation cores are available
            parallelizationThreshold = weightOfNodeSubTree[tree.getRoot().getNr()] * minimalProportionForParallelization;
        }
    }

    @Override
    public void init(PrintStream out) {
    }

    @Override
    public void log(long sample, PrintStream out) {
    }

    @Override
    public void close(PrintStream out) {
    }
}