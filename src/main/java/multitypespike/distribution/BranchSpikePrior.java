package multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.util.InputUtil;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
import beast.base.spec.type.Simplex;
import beast.base.util.Randomizer;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.special.Gamma;

import java.util.*;
import java.util.concurrent.CompletionException;
import java.util.concurrent.Executor;
import java.util.concurrent.ForkJoinPool;



@Description("Calculates the prior probability of branch spikes, accounting for the expected number of hidden speciation events " +
        "derived from the birth-death-migration process. It models the total spike size on a branch as a sum of Gamma-distributed events")
public class BranchSpikePrior extends Distribution {

    final public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-Prime parameterization (default: that of bdmDistr).",
            Input.Validate.OPTIONAL);

    final public Input<Tree> treeInput = new Input<>("tree", "tree input.", Input.Validate.REQUIRED);

    final public Input<RealVector<? extends PositiveReal>> spikeShapeInput = new Input<>("spikeShape", "shape parameter for the " +
            "gamma distribution of the spikes.", Input.Validate.REQUIRED);

    final public Input<RealVectorParam<? extends NonNegativeReal>> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree.",
            Input.Validate.REQUIRED);

    final public Input<RealScalar<? extends NonNegativeReal>> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "Difference in time between the final sample and the end of the BD process " +
            "(default: that of bdmDistr, or 0).");

    final public Input<Simplex> startTypePriorProbsInput = new Input<>("startTypePriorProbs",
            "Prior probabilities of the type at the start of the process (default: those of bdmDistr).");

    final public Input<BirthDeathMigrationDistribution> bdmDistrInput = new Input<>("bdmDistr",
            "BDMM-Prime tree prior; supplies the parameterization and start-type probabilities unless given " +
            "explicitly. Required for multi-type analyses.", Input.Validate.OPTIONAL);

    public Input<Boolean> useAnalyticalSingleTypeSolutionInput = new Input<>("useAnalyticalSingleTypeSolution",
            "Use the analytical branch spike prior when the model has only one type.",
            true);

    public Input<Boolean> initializeSpikesInput = new Input<>("initializeSpikes",
            "Initialize spike values by sampling from the BranchSpikePrior distribution (default true).",
            true);

    public Input<Double> relativeToleranceInput = new Input<>("relTolerance",
            "Relative tolerance for multi-type hidden events numerical integration.",
            1e-6);

    public Input<Double> absoluteToleranceInput = new Input<>("absTolerance",
            "Absolute tolerance for multi-type hidden events numerical integration.",
            1e-6);

    public Input<Boolean> parallelizeInput = new Input<>(
            "parallelize","Whether or not to parallelized the calculation of subtree likelihoods. " +
            "(Default true).",
            true);

    /* If a large number a cores is available (more than 8 or 10) the calculation speed can be increased by diminishing
    the parallelization factor. On the other hand, if only 2-4 cores are available, a slightly higher value (1/5 to 1/8)
    can be beneficial to the calculation speed. */
    public Input<Double> minimalProportionForParallelizationInput = new Input<>(
            "parallelizationFactor", "The minimal relative size a child subtree must have to start parallel calculations. " +
            "Default adapts to available cores: 1/20 (≥16 cores), 1/15 (8-15 cores), 1/10 (4-7 cores), 1/8 (2-3 cores).",
            getDefaultParallelizationFactor()
    );

    private Parameterization parameterization;
    private Simplex startTypePriorProbs;
    private RealScalar<? extends NonNegativeReal> finalSampleOffsetScalar;
    private boolean initialised = false;
    private double[] intervalEndTimes, A, B, weightOfNodeSubTree;
    private double[] expectedHiddenEvents, piVals, storedExpectedHiddenEvents, storedPiVals;
    private double[] nodePiVals, storedNodePiVals;
    private double lambda_i, mu_i, psi_i, t_i, A_i, B_i, finalSampleOffset;
    private boolean spikesInitialised = false;
    private boolean isParallelizedCalculation;
    private Executor pool = null;
    private boolean hiddenEventsCached = false, requiresReintegration;
    private boolean storedHiddenEventsCached = false;
    private double relTol, absTol;
    public int nodeCount, nTypes;
    public double minimalProportionForParallelization;

    @Override
    public void initAndValidate() {
        BirthDeathMigrationDistribution bdm = bdmDistrInput.get();
        parameterization = parameterizationInput.get() != null ? parameterizationInput.get()
                : bdm != null ? bdm.parameterizationInput.get() : null;
        startTypePriorProbs = startTypePriorProbsInput.get() != null ? startTypePriorProbsInput.get()
                : bdm != null ? bdm.startTypePriorProbsInput.get() : null;
        finalSampleOffsetScalar = finalSampleOffsetInput.get() != null ? finalSampleOffsetInput.get()
                : bdm != null ? bdm.finalSampleOffsetInput.get()
                : new RealScalarParam<>(0.0, NonNegativeReal.INSTANCE);

        // BEAUti may create this prior before the BDMM-Prime tree prior is linked
        if (parameterization == null) return;
        initialised = true;

        nTypes = parameterization.getNTypes();
        nodeCount = treeInput.get().getNodeCount();
        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetScalar.get();
        requiresReintegration = true;

        expectedHiddenEvents = new double[nodeCount * nTypes];
        piVals = new double[nodeCount * nTypes];
        storedExpectedHiddenEvents = new double[nodeCount * nTypes];
        storedPiVals = new double[nodeCount * nTypes];
        nodePiVals = new double[nodeCount * nTypes];
        storedNodePiVals = new double[nodeCount * nTypes];

        weightOfNodeSubTree = new double[treeInput.get().getLeafNodeCount() * 2];

        relTol = relativeToleranceInput.get();
        absTol = absoluteToleranceInput.get();

        if (nTypes != 1 || !useAnalyticalSingleTypeSolutionInput.get()) {
            if (startTypePriorProbs == null) {
                throw new IllegalArgumentException("'startTypePriorProbs' must be specified for multi-type analyses.");
            }

            if (bdmDistrInput.get() == null) {
                throw new IllegalArgumentException("BirthDeathMigrationDistribution,'bdmDistr', must be specified for multi-type analyses.");
            }

            // The spike prior reads the tree prior's stored p0/ge integration results, so switch them on
            BirthDeathMigrationDistribution bdmDistr = bdmDistrInput.get();
            if (!bdmDistr.saveIntegrationResultsInput.get()) {
                bdmDistr.saveIntegrationResultsInput.setValue(true, bdmDistr);
                bdmDistr.initAndValidate();
            }

            isParallelizedCalculation = parallelizeInput.get();
            minimalProportionForParallelization = minimalProportionForParallelizationInput.get();
        }

        // Spike shape dimension check
        int spikeShapeDim = spikeShapeInput.get().size();
        if (nTypes == 1 && spikeShapeDim > 1) {
            throw new IllegalArgumentException("Single-type model requires exactly one spikeShape parameter.");
        }
        if (nTypes > 1 && spikeShapeDim != 1 && spikeShapeDim != nTypes) {
            throw new IllegalArgumentException("For multi-type models, 'spikeShape' must have dimension 1 (shared) or nTypes (" + nTypes + ").");
        }

        if (nTypes == 1) {
            A = new double[parameterization.getTotalIntervalCount()];
            B = new double[parameterization.getTotalIntervalCount()];
            computeConstants(A, B);

            spikesInput.get().setDimension(nodeCount);

        } else {
            spikesInput.get().setDimension(nodeCount * nTypes);
        }

        if (isParallelizedCalculation) {
            pool = ForkJoinPool.commonPool();
        }
    }

    private void ensureInitialised() {
        if (initialised) return;
        initAndValidate();
        if (!initialised)
            throw new IllegalArgumentException("BranchSpikePrior '" + getID() + "' needs the BDMM-Prime tree prior: " +
                    "set its bdmDistr input (in BEAUti, select BDMM-Prime as the tree prior of this partition).");
    }

    private void initialiseSpikes() {

        // Initialise spike values by sampling from the spike prior distribution
        if (nTypes > 1) {
            // p0/ge integration results must be current, whatever the order of the distributions
            bdmDistrInput.get().calculateLogP();
            sampleMultiTypeSpikes();
        } else {
            sampleSingleTypeSpikes();
        }
    }


    private void computeConstants(double[] A, double[] B) {

        for (int i = parameterization.getTotalIntervalCount() - 1; i >= 0; i--) {

            double p_i_prev;

            if (i + 1 < parameterization.getTotalIntervalCount()) {
                p_i_prev = get_p_i(parameterization.getBirthRates()[i + 1][0],
                        parameterization.getDeathRates()[i + 1][0],
                        parameterization.getSamplingRates()[i + 1][0],
                        A[i + 1], B[i + 1],
                        parameterization.getIntervalEndTimes()[i + 1],
                        parameterization.getIntervalEndTimes()[i]);
            } else {
                p_i_prev = 1.0;
            }

            double rho_i = parameterization.getRhoValues()[i][0];
            double lambda_i = parameterization.getBirthRates()[i][0];
            double mu_i = parameterization.getDeathRates()[i][0];
            double psi_i = parameterization.getSamplingRates()[i][0];

            A[i] = Math.sqrt((lambda_i - mu_i - psi_i) * (lambda_i - mu_i - psi_i) + 4 * lambda_i * psi_i);
            B[i] = ((1 - 2 * (1 - rho_i) * p_i_prev) * lambda_i + mu_i + psi_i) / A[i];
        }
    }


    private double get_p_i(double lambda, double mu, double psi, double A, double B, double t_i, double t) {

        if (lambda > 0.0) {
            double v = Math.exp(A * (t_i - t)) * (1 + B);
            return (lambda + mu + psi - A * (v - (1 - B)) / (v + (1 - B)))
                    / (2 * lambda);
        } else {
            // The limit of p_i as lambda -> 0
            return 0.5;
        }
    }

    private void updateParametersForInterval(int i) {
        // update parameters for interval index i
        lambda_i = parameterization.getBirthRates()[i][0];
        mu_i = parameterization.getDeathRates()[i][0];
        psi_i = parameterization.getSamplingRates()[i][0];
        t_i = parameterization.getIntervalEndTimes()[i];

        A_i = A[i];
        B_i = B[i];
    }


    /**
     * Single type expected number of hidden events for interval (t0,t1)
     */
    private double integral_2lambda_i_p_i(double t_0, double t_1) {
        double t0 = t_i - t_0;
        double t1 = t_i - t_1;

        return ((t0 - t1) * (mu_i + psi_i + lambda_i + A_i) + 2.0 * Math
                .log(((-B_i - 1) * Math.exp(A_i * t1) + B_i - 1)
                        / ((-B_i - 1) * Math.exp(A_i * t0) + B_i - 1)));
    }


    /**
     * Single type expected number of hidden events for branch
     */
    public double getExpNrHiddenEventsForBranch(Node node) {
        if (node.isRoot() || node.isDirectAncestor()) return 0;

        double expNrHiddenEvents = 0;
        int nodeIndex = parameterization.getNodeIntervalIndex(node, finalSampleOffset);
        int parentIndex = parameterization.getNodeIntervalIndex(node.getParent(), finalSampleOffset);
        double t0 = parameterization.getNodeTime(node.getParent(), finalSampleOffset);
        double T = parameterization.getNodeTime(node, finalSampleOffset);
        updateParametersForInterval(parentIndex);

        if (nodeIndex == parentIndex) return integral_2lambda_i_p_i(t0, T);

        for (int k = parentIndex; k < nodeIndex; k++) {
            if (k > parentIndex) updateParametersForInterval(k);
            double t1 = intervalEndTimes[k];
            expNrHiddenEvents += integral_2lambda_i_p_i(t0, t1);
            t0 = t1;
        }

        updateParametersForInterval(nodeIndex);
        expNrHiddenEvents += integral_2lambda_i_p_i(t0, T);

        return expNrHiddenEvents;
    }

    @Override
    public double calculateLogP() {
        ensureInitialised();

        if (!spikesInitialised) {
            if (initializeSpikesInput.get()) initialiseSpikes();
            spikesInitialised = true;
        }

        if (!useAnalyticalSingleTypeSolutionInput.get() && nTypes == 1) return multiTypeCalculateLogP();
        else if (nTypes == 1) return singleTypeCalculateLogP();
        else {
            try {
                return multiTypeCalculateLogP();
            } catch (MaxCountExceededException ex) {
                // Catches it if thrown directly on the main thread
                 Log.warning("Integration error encountered in prior calculation (sync)");
                return Double.NEGATIVE_INFINITY;
            } catch (CompletionException ex) {
                // Catches it if thrown in the CompletableFuture threads.
                if (ex.getCause() instanceof MaxCountExceededException) {
                     Log.warning("Integration error encountered in prior calculation (async)");
                    return Double.NEGATIVE_INFINITY;
                } else if (ex.getCause() instanceof OutOfMemoryError) {
                    return handleOutOfMemory((OutOfMemoryError) ex.getCause());
                } else {
                    // If it was a different multithreading crash (e.g. NullPointer), re-throw it
                    throw ex;
                }
            } catch (OutOfMemoryError ex) {
                // Catches it if thrown directly on the main thread (no parallelisation, or the
                // failure happened outside the CompletableFuture-managed part of the computation).
                return handleOutOfMemory(ex);
            }
        }
    }

    /**
     * Handles integration OOM errors by rejecting the proposal (logP = -Infinity).
     */
    private double handleOutOfMemory(OutOfMemoryError ex) {
        Runtime rt = Runtime.getRuntime();
        long usedMB = (rt.totalMemory() - rt.freeMemory()) / (1024 * 1024);
        long maxMB = rt.maxMemory() / (1024 * 1024);
        Log.warning("Out of memory during multi-type hidden events integration " +
                "(nTypes=" + nTypes + ", nodeCount=" + nodeCount + ", heap ~" + usedMB + "/" + maxMB +
                " MB used) - treating this proposal as rejected (logP = -Infinity).");

        // Explicit GC call to attempt heap recovery before the next MCMC proposal.
        System.gc();

        return Double.NEGATIVE_INFINITY;
    }


    public double singleTypeCalculateLogP() {
        logP = 0.0;
        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetScalar.get();

        // Check spikeShape is positive
        double spikeShape = spikeShapeInput.get().get(0);
        if (spikeShape <= 0) {
            return Double.NEGATIVE_INFINITY;
        }

        computeConstants(A, B);

        double[] logProbs = new double[2];

        // Loop over all nodes in the tree
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);
            double branchSpike = spikesInput.get().get(nodeNr);

            // Handle origin branch and sampled ancestor branches
            if (node.isRoot() || node.isDirectAncestor()) {
                expectedHiddenEvents[nodeNr] = 0.0;
                // Spikes = 0 for origin branch and sampled ancestor branches.
                // Use a Gamma(1,1) = Exp(1) pseudo-prior for latent spikes when they are not
                // included in the model. This facilitates transitions between models of
                // different dimensions.

                // Log density of Exp(1) is -x
                logP -= branchSpike;
                continue;
            }

            // Compute expected number of hidden speciation events for this branch
            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);
            expectedHiddenEvents[nodeNr] = expNrHiddenEvents;

            // The observed speciation event adds a spike unless the parent is a sampled ancestor
            branchSpikeLogProbs(branchSpike, expNrHiddenEvents, spikeShape, logProbs);
            logP += node.getParent().isFake() ? logProbs[0] : logProbs[1];
        }

        // Numerical issue
        if (logP == Double.POSITIVE_INFINITY) logP = Double.NEGATIVE_INFINITY;

        return logP;
    }


    public double multiTypeCalculateLogP() {
        logP = 0.0;
        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetScalar.get();

        for (int i = 0; i < spikeShapeInput.get().size(); i++) {
            if (spikeShapeInput.get().get(i) <= 0) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        if (requiresReintegration) {
            MultiTypeHiddenEventsIntegrator hiddenEventsIntegrator = new MultiTypeHiddenEventsIntegrator(
                    parameterization, treeInput.get(), bdmDistrInput.get().getIntegrationResults(),
                    absTol, relTol, false,
                    isParallelizedCalculation, pool, weightOfNodeSubTree, minimalProportionForParallelization
            );
            hiddenEventsIntegrator.integrateHiddenEvents(
                    startTypeProbabilities(), parameterization, finalSampleOffset
            );

            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                Node node = treeInput.get().getNode(nodeNr);
                // A sampled ancestor has the type of its (fake) parent node
                double[] nodePi = hiddenEventsIntegrator.getPiAtNode(
                        node.isDirectAncestor() ? node.getParent().getNr() : nodeNr);
                for (int i = 0; i < nTypes; i++) {
                    nodePiVals[nodeNr * nTypes + i] = Math.min(Math.max(nodePi[i], 0.0), 1.0);
                }

                if (node.isRoot() || node.isDirectAncestor()) {
                    for (int i = 0; i < nTypes; i++) {
                        expectedHiddenEvents[nodeNr * nTypes + i] = 0.0;
                    }
                    continue;
                }
                double[] expHidden = hiddenEventsIntegrator.getExpNrHiddenEventsForNode(nodeNr);
                // π is evaluated at the parent node: the observed speciation event is at the
                // start (past-most point) of the branch, which is the parent node.
                int parentNr = node.getParent().getNr();
                double[] pi = hiddenEventsIntegrator.getPiAtNode(parentNr);
                System.arraycopy(expHidden, 0, expectedHiddenEvents, nodeNr * nTypes, nTypes);
                for (int i = 0; i < nTypes; i++) {
                    piVals[nodeNr * nTypes + i] = Math.min(Math.max(pi[i], 0.0), 1.0);
                }
            }
            hiddenEventsCached = true;
        }

        double[] logP0 = new double[nTypes];
        double[] logP1 = new double[nTypes];
        double[] logProbs = new double[2];

        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);

            if (node.isRoot() || node.isDirectAncestor()) {
                // Spikes = 0 for origin branch and sampled ancestor branches.
                // Use a Gamma(1,1) = Exp(1) pseudo-prior for latent spikes when they are not
                // included in the model. This facilitates transitions between models of
                // different dimensions.

                for (int i = 0; i < nTypes; i++) {
                    double branchSpike = spikesInput.get().get(nodeNr * nTypes + i);

                    // Log density of Exp(1) is -x
                    logP -= branchSpike;
                }
                continue;
            }

            // Calculate P0 (no observed event) and P1 (1 observed event) for all types
            for (int i = 0; i < nTypes; i++) {
                double expNrHiddenEvents = expectedHiddenEvents[nodeNr * nTypes + i];
                if (expNrHiddenEvents > 1000.0) return Double.NEGATIVE_INFINITY;

                branchSpikeLogProbs(spikesInput.get().get(nodeNr * nTypes + i),
                        expNrHiddenEvents, getSpikeShape(i), logProbs);
                logP0[i] = logProbs[0];
                logP1[i] = logProbs[1];
            }

            logP += nodeLogPrior(logP0, logP1, piVals, nodeNr * nTypes, node.getParent().isFake());
        }

        // Numerical issue
        if (logP == Double.POSITIVE_INFINITY) logP = Double.NEGATIVE_INFINITY;
        return logP;
    }


    /**
     * Log density of a spike made up of nEvents Gamma(spikeShape, rate spikeShape) increments.
     */
    public static double logSpikeDensity(double spike, int nEvents, double spikeShape) {
        double alpha = spikeShape * nEvents;
        return alpha * Math.log(spikeShape) + (alpha - 1.0) * Math.log(spike) - spikeShape * spike
                - Gamma.logGamma(alpha);
    }

    private static final double LOG_TAIL_TOLERANCE = Math.log(1e-12);
    private static final int MAX_HIDDEN_EVENTS = 100000;

    /**
     * Log probabilities of the spike of one type on a branch, summed over the Poisson(expNrHiddenEvents)
     * number of hidden events: logProbs[0] if the observed speciation event is not of this type (or the
     * parent is a sampled ancestor), logProbs[1] if it is. The summands are log-concave in the number
     * of hidden events, so the sum stops once they decrease and the remaining tail is negligible.
     */
    public static void branchSpikeLogProbs(double spike, double expNrHiddenEvents, double spikeShape,
                                           double[] logProbs) {
        if (spike == 0.0) {
            logProbs[0] = -expNrHiddenEvents;
            logProbs[1] = Double.NEGATIVE_INFINITY;
            return;
        }
        if (!(spike > 0.0) || !(spikeShape > 0.0)) {
            logProbs[0] = Double.NEGATIVE_INFINITY;
            logProbs[1] = Double.NEGATIVE_INFINITY;
            return;
        }
        if (!(expNrHiddenEvents > 0.0)) {
            logProbs[0] = Double.NEGATIVE_INFINITY;
            logProbs[1] = logSpikeDensity(spike, 1, spikeShape);
            return;
        }

        double logMu = Math.log(expNrHiddenEvents);
        double logPk = -expNrHiddenEvents;           // log Poisson probability of k hidden events
        double logDensityK = Double.NEGATIVE_INFINITY; // log density of k spike increments
        double sum0 = Double.NEGATIVE_INFINITY, sum1 = Double.NEGATIVE_INFINITY;
        double prev0 = Double.NEGATIVE_INFINITY, prev1 = Double.NEGATIVE_INFINITY;

        for (int k = 0; k < MAX_HIDDEN_EVENTS; k++) {
            double logDensityNext = logSpikeDensity(spike, k + 1, spikeShape);
            double term0 = logPk + logDensityK;
            double term1 = logPk + logDensityNext;
            sum0 = logAdd(sum0, term0);
            sum1 = logAdd(sum1, term1);

            if (k >= expNrHiddenEvents
                    && tailIsNegligible(term0, prev0, sum0)
                    && tailIsNegligible(term1, prev1, sum1))
                break;

            prev0 = term0;
            prev1 = term1;
            logDensityK = logDensityNext;
            logPk += logMu - Math.log(k + 1);
        }

        logProbs[0] = sum0;
        logProbs[1] = sum1;
    }

    // For a decreasing log-concave sequence the tail after the current term is at most term * r / (1 - r)
    private static boolean tailIsNegligible(double term, double prevTerm, double logSum) {
        if (term == Double.NEGATIVE_INFINITY) return true;
        if (!(term < prevTerm)) return false;
        double logRatio = term - prevTerm;
        double logTail = term + logRatio - Math.log(-Math.expm1(logRatio));
        return logTail < logSum + LOG_TAIL_TOLERANCE;
    }

    public static double logAdd(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) return b;
        if (b == Double.NEGATIVE_INFINITY) return a;
        double max = Math.max(a, b);
        return max + Math.log1p(Math.exp(-Math.abs(a - b)));
    }

    /**
     * Joint log prior of the spikes on one branch across types, given the per-type log probabilities
     * without (logP0) and with (logP1) the observed speciation event. The observed event belongs to
     * exactly one type, with probabilities pi[piOffset + i]; sampled ancestors are not speciation events.
     */
    public static double nodeLogPrior(double[] logP0, double[] logP1, double[] pi, int piOffset,
                                      boolean hasFakeParent) {
        int nTypes = logP0.length;
        double logP = 0.0;

        if (hasFakeParent) {
            for (int i = 0; i < nTypes; i++) logP += logP0[i];
            return logP;
        }

        double maxLogTerm = Double.NEGATIVE_INFINITY;
        double[] logTerms = new double[nTypes];
        for (int i = 0; i < nTypes; i++) {
            double p = pi[piOffset + i];
            logTerms[i] = Double.NEGATIVE_INFINITY;
            if (p > 0) {
                double term = Math.log(p) + logP1[i];
                for (int j = 0; j < nTypes; j++) {
                    if (j != i) term += logP0[j];
                }
                logTerms[i] = term;
                if (term > maxLogTerm) maxLogTerm = term;
            }
        }
        if (maxLogTerm == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;

        double sumExp = 0.0;
        for (int i = 0; i < nTypes; i++) {
            if (logTerms[i] > Double.NEGATIVE_INFINITY) sumExp += Math.exp(logTerms[i] - maxLogTerm);
        }
        return maxLogTerm + Math.log(sumExp);
    }

    private double[] startTypeProbabilities() {
        double[] probs = new double[startTypePriorProbs.size()];
        for (int i = 0; i < probs.length; i++) probs[i] = startTypePriorProbs.get(i);
        return probs;
    }

    @Override
    public List<String> getArguments() {
        List<String> args = new ArrayList<>();
        args.add(spikesInput.get().getID());
        return args;
    }

    @Override
    public List<String> getConditions() {
        List<String> conds = new ArrayList<>();
        if (treeInput.get() != null) conds.add(treeInput.get().getID());
        if (parameterization != null) conds.add(parameterization.getID());
        if (spikeShapeInput.get() instanceof BEASTInterface) conds.add(((BEASTInterface) spikeShapeInput.get()).getID());
        if (startTypePriorProbs instanceof BEASTInterface) conds.add(((BEASTInterface) startTypePriorProbs).getID());
        if (bdmDistrInput.get() != null) conds.add(bdmDistrInput.get().getID());
        return conds;
    }


    @Override
    public void sample(State state, Random random) {

        if (sampledFlag) return;
        sampledFlag = true;
        ensureInitialised();
        // Cause conditional parameters to be sampled
        sampleConditions(state, random);

        // Single-type case
        if (nTypes == 1) {

            sampleSingleTypeSpikes();

            // Multi-type case
        } else {

            // Call calculate LogP to get p0ge integration results
            bdmDistrInput.get().calculateLogP();

            sampleMultiTypeSpikes();

        }
    }

    private void sampleSingleTypeSpikes() {

        double spikeShape = spikeShapeInput.get().get(0);
        spikesInput.get().setDimension(nodeCount);

        if (spikeShape <= 0) {
            throw new IllegalArgumentException("Cannot sample spikes because spikeShape is non-positive " + spikeShape);
        }

        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {

            Node node = treeInput.get().getNode(nodeNr);

            // Handle origin branch and sampled ancestor branch
            if (node.isRoot() || node.isDirectAncestor()) {
                spikesInput.get().set(nodeNr, 0.0);
                continue;
            }

            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);
            int nHiddenEvents = (int) Randomizer.nextPoisson(expNrHiddenEvents);
            int nSpikes = node.getParent().isFake() ? nHiddenEvents : nHiddenEvents + 1;
            double alpha = spikeShape * nSpikes;

            // Sample spike from Gamma distribution if nSpikes > 0
            // Uses spikeShape instead of 1/spikeShape due to different parameterisation of the Gamma distribution
            double spike = (nSpikes == 0) ? 0.0 : Randomizer.nextGamma(alpha, spikeShape);
            spikesInput.get().set(nodeNr, spike);

        }
    }

    private void sampleMultiTypeSpikes() {

        spikesInput.get().setDimension(nodeCount * nTypes);

        for (int i = 0; i < nTypes; i++) {
            if (getSpikeShape(i) <= 0) {
                throw new IllegalArgumentException("Cannot sample spikes because spikeShape is non-positive " + getSpikeShape(i));
            }
        }

        MultiTypeHiddenEventsIntegrator hiddenEventsIntegrator = new MultiTypeHiddenEventsIntegrator(
                parameterization, treeInput.get(), bdmDistrInput.get().getIntegrationResults(),
                1e-10, 1e-10, false, isParallelizedCalculation, pool,
                weightOfNodeSubTree, minimalProportionForParallelization
        );
        hiddenEventsIntegrator.integrateHiddenEvents(
                startTypeProbabilities(), parameterization, finalSampleOffset
        );

        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {

            Node node = treeInput.get().getNode(nodeNr);

            if (node.isRoot() || node.isDirectAncestor()) {
                // Zero spikes for root and direct ancestors
                for (int i = 0; i < nTypes; i++) {
                    spikesInput.get().set(nodeNr * nTypes + i, 0.0);
                }
                continue;
            }

            // Compute expected number of hidden speciation events for this branch for each types
            double[] expNrHiddenEventsArray = hiddenEventsIntegrator.getExpNrHiddenEventsForNode(nodeNr);

            // Compute π at time of the observed speciation event of the node, π(t₀)
            Node parent = node.getParent();
            int parentNr = parent.getNr();
            double[] piArray = hiddenEventsIntegrator.getPiAtNode(parentNr);

            // Sample the specific type of the observed speciation event from Categorical(piArray) distribution
            int obsEventType = -1;
            if (!parent.isFake()) {
                double u = Randomizer.nextDouble();
                double cumsum = 0.0;
                for (int i = 0; i < nTypes; i++) {
                    cumsum += Math.max(piArray[i], 0.0);
                    if (u <= cumsum) {
                        obsEventType = i;
                        break;
                    }
                }
                // Guard against piArray not summing to exactly 1 due to numerical integration error
                if (obsEventType == -1) obsEventType = nTypes - 1;
            }

            for (int i = 0; i < nTypes; i++) {

                double expNrHiddenEvents = expNrHiddenEventsArray[i];

                int obsEvent = (i == obsEventType) ? 1 : 0;

                int nHiddenEvents = (int) Randomizer.nextPoisson(expNrHiddenEvents);
                int nSpikes = node.getParent().isFake() ? nHiddenEvents : nHiddenEvents + obsEvent;
                double spikeShape = getSpikeShape(i);
                double alpha = spikeShape * nSpikes;

                double spike = (nSpikes == 0) ? 0.0 : Randomizer.nextGamma(alpha, spikeShape);
                spikesInput.get().set(nodeNr * nTypes + i, spike);

            }
        }
    }


    private static double getDefaultParallelizationFactor() {
        int cores = Runtime.getRuntime().availableProcessors();

        if (cores >= 16) return 1.0 / 20;   // 0.05
        if (cores >= 8)  return 1.0 / 15;   // ~0.067
        if (cores >= 4)  return 1.0 / 10;   // 0.1
        return 1.0 / 8;                      // 0.125 for 2 cores
    }


    /**
     * Methods for passing precomputed values to loggers
     */
    public double getExpectedHiddenEvents(int nodeNr){
        return expectedHiddenEvents[nodeNr];
    }

    public double getExpectedHiddenEvents(int nodeNr, int type) {
        if (nTypes == 1) return expectedHiddenEvents[nodeNr];
        return expectedHiddenEvents[nodeNr * nTypes + type];
    }

    /** Type probabilities at the parent node, i.e. of the observed speciation event at the start of the branch. */
    public double getPiVals(int nodeNr, int type) {
        return piVals[nodeNr * nTypes + type];
    }

    /** Type probabilities of the lineage at the node itself. */
    public double getNodeTypeProbability(int nodeNr, int type) {
        return nodePiVals[nodeNr * nTypes + type];
    }

    public double getSpikeShape(int type) {
        int spikeShapeDim = spikeShapeInput.get().size();
        if (nTypes == 1 || spikeShapeDim == 1) return spikeShapeInput.get().get(0);
        else return spikeShapeInput.get().get(type);
    }

    @Override
    protected boolean requiresRecalculation() {
        // If only spikes and spikeShape are dirty, then no need to reintegrate expNrHiddenEvents
        requiresReintegration = (!hiddenEventsCached ||
                InputUtil.isDirty(treeInput) ||
                InputUtil.isDirty(parameterizationInput) ||
                InputUtil.isDirty(bdmDistrInput) ||
                InputUtil.isDirty(startTypePriorProbsInput) ||
                InputUtil.isDirty(finalSampleOffsetInput)
        );

        return InputUtil.isDirty(spikesInput) ||
                InputUtil.isDirty(spikeShapeInput) ||
                InputUtil.isDirty(treeInput) ||
                InputUtil.isDirty(parameterizationInput) ||
                InputUtil.isDirty(bdmDistrInput) ||
                InputUtil.isDirty(startTypePriorProbsInput) ||
                InputUtil.isDirty(finalSampleOffsetInput);
    }

    @Override
    public void store() {
        System.arraycopy(expectedHiddenEvents, 0, storedExpectedHiddenEvents, 0, expectedHiddenEvents.length);
        System.arraycopy(piVals, 0, storedPiVals, 0, piVals.length);
        System.arraycopy(nodePiVals, 0, storedNodePiVals, 0, nodePiVals.length);
        storedHiddenEventsCached = hiddenEventsCached;
        super.store();
    }

    @Override
    public void restore() {
        double[] tmpExp = storedExpectedHiddenEvents;
        storedExpectedHiddenEvents = expectedHiddenEvents;
        expectedHiddenEvents = tmpExp;

        double[] tmpPi = storedPiVals;
        storedPiVals = piVals;
        piVals = tmpPi;

        double[] tmpNodePi = storedNodePiVals;
        storedNodePiVals = nodePiVals;
        nodePiVals = tmpNodePi;

        hiddenEventsCached = storedHiddenEventsCached;
        super.restore();
    }

}