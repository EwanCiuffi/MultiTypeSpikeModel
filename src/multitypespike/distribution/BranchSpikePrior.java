package multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math3.exception.MaxCountExceededException;

import java.util.*;
import java.util.concurrent.Executor;
import java.util.concurrent.ForkJoinPool;



@Description("Calculates the prior probability of branch spikes, accounting for the expected number of hidden speciation events " +
        "derived from the birth-death-migration process. It models the total spike size on a branch as a sum of Gamma-distributed events")
public class BranchSpikePrior extends Distribution {

    final public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-prime parameterization object (see BDMM-prime package for available parameterizations).",
            Input.Validate.REQUIRED);

    final public Input<Tree> treeInput = new Input<>("tree", "tree input.", Input.Validate.REQUIRED);

    final public Input<RealParameter> spikeShapeInput = new Input<>("spikeShape", "shape parameter for the " +
            "gamma distribution of the spikes.", Input.Validate.REQUIRED);

    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree.",
            Input.Validate.REQUIRED);

    final public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "if provided, the difference in time between the final sample and the end of the BD process.",
            new RealParameter("0.0"));

    final public Input<RealParameter> startTypePriorProbsInput = new Input<>("startTypePriorProbs",
            "The prior probabilities for the initial individual type.",
            Input.Validate.OPTIONAL);

    final public Input<BirthDeathMigrationDistribution> bdmDistrInput = new Input<>("bdmDistr",
            "Birth-death-migration model distribution.", Input.Validate.OPTIONAL);

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
    private double[] intervalEndTimes, A, B, weightOfNodeSubTree;
    private double[] expectedHiddenEvents, piVals, storedExpectedHiddenEvents, storedPiVals;
    private double lambda_i, mu_i, psi_i, t_i, A_i, B_i, finalSampleOffset;
    private boolean spikesInitialised = false;
    private static boolean isParallelizedCalculation;
    private Executor pool = null;
    private boolean hiddenEventsCached = false, requiresReintegration = true;
    private boolean storedHiddenEventsCached = false;
    private static double relTol, absTol;
    public int nodeCount, nTypes;
    public double minimalProportionForParallelization;

    // Threshold for treating a spike as numerically zero
    private static final double spikeZeroTol = 1e-9;

    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        nodeCount = treeInput.get().getNodeCount();
        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetInput.get().getArrayValue(0);

        expectedHiddenEvents = new double[nodeCount * nTypes];
        piVals = new double[nodeCount * nTypes];
        storedExpectedHiddenEvents = new double[nodeCount * nTypes];
        storedPiVals = new double[nodeCount * nTypes];

        weightOfNodeSubTree = new double[treeInput.get().getLeafNodeCount() * 2];

        relTol = relativeToleranceInput.get();
        absTol = absoluteToleranceInput.get();

        if (nTypes != 1) {
            if (startTypePriorProbsInput.get() == null) {
                throw new IllegalArgumentException("'startTypePriorProbs' must be specified for multi-type analyses.");
            }

            if (bdmDistrInput.get() == null) {
                throw new IllegalArgumentException("BirthDeathMigrationDistribution,'bdmDistr', must be specified for multi-type analyses.");
            }

            isParallelizedCalculation = parallelizeInput.get();
            minimalProportionForParallelization = minimalProportionForParallelizationInput.get();
        }

        // Spike shape dimension check
        int spikeShapeDim = spikeShapeInput.get().getDimension();
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

        if (!initializeSpikesInput.get()) {
            // Ensure spike values are initialised to positive values
            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                for (int i = 0; i < nTypes; i++) {
                    int index = nodeNr * nTypes + i;
                    if (spikesInput.get().getValue(index) == 0) spikesInput.get().setValue(index, 1e-10);
                }
            }
        }
    }

    private void initialiseSpikes() {

        // Initialise spike values by sampling from the spike prior distribution
        if (nTypes > 1) sampleMultiTypeSpikes();
        else sampleSingleTypeSpikes();

        // Ensure spike values are initialised to positive values
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            for (int i = 0; i < nTypes; i++) {
                int index = nodeNr * nTypes + i;
                if (spikesInput.get().getValue(index) == 0) spikesInput.get().setValue(index, 1e-10);
            }
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

        if (!spikesInitialised) {
            if (initializeSpikesInput.get()) initialiseSpikes();
            spikesInitialised = true;
        }

        if (!useAnalyticalSingleTypeSolutionInput.get() && nTypes == 1) return multiTypeCalculateLogP();
        else if (nTypes == 1) return singleTypeCalculateLogP();
        else try {return multiTypeCalculateLogP();
            } catch(MaxCountExceededException ex) {
                System.err.println("Warning: integration error encountered in prior calculation");
                return Double.NEGATIVE_INFINITY;
            }
    }


    // If there are too many hidden events on a branch (e.g. during mixing) then the gamma distribution shape is large,
    // which causes instabilities
    final double MAX_CUM_SUM = 0.999;

    public double singleTypeCalculateLogP() {
        logP = 0.0;
        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetInput.get().getArrayValue(0);

        // Check spikeShape is positive
        double spikeShape = spikeShapeInput.get().getArrayValue(0);
        if (spikeShape <= 0) {
            return Double.NEGATIVE_INFINITY;
        }

        computeConstants(A, B);

        // Reuse a single GammaDistributionImpl across all nodes
        GammaDistributionImpl gamma = new GammaDistributionImpl(spikeShape, 1.0 / spikeShape);

        // Loop over all nodes in the tree
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);
            double branchSpike = spikesInput.get().getValue(nodeNr);

            // Handle origin branch and sampled ancestor branches
            if (node.isRoot() || node.isDirectAncestor()) {
                expectedHiddenEvents[nodeNr] = 0.0;
                // Spikes = 0 for origin branch and sampled ancestor branches.
                // Use a Gamma(2.0,0.5) pseudo-prior for latent spikes when they are not
                // included in the model. This facilitates transitions between models of
                // different dimensions.

                gamma.setAlpha(2.0);
                gamma.setBeta(0.5);
                logP += gamma.logDensity(branchSpike);
                continue;
            }

            // Compute expected number of hidden speciation events for this branch
            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);
            expectedHiddenEvents[nodeNr] = expNrHiddenEvents;

            boolean isZeroSpike = (branchSpike < spikeZeroTol);

            if (expNrHiddenEvents > 0) {
                // Integrate over all possible spike values
                double branchP = 0.0;
                double cumsum = 0.0;
                int k = 0;
                double logFactorialK = 0.0;

                while (cumsum < MAX_CUM_SUM) {
                    // Probability of k hidden events P(k) under a Poisson(mu)
                    double logpk = -expNrHiddenEvents + k * Math.log(expNrHiddenEvents) - logFactorialK;
                    cumsum += Math.exp(logpk);

                    // Number of spikes is k + 1 unless parent of the node is a sampled ancestor (fake), in which case it is k
                    int nSpikes = node.getParent().isFake() ? k : k + 1;

                    if (nSpikes == 0) {
                        // Valid zero spike
                        if (isZeroSpike) branchP += Math.exp(logpk);

                    } else {
                        // Compute log-probability of observed spike under Gamma distribution
                        gamma.setAlpha(spikeShape * nSpikes);
                        gamma.setBeta(1.0 / spikeShape);
                        double gammaLogP = gamma.logDensity(branchSpike);
                        if (!isZeroSpike && Double.isFinite(gammaLogP)) {
                            branchP += Math.exp(logpk + gammaLogP);
                        }
                    }
                    k++;
                    logFactorialK += Math.log(k);
                }
                // Add log-likelihood for this branch
                logP += Math.log(branchP);

            } else {
                if (!node.getParent().isFake()) {
                    gamma.setAlpha(spikeShape);
                    gamma.setBeta(1.0 / spikeShape);
                    logP += gamma.logDensity(branchSpike);
                } else if (!isZeroSpike) {
                    logP += Double.NEGATIVE_INFINITY;
                }
            }
        }

        // Numerical issue
        if (logP == Double.POSITIVE_INFINITY) logP = Double.NEGATIVE_INFINITY;

        return logP;
    }


    public double multiTypeCalculateLogP() {
        logP = 0.0;
        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetInput.get().getArrayValue(0);

        for (int i = 0; i < spikeShapeInput.get().getDimension(); i++) {
            if (spikeShapeInput.get().getArrayValue(i) <= 0) {
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
                    startTypePriorProbsInput.get().getDoubleValues(), parameterization, finalSampleOffset
            );

            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                Node node = treeInput.get().getNode(nodeNr);
                if (node.isRoot() || node.isDirectAncestor()) {
                    // The root/SA pseudo-prior contribution to logP is intentionally NOT
                    // applied here. It is applied below in the main per-node loop.
                    for (int i = 0; i < nTypes; i++) {
                        expectedHiddenEvents[nodeNr * nTypes + i] = 0.0;
                    }
                    continue;
                }
                double[] expHidden = hiddenEventsIntegrator.getExpNrHiddenEventsForNode(nodeNr);
                // π is evaluated at the parent node: the observed speciation event is at the
                // start of the branch, which is the parent node.
                int parentNr = node.getParent().getNr();
                double[] pi = hiddenEventsIntegrator.getPiAtNode(parentNr);
                System.arraycopy(expHidden, 0, expectedHiddenEvents, nodeNr * nTypes, nTypes);
                for (int i = 0; i < nTypes; i++) {
                    piVals[nodeNr * nTypes + i] = Math.min(Math.max(pi[i], 0.0), 1.0);
                }
            }
            hiddenEventsCached = true;
        }

        // Reuse a single GammaDistributionImpl across all nodes
        GammaDistributionImpl gamma = new GammaDistributionImpl(1.0, 1.0);
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);

            if (node.isRoot() || node.isDirectAncestor()) {
                // Spikes = 0 for origin branch and sampled ancestor branches.
                // Use a Gamma(2.0,0.5) pseudo-prior for latent spikes when they are not
                // included in the model. This facilitates transitions between models of
                // different dimensions.

                for (int i = 0; i < nTypes; i++) {
                    double branchSpike = spikesInput.get().getValue(nodeNr * nTypes + i);

                    gamma.setAlpha(2);
                    gamma.setBeta(0.5);
                    logP += gamma.logDensity(branchSpike);
                }
                continue;
            }

            Node parent = node.getParent();
            boolean isFakeParent = parent.isFake();

            for (int i = 0; i < nTypes; i++) {
                double branchSpike = spikesInput.get().getValue(nodeNr * nTypes + i);
                double expNrHiddenEvents = expectedHiddenEvents[nodeNr * nTypes + i];
                double pi = piVals[nodeNr * nTypes + i];

                double logPi = (pi > 0) ? Math.log(pi) : Double.NEGATIVE_INFINITY;
                double logOneMinusPi = (pi < 1.0) ? Math.log(1.0 - pi) : Double.NEGATIVE_INFINITY;

                boolean isZeroSpike = (branchSpike < spikeZeroTol);

                double spikeShape = getSpikeShape(i);

                if (expNrHiddenEvents > 0) {
                    double branchP = 0.0;
                    double cumsum = 0.0;
                    int k = 0;
                    double logFactorialK = 0.0;

                    while (cumsum < MAX_CUM_SUM) {
                        double logpk = -expNrHiddenEvents + k * Math.log(expNrHiddenEvents) - logFactorialK;
                        double pk = Math.exp(logpk);
                        cumsum += pk;

                        // obsEvent == 0: the observed speciation is NOT of type i
                        {
                            int nSpikes = k; // k + 0
                            double logPobs0 = logOneMinusPi;
                            if (nSpikes == 0) {
                                if (isZeroSpike && logPobs0 > Double.NEGATIVE_INFINITY)
                                    branchP += pk * Math.exp(logPobs0);
                            } else if (!isZeroSpike && logPobs0 > Double.NEGATIVE_INFINITY) {
                                gamma.setAlpha(spikeShape * nSpikes);
                                gamma.setBeta(1.0 / spikeShape);
                                double gammaLogP = gamma.logDensity(branchSpike);
                                if (Double.isFinite(gammaLogP))
                                    branchP += pk * Math.exp(gammaLogP + logPobs0);
                            }
                        }

                        // obsEvent == 1: the observed speciation IS of type i
                        {
                            int nSpikes = isFakeParent ? k : k + 1;
                            double logPobs1 = logPi;
                            if (nSpikes == 0) {
                                if (isZeroSpike && logPobs1 > Double.NEGATIVE_INFINITY)
                                    branchP += pk * Math.exp(logPobs1);
                            } else if (!isZeroSpike && logPobs1 > Double.NEGATIVE_INFINITY) {
                                gamma.setAlpha(spikeShape * nSpikes);
                                gamma.setBeta(1.0 / spikeShape);
                                double gammaLogP = gamma.logDensity(branchSpike);
                                if (Double.isFinite(gammaLogP))
                                    branchP += pk * Math.exp(gammaLogP + logPobs1);
                            }
                        }

                        k++;
                        logFactorialK += Math.log(k);
                    }
                    logP += Math.log(branchP);
                } else {
                    // No hidden events expected
                    if (!isFakeParent) {
                        double branchP = 0.0;

                        // obsEvent == 0: zero spike, weight = (1 - pi)
                        if (isZeroSpike && logOneMinusPi > Double.NEGATIVE_INFINITY)
                            branchP += Math.exp(logOneMinusPi);

                        // obsEvent == 1: nonzero spike, weight = pi
                        if (!isZeroSpike && logPi > Double.NEGATIVE_INFINITY) {
                            gamma.setAlpha(spikeShape);
                            gamma.setBeta(1.0 / spikeShape);
                            double gammaLogP = gamma.logDensity(branchSpike);
                            if (Double.isFinite(gammaLogP))
                                branchP += Math.exp(gammaLogP + logPi);
                        }

                        logP += Math.log(branchP);
                    } else if (!isZeroSpike) {
                        logP += Double.NEGATIVE_INFINITY;
                    }
                }
            }
        }

        // Numerical issue
        if (logP == Double.POSITIVE_INFINITY) logP = Double.NEGATIVE_INFINITY;
        return logP;
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
        if (parameterizationInput.get() != null) conds.add(parameterizationInput.get().getID());
        if (spikeShapeInput.get() != null) conds.add(spikeShapeInput.get().getID());
        return conds;
    }


    @Override
    public void sample(State state, Random random) {

        if (sampledFlag) return;
        sampledFlag = true;
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

        double spikeShape = spikeShapeInput.get().getValue();
        spikesInput.get().setDimension(nodeCount);

        if (spikeShape <= 0) {
            throw new IllegalArgumentException("Cannot sample spikes because spikeShape is non-positive " + spikeShape);
        }

        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {

            Node node = treeInput.get().getNode(nodeNr);

            // Handle origin branch and sampled ancestor branch
            if (node.isRoot() || node.isDirectAncestor()) {
                spikesInput.get().setValue(nodeNr, 0.0);
                continue;
            }

            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);
            int nHiddenEvents = (int) Randomizer.nextPoisson(expNrHiddenEvents);
            int nSpikes = node.getParent().isFake() ? nHiddenEvents : nHiddenEvents + 1;
            double alpha = spikeShape * nSpikes;

            // Sample spike from Gamma distribution if nSpikes > 0
            // Uses spikeShape instead of 1/spikeShape due to different parameterisation of the Gamma distribution
            double spike = (nSpikes == 0) ? 0.0 : Randomizer.nextGamma(alpha, spikeShape);
            spikesInput.get().setValue(nodeNr, spike);

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
                startTypePriorProbsInput.get().getDoubleValues(), parameterization, finalSampleOffset
        );

        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {

            Node node = treeInput.get().getNode(nodeNr);

            if (node.isRoot() || node.isDirectAncestor()) {
                // Zero spikes for root and direct ancestors
                for (int i = 0; i < nTypes; i++) {
                    spikesInput.get().setValue(nodeNr * nTypes + i, 0.0);
                }
                continue;
            }

            // Compute expected number of hidden speciation events for this branch for each types
            double[] expNrHiddenEventsArray = hiddenEventsIntegrator.getExpNrHiddenEventsForNode(nodeNr);

            // Compute π at time of the observed speciation event of the node, π(t₀)
            Node parent = node.getParent();
            int parentNr = parent.getNr();
            double[] piArray = hiddenEventsIntegrator.getPiAtNode(parentNr);


            for (int i = 0; i < nTypes; i++) {

                double expNrHiddenEvents = expNrHiddenEventsArray[i];
                double pi = Math.max(piArray[i], 0.0);

                int obsEvent = (Randomizer.nextDouble() < pi) ? 1 : 0;

                int nHiddenEvents = (int) Randomizer.nextPoisson(expNrHiddenEvents);
                int nSpikes = node.getParent().isFake() ? nHiddenEvents : nHiddenEvents + obsEvent;
                double spikeShape = getSpikeShape(i);
                double alpha = spikeShape * nSpikes;

                double spike = (nSpikes == 0) ? 0.0 : Randomizer.nextGamma(alpha, spikeShape);
                spikesInput.get().setValue(nodeNr * nTypes + i, spike);

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

    public double getPiVals(int nodeNr, int type) {
        return piVals[nodeNr * nTypes + type];
    }

    public double getSpikeShape(int type) {
        int spikeShapeDim = spikeShapeInput.get().getDimension();
        if (nTypes == 1 || spikeShapeDim == 1) return spikeShapeInput.get().getArrayValue(0);
        else return spikeShapeInput.get().getArrayValue(type);
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

        hiddenEventsCached = storedHiddenEventsCached;
        super.restore();
    }

}