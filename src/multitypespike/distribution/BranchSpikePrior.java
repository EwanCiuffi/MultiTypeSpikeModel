package multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.*;


// Based on <GammaSpikeModel>  Copyright (C) <2025>  <Jordan Douglas>


@Description("Prior distribution on the spike size of a branch, " +
        "dependent on the number of hidden events along that branch")
public class BranchSpikePrior extends Distribution {

    final public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-prime parameterization object (see BDMM-prime package for available parameterizations)",
            Input.Validate.REQUIRED);

    final public Input<Tree> treeInput = new Input<>("tree", "tree input", Input.Validate.REQUIRED);

    final public Input<RealParameter> spikeShapeInput = new Input<>("spikeShape", "shape parameter for the " +
            "gamma distribution of the spikes", Input.Validate.REQUIRED);

    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree",
            Input.Validate.REQUIRED);

    final public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "if provided, the difference in time between the final sample and the end of the BD process",
            new RealParameter("0.0"));

    final public Input<RealParameter> startTypePriorProbsInput = new Input<>("startTypePriorProbs",
            "The prior probabilities for the initial individual type",
            Input.Validate.OPTIONAL);

    final public Input<BirthDeathMigrationDistribution> bdmDistrInput = new Input<>("bdmDistr",
            "Birth-death-migration model distribution.", Input.Validate.OPTIONAL);


    private Parameterization parameterization;
    private double[] intervalEndTimes, A, B, expectedHiddenEvents;
    private double lambda_i, mu_i, psi_i, t_i, A_i, B_i, finalSampleOffset;
    public int nodeCount, nTypes, spikeShapeDim;

    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        nodeCount = treeInput.get().getNodeCount();

        // Model not implemented for non-zero births among demes. Rates among demes parameter should not be provided.

        if (nTypes != 1) {
            if (startTypePriorProbsInput.get() == null) {
                throw new IllegalArgumentException("'startTypePriorProbs' must be specified for multi-type analyses.");
            }

            if (bdmDistrInput.get() == null) {
                throw new IllegalArgumentException("BirthDeathMigrationDistribution,'bdmDistr', must be specified for multi-type analyses.");
            }
        }

        // Spike shape dimension checks
        spikeShapeDim = spikeShapeInput.get().getDimension();
        if (nTypes == 1 && spikeShapeDim > 1) {
            throw new IllegalArgumentException("Single-type model requires exactly one spikeShape parameter.");
        }
        if (nTypes > 1 && spikeShapeDim != 1 && spikeShapeDim != nTypes) {
            throw new IllegalArgumentException("For multi-type models, 'spikeShape' must have dimension 1 (shared) or nTypes (" + nTypes + ").");
        }

        A = new double[parameterization.getTotalIntervalCount()];
        B = new double[parameterization.getTotalIntervalCount()];

        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetInput.get().getArrayValue(0);
        computeConstants(A, B);

        if(nTypes==1) {
            expectedHiddenEvents = new double[nodeCount];
        } else {
            expectedHiddenEvents = new double[nTypes * nodeCount];
        }

        if (nTypes == 1) {
            spikesInput.get().setDimension(nodeCount);
        } else {
            spikesInput.get().setDimension(nodeCount * nTypes);
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

    public double getExpectedHiddenEvents(int nodeNr){
        return expectedHiddenEvents[nodeNr];
    }

    public double getExpectedHiddenEvents(int nodeNr, int type) {
        if (nTypes == 1) return expectedHiddenEvents[nodeNr];
        return expectedHiddenEvents[nodeNr * nTypes + type];
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

    // Single type expected number of hidden events for interval (t0,t1)
    private double integral_2lambda_i_p_i(double t_0, double t_1) {
        double t0 = t_i - t_0;
        double t1 = t_i - t_1;

        return ((t0 - t1) * (mu_i + psi_i + lambda_i + A_i) + 2.0 * Math
                .log(((-B_i - 1) * Math.exp(A_i * t1) + B_i - 1)
                        / ((-B_i - 1) * Math.exp(A_i * t0) + B_i - 1)));
    }


    // Single type expected number of hidden events for branch
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
        if (nTypes == 1) return singleTypeCalculateLogP();
        else return multiTypeCalculateLogP();
    }

    // If there are too many hidden events on a branch (e.g. during mixing) then the gamma distribution shape is large, which causes
    // instabilities
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

        // Loop over all nodes in the tree
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);
            double branchSpike = spikesInput.get().getValue(nodeNr);

            // Handle origin branch and sampled ancestor branches
            if (node.isRoot() || node.isDirectAncestor()) {
                expectedHiddenEvents[nodeNr] = 0.0;
                // Scaled spikes = 0 for origin branch and sampled ancestor branches
                // Set a pseudo-prior for spikes when they are not included in the model
                // This facilitates transitions between models of different dimensions
                GammaDistribution gamma = new GammaDistributionImpl(spikeShape, 1 / spikeShape);
                logP += gamma.logDensity(branchSpike);
                continue;
            }

            // Compute expected number of hidden speciation events for this branch
            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);
            expectedHiddenEvents[nodeNr] = expNrHiddenEvents;

            if (expNrHiddenEvents > 0) {
                // Integrate over all possible spike values
                double branchP = 0.0;
                double cumsum = 0.0;
                int k = 0;

                while (cumsum < MAX_CUM_SUM) {
                    // Probability of observing k hidden events P(k) under a Poisson(mu)
                    double logpk = -expNrHiddenEvents + k * Math.log(expNrHiddenEvents);
                    for (int n = 2; n <= k; n++) logpk -= Math.log(n); // - log(k!)
                    cumsum += Math.exp(logpk);

                    // Number of spikes is k + 1 unless parent of the node is a sampled ancestor (fake), in which case it is k
                    int nSpikes = node.getParent().isFake() ? k : k + 1;

                    if (nSpikes == 0) {
                        // Valid zero spike
                        if (branchSpike == 0) branchP += Math.exp(logpk);

                    } else {
                        // Compute log-probability of observed spike under Gamma distribution
                        GammaDistribution gamma = new GammaDistributionImpl(
                                spikeShape * nSpikes, 1/spikeShape);
                        double gammaLogP = gamma.logDensity(branchSpike);
                        if (branchSpike != 0 && Double.isFinite(gammaLogP)) {
                            branchP += Math.exp(logpk + gammaLogP);
                        }
                    }
                    k++;
                }
                // Add log-likelihood for this branch
                logP += Math.log(branchP);

            } else {
                throw new RuntimeException("Expected number of hidden events is zero for non-sampled ancestor branch");
            }
        }

        // Numerical issue
        if (logP == Double.POSITIVE_INFINITY) logP = Double.NEGATIVE_INFINITY;

        return logP;
    }


    public double[] getMultiTypeExpForBranch(Node node, ContinuousOutputModel piValuesForBranch) {
        MultiTypeHiddenEvents multiTypeODEs = new MultiTypeHiddenEvents(node.getNr(),
                parameterization, bdmDistrInput.get().getIntegrationResults(),
                piValuesForBranch, 1e-8, 1e-8
                );
        multiTypeODEs.integrateForBranch(node, parameterization, finalSampleOffset);
        return multiTypeODEs.getExpectedHiddenEvents();
    }


    public double multiTypeCalculateLogP() {
        logP = 0.0;
        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetInput.get().getArrayValue(0);
        double[] spikeShapeArray = spikeShapeInput.get().getDoubleValues();

        // Check spikeShape is positive
        if (Arrays.stream(spikeShapeArray).anyMatch(x -> x <= 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        computeConstants(A, B);

        PiSystem piSystem = new PiSystem(parameterization, treeInput.get(),
                bdmDistrInput.get().getIntegrationResults(),1e-8, 1e-8);
        piSystem.integratePi(startTypePriorProbsInput.get().getDoubleValues(), parameterization, finalSampleOffset);

        // Loop over all nodes in the tree
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {

            Node node = treeInput.get().getNode(nodeNr);

            // Handle origin branch and sampled ancestor branches
            if (node.isRoot() || node.isDirectAncestor()) {
                Arrays.fill(expectedHiddenEvents, nodeNr * nTypes, nodeNr * nTypes + nTypes, 0.0);

                // Scaled spikes = 0 for origin branch and sampled ancestor branches
                // Set a pseudo-prior for spikes when they are not included in the model
                // This facilitates transitions between models of different dimensions
                for (int i = 0; i < nTypes; i++) {
                    double spikeShape = getSpikeShape(spikeShapeArray, i);
                    GammaDistribution gamma = new GammaDistributionImpl(spikeShape, 1 / spikeShape);
                    logP += gamma.logDensity(spikesInput.get().getValue(nodeNr + nodeCount * i));
                }
                continue;
            }

            ContinuousOutputModel piCom = piSystem.getIntegrationResultsForNode(node.getNr());

            // Compute expected number of hidden speciation events for this branch for all types
            double[] expNrHiddenEventsArray = getMultiTypeExpForBranch(node, piCom);
            System.arraycopy(expNrHiddenEventsArray, 0, expectedHiddenEvents, nodeNr * nTypes, nTypes);

            // Get pi at time of the observed speciation event of the node
            piCom.setInterpolatedTime(parameterization.getNodeTime(node.getParent(), finalSampleOffset));
            double[] piArray = piCom.getInterpolatedState();

            // Integrate over all possible spike values
            for (int i = 0; i < nTypes; i++) {
                double branchSpike = spikesInput.get().getValue(nodeNr + nodeCount * i);
                double spikeShape = getSpikeShape(spikeShapeArray, i);
                double expNrHiddenEvents = expNrHiddenEventsArray[i];

                if (expNrHiddenEvents > 0) {
                    double branchP = 0.0;
                    double cumsum = 0.0;
                    int k = 0;

                    while (cumsum < MAX_CUM_SUM) {
                        // Probability of observing k hidden events P(k) of type i under a Poisson(mu)
                        double logpk = -expNrHiddenEvents + k * Math.log(expNrHiddenEvents);
                        for (int n = 2; n <= k; n++) logpk -= Math.log(n); // - log(k!)
                        cumsum += Math.exp(logpk);

                        // Number of spikes is k + π(t₀) unless parent of the node is a sampled ancestor (fake), in which case it is k
                        double nSpikes = node.getParent().isFake() ? k : k + piArray[i];

                        if (nSpikes == 0) {
                            // Valid zero spike
                            if (branchSpike == 0) branchP += Math.exp(logpk);

                        } else {
                            // Compute log-probability of observed spike under Gamma distribution
                            GammaDistribution gamma = new GammaDistributionImpl(
                                    spikeShape * nSpikes, 1 / spikeShape);
                            double gammaLogP = gamma.logDensity(branchSpike);
                            if (branchSpike != 0 && Double.isFinite(gammaLogP)) {
                                branchP += Math.exp(logpk + gammaLogP);
                            }
                        }
                        k++;
                    }
                    // Add log-likelihood for this branch
                    logP += Math.log(branchP);

                } else {
                    throw new RuntimeException("Expected number of hidden events is zero for non-sampled ancestor branch");
                }
            }
        }
        // Numerical issue
        if (logP == Double.POSITIVE_INFINITY) logP = Double.NEGATIVE_INFINITY;

        return logP;
    }


    private double getSpikeShape(double[] spikeShapeArray, int type) {
        if (nTypes == 1 || spikeShapeDim == 1) return spikeShapeArray[0];
        else return spikeShapeArray[type];
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

    // Samples spikes for single-type models
    @Override
    public void sample(State state, Random random) {

        if (sampledFlag) return;
        sampledFlag = true;

        if (nTypes > 1) {
            throw new UnsupportedOperationException("Sampling not implemented for multi-type spikes yet");
        }

        // Cause conditional parameters to be sampled
        sampleConditions(state, random);

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


    @Override
    protected boolean requiresRecalculation() {
        return InputUtil.isDirty(spikesInput) ||
                InputUtil.isDirty(spikeShapeInput) ||
                InputUtil.isDirty(treeInput) ||
                InputUtil.isDirty(parameterizationInput) ||
                InputUtil.isDirty(finalSampleOffsetInput);
    }


    public static void main(String[] args) {

        String newick = "(t1[&state=0] : 1.0, t2[&state=1] : 1.0);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealParameter origin = new RealParameter("2.0");
        RealParameter startTypePriorProbs = new RealParameter("0.5 0.5");


        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("3.0 3.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.5"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.3"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "rhoSampling", new TimedParameter(
                        origin,
                        new RealParameter("0.2"), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();

        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", startTypePriorProbs,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true
        );

        density.calculateLogP();  // Calculate LogP to call integration method

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density
        );

        System.out.println(bsp.multiTypeCalculateLogP());

    }
}

