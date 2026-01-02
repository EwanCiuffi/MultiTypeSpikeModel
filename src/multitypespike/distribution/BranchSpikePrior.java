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
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.special.Gamma;

import java.util.*;


/**
 * Based on <GammaSpikeModel>  Copyright (C) <2025>  <Jordan Douglas>
 */


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

    public Input<Boolean> useAnalyticalSingleTypeSolutionInput = new Input<>("useAnalyticalSingleTypeSolution",
            "Use the analytical branch spike prior when the model has only one type.",
            true);


    private Parameterization parameterization;
    private double[] intervalEndTimes, A, B, expectedHiddenEvents, piVals;
    private double lambda_i, mu_i, psi_i, t_i, A_i, B_i, finalSampleOffset;
    public int nodeCount, nTypes, spikeShapeDim;


    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        nodeCount = treeInput.get().getNodeCount();

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
            piVals = new double[nTypes * nodeCount];
        }

        if (nTypes == 1) {
            spikesInput.get().setDimension(nodeCount);
        } else {
            spikesInput.get().setDimension(nodeCount * nTypes);
        }

        // Ensure spike values are positive for all root and direct-ancestor nodes
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);

            if (!node.isDirectAncestor() && !node.isRoot()) {
                continue;
            }

            for (int i = 0; i < nTypes; i++) {
                int index = nodeNr * nTypes + i;
                if (spikesInput.get().getValue(index)==0) spikesInput.get().setValue(index, 0.5);
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
     *
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
        if (!useAnalyticalSingleTypeSolutionInput.get() && nTypes == 1) return multiTypeCalculateLogP();
        else if (nTypes == 1) return singleTypeCalculateLogP();
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
        GammaDistribution gamma;

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
                gamma = new GammaDistributionImpl(spikeShape, 1 / spikeShape);
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
                    // Probability of k hidden events P(k) under a Poisson(mu)
                    double logpk = -expNrHiddenEvents + k * Math.log(expNrHiddenEvents) - Gamma.logGamma(k + 1);
                    cumsum += Math.exp(logpk);

                    // Number of spikes is k + 1 unless parent of the node is a sampled ancestor (fake), in which case it is k
                    int nSpikes = node.getParent().isFake() ? k : k + 1;

                    if (nSpikes == 0) {
                        // Valid zero spike
                        if (branchSpike == 0) branchP += Math.exp(logpk);

                    } else {
                        // Compute log-probability of observed spike under Gamma distribution
                        gamma = new GammaDistributionImpl(
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
                if (!node.getParent().isFake()) {
                    gamma = new GammaDistributionImpl(spikeShape, 1 / spikeShape);
                    logP += gamma.logDensity(branchSpike);
                } else if (branchSpike != 0.0) {
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
        double[] spikeShapeArray = spikeShapeInput.get().getDoubleValues();
        GammaDistribution gamma;

        // Check spikeShape is positive
        if (Arrays.stream(spikeShapeArray).anyMatch(x -> x <= 0.0)) {
            return Double.NEGATIVE_INFINITY;
        }

        // Integrate hidden events along each branch of the tree
        MultiTypeHiddenEventsIntegrator hiddenEventsIntegrator = new MultiTypeHiddenEventsIntegrator(
                parameterization,treeInput.get(), bdmDistrInput.get().getIntegrationResults(),
                1e-8, 1e-8, false
        );
        hiddenEventsIntegrator.integrateHiddenEvents(
                startTypePriorProbsInput.get().getDoubleValues(), parameterization, finalSampleOffset
        );

        // Loop over all nodes in the tree
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {

            Node node = treeInput.get().getNode(nodeNr);


            // Handle origin branch and sampled ancestor branches
            if (node.isRoot() || node.isDirectAncestor()) {
                // Scaled spikes = 0 for origin branch and sampled ancestor branches
                // Set a pseudo-prior for spikes when they are not included in the model
                // This facilitates transitions between models of different dimensions
                for (int i = 0; i < nTypes; i++) {
                    expectedHiddenEvents[nodeNr * nTypes + i] = 0;
                    piVals[nodeNr * nTypes + i] = 0;

                    double branchSpike = spikesInput.get().getValue(nodeNr * nTypes + i);
                    double spikeShape = getSpikeShape(spikeShapeArray, i);
                    gamma = new GammaDistributionImpl(spikeShape, 1 / spikeShape);
                    logP += gamma.logDensity(branchSpike);
                }
                continue;
            }

            // Compute expected number of hidden speciation events for this branch for all types
            double[] expNrHiddenEventsArray = hiddenEventsIntegrator.getExpNrHiddenEventsForNode(nodeNr);

            // Compute π at time of the observed speciation event of the node, π(t₀)
            double[] piArray = hiddenEventsIntegrator.getPiAtNode(nodeNr);

            // Integrate over types
            for (int i = 0; i < nTypes; i++) {
                double branchSpike = spikesInput.get().getValue(nodeNr * nTypes + i);
                double spikeShape = getSpikeShape(spikeShapeArray, i);
                double expNrHiddenEvents = expNrHiddenEventsArray[i];
                double pi = Math.min(Math.max(piArray[i], 0.0), 1.0);

                expectedHiddenEvents[nodeNr * nTypes + i] = expNrHiddenEvents;
                piVals[nodeNr * nTypes + i] = pi;

                if (expNrHiddenEvents > 0) {
                    double branchP = 0.0;
                    double cumsum = 0.0;
                    int k = 0;

                    while (cumsum < MAX_CUM_SUM) {
                        // Probability of k hidden events P(k) of type i under a Poisson(mu)
                        double logpk = -expNrHiddenEvents + k*Math.log(expNrHiddenEvents) - Gamma.logGamma(k+1);
                        cumsum += Math.exp(logpk);

                        // Integrate over type of the observed speciation event
                        for (int obsEvent = 0; obsEvent <= 1; obsEvent++) {
                            double pObs = obsEvent == 1 ? pi : (1 - pi);

                            // Number of spikes is k + obsEvent unless parent of the node is a sampled ancestor (fake), in which case it is k
                            int nSpikes = node.getParent().isFake() ? k : k + obsEvent;

                            if (nSpikes == 0) {

                                // Valid zero spike
                                if (branchSpike == 0) branchP += Math.exp(logpk + Math.log(pObs));

                            } else {

                                // Compute log-probability of observed spike under Gamma distribution
                                gamma = new GammaDistributionImpl(spikeShape * nSpikes, 1 / spikeShape);
                                double gammaLogP = gamma.logDensity(branchSpike);

                                if (branchSpike != 0 && Double.isFinite(gammaLogP)) {
                                    branchP += Math.exp(logpk + gammaLogP + Math.log(pObs));
                                }
                            }
                        }
                        k++;
                    }

                    // Add log-likelihood for this branch
                    logP += Math.log(branchP);

                } else {
                    if (!node.getParent().isFake()) {
                        double branchP = 0.0;
                        for (int obsEvent = 0; obsEvent <= 1; obsEvent++) {
                            double pObs = obsEvent == 1 ? pi : (1 - pi);

                            // Valid zero spike
                            if (branchSpike == 0 && obsEvent == 0) {
                                branchP += pObs;
                                continue;
                            }

                            // Valid non-zero spike
                            if (obsEvent == 1 && branchSpike != 0) {
                                gamma = new GammaDistributionImpl(spikeShape, 1 / spikeShape);
                                branchP += Math.exp(gamma.logDensity(branchSpike) + Math.log(pObs));
                            }
                        }
                        logP +=  Math.log(branchP);
                    } else if (branchSpike != 0.0) {
                        logP += Double.NEGATIVE_INFINITY;
                    }
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


    @Override
    public void sample(State state, Random random) {

        if (sampledFlag) return;
        sampledFlag = true;
        // Cause conditional parameters to be sampled
        sampleConditions(state, random);

        // Single-type case
        if (nTypes == 1) {

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

        // Multi-type case
        } else {
            // Call calculate LogP to get p0ge integration results
            bdmDistrInput.get().calculateLogP();

            double[] spikeShapeArray = spikeShapeInput.get().getDoubleValues();
            spikesInput.get().setDimension(nodeCount * nTypes);

            MultiTypeHiddenEventsIntegrator hiddenEventsIntegrator = new MultiTypeHiddenEventsIntegrator(
                    parameterization, treeInput.get(), bdmDistrInput.get().getIntegrationResults(),
                    1e-8, 1e-8, false
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
                double[] piArray = hiddenEventsIntegrator.getPiAtNode(nodeNr);


                for (int i = 0; i < nTypes; i++) {

                    double expNrHiddenEvents = expNrHiddenEventsArray[i];
                    double spikeShape = getSpikeShape(spikeShapeArray, i);
                    double pi = Math.max(piArray[i], 0.0);

                    int obsEvent = (Randomizer.nextDouble() < pi) ? 1 : 0;

                    int nHiddenEvents = (int) Randomizer.nextPoisson(expNrHiddenEvents);
                    int nSpikes = node.getParent().isFake() ? nHiddenEvents : nHiddenEvents + obsEvent;
                    double alpha = spikeShape * nSpikes;

                    double spike = (nSpikes == 0) ? 0.0 : Randomizer.nextGamma(alpha, spikeShape);
                    spikesInput.get().setValue(nodeNr * nTypes + i, spike);

                }
            }
        }
    }


    /**
     * Methods for passing expectedNrHiddenEvents, pi values to loggers
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




    @Override
    protected boolean requiresRecalculation() {
        return InputUtil.isDirty(spikesInput) ||
                InputUtil.isDirty(spikeShapeInput) ||
                InputUtil.isDirty(treeInput) ||
                InputUtil.isDirty(parameterizationInput) ||
                InputUtil.isDirty(bdmDistrInput) ||
                InputUtil.isDirty(startTypePriorProbsInput) ||
                InputUtil.isDirty(finalSampleOffsetInput);
    }

}

