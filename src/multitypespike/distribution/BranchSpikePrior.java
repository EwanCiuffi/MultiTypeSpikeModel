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
import org.apache.commons.math.distribution.ExponentialDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import java.util.*;


// Based on <GammaSpikeModel>  Copyright (C) <2025>  <Jordan Douglas>


@Description("Prior distribution on the spike size of a branch, " +
        "dependent on the number of hidden events along that branch")
public class BranchSpikePrior extends Distribution {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-prime parameterization object (see BDMM-prime package for available parameterizations)",
            Input.Validate.REQUIRED);

    final public Input<Tree> treeInput = new Input<>("tree", "tree input", Input.Validate.REQUIRED);

    final public Input<RealParameter> spikeShapeInput = new Input<>("spikeShape", "shape parameter for the " +
            "gamma distribution of the spikes", Input.Validate.REQUIRED);

    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree",
            Input.Validate.REQUIRED);

    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "if provided, the difference in time between the final sample and the end of the BD process",
            new RealParameter("0.0"));

    public Input<RealParameter> startTypePriorProbsInput = new Input<>("startTypePriorProbs",
            "The prior probabilities for the initial individual type",
            Input.Validate.OPTIONAL);

    public Input<BirthDeathMigrationDistribution> bdmDistrInput = new Input<>("bdmDistr",
            "Birth-death-migration model distribution.", Input.Validate.OPTIONAL);


    public Parameterization parameterization;
    public int nTypes;
    private double[] intervalEndTimes, A, B;
    private double lambda_i, mu_i, psi_i, t_i, A_i, B_i, finalSampleOffset;
    private double[][] lambda_ij;
    public RealParameter expectedHiddenEvents = new RealParameter("0.0");
    private PiState piState;


    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        if (nTypes != 1) {
//            throw new RuntimeException("Error: Not implemented for models with >1 type yet");
            if (startTypePriorProbsInput == null || startTypePriorProbsInput.get() == null) {
                throw new IllegalArgumentException("'startTypePriorProbs' must be specified for multi-type analyses.");
            }

            if (bdmDistrInput == null || bdmDistrInput.get() == null) {
                throw new IllegalArgumentException("BirthDeathMigrationDistribution,'bdmDistr', must be specified for multi-type analyses.");
            }
        }

        A = new double[parameterization.getTotalIntervalCount()];
        B = new double[parameterization.getTotalIntervalCount()];

        computeConstants(A, B);
        spikesInput.get().setDimension(treeInput.get().getNodeCount());
        expectedHiddenEvents.setDimension(treeInput.get().getNodeCount());
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
        lambda_ij = parameterization.getCrossBirthRates()[i];
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

        for (int k = parentIndex; k <= nodeIndex - 1; k++) {
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
        computeConstants(A, B);

        // Check spikeShape is positive
        double spikeShape = spikeShapeInput.get().getValue();
        if (spikeShape <= 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // Loop over all nodes in the tree
        for (int nodeNr = 0; nodeNr < treeInput.get().getNodeCount(); nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);
            double branchSpike = spikesInput.get().getValue(nodeNr);

            // Handle origin branch and sampled ancestor branches
            if (node.isRoot() || node.isDirectAncestor()) {
                expectedHiddenEvents.setValue(nodeNr, 0.0);
                // Scaled spikes = 0 for origin branch and sampled ancestor branches
                // Set a pseudo-prior for spikes when they are not included in the model
                // This facilitates transitions between models of different dimensions
                GammaDistribution gamma = new GammaDistributionImpl(spikeShape, 1 / spikeShape);
                logP += gamma.logDensity(branchSpike);
                continue;
            }

            // Compute expected number of hidden speciation events for this branch
            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);
            expectedHiddenEvents.setValue(nodeNr, expNrHiddenEvents);

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
                        GammaDistribution gamma = new GammaDistributionImpl(spikeShape * nSpikes, 1/spikeShape);
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




//    public double typedExpNrHiddenEventsForInterval(int i, int nodeNr, double t0, double t1) {
//        int nSteps = 20;
//        double dt = (t1 - t0) / nSteps;
//        double expNrHiddenEvents = 0.0;
//
//        ContinuousOutputModel p0geCom = bdmmInput.get().getCOMs()[nodeNr];
//
//        for (int step = 0; step < nSteps; step++) {
//            double tLeft = t0 + step * dt;
//            double tRight = tLeft + dt;
//            double tMid = (tLeft + tRight) / 2.0;
//
//            p0geCom.setInterpolatedTime(tMid);
//            double[] p0 = p0geCom.getInterpolatedState();
//            double pi_i = piState.getPiAtTime(tMid)[i];
//
//            expNrHiddenEvents += 2 * pi_i * lambda_i * p0[i] * dt;
//
//            for (int j = 0; j < nTypes; j++) {
//                if(i==j) continue;
//                expNrHiddenEvents += pi_i * lambda_ij[i][j] * p0[j] * dt;
//            }
//        }
//
//        return expNrHiddenEvents;
//    }


    public double typedExpNrHiddenEventsForInterval(int i, int nodeNr, double t0, double t1) {

        TypedHiddenEventsODE ode = new TypedHiddenEventsODE(i, nodeNr,
                lambda_i, lambda_ij, nTypes, piState, bdmDistrInput.get().getIntegrationResults());

        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                1e-6,     // min step
                1.0,      // max step
                1e-8,    // abs tolerance
                1e-8     // rel tolerance
        );

        double[] y = new double[] {0.0}; // initial value
        integrator.integrate(ode, t0, y, t1, y);

        return y[0];
    }


    public double getTypedExpNrHiddenEventsForBranch(Node node, int i) {
        if (node.isRoot() || node.isDirectAncestor()) return 0;

        double expNrHiddenEvents = 0;
        int nodeNr = node.getNr();
        int nodeIndex = parameterization.getNodeIntervalIndex(node, finalSampleOffset);
        int parentIndex = parameterization.getNodeIntervalIndex(node.getParent(), finalSampleOffset);
        double t0 = parameterization.getNodeTime(node.getParent(), finalSampleOffset);
        double T = parameterization.getNodeTime(node, finalSampleOffset);
        updateParametersForInterval(parentIndex);

        if (nodeIndex == parentIndex) return typedExpNrHiddenEventsForInterval(i, nodeNr, t0, T);

        for (int k = parentIndex; k <= nodeIndex - 1; k++) {
            if (k > parentIndex) updateParametersForInterval(k);
            double t1 = intervalEndTimes[k];
            expNrHiddenEvents += typedExpNrHiddenEventsForInterval(i, nodeNr, t0, t1);
            t0 = t1;
        }

        updateParametersForInterval(nodeIndex);
        expNrHiddenEvents += typedExpNrHiddenEventsForInterval(i, nodeNr, t0, T);

        return expNrHiddenEvents;
    }


    private void integratePi(){
        this.piState = new PiState(parameterization.getNTypes());
        PiSystem piSystem = new PiSystem(
                parameterization,
                treeInput.get(),
                bdmDistrInput.get().getIntegrationResults(),
                1e-8,
                1e-8
        );
        piSystem.integratePi(treeInput.get(), piState, startTypePriorProbsInput.get().getDoubleValues(),
        parameterization, finalSampleOffset);
    }

    // To be updated //
    public double multiTypeCalculateLogP(){
        logP = 0.0;

        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetInput.get().getArrayValue(0);
        computeConstants(A, B);

        this.integratePi();

        // Check spikeShape is positive
        double spikeShape = spikeShapeInput.get().getValue();
        if (spikeShape <= 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // Loop over all nodes in the tree
        for (int nodeNr = 0; nodeNr < treeInput.get().getNodeCount(); nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);
            double branchSpike = spikesInput.get().getValue(nodeNr);

            // Handle origin branch and sample ancestor branches
            if (node.isRoot() || node.isDirectAncestor()) {
                expectedHiddenEvents.setValue(nodeNr, 0.0);

                // Scaled spikes = 0 for origin branch and sampled ancestor branches
                // Set a pseudo-prior for spikes when they are not included in the model
                // This facilitates transitions between models of different dimensions

//                GammaDistribution gamma = new GammaDistributionImpl(spikeShape, 1.0 / spikeShape);
//                logP += gamma.logDensity(branchSpike);

                ExponentialDistribution exponential = new ExponentialDistributionImpl(0.005);
                logP += exponential.logDensity(branchSpike);
                continue;
            }


            for (int i = 0; nodeNr < nTypes; i++) {

                // Compute expected number of hidden speciation events for this branch
                double expNrHiddenEventsForTypeI = getTypedExpNrHiddenEventsForBranch(node, i);
//                expectedHiddenEvents.setValue(nodeNr, expNrHiddenEventsForTypeI);

                if (expNrHiddenEventsForTypeI > 0) {
                    // Integrate over all possible spike values
                    double branchP = 0.0;
                    double cumsum = 0.0;
                    int k = 0;

                    while (cumsum < MAX_CUM_SUM) {
                        // Probability of observing k hidden events P(k) under a Poisson(mu)
                        double logpk = -expNrHiddenEventsForTypeI + k * Math.log(expNrHiddenEventsForTypeI);
                        for (int n = 2; n <= k; n++) logpk -= Math.log(n); // - log(k!)
                        double pk = Math.exp(logpk);
                        cumsum += pk;

                        // Number of spikes is k+1 unless parent of the node is a sampled ancestor (fake), in which case it is k
                        int nSpikes = node.getParent().isFake() ? k : k + 1;
                        double alpha = spikeShape * nSpikes;


                        double gammaLogP;
                        if (nSpikes == 0) {
                            // Spike must be zero if no speciation events (should not be needed as sampled ancestors already dealt with)
                            gammaLogP = (branchSpike == 0) ? 0.0 : Double.NEGATIVE_INFINITY;
                        } else {
                            // Compute log-probability of observed spike under Gamma distribution
                            GammaDistribution gamma = new GammaDistributionImpl(alpha, 1.0 / spikeShape);
                            gammaLogP = gamma.logDensity(branchSpike);
                        }

                        // Only add valid probabilities to the branch likelihood
                        if (Double.isFinite(gammaLogP)) {
                            branchP += Math.exp(logpk + gammaLogP);
                        }

                        k++;
                    }

                    // Add log-likelihood for this branch
                    logP += (branchP == 0.0) ? Double.NEGATIVE_INFINITY : Math.log(branchP);

                } else {
                    // Consistency check: If expectedNrHiddenEvents is zero, spike must be zero
                    System.out.println("Expected number of hidden events is zero for non-sampled ancestor branch");
                    logP += (branchSpike <= 0.0) ? 0.0 : Double.NEGATIVE_INFINITY;
                }
            }
            // Numerical safeguard
            if (Double.isInfinite(logP)) {
                logP = Double.NEGATIVE_INFINITY;
            }
        }

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
        conds.add(spikeShapeInput.get().getID());
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

        int dimension = treeInput.get().getNodeCount();
        double spikeShape = spikeShapeInput.get().getValue();
        spikesInput.get().setDimension(dimension);

        if (spikeShape <= 0) {
            throw new IllegalArgumentException("Cannot sample spikes because spikeShape is non-positive " + spikeShape);
        }

        for (int nodeNr = 0; nodeNr < dimension; nodeNr++) {

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
            // Uses spikeShape instead of 1/spikeShape due to different Gamma parameterisation
            double spike = (nSpikes == 0) ? 0.0: Randomizer.nextGamma(alpha, spikeShape);

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


}

