package multitypespike.distribution;

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

    final public Input<RealParameter> spikesInput = new Input<>("spikes", "vector of spike amplitudes for each branch",
            Input.Validate.REQUIRED);

    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "if provided, the difference in time between the final sample and the end of the BD process",
            new RealParameter("0.0"));



    public Parameterization parameterization;
    int nTypes;
    double[] intervalEndTimes, A, B;
    private double lambda_i, mu_i, psi_i, t_i, A_i, B_i, finalSampleOffset;

    
    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        if (nTypes != 1) {
            throw new RuntimeException("Error: Not implemented for models with >1 type yet");
        }

        A = new double[parameterization.getTotalIntervalCount()];
        B = new double[parameterization.getTotalIntervalCount()];

        if (spikesInput.get().getDimension() != treeInput.get().getNodeCount()){
            throw new RuntimeException("Error: Dimension of spikesInput, " + spikesInput.get().getDimension() +
                    " does not match the number of nodes, " + (treeInput.get().getNodeCount()) );
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



    private double integral_2lambda_i_p_i(double t_0, double t_1) {
        double t0 = t_i - t_0;
        double t1 = t_i - t_1;

        return ((t0 - t1) * (mu_i + psi_i + lambda_i + A_i) + 2.0 * Math
                .log(((-B_i - 1) * Math.exp(A_i * t1) + B_i - 1)
                        / ((-B_i - 1) * Math.exp(A_i * t0) + B_i - 1)));
    }



    private void updateParametersforInterval(int i) {
        lambda_i = parameterization.getBirthRates()[i][0];
        mu_i = parameterization.getDeathRates()[i][0];
        psi_i = parameterization.getSamplingRates()[i][0];
        t_i = parameterization.getIntervalEndTimes()[i];

//        System.out.println("i = " + i);
//        System.out.println("int end times = " + parameterization.getIntervalEndTimes()[i]);
//        System.out.println("t_i = " + t_i);

        A_i = A[i];
        B_i = B[i];
    }



    public double getExpNrHiddenEventsForBranch(Node node) {
        if (node.isRoot()) return 0;
//        System.out.println("node Nr = " + node.getNr());
        // Forward in time calculation i.e. Time(parentNode) < Time(node)
        double expNrHiddenEvents = 0;
        int nodeIndex = parameterization.getNodeIntervalIndex(node, finalSampleOffset);
        int parentIndex = parameterization.getNodeIntervalIndex(node.getParent(), finalSampleOffset);
        double t0 = parameterization.getNodeTime(node.getParent(), finalSampleOffset);
        double T = parameterization.getNodeTime(node, finalSampleOffset);
        updateParametersforInterval(parentIndex);

        if (nodeIndex == parentIndex) return(integral_2lambda_i_p_i(t0, T));

        for (int i = parentIndex; i <= nodeIndex - 1; i++) {
            if (i > parentIndex) updateParametersforInterval(i);

            double t1 = intervalEndTimes[i];
            expNrHiddenEvents += integral_2lambda_i_p_i(t0, t1);
            t0 = t1;
        }

        updateParametersforInterval(nodeIndex);
        expNrHiddenEvents += integral_2lambda_i_p_i(t0, T);

        return (expNrHiddenEvents);
    }


    // If there are too many hidden events on a branch (e.g. during mixing) then the gamma distribution shape is large, which causes
    // instabilities
    final double MAX_CUM_SUM = 0.999;

    public double calculateLogP() {
        logP = 0.0;

        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetInput.get().getArrayValue(0);
        computeConstants(A, B);

        // Check spikeShape is positive
        double spikeShape = spikeShapeInput.get().getValue();
        if (spikeShape <= 0) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        // Calculate density of the spike size of each branch, assuming that each spike is drawn from a
        // Gamma(alpha, beta) distribution, integrating across all possible numbers of hidden speciation events.
        for (int nodeNr = 0; nodeNr < treeInput.get().getNodeCount(); nodeNr++) {
//        for (int nodeNr = 0; nodeNr < treeInput.get().getNodeCount()-1; nodeNr++) {

            Node node = treeInput.get().getNode(nodeNr);
            double branchSpike = spikesInput.get().getValue(nodeNr);

            // skip root node
            if (node.isRoot()) {
//                GammaDistribution gamma = new GammaDistributionImpl(spikeShape, 1 / spikeShape);
//                logP += gamma.logDensity(branchSpike);
                continue;
            }

            // Integrate over all possible spike amplitude values
            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);

            if (expNrHiddenEvents > 0 ) {

                double branchP = 0;
                int k = 0;
                double cumsum = 0;

                while (cumsum < MAX_CUM_SUM) {

                    // Probability of observing k hidden events P(k) under a Poisson(mu),
                    // where mu = Expected number of hidden events
                    double logpk = -expNrHiddenEvents + k * Math.log(expNrHiddenEvents);
                    for (int i = 2; i <= k; i++) logpk -= Math.log(i); // - log(k!)

                    double pk = Math.exp(logpk);
                    cumsum += pk;

                    double alpha;

                    // If parent of node is 'fake' (parent of sampled ancestor),
                    // then do not add + 1 to alpha, since this is not a true observed speciation event
                    if (!node.getParent().isFake()) {
                        alpha = spikeShape * (k + 1);
                    } else {
                        alpha = spikeShape * k;
                    }

                    double gammaLogP;

                    // Case where alpha = 0 (no spike expected)
                    if (alpha == 0) {
                        if (branchSpike != 0) {
                            gammaLogP = Double.NEGATIVE_INFINITY; // invalid
                        } else {
                            gammaLogP = 0.0;
                        }

                    } else {
                        // Compute the log-probability of the observed spike under the Gamma distribution
                        GammaDistribution gamma = new GammaDistributionImpl(alpha, 1 / spikeShape);
                        gammaLogP = gamma.logDensity(branchSpike);
                    }


                    // Only add valid probabilities to the branch likelihood
                    if (gammaLogP != Double.NEGATIVE_INFINITY && !Double.isNaN(gammaLogP)) {
                        branchP += Math.exp(logpk + gammaLogP);
                    }

                    // Integrate across all possible numbers of hidden events
                    k++;

                }
                logP += Math.log(branchP);

            }
            else {
                // Case when the branch length is 0 (e.g. sampled ancestor), spike is zero
                if (branchSpike != 0) logP += Double.NEGATIVE_INFINITY;
            }

            // Numerical issue
            if (logP == Double.POSITIVE_INFINITY) {
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

        for (int nodeNr = 0; nodeNr < dimension; nodeNr++) {

            Node node = treeInput.get().getNode(nodeNr);

            // Skip root node
            if (node.isRoot()) {
                spikesInput.get().setValue(nodeNr, 0.0);
                continue;
            }

            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);
            long NrHiddenEvents = Randomizer.nextPoisson(expNrHiddenEvents);

            // Determine alpha based on whether parent is a true observed speciation event
            double alpha;
            if (!node.getParent().isFake()) {
                alpha = spikeShape * (NrHiddenEvents + 1);
            } else {
                alpha = spikeShape * NrHiddenEvents;
            }

            double spike;
            if (alpha == 0) {
                // No spike expected
                spike = 0.0;
            } else {
                // Sample spike from Gamma distribution
                spike = Randomizer.nextGamma(alpha, spikeShape);
            }

            spikesInput.get().setValue(nodeNr, spike);
        }
    }





    @Override
    protected boolean requiresRecalculation() {

        return true;
//        return super.requiresRecalculation() ||
//                InputUtil.isDirty(spikesInput) ||
//                InputUtil.isDirty(spikeShapeInput);
    }

}

