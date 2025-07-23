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
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

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


    public Parameterization parameterization;
    int nTypes;
    double[] intervalEndTimes, A, B;
    public RealParameter expectedHiddenEvents = new RealParameter("0.0");
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

        if (spikesInput.get().getDimension() != treeInput.get().getNodeCount()) {
            spikesInput.get().setDimension(treeInput.get().getNodeCount());
        }

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

        A_i = A[i];
        B_i = B[i];
    }



    public double getExpNrHiddenEventsForBranch(Node node) {
        if (node.isRoot() || node.isDirectAncestor()) return 0;

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

    @Override
    public double calculateLogP() {
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
                    for (int i = 2; i <= k; i++) logpk -= Math.log(i); // - log(k!)
                    double pk = Math.exp(logpk);
                    cumsum += pk;

                    // Number of spikes is k+1 unless parent of the node is a sampled ancestor (fake), in which case it's just k
                    int nSpikes = node.getParent().isFake() ? k : k + 1;
                    double alpha = spikeShape * nSpikes;


                    double gammaLogP;
                    if (alpha == 0) {
                        // Spike must be zero if no speciation events
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
                // If branch length is 0, spike must be zero
                logP += (branchSpike == 0.0) ? 0.0 : Double.NEGATIVE_INFINITY;
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

        for (int nodeNr = 0; nodeNr < dimension; nodeNr++) {

            Node node = treeInput.get().getNode(nodeNr);

            // Handle origin branch and sampled ancestor branch
            if (node.isRoot() || node.isDirectAncestor()) {
                spikesInput.get().setValue(nodeNr, 0.0);
                continue;
            }

            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);
            long NrHiddenEvents = Randomizer.nextPoisson(expNrHiddenEvents);

            // Determine alpha based on whether node is a true observed speciation event
            double alpha;
            alpha = (node.isFake()) ? spikeShape*NrHiddenEvents : spikeShape*(NrHiddenEvents + 1);

            // Sample spike from Gamma distribution if alpha > 0
            double spike;
            spike = (alpha == 0) ? 0.0: Randomizer.nextGamma(alpha, spikeShape);

            spikesInput.get().setValue(nodeNr, spike);
        }
    }



    @Override
    protected boolean requiresRecalculation() {
        if( InputUtil.isDirty(spikesInput) ||
            InputUtil.isDirty(spikeShapeInput) ||
            InputUtil.isDirty(treeInput) ||
            InputUtil.isDirty(parameterizationInput) ||
            InputUtil.isDirty(finalSampleOffsetInput)){

            return true;
        }

        return false;

    }



}

