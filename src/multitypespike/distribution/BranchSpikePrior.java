package mutitypespike.distribution;

import bdmmprime.parameterization.*;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import java.util.*;

// Based on <GammaSpikeModel>  Copyright (C) <2024>  <Jordan Douglas>


@Description("Prior distribution on the spike size of a branch, " +
        "dependent on the number of hidden events along that branch")
public class BranchSpikePrior extends Distribution {
    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-prime parameterization object (see BDMM-prime package for available parameterizations)",
            Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "tree input", Input.Validate.REQUIRED);
    final public Input<RealParameter> gammaShapeInput = new Input<>("gammaShape", "shape hyper-parameter for the " +
            "gamma distribution of the spikes", Input.Validate.REQUIRED);
    final public Input<RealParameter> spikesInput = new Input<>("spikes", "vector of spike amplitudes for each branch",
            Input.Validate.REQUIRED);
    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "If provided, the difference in time between the final sample and the end of the BD process.",
            new RealParameter("0.0"));

//  For model selection: to be implemented.
    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator",
            "indicator for presence/absence of spikes", Input.Validate.OPTIONAL);


    public Parameterization parameterization;
    int nTypes;
    double[] intervalEndTimes, A, B;
    private double lambda_i, mu_i, psi_i, t_i, A_i, B_i, finalSampleOffset;
    RealParameter gammaShape;

    
    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        if (nTypes != 1) {
            throw new RuntimeException("Error: Not implemented for models with >1 type yet");
        }

        gammaShape = gammaShapeInput.get();
        intervalEndTimes = parameterization.getIntervalEndTimes();
        finalSampleOffset = finalSampleOffsetInput.get().getArrayValue(0);

        A = new double[parameterization.getTotalIntervalCount()];
        B = new double[parameterization.getTotalIntervalCount()];
        computeConstants(A,B);

        if (spikesInput.get().getDimension() != treeInput.get().getInternalNodeCount() ){
            throw new RuntimeException("Error: Dimension of spikesInput does not match the number of internal nodes");
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



    private double integral_p_i(double t_0, double t_1) {
        double t0 = t_i - t_0;
        double t1 = t_i - t_1;

        return (1.0 / (2.0 * lambda_i)) * ((t0 - t1) * (mu_i + psi_i + lambda_i + A_i) + 2.0 * Math
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
        System.out.println("i=" + i + ", A_i=" + A_i + ", B_i=" + B_i);
        
    }



    public double getExpNrHiddenEventsForBranch(Node node) {
        // Forward in time calculation i.e. Time(parentNode) < Time(node)
        double expNrHiddenEvents = 0;
        int nodeIndex = parameterization.getNodeIntervalIndex(node, finalSampleOffset);
        int parentIndex = parameterization.getNodeIntervalIndex(node.getParent(), finalSampleOffset);
        double t0 = parameterization.getNodeTime(node.getParent(), finalSampleOffset);
        double T = parameterization.getNodeTime(node, finalSampleOffset);
        updateParametersforInterval(parentIndex);

        if (nodeIndex == parentIndex) return(2*lambda_i*integral_p_i(t0, T));

        for (int i = parentIndex; i <= nodeIndex - 1; i++) {
            if (i > parentIndex) updateParametersforInterval(i);

            double t1 = intervalEndTimes[i];
            expNrHiddenEvents += 2*lambda_i*integral_p_i(t0, t1);
            t0 = t1;
        }

        updateParametersforInterval(nodeIndex);
        expNrHiddenEvents += 2*lambda_i*integral_p_i(t0, T);

        return (expNrHiddenEvents);
    }


    // If there are too many hidden events on a branch (e.g. during mixing) then the gamma distribution shape is large, which causes
    // instabilities
    final double MAX_CUM_SUM = 0.999;

    public double calculateLogP() {
        logP = 0;

        // Calculate density of the spike size of each branch, assuming that each spike is drawn from a
        // Gamma(alpha, beta) distribution, integrating across all possible numbers of hidden speciation events.
        for (int nodeNr = 0; nodeNr < treeInput.get().getInternalNodeCount(); nodeNr++) {

            double branchSpike = spikesInput.get().getValue(nodeNr);

            // Integrate over all possible spike amplitude values
            Node node = treeInput.get().getNode(nodeNr);
            double expNrHiddenEvents = getExpNrHiddenEventsForBranch(node);

            if (expNrHiddenEvents > 0) {

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

                    double alpha = gammaShape.getValue() * (k + 1);
                    double beta = 1 / gammaShape.getValue();

                    GammaDistribution gamma = new GammaDistributionImpl(alpha, beta);
                    double gammaLogP = gamma.logDensity(branchSpike);

                    if (gammaLogP == Double.NEGATIVE_INFINITY || Double.isNaN(gammaLogP)) {
                        branchP += 0;
                    } else {
                        branchP += Math.exp(logpk + gammaLogP);
                    }

                    // Integrate across all possible numbers of hidden events
                    k++;

                }

                logP += Math.log(branchP);

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
        conds.add(gammaShapeInput.get().getID());
        if (treeInput.get() != null) conds.add(treeInput.get().getID());
        if (parameterizationInput.get() != null) conds.add(parameterizationInput.get().getID());
        if (gammaShapeInput.get() != null) conds.add(gammaShapeInput.get().getID());
        return conds;
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation() ||
                InputUtil.isDirty(spikesInput) ||
                InputUtil.isDirty(gammaShapeInput);
    }





    /*Testing*/
    public static void main(String[] args) {
     String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";
        TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
        Tree myTree = treeParser;
        //Node[] node = myTree.getNodesAsArray();

        RealParameter originParam = new RealParameter("2.0");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        new RealParameter("0.25"),
                        new RealParameter("0.5 1.0"), 1),
                "deathRate", new SkylineVectorParameter(
                        new RealParameter("0.25"),
                        new RealParameter("0.3 0.5"), 1),
                "samplingRate", new SkylineVectorParameter(
                        new RealParameter("0.25"),
                        new RealParameter("0.1 0.25"), 1),
                "removalProb", new SkylineVectorParameter(
                        new RealParameter("0.25"),
                        new RealParameter("1.0 1.0"), 1),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("0.0")));
//        parameterization.initByName(
//                "typeSet", new TypeSet(1),
//                "processLength", originParam,
//                "birthRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.75"), 1),
//                "deathRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.3"), 1),
//                "samplingRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.1"), 1),
//                "removalProb", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0"), 1),
//                "rhoSampling", new TimedParameter(
//                        originParam,
//                        new RealParameter("1.0"))
//        );


        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization, "tree", myTree, "gammaShape", "1.0", "spikes", "1.0 0.5 0.1");

        double[] intervalEndTimes = parameterization.getIntervalEndTimes();
        Node node = myTree.getNode(5);

        System.out.println(bsp.getExpNrHiddenEventsForBranch(node));

    }
}

