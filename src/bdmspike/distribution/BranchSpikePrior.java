package bdmspike.distribution;

import bdmmprime.parameterization.*;
import beast.base.core.Description;
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

// Elements borrowed from <GammaSpikeModel>  Copyright (C) <2024>  <Jordan Douglas>

@Description("Prior distribution on the spike size of a branch, " +
        "dependent on the number of hidden events along that branch")
public class BranchSpikePrior extends Distribution {
    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "tree input", Input.Validate.REQUIRED);
    final public Input<RealParameter> gammaShapeInput = new Input<>("gammaShape", "shape hyper-parameter for the " +
            "gamma distribution of the spikes", Input.Validate.REQUIRED);
    final public Input<RealParameter> spikesInput = new Input<>("spikes", "vector of spike amplitudes for each branch",
            Input.Validate.REQUIRED);

    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator",
            "indicator for presence/absence of a spike", Input.Validate.OPTIONAL);

    public Parameterization parameterization;
    int nTypes;
    double[] birthRateChangeTimes, deathRateChangeTimes, samplingRateChangeTimes, rhoSamplingTimes, intervalEndTimes;
    double[][] birthRates, deathRates, samplingRates, rhoValues;

    RealParameter gammaShape;

    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        nTypes = parameterization.getNTypes();
        if (nTypes != 1) {
            throw new RuntimeException("Error: Not implemented for models with >1 type yet.");
        }

        birthRateChangeTimes = parameterization.getBirthRateChangeTimes();
        deathRateChangeTimes = parameterization.getDeathRateChangeTimes();
        samplingRateChangeTimes = parameterization.getSamplingRateChangeTimes();
        rhoSamplingTimes = parameterization.getRhoSamplingTimes();
        birthRates = parameterization.getBirthRates();
        deathRates = parameterization.getDeathRates();
        samplingRates = parameterization.getSamplingRates();
        rhoValues = parameterization.getRhoValues();
        gammaShape = gammaShapeInput.get();
        intervalEndTimes = parameterization.getIntervalEndTimes();

        if (spikesInput.get().getDimension() != treeInput.get().getInternalNodeCount() ){
            throw new RuntimeException("Error: Dimension of spikesInput does not match the number of internal nodes");
        }
    }


    // Calculate the expected number of hidden events along a branch within a given time interval
    private double getC1(double lambda, double mu, double psi, double rho) {
        return Math.sqrt(Math.abs(Math.pow(lambda - mu - psi, 2) + 4 * lambda * psi));
    }

    private double getC2(double lambda, double mu, double psi, double rho) {
        double c1 = getC1(lambda, mu, psi, rho);
        return -(lambda - mu - 2 * lambda * rho - psi) / c1;
    }

    // Equation 1 from
    // Stadler, Tanja. "Sampling-through-time in birth–death trees."
    // Journal of theoretical biology 267.3 (2010): 396-404.
    private double getP0(double height, double lambda, double mu, double psi, double rho,
                         double c1, double c2) throws Exception {
        double a = lambda + mu + psi;
        double exp = Math.exp(-c1 * height) * (1 - c2);
        double b = 1 + c2;
        double top = a + c1 * (exp - b) / (exp + b);
        double bottom = 2 * lambda;
        double result = top / bottom;
        if (Double.isNaN(result) || result <= 0) {
            //Log.warning("c1=" + c1 + " c2=" + c2 + " lambda =" + lambda + " mu=" + mu + " psi=" + psi);
            throw new Exception("Numerical error: p0 is non-positive " + top + "/" + bottom + "=" + result);
        }
        return result;
    }

    // Equation 2 from
    // Stadler, Tanja. "Sampling-through-time in birth–death trees."
    // Journal of theoretical biology 267.3 (2010): 396-404.
    private double getP1(double height, double lambda, double mu, double psi, double rho,
                         double c1, double c2) throws Exception {
        double top = 4 * rho;
        double bottom = 2 * (1 - c2 * c2) + Math.exp(-c1 * height) * (1 - c2) * (1 - c2) + Math.exp(c1 * height) * (1 + c2) * (1 + c2);
        double result = top / bottom;
        if (Double.isNaN(result) || result <= 0) {
            //Log.warning("c1=" + c1 + " c2=" + c2 + " lambda =" + lambda + " mu=" + mu + " psi=" + psi);
            throw new Exception("Numerical error: p1 is non-positive " + top + "/" + bottom + "=" + result);
        }
        return result;
    }


    // Equation 2 from Bokma,Folmer,van den Brink, Valentijn, Stadler, Tanja. "Unexpectedly many extinct hominins."
    // Evolution, Volume 66, Issue 9, 1 September 2012, Pages 2969–2974. https://doi.org/10.1111/j.1558-5646.2012.01660.x
    public double getExpNrHiddenEventsForInterval(double t0, double t1) {
        int index = parameterization.getIntervalIndex(t0);

        double mu = birthRates[index][0];
        double lambda = deathRates[index][0];
        double psi = samplingRates[index][0];
        double rho = rhoValues[index][0];


        double c1 = getC1(lambda, mu, psi, rho);
        double c2 = getC2(lambda, mu, psi, rho);
        double f = ((c2 - 1) * Math.exp(-c1 * t0) - c2 - 1) / ((c2 - 1) * Math.exp(-c1 * t1) - c2 - 1);

        return (t1 - t0) * (lambda + mu + psi - c1) + 2 * Math.log(f);
    }


    public double getExpNrHiddenEventsForBranch(double branchStartTime, double branchEndTime) {

        double ExpNrHiddenEvents = 0.0;

        // Handle the first interval separately if needed, (due to intervalEndTimes not starting at 0)
        if (intervalEndTimes[0] > branchStartTime) {
            double t0 = branchStartTime;
            double t1 = intervalEndTimes[0];
            if (t1 > branchStartTime && t0 < branchEndTime) {
                t0 = Math.max(t0, branchStartTime);
                t1 = Math.min(t1, branchEndTime);
                ExpNrHiddenEvents += getExpNrHiddenEventsForInterval(t0, t1);
            }
        }
        // Loop through the rest of the intervals
        for (int i = 0; i < intervalEndTimes.length - 1; i++) {
            double t0 = intervalEndTimes[i];
            double t1 = intervalEndTimes[i + 1];

            // Check if the interval is within the branch time span
            if (t1 <= branchStartTime || t0 >= branchEndTime) {
                continue; // Skip intervals outside the branch time span
            }

            // Adjust the interval to fit within the branch time span
            t0 = Math.max(t0, branchStartTime);
            t1 = Math.min(t1, branchEndTime);

            ExpNrHiddenEvents += getExpNrHiddenEventsForInterval(t0, t1);
        }
        return ExpNrHiddenEvents;
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
            double branchStartTime = node.getHeight();
            double branchEndTime = node.isRoot() ? branchStartTime : node.getParent().getHeight();
            double ExpNrHiddenEvents = getExpNrHiddenEventsForBranch(branchStartTime, branchEndTime);

            if (ExpNrHiddenEvents > 0) {

                double branchP = 0;
                int k = 0;
                double cumsum = 0;
                while (cumsum < MAX_CUM_SUM) {

                    // Probability of observing k hidden events P(k) under a Poisson(mu),
                    // where mu = Expected number of hidden events
                    double logpk = -ExpNrHiddenEvents + k * Math.log(ExpNrHiddenEvents);
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
        //conds.add(meanInput.get().getID());
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
                        new RealParameter("4.0 5.0"), 1),
                "deathRate", new SkylineVectorParameter(
                        new RealParameter("0.5"),
                        new RealParameter("3.0 2.0"), 1),
                "samplingRate", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("1.5 2.5"), 1),
                "removalProb", new SkylineVectorParameter(
                        new RealParameter("1.5"),
                        new RealParameter("1.0 3.0"), 1),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("0.0")));
//        parameterization.initByName(
//                "typeSet", new TypeSet(1),
//                "processLength", originParam,
//                "birthRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("4.0"), 1),
//                "deathRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("3.0"), 1),
//                "samplingRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("1.5"), 1),
//                "removalProb", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("1.0"), 1),
//                "rhoSampling", new TimedParameter(
//                        originParam,
//                        new RealParameter("0.0"))
//        );

//        int nTypes;
//        double[] birthRateChangeTimes, deathRateChangeTimes, samplingRateChangeTimes, rhoSamplingTimes, intervalEndTimes;
//        double[][] birthRates, deathRates, samplingRates, rhoValues;

//        birthRateChangeTimes = parameterization.getBirthRateChangeTimes();
//        deathRateChangeTimes = parameterization.getDeathRateChangeTimes();
//        samplingRateChangeTimes = parameterization.getSamplingRateChangeTimes();
//        rhoSamplingTimes = parameterization.getRhoSamplingTimes();
//        birthRates = parameterization.getBirthRates();
//        deathRates = parameterization.getDeathRates();
//        samplingRates = parameterization.getSamplingRates();
//        rhoValues = parameterization.getRhoValues();
//        intervalEndTimes = parameterization.getIntervalEndTimes();


        //System.out.println(parameterization.getNTypes());
        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization, "tree", myTree, "gammaShape", "1.0", "spikes", "1.0 1.0 1.0");
//        System.out.println(bsp.getExpNrHiddenEventsForInterval(0.4, 0.5) + bsp.getExpNrHiddenEventsForInterval(0.5, 1.0)
//                + bsp.getExpNrHiddenEventsForInterval(1.0, 1.5) + bsp.getExpNrHiddenEventsForInterval(1.5, 1.79)
//        );
        Node node = myTree.getNode(2);
        System.out.println(myTree.getInternalNodeCount());
        //System.out.println(node.getHeight() + "," + node.getParent().getHeight());

        //System.out.println(bsp.getExpNrHiddenEventsForBranch(node.getHeight(), node.getParent().getHeight()));

        System.out.println(bsp.calculateLogP());
        //System.out.println(Arrays.toString(parameterization.getIntervalEndTimes()));

    }

}

