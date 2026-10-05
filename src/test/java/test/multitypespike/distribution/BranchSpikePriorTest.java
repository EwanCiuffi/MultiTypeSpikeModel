package test.multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.tree.*;
import beast.base.spec.evolution.tree.coalescent.ConstantPopulation;
import beast.base.spec.evolution.tree.coalescent.RandomTree;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.SimplexParam;
import test.multitypespike.Params;
import beast.base.util.Randomizer;
import multitypespike.distribution.BranchSpikePrior;
import multitypespike.logger.HiddenEventsLogger;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class BranchSpikePriorTest {


    @Test
    public void BSPMultiTypeSingleTypeComparisonTest() {

        String newick = "(t1[&state=1]:1.5, t2[&state=1]:0.5);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("2.5");
        SimplexParam startTypePriorProbs = Params.simplex("1.0");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, Params.real("2.0"), 1),
                "deathRate", new SkylineVectorParameter(null, Params.real("1.0"), 1),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.5"), 1),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0"), 1)
        );

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
        density.calculateLogP();

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", false,
                "bdmDistr", density,
                "initializeSpikes", false,
                "parallelize", false
        );

        System.out.println("Before multiTypeResult"); double multiTypeResult = bsp.calculateLogP(); System.out.println("After multiTypeResult");

        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", true,
                "bdmDistr", density,
                "initializeSpikes", false
        );
        double singleTypeResult = bsp.calculateLogP();

        double tolerance = 1e-3;

        assertEquals(singleTypeResult, multiTypeResult, tolerance, "Mismatch:"  + "single-type result = " + singleTypeResult
                        + " does not match multi-type result = " +  multiTypeResult);

    }


    @Test
    public void BSPMultiTypeSingleTypeComparisonSkylineTest() {

        String newick = "(t1[&state=1]:1.5, t2[&state=1]:0.5);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("2.5");
        SimplexParam startTypePriorProbs = Params.simplex("1.0");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(Params.real("1.0"), Params.real("2.0 0.5"), 1),
                "deathRate", new SkylineVectorParameter(Params.real("1.5"), Params.real("1.0 0.2"), 1),
                "samplingRate", new SkylineVectorParameter(Params.real("2.0"), Params.real("0.5 1.8"), 1),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0"), 1)
        );

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
        density.calculateLogP();

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", false,
                "bdmDistr", density,
                "initializeSpikes", false
        );

        System.out.println("Before multiTypeResult"); double multiTypeResult = bsp.calculateLogP(); System.out.println("After multiTypeResult");

        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", true,
                "bdmDistr", density,
                "initializeSpikes", false
        );
        double singleTypeResult = bsp.calculateLogP();

        double tolerance = 1e-3;

        assertEquals(singleTypeResult, multiTypeResult, tolerance, "Mismatch:"  + "single-type result = " + singleTypeResult
                        + " does not match multi-type result = " +  multiTypeResult);

    }

    /**
     * Single-type test for expected number of hidden events on a branch with no rate shifts
     * Compares with simulations of birth-death trajectories performed in R
     */
    @Test
    public void singleTypeExpectedHiddenEventsTest() {
        String newick = "(t1:1.0, t2:1.0);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("2.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, Params.real("3.0"), 1),
                "deathRate", new SkylineVectorParameter(null, Params.real("0.5"), 1),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.0"), 1),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0"), 1),
                "rhoSampling", new TimedParameter(Params.asVector(origin), Params.real("0.2"))
        );

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true
        );

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
                "bdmDistr", density
        );

        // Get expected hidden events for node 0 (t1)
        Node node = tree.getNode(0);
        double expectedHiddenEvents = bsp.getExpNrHiddenEventsForBranch(node);

        // Compare to R simulation result (hiddenEventsSim.R)
        double expectedFromSimulation = 3.39054;
        double tolerance = 5e-3;

        assertEquals(expectedFromSimulation, expectedHiddenEvents, tolerance, "Expected hidden events should match simulation");
    }

    /**
     * Single-type test for expected number of hidden events on a branch with no rate shifts
     * Compares with the output of GammaSpike Model - 01/03/2025
     */
    @Test
    public void noRateShiftsSingleTypeCaseTest() {

        String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";
        TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
        Tree tree = treeParser;

        RealScalarParam<NonNegativeReal> originParam = Params.scalar("2.0");
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.75"), 1),
                "deathRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.3"), 1),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.1"), 1),
                "removalProb", new SkylineVectorParameter(
                        null,
                        Params.real("0"), 1),
                "rhoSampling", new TimedParameter(Params.asVector(originParam),
                        Params.real("1.0"))
        );
        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization, "tree", tree, "spikeShape", Params.positive("1.0"), "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"));
        Node node = tree.getNode(5);

        /*
        Stubs stub = new Stubs();
        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();
        StumpedTreePrior stp = new StumpedTreePrior();
        stp.initByName("lambda", "0.75", "r0", "2.5", "samplingProportion", "0.25", "tree", tree);
        stub.initByName("tree", tree, "prior", stp);
        gamma_bsp.initByName("spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"), "shape", "1.0", "stubs", stub, "tree", tree);

        stp.getMeanStubNumber(node.getHeight(), node.getParent().getHeight())  = 0.186082405828208.
        gamma_bsp.calculateLogP() = -3.4385926164063445.
        */

        assertEquals(0.186082405828208, bsp.getExpNrHiddenEventsForBranch(node), 1e-5);
//        assertEquals(-3.4385926164063445, bsp.calculateLogP(), 1e-5); // MTS and GS logPs do not match due to differences in handling of the root and sampled ancestors

    }


    /**
     * Test to check the logic of the branch traversal process used for calculating the
     * expected number of hidden events along a branch (by computing the lengths of branches instead of
     * the expected number of hidden events)
     */
    @Test
    public void branchTraversalTest() {
        RealScalarParam<NonNegativeReal> originParam = Params.scalar("2.0");
        Parameterization parameterization = new CanonicalParameterization();

        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        Params.real("0.1 0.2 0.3 0.4 0.5 0.6"),
                        Params.real("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "deathRate", new SkylineVectorParameter(
                        Params.real("0.7 0.8 0.9 1.0 1.1 1.2"),
                        Params.real("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "samplingRate", new SkylineVectorParameter(
                        Params.real("1.3 1.4 1.5 1.6 1.7 1.8"),
                        Params.real("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "removalProb", new SkylineVectorParameter(
                        Params.real("1.9 2.0 2.1 2.2 2.5 3.0"),
                        Params.real("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "rhoSampling", new TimedParameter(Params.asVector(originParam),
                        Params.real("0.0")));

        // Code to generate random trees for testing

        List<Sequence> seqList = new ArrayList<Sequence>();
        int Nleaves = 10;
        int Niter = 10000;

        for (int i = 0; i < Nleaves; i++) {
            String taxonID = "t " + i;
            seqList.add(new Sequence(taxonID, "?"));
        }
        Alignment alignment = new Alignment(seqList, "nucleotide");

        // Population model
        ConstantPopulation populationModel = new ConstantPopulation();
        populationModel.initByName("popSize", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));

        // Create RandomTree
        RandomTree randomTree = new RandomTree();
        for (int k = 0; k < Niter; k++) {
            randomTree.initByName("taxa", alignment, "populationModel", populationModel);

            for (int n = 1; n < randomTree.getNodeCount() - 1; n++) {
                Node node = randomTree.getNode(n);

                double trueBranchTime = parameterization.getNodeTime(node, 0)
                        - parameterization.getNodeTime(node.getParent(), 0);

                double[] intervalEndTimes = parameterization.getIntervalEndTimes();
                double branchTime = 0;
                int nodeIndex = parameterization.getNodeIntervalIndex(node, 0);
                int parentIndex = parameterization.getNodeIntervalIndex(node.getParent(), 0);
                double t0 = parameterization.getNodeTime(node.getParent(), 0);
                double T = parameterization.getNodeTime(node, 0);

                if (nodeIndex == parentIndex) branchTime = (T - t0);
                else {
                    for (int i = parentIndex; i <= nodeIndex - 1; i++) {

                        double t1 = intervalEndTimes[i];
                        branchTime += (t1 - t0);
                        t0 = t1;
                    }

                    branchTime += (T - t0);
                }

                assertEquals(branchTime, trueBranchTime, 1e-9);
            }
        }
    }

    /**
     * Test to validate the calculation of the single-type branch spike density.
     * We override getExpNrHiddenEventsForBranch to condition on known expected hidden events,
     * allowing direct comparison against exact values calculated in R.
     */
    @Test
    public void exactSingleTypeDensityTest() throws Exception {
        // Simple 2-taxon tree
        String newick = "(t1:1.0, t2:1.0);";
        Tree tree = new TreeParser(newick, false, false, true, 0);

        BranchSpikePrior bsp = new BranchSpikePrior() {
            @Override
            public double getExpNrHiddenEventsForBranch(Node node) {
                // Hardcode expectations
                if (node.getNr() == 0) return 0.5; // t1
                if (node.getNr() == 1) return 1.2; // t2
                return 0.0; // root/internal nodes
            }
        };

        // Because we overrode getExpNrHiddenEventsForBranch, the actual rates here do not matter.
        Parameterization param = new CanonicalParameterization();
        param.initByName("processLength", Params.scalar("2.0"),
                "birthRate", new SkylineVectorParameter(null, Params.real("2.0"), 1),
                "deathRate", new SkylineVectorParameter(null, Params.real("1.0"), 1),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.5"), 1),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0"), 1)
        );


        // Initialise the prior
        bsp.initByName(
                "parameterization", param,
                "tree", tree,
                "spikeShape", Params.positive("2.0"),
                "spikes", Params.real("0.8 2.5 0.5"), // 0.8 for node 0, 2.5 for node 1, 0.5 for root (pseudo-prior)
                "initializeSpikes", false
        );

        double logP = bsp.calculateLogP();

        // This expected value comes directly from the total_logP output of generate_true_logdensity.R
        double expectedLogP = -2.76114731;

        double tolerance = 1e-5;
        assertEquals(expectedLogP, logP, tolerance, "logP does not match exact R calculation.");
    }

    /**
     * Test to validate the calculation of the multi-type branch spike density.
     * We override getExpNrHiddenEventsForBranch to condition on known expected hidden events,
     * allowing direct comparison against exact values calculated in R.
     */
    @Test
    public void exactMultiTypeDensityTest() throws Exception {
        String newick = "(t1[&state=0]:1.0, t2[&state=1]:1.0);";
        Tree tree = new TreeParser(newick, false, false, true, 0);

        BranchSpikePrior bsp = new BranchSpikePrior();

        // Parameterization for 2 types
        Parameterization param = new CanonicalParameterization();
        param.initByName(
                "typeSet", new TypeSet(2),
                "processLength", Params.scalar("2.0"),
                "birthRate", new SkylineVectorParameter(null, Params.real("2.0 2.0"), 2),
                "deathRate", new SkylineVectorParameter(null, Params.real("1.0 1.0"), 2),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.5 0.5"), 2),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0 0.0"), 2)
        );

        // Dummy BDM distribution to satisfy initAndValidate() requirements
        BirthDeathMigrationDistribution bdm = new BirthDeathMigrationDistribution();
        bdm.initByName("parameterization", param, "tree", tree, "typeLabel", "state",
                "startTypePriorProbs", Params.simplex("0.5 0.5"), "storeIntegrationResults", true);

        // Initialise the prior
        bsp.initByName(
                "parameterization", param,
                "tree", tree,
                "bdmDistr", bdm,
                "startTypePriorProbs", Params.simplex("0.5 0.5"),
                "spikeShape", Params.positive("2.0"),
                // Spikes dimension = 3 nodes * 2 types = 6 values
                // Node0_Type0, Node0_Type1, Node1_Type0, Node1_Type1, Root_Type0, Root_Type1
                "spikes", Params.real("0.8 0.3 2.5 1.1 0.5 0.5"),
                "initializeSpikes", false,
                "useAnalyticalSingleTypeSolution", false // Force multi-type execution
        );

        // Force the recalculation flag to false so the integrator is completely skipped
        Field reqReintField = BranchSpikePrior.class.getDeclaredField("requiresReintegration");
        reqReintField.setAccessible(true);
        reqReintField.set(bsp, false);

        // Inject expected hidden events
        Field expHiddenEventsField = BranchSpikePrior.class.getDeclaredField("expectedHiddenEvents");
        expHiddenEventsField.setAccessible(true);
        double[] expHidden = (double[]) expHiddenEventsField.get(bsp);

        // Node 0 (t1)
        expHidden[0] = 0.5; // Type 0
        expHidden[1] = 0.2; // Type 1
        // Node 1 (t2)
        expHidden[2] = 1.2; // Type 0
        expHidden[3] = 0.9; // Type 1
        // Node 2 (Root) stays 0.0

        // Inject pi values
        Field piValsField = BranchSpikePrior.class.getDeclaredField("piVals");
        piValsField.setAccessible(true);
        double[] piVals = (double[]) piValsField.get(bsp);

        // Node 0 (t1)
        piVals[0] = 0.8; // Type 0
        piVals[1] = 0.6; // Type 1
        // Node 1 (t2)
        piVals[2] = 0.7; // Type 0
        piVals[3] = 0.4; // Type 1
        // Node 2 (Root) stays 0.0

        double logP = bsp.calculateLogP();

        // This expected value comes directly from the total_logP output of generate_true_multitype_logdensity.R
        double expectedLogP = -6.07916977;

        double tolerance = 1e-5;
        assertEquals(expectedLogP, logP, tolerance, "Multi-type logP does not match exact R calculation.");
    }

    /**
     * Test to verify that the multi-type model correctly treats the observed speciation
     * event as mutually exclusive across types.
     * * In the generative model, an observed speciation event belongs to exactly one type.
     * Therefore, it is impossible for all types on a branch to simultaneously have
     * 0 observed events (which would result in a 0.0 spike if there are no hidden events).
     * * This test has 0.0 spikes for all types on all branches.
     * Under an incorrect assumption of statistical independence between types,
     * the model would allow the observed event to "vanish" (evaluate to 0 for all types)
     * and return a finite logP.
     * Under the correct mutually exclusive joint probability model, it must return NEGATIVE_INFINITY.
     */
    @Test
    public void multiTypeMutuallyExclusiveObservedEventTest() {
        String newick = "(t1[&state=0]:1.0, t2[&state=1]:1.0);";
        Tree tree = new TreeParser(newick, false, false, true, 0);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", Params.scalar("2.0"),
                "birthRate", new SkylineVectorParameter(null, Params.real("2.0 2.0"), 2),
                "deathRate", new SkylineVectorParameter(null, Params.real("1.0 1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(null, Params.real("0.5 0.5"), 2),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.5 0.5"), 2),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0 0.0"), 2)
        );

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", Params.simplex("0.5 0.5"),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "storeIntegrationResults", true
        );
        density.calculateLogP();

        // Scenario 1: All zero spikes
        BranchSpikePrior bspZeroes = new BranchSpikePrior();
        bspZeroes.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("0.0 0.0 0.0 0.0 0.0 0.0"),
                "startTypePriorProbs", Params.simplex("0.5 0.5"),
                "bdmDistr", density,
                "initializeSpikes", false,
                "useAnalyticalSingleTypeSolution", false
        );

        double logPZeroes = bspZeroes.calculateLogP();

        assertEquals(Double.NEGATIVE_INFINITY, logPZeroes, 1e-10, "An observed speciation event must occur on exactly one type. " +
                        "A configuration with 0.0 spikes across all types should be impossible.");

        // Scenario 2: Exactly one type has a spike per branch
        BranchSpikePrior bspValid = new BranchSpikePrior();
        bspValid.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.0 0.0 1.5 0.5 0.0"),
                "startTypePriorProbs", Params.simplex("0.5 0.5"),
                "bdmDistr", density,
                "initializeSpikes", false,
                "useAnalyticalSingleTypeSolution", false
        );

        double logPValid = bspValid.calculateLogP();

        // Assert that assigning the event to exactly one type yields a valid, finite logP
        assertFalse(Double.isInfinite(logPValid) || Double.isNaN(logPValid), "Assigning the observed speciation event to exactly one type should yield a valid finite logP.");
    }

    private static double bruteForceLogProb(double spike, double mu, double shape, int nObs) {
        double sum = 0.0;
        for (int k = 0; k < 2000; k++) {
            int n = k + nObs;
            double logPk = -mu + k * Math.log(mu) - org.apache.commons.math3.special.Gamma.logGamma(k + 1);
            double logDensity = spike == 0.0 ? (n == 0 ? 0.0 : Double.NEGATIVE_INFINITY)
                    : (n == 0 ? Double.NEGATIVE_INFINITY : BranchSpikePrior.logSpikeDensity(spike, n, shape));
            sum += Math.exp(logPk + logDensity);
        }
        return Math.log(sum);
    }

    // Sum over hidden events must match brute force, including large spikes on branches with few expected events
    @Test
    public void branchSpikeLogProbsTest() {
        double[][] cases = {{1.0, 0.05, 7.4}, {3.0, 0.05, 7.4}, {5.0, 0.3, 7.4}, {5.0, 1.0, 7.4},
                {0.2, 20.0, 2.0}, {40.0, 3.0, 1.5}, {0.0, 0.7, 5.0}};
        double[] logProbs = new double[2];
        for (double[] c : cases) {
            BranchSpikePrior.branchSpikeLogProbs(c[0], c[1], c[2], logProbs);
            assertEquals(bruteForceLogProb(c[0], c[1], c[2], 0), logProbs[0], 1e-9);
            assertEquals(bruteForceLogProb(c[0], c[1], c[2], 1), logProbs[1], 1e-9);
        }

        // No hidden events expected
        BranchSpikePrior.branchSpikeLogProbs(2.0, 0.0, 3.0, logProbs);
        assertEquals(Double.NEGATIVE_INFINITY, logProbs[0]);
        assertEquals(BranchSpikePrior.logSpikeDensity(2.0, 1, 3.0), logProbs[1], 1e-12);
        BranchSpikePrior.branchSpikeLogProbs(0.0, 0.0, 3.0, logProbs);
        assertEquals(0.0, logProbs[0], 1e-12);
        assertEquals(Double.NEGATIVE_INFINITY, logProbs[1]);
    }

    // Sampled hidden events must follow their conditional distribution given the spike
    @Test
    public void hiddenEventSamplerTest() {
        Randomizer.setSeed(7);
        double spike = 3.0, mu = 0.8, shape = 4.0;
        double[] logProbs = new double[2];
        BranchSpikePrior.branchSpikeLogProbs(spike, mu, shape, logProbs);

        for (int nObs = 0; nObs <= 1; nObs++) {
            double exactMean = 0.0;
            double logTotal = logProbs[nObs];
            for (int k = 0; k < 200; k++) {
                double logPk = -mu + k * Math.log(mu) - org.apache.commons.math3.special.Gamma.logGamma(k + 1);
                if (k + nObs > 0)
                    exactMean += k * Math.exp(logPk + BranchSpikePrior.logSpikeDensity(spike, k + nObs, shape) - logTotal);
            }

            int n = 200000;
            double sampleMean = 0.0;
            for (int i = 0; i < n; i++)
                sampleMean += HiddenEventsLogger.sampleHiddenEventCount(spike, mu, shape, nObs, logTotal);
            sampleMean /= n;

            assertEquals(exactMean, sampleMean, 0.01, "Hidden event sampler mean (nObs = " + nObs + ")");
        }
    }


    /** A tree prior left on its defaults is switched to storing its integration results. */
    @Test
    public void storedResultsEnabledAutomaticallyTest() {
        Tree tree = new TreeParser("(t1[&state=1]:1.5, t2[&state=0]:0.5);", false, false, true, 0);
        Parameterization param = new CanonicalParameterization();
        param.initByName(
                "typeSet", new TypeSet(2),
                "processLength", Params.scalar("2.5"),
                "birthRate", new SkylineVectorParameter(null, Params.real("2.0 1.5"), 2),
                "deathRate", new SkylineVectorParameter(null, Params.real("1.0 0.5"), 2),
                "migrationRate", new SkylineMatrixParameter(null, Params.real("0.2 0.3"), 2),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.5 0.5"), 2),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0"), 2));

        BirthDeathMigrationDistribution bdm = new BirthDeathMigrationDistribution();
        bdm.initByName("parameterization", param, "tree", tree, "typeLabel", "state",
                "startTypePriorProbs", Params.simplex("0.5 0.5"));

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("tree", tree, "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
                "bdmDistr", bdm, "initializeSpikes", false, "parallelize", false);

        bdm.calculateLogP();
        assertTrue(bdm.saveIntegrationResultsInput.get());
        assertTrue(Double.isFinite(bsp.calculateLogP()));
    }
}
