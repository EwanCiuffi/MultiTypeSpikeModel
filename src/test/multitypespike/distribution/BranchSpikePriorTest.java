package test.multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.tree.*;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.coalescent.RandomTree;
import beast.base.inference.parameter.RealParameter;
import multitypespike.distribution.BranchSpikePrior;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertEquals;

public class BranchSpikePriorTest {


    @Test
    public void BSPMultiTypeSingleTypeComparisonTest() {

        String newick = "(t1[&state=1]:1.5, t2[&state=1]:0.5);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealParameter origin = new RealParameter("2.5");
        RealParameter startTypePriorProbs = new RealParameter("1.0");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, new RealParameter("2.0"), 1),
                "deathRate", new SkylineVectorParameter(null, new RealParameter("1.0"), 1),
                "samplingRate", new SkylineVectorParameter(null, new RealParameter("0.5"), 1),
                "removalProb", new SkylineVectorParameter(null, new RealParameter("0.0"), 1)
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
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", false,
                "bdmDistr", density,
                "initializeSpikes", false
        );

        double multiTypeResult = bsp.calculateLogP();

        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", true,
                "bdmDistr", density,
                "initializeSpikes", false
        );
        double singleTypeResult = bsp.calculateLogP();

        double tolerance = 1e-3;

        assertEquals("Mismatch:"  + "single-type result = " + singleTypeResult
                        + " does not match multi-type result = " +  multiTypeResult,
                singleTypeResult, multiTypeResult, tolerance);

    }

    @Test
    public void BSPMultiTypeSingleTypeComparisonSkylineTest() {

        String newick = "(t1[&state=1]:1.5, t2[&state=1]:0.5);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealParameter origin = new RealParameter("2.5");
        RealParameter startTypePriorProbs = new RealParameter("1.0");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(new RealParameter("1.0"), new RealParameter("2.0 0.5"), 1),
                "deathRate", new SkylineVectorParameter(new RealParameter("1.5"), new RealParameter("1.0 0.2"), 1),
                "samplingRate", new SkylineVectorParameter(new RealParameter("2.0"), new RealParameter("0.5 1.8"), 1),
                "removalProb", new SkylineVectorParameter(null, new RealParameter("0.0"), 1)
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
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", false,
                "bdmDistr", density,
                "initializeSpikes", false
        );

        double multiTypeResult = bsp.calculateLogP();

        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", true,
                "bdmDistr", density,
                "initializeSpikes", false
        );
        double singleTypeResult = bsp.calculateLogP();

        double tolerance = 1e-3;

        assertEquals("Mismatch:"  + "single-type result = " + singleTypeResult
                        + " does not match multi-type result = " +  multiTypeResult,
                singleTypeResult, multiTypeResult, tolerance);

    }

    /**
     * Single-type test for expected number of hidden events on a branch with no rate shifts
     * Compares with simulations of birth-death trajectories performed in R
     */
    @Test
    public void singleTypeExpectedHiddenEventsTest() {
        String newick = "(t1:1.0, t2:1.0);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealParameter origin = new RealParameter("2.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, new RealParameter("3.0"), 1),
                "deathRate", new SkylineVectorParameter(null, new RealParameter("0.5"), 1),
                "samplingRate", new SkylineVectorParameter(null, new RealParameter("0.0"), 1),
                "removalProb", new SkylineVectorParameter(null, new RealParameter("0.0"), 1),
                "rhoSampling", new TimedParameter(origin, new RealParameter("0.2"))
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
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "bdmDistr", density
        );

        // Get expected hidden events for node 0 (t1)
        Node node = tree.getNode(0);
        double expectedHiddenEvents = bsp.getExpNrHiddenEventsForBranch(node);

        // Compare to R simulation result (hiddenEventsSim.R)
        double expectedFromSimulation = 3.39054;
        double tolerance = 5e-3;

        assertEquals("Expected hidden events should match simulation", expectedFromSimulation, expectedHiddenEvents, tolerance);
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

        RealParameter originParam = new RealParameter("2.0");
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.75"), 1),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.3"), 1),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.1"), 1),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0"), 1),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("1.0"))
        );
        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization, "tree", tree, "spikeShape", "1.0", "spikes", "1.0 0.5 0.1 0.2 0.7 0.1");
        Node node = tree.getNode(5);

        /*
        Stubs stub = new Stubs();
        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();
        StumpedTreePrior stp = new StumpedTreePrior();
        stp.initByName("lambda", "0.75", "r0", "2.5", "samplingProportion", "0.25", "tree", tree);
        stub.initByName("tree", tree, "prior", stp);
        gamma_bsp.initByName("spikes", "1.0 0.5 0.1 0.2 0.7 0.1", "shape", "1.0", "stubs", stub, "tree", tree);

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
        RealParameter originParam = new RealParameter("2.0");
        Parameterization parameterization = new CanonicalParameterization();

        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        new RealParameter("0.1 0.2 0.3 0.4 0.5 0.6"),
                        new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "deathRate", new SkylineVectorParameter(
                        new RealParameter("0.7 0.8 0.9 1.0 1.1 1.2"),
                        new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "samplingRate", new SkylineVectorParameter(
                        new RealParameter("1.3 1.4 1.5 1.6 1.7 1.8"),
                        new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "removalProb", new SkylineVectorParameter(
                        new RealParameter("1.9 2.0 2.1 2.2 2.5 3.0"),
                        new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("0.0")));

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
        populationModel.initByName("popSize", new RealParameter("1.0"));

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




//    public static void main(String[] args) {
//
//        String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";
//        TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
//        Tree tree = treeParser;
//
//        RealParameter originParam = new RealParameter("2.0");
//
//        Parameterization parameterization = new CanonicalParameterization();
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
//
//        Node node = tree.getNode(5);
//        BranchSpikePrior bsp = new BranchSpikePrior();
//        bsp.initByName("parameterization", parameterization, "tree", tree, "spikeShape", "1.0", "spikes", "1.0 0.5 0.1 0.2 0.7 0.1");
//        System.out.println(bsp.getExpNrHiddenEventsForBranch(node));
//        System.out.println(bsp.calculateLogP());

//        System.out.println(parameterization.getNodeTime(node,0));
//        System.out.println(parameterization.getNodeTime(node.getParent(),0));
//
//        System.out.println(node.getParent().getHeight());
//        System.out.println(node.getHeight());



//        Stubs stub = new Stubs();
//        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();
//        StumpedTreePrior stp = new StumpedTreePrior();
//        stp.initByName("lambda", "0.75", "r0", "2.5", "samplingProportion", "0.25", "tree", tree);
//        stub.initByName("tree", tree, "prior", stp);
//        gamma_bsp.initByName("spikes", "1.0 0.5 0.1 0.2 0.7 0.1", "shape", "1.0", "stubs", stub, "tree", tree);
//        System.out.println(stp.getMeanStubNumber(node.getHeight(), node.getParent().getHeight()));
//        System.out.println(gamma_bsp.calculateLogP());


//}


}

