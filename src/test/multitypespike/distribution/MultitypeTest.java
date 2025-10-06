package test.multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.tree.*;
import beast.base.inference.parameter.RealParameter;
import multitypespike.distribution.BranchSpikePrior;
import multitypespike.distribution.MultiTypeHiddenEventsODEs;
import multitypespike.distribution.PiSystem;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.junit.Test;


import static junit.framework.Assert.assertEquals;


public class MultitypeTest {

//    public static void main(String[] args) {
//
//        String newick = "(t1[&state=1] : 1.5, t2[&state=1] : 0.5);";
//        Tree tree = new TreeParser(newick, false,false,true,0);
//        RealParameter origin = new RealParameter("2.5");
//        RealParameter startTypePriorProbs = new RealParameter("1.0");
//
//        Parameterization parameterization = new CanonicalParameterization();
//        parameterization.initByName(
//                "typeSet", new TypeSet(1),
//                "processLength", origin,
//                "birthRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("2.0"), 1),
//                "deathRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("1.0"), 1),
//                "birthRateAmongDemes", new SkylineMatrixParameter(
//                        null,
//                        new RealParameter("0.0"), 1),
//                "migrationRate", new SkylineMatrixParameter(
//                        null,
//                        new RealParameter("0"), 1),
//                "samplingRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.5"), 1),
//                "removalProb", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.0"), 1));
//
//        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
//
//        density.initByName(
//                "parameterization", parameterization,
//                "startTypePriorProbs", startTypePriorProbs,
//                "conditionOnSurvival", false,
//                "tree", tree,
//                "typeLabel", "state",
//                "parallelize", false,
//                "useAnalyticalSingleTypeSolution", false,
//                "storeIntegrationResults", true);
//
//
//        density.calculateLogP();  // Calculate LogP to call integration method
//        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();
//        // nodeNr 0 = t1
//        // nodeNr 1 = t2
//        // nodeNr 2 = root
//
//        int nodeNr = 0;
//        Node node = tree.getNode(nodeNr);
////        double nodeTime = parameterization.getNodeTime(node, 0);
////        double end =  node.isRoot() ? parameterization.getNodeTime(node, 0): parameterization.getNodeTime(node.getParent(), 0);
//
////        System.out.println("nodeTime = "+ nodeTime);
////        System.out.println("endTime = "+ end);
//
////        int steps = 50;
////        ContinuousOutputModel model = p0geResults[nodeNr];
////        System.out.println(model);
//
//        PiSystem piSystem = new PiSystem(parameterization, tree, p0geResults, 1e-100, 1e-7);
//        piSystem.integratePi(startTypePriorProbs.getDoubleValues(), parameterization, 0.0);
//
//        BranchSpikePrior bsp = new BranchSpikePrior();
//        bsp.initByName("parameterization", parameterization,
//                "tree", tree,
//                "spikeShape", "1.0",
//                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
//                "startTypePriorProbs", startTypePriorProbs,
//                "bdmDistr", density);
//
//
//        MultiTypeHiddenEventsODEs multitypeHiddenEvents = new MultiTypeHiddenEventsODEs(nodeNr, parameterization,
//                p0geResults, piSystem.getIntegrationResultsForNode(nodeNr) ,1e-8,1e-8);
//        multitypeHiddenEvents.integrateForBranch(node, parameterization, 0);
//        System.out.println("multi-type method = " + multitypeHiddenEvents.getExpectedHiddenEvents()[0]);
//
//        System.out.println("single-type method = " + bsp.getExpNrHiddenEventsForBranch(node));
//
//
//    }

    @Test
    public void piSystemTest() {

        String newick = "(t5[&type=0]:5.7,((t1[&type=0]:1,t2[&type=0]:2):1,(t3[&type=1]:3,t4[&type=1]:4):0.5):1.3):0.0;";
        Tree tree = new TreeParser(newick, false, false, true, 0);

        RealParameter origin = new RealParameter("6.0");
        RealParameter startTypePriorProbs = new RealParameter("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, new RealParameter("1.2 1.2"), 2),
                "deathRate", new SkylineVectorParameter(null, new RealParameter("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(null, new RealParameter("0.1 0.1"), 2),
                "samplingRate", new SkylineVectorParameter(null, new RealParameter("0.1"), 2),
                "removalProb", new SkylineVectorParameter(null, new RealParameter("1.0"), 2)
        );

        // Integrate p0ge system ===
        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", startTypePriorProbs,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true
        );

        density.calculateLogP();
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();
        PiSystem piSystem = new PiSystem(parameterization, tree, p0geResults, 1e-100, 1e-8);
        piSystem.integratePi(startTypePriorProbs.getDoubleValues(), parameterization, 0.0);


        //  Empirically computed pi values from BDMM-Prime stochastic mapping method
        double[][] expected = new double[][] {
                {0.9494, 0.0506}, // node 5
                {0.4749, 0.5251}, // node 6
                {0.6949, 0.3051}, // node 7
        };
        int[] nodeNumbers = new int[] {5, 6, 7};

        // Compare interpolated pi values
        double tolerance = 0.01; // numerical tolerance

        for (int i = 0; i < nodeNumbers.length; i++) {
            int nodeNr = nodeNumbers[i];
            Node node = tree.getNode(nodeNr);
            ContinuousOutputModel com = piSystem.integrationResults[nodeNr];

            com.setInterpolatedTime(parameterization.getNodeTime(node, 0.0));
            double[] pi = com.getInterpolatedState();

            System.out.printf("Node %d -> π₀=%.4f π₁=%.4f (expected %.4f %.4f)%n",
                    nodeNr, pi[0], pi[1], expected[i][0], expected[i][1]);

            assertEquals("π₀ mismatch at node " + nodeNr, expected[i][0], pi[0], tolerance);
            assertEquals("π₁ mismatch at node " + nodeNr, expected[i][1], pi[1], tolerance);
        }
    }

    @Test
    public void hiddenEventsODEAnalyticalSingleTypeTest() {

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
                "birthRateAmongDemes", new SkylineMatrixParameter(null, new RealParameter("0.0"), 1),
                "migrationRate", new SkylineMatrixParameter(null, new RealParameter("0"), 1),
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
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        PiSystem piSystem = new PiSystem(parameterization, tree, p0geResults, 1e-100, 1e-7);
        piSystem.integratePi(startTypePriorProbs.getDoubleValues(), parameterization, 0.0);

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density
        );

        double tolerance = 1e-6;
        for (int nodeNr = 0; nodeNr <= 1; nodeNr++) {
            Node node = tree.getNode(nodeNr);

            MultiTypeHiddenEventsODEs multitypeHiddenEvents = new MultiTypeHiddenEventsODEs(
                    nodeNr, parameterization, p0geResults,
                    piSystem.getIntegrationResultsForNode(nodeNr), 1e-8, 1e-8
            );
            multitypeHiddenEvents.integrateForBranch(node, parameterization, 0);
            double multiTypeResult = multitypeHiddenEvents.getExpectedHiddenEvents()[0];
            double singleTypeResult = bsp.getExpNrHiddenEventsForBranch(node);

            System.out.printf("Node %d: multi-type = %.10f, single-type = %.10f%n",
                    nodeNr, multiTypeResult, singleTypeResult);

            assertEquals("Mismatch at node " + nodeNr + "single-type result = " + singleTypeResult
                    + " does not match multi-type result = " +  multiTypeResult,
                    singleTypeResult, multiTypeResult, tolerance);

        }
    }

    @Test
    public void hiddenEventsODESingleTypeSkylineTest() {

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
                "birthRateAmongDemes", new SkylineMatrixParameter(null, new RealParameter("0.0"), 1),
                "migrationRate", new SkylineMatrixParameter(null, new RealParameter("0.0"), 1),
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
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        PiSystem piSystem = new PiSystem(parameterization, tree, p0geResults, 1e-100, 1e-7);
        piSystem.integratePi(startTypePriorProbs.getDoubleValues(), parameterization, 0.0);

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density
        );

        double tolerance = 1e-6;
        for (int nodeNr = 0; nodeNr <= 1; nodeNr++) {
            Node node = tree.getNode(nodeNr);

            MultiTypeHiddenEventsODEs multitypeHiddenEvents = new MultiTypeHiddenEventsODEs(
                    nodeNr, parameterization, p0geResults,
                    piSystem.getIntegrationResultsForNode(nodeNr), 1e-8, 1e-8
            );
            multitypeHiddenEvents.integrateForBranch(node, parameterization, 0);
            double multiTypeResult = multitypeHiddenEvents.getExpectedHiddenEvents()[0];
            double singleTypeResult = bsp.getExpNrHiddenEventsForBranch(node);

            System.out.printf("Node %d: multi-type = %.10f, single-type = %.10f%n",
                    nodeNr, multiTypeResult, singleTypeResult);

            assertEquals("Mismatch at node " + nodeNr + "single-type result = " + singleTypeResult
                            + " does not match multi-type result = " +  multiTypeResult,
                    singleTypeResult, multiTypeResult, tolerance);

        }
    }
}
