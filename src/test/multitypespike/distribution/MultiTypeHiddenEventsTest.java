package test.multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.tree.*;
import beast.base.inference.parameter.RealParameter;
import multitypespike.distribution.BranchSpikePrior;
import multitypespike.distribution.MultiTypeHiddenEventsIntegrator;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.junit.Test;

import java.util.concurrent.Executor;
import java.util.concurrent.ForkJoinPool;

import static junit.framework.Assert.assertEquals;

public class MultiTypeHiddenEventsTest {

    // Test multi-type hidden events expectation calculation for non-zero birth rate among demes
    @Test
    public void multiTypeEventsBirthAmongDemesTest1() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealParameter origin = new RealParameter("1.0");
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
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1.5 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "rhoSampling", new TimedParameter(
                        origin,
                        new RealParameter("0.2 0.0"), 2));

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
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        int nodeNr = 0;

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density);

        Executor pool = ForkJoinPool.commonPool();

        double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

        MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                new MultiTypeHiddenEventsIntegrator(
                        parameterization, tree, p0geResults,
                        1e-6, 1e-6,
                        false, true, pool, weightOfNodeSubTree,0.1
                );

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getDoubleValues(), parameterization,0.0, 1.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 3.39;
        double sim_hiddenEvents_1 = 0.0;
        double tolerance = 1e-2;

        assertEquals("Type 0 hidden events mismatch", sim_hiddenEvents_0, hiddenEvents[0], tolerance);
        assertEquals("Type 1 hidden events mismatch", sim_hiddenEvents_1, hiddenEvents[1], tolerance);
    }


    // Test multi-type hidden events expectation calculation for non-zero birth rate among demes
    @Test
    public void multiTypeEventsBirthAmongDemesTest2() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealParameter origin = new RealParameter("1.0");
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
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1.5 1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.4"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "rhoSampling", new TimedParameter(
                        origin,
                        new RealParameter("0.2 0.0"), 2));

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
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        int nodeNr = 0;


        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density);

        Executor pool = ForkJoinPool.commonPool();
        double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

        MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                new MultiTypeHiddenEventsIntegrator(
                        parameterization, tree, p0geResults,
                        1e-6, 1e-6,
                        false, true, pool, weightOfNodeSubTree, 0.1
                );

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getDoubleValues(), parameterization,0.0, 1.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 2.02553;
        double sim_hiddenEvents_1 = 1.4639;
        double tolerance = 1e-2;

        assertEquals("Type 0 hidden events mismatch", sim_hiddenEvents_0, hiddenEvents[0], tolerance);
        assertEquals("Type 1 hidden events mismatch", sim_hiddenEvents_1, hiddenEvents[1], tolerance);
    }


    // Test multi-type hidden events expectation calculation for the zero migration case
    @Test
    public void multiTypeEventsNoMigrationTest() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
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
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "rhoSampling", new TimedParameter(
                        origin,
                        new RealParameter("0.2 0.0"), 2));

        // Integrate p0ge system
        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", startTypePriorProbs,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true);


        density.calculateLogP();
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        int nodeNr = 0;

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density);

        Executor pool = ForkJoinPool.commonPool();
        double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

        MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                new MultiTypeHiddenEventsIntegrator(
                        parameterization, tree, p0geResults,
                        1e-6, 1e-6,
                        false, true, pool, weightOfNodeSubTree, 0.1
                );

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getDoubleValues(), parameterization,1.0, 2.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 3.3935;
        double sim_hiddenEvents_1 = 0.0;
        double tolerance = 1e-2;

        assertEquals("Type 0 hidden events mismatch", sim_hiddenEvents_0, hiddenEvents[0], tolerance);
        assertEquals("Type 1 hidden events mismatch", sim_hiddenEvents_1, hiddenEvents[1], tolerance);
    }


    // Test multi-type hidden events expectation calculation for equal migration rates between types
    @Test
    public void multiTypeEventsTest1() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealParameter origin = new RealParameter("1.0");
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
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.4 0.4"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "rhoSampling", new TimedParameter(
                        origin,
                        new RealParameter("0.2 0.0"), 2));

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
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        int nodeNr = 0;

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density);

        Executor pool = ForkJoinPool.commonPool();
        double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

        MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                new MultiTypeHiddenEventsIntegrator(
                        parameterization, tree, p0geResults,
                        1e-6, 1e-6,
                        false, true, pool,
                        weightOfNodeSubTree, 0.1
                );

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getDoubleValues(), parameterization,0.0, 1.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 2.56971;
        double sim_hiddenEvents_1 = 1.67156;
        double tolerance = 1e-2;

        assertEquals("Type 0 hidden events mismatch", sim_hiddenEvents_0, hiddenEvents[0], tolerance);
        assertEquals("Type 1 hidden events mismatch", sim_hiddenEvents_1, hiddenEvents[1], tolerance);
    }


    // Test multi-type hidden events expectation calculation for unequal migration rates between types
    @Test
    public void multiTypeEventsTest2() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealParameter origin = new RealParameter("1.0");
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
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.9 0.23"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "rhoSampling", new TimedParameter(
                        origin,
                        new RealParameter("0.2 0.0"), 2));

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
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        int nodeNr = 0;

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density);

        Executor pool = ForkJoinPool.commonPool();
        double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

        MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                new MultiTypeHiddenEventsIntegrator(
                        parameterization, tree, p0geResults,
                        1e-6, 1e-6,
                        false, true, pool, weightOfNodeSubTree, 0.1
                );

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getDoubleValues(), parameterization,0.0, 1.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 2.78364;
        double sim_hiddenEvents_1 = 1.80522;
        double tolerance = 1e-2;

        assertEquals("Type 0 hidden events mismatch", sim_hiddenEvents_0, hiddenEvents[0], tolerance);
        assertEquals("Type 1 hidden events mismatch", sim_hiddenEvents_1, hiddenEvents[1], tolerance);
    }


    @Test
    public void piIntegrationTest() {

        String newick = "(t5[&type=0]:5.7,((t1[&type=0]:1,t2[&type=0]:2):1,"
                + "(t3[&type=1]:3,t4[&type=1]:4):0.5):1.3):0.0;";
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

        // Integrate p0/ge system
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

        Executor pool = ForkJoinPool.commonPool();
        double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

        MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                new MultiTypeHiddenEventsIntegrator(
                        parameterization, tree, p0geResults,
                        1e-6, 1e-6,
                        true  , true, pool, weightOfNodeSubTree, 0.1
                        // Store π trajectories for testing
                );

        multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getDoubleValues(),
                parameterization, 0.0);


        // Empirical π expectations from BDMM-Prime stochastic mapping
        double[][] expected = new double[][] {
                {0.9494, 0.0506}, // node 5
                {0.4749, 0.5251}, // node 6
                {0.6949, 0.3051}, // node 7
        };
        int[] nodeNumbers = new int[] {5, 6, 7};

        // Compare trajectories at internal nodes
        double tolerance = 0.01;

        for (int i = 0; i < nodeNumbers.length; i++) {
            int nodeNr = nodeNumbers[i];
            Node node = tree.getNode(nodeNr);

            ContinuousOutputModel com = multitypeHiddenEvents.getPiIntegrationResultsForNode(nodeNr);

            com.setInterpolatedTime(parameterization.getNodeTime(node, 0.0));
            double[] pi = com.getInterpolatedState();

            System.out.printf(
                    "Node %d -> π₀=%.4f π₁=%.4f (expected %.4f %.4f)%n",
                    nodeNr, pi[0], pi[1], expected[i][0], expected[i][1]
            );

            assertEquals("π₀ mismatch at node " + nodeNr, expected[i][0], pi[0], tolerance);
            assertEquals("π₁ mismatch at node " + nodeNr, expected[i][1], pi[1], tolerance);
        }
    }


    @Test
    public void piIntegrationTest2() {

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
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0"), 2),
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
                "storeIntegrationResults", true);


        density.calculateLogP();  // Calculate LogP to call integration method
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();
        // nodeNr 0 = t1
        // nodeNr 1 = t2
        // nodeNr 2 = root

        double[] totalTimeInType = new double[2];  // For type 0 and type 1

        for (int nodeNr = 0; nodeNr < 2; nodeNr++) {
            Node node = tree.getNode(nodeNr);
            if (node.isRoot()) continue;

            double nodeTime = parameterization.getNodeTime(node, 0);
            double parentTime = parameterization.getNodeTime(node.getParent(), 0);
            double branchLength = nodeTime - parentTime;
            int steps = 200;
            double stepSize = branchLength / steps;

            Executor pool = ForkJoinPool.commonPool();
            double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

            MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                    new MultiTypeHiddenEventsIntegrator(
                            parameterization, tree, p0geResults,
                            1e-6, 1e-6,
                            true, true, pool, weightOfNodeSubTree, 0.1
                    );

            multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getDoubleValues(),
                    parameterization, 0.0);
            ContinuousOutputModel model = multitypeHiddenEvents.getPiIntegrationResultsForNode(nodeNr);

            for (int i = 0; i < steps; i++) {
                double t = parentTime + (i + 0.5) * stepSize;  // Midpoint of interval
                model.setInterpolatedTime(t);
                double[] state = model.getInterpolatedState();  // [π0, π1, ge0, ge1]

                totalTimeInType[0] += state[0] * stepSize;
                totalTimeInType[1] += state[1] * stepSize;
            }
        }
        // Tree type statistics from BDMM-Prime stochastic mapping method (mapper2.xml)
        double typedTreeLength_0 = 0.93;
        double typedTreeLength_1 = 1.07;
        double tolerance = 5e-3;

        assertEquals("Type 0 edge length mismatch", typedTreeLength_0, totalTimeInType[0], tolerance);
        assertEquals("Type 1 edge length mismatch", typedTreeLength_1, totalTimeInType[1], tolerance);
    }


    @Test
    public void piIntegrationBirthAmongDemesTest() {

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
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1.5 2.0"), 2),
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
                "storeIntegrationResults", true);


        density.calculateLogP();  // Calculate LogP to call integration method
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();
        // nodeNr 0 = t1
        // nodeNr 1 = t2
        // nodeNr 2 = root

        double[] totalTimeInType = new double[2];  // For type 0 and type 1

        for (int nodeNr = 0; nodeNr < 2; nodeNr++) {
            Node node = tree.getNode(nodeNr);
            if (node.isRoot()) continue;

            double nodeTime = parameterization.getNodeTime(node, 0);
            double parentTime = parameterization.getNodeTime(node.getParent(), 0);
            double branchLength = nodeTime - parentTime;
            int steps = 2000;
            double stepSize = branchLength / steps;

            Executor pool = ForkJoinPool.commonPool();
            double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

            MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                    new MultiTypeHiddenEventsIntegrator(
                        parameterization, tree, p0geResults,
                            1e-6, 1e-6,
                        true, true, pool, weightOfNodeSubTree, 0.1
                    );

            multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getDoubleValues(),
                    parameterization, 0.0);
            ContinuousOutputModel model = multitypeHiddenEvents.getPiIntegrationResultsForNode(nodeNr);

            for (int i = 0; i < steps; i++) {
                double t = parentTime + (i + 0.5) * stepSize;  // Midpoint of interval
                model.setInterpolatedTime(t);
                double[] state = model.getInterpolatedState();  // [π0, π1, ge0, ge1]

                totalTimeInType[0] += state[0] * stepSize;
                totalTimeInType[1] += state[1] * stepSize;
            }
        }
        // Tree type statistics from BDMM-Prime stochastic mapping method (mapper2.xml)
        double typedTreeLength_0 = 1.0649;
        double typedTreeLength_1 = 0.9351;
        double tolerance = 1e-2;

        System.out.println("Type 0 edge length = " +  totalTimeInType[0]);
        System.out.println("Type 1 edge length = " +  totalTimeInType[1]);

        assertEquals("Type 0 edge length mismatch", typedTreeLength_0, totalTimeInType[0], tolerance);
        assertEquals("Type 1 edge length mismatch", typedTreeLength_1, totalTimeInType[1], tolerance);
    }


    @Test
    public void piIntegrationSkylineTest() {

    String newick = "(t1[&state=0] : 1.0, t2[&state=1] : 1.0);";
    Tree tree = new TreeParser(newick, false, false, true, 0);
    RealParameter origin = new RealParameter("1.0001");
    RealParameter startTypePriorProbs = new RealParameter("0.5 0.5");

    Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
            "birthRate", new SkylineVectorParameter(
                        new RealParameter("0.33 0.66"),
                        new RealParameter("3.0 2.4 0.5 1.5 1.0 2.0"), 2),
            "deathRate", new SkylineVectorParameter(
                        new RealParameter("0.1 0.2"),
                                new RealParameter("1.0 0.7 0.1 0.2 0.4 0.7"), 2),
            "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                                new RealParameter("0.0 0.0"), 2),
            "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("0.7 0.9"),
                                new RealParameter("0.2 0.3 1.0 0.7 0.6 0.2"), 2),
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
        "storeIntegrationResults", true);

    density.calculateLogP();  // Calculate LogP to call integration method

    ContinuousOutputModel[] p0geResults = density.getIntegrationResults();
    // nodeNr 0 = t1
    // nodeNr 1 = t2
    // nodeNr 2 = root

    double[] totalTimeInType = new double[2];  // For type 0 and type 1

    for (int nodeNr = 0; nodeNr < 2; nodeNr++) {
        Node node = tree.getNode(nodeNr);
        if (node.isRoot()) continue;

        double nodeTime = parameterization.getNodeTime(node, 0);
        double parentTime = parameterization.getNodeTime(node.getParent(), 0);
        double branchLength = nodeTime - parentTime;
        int steps = 200;
        double stepSize = branchLength / steps;

        Executor pool = ForkJoinPool.commonPool();
        double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

        MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                new MultiTypeHiddenEventsIntegrator(
                        parameterization, tree, p0geResults,
                        1e-6, 1e-6,
                        true, true, pool, weightOfNodeSubTree, 0.1
                );

        multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getDoubleValues(),
                parameterization, 0.0);
        ContinuousOutputModel model = multitypeHiddenEvents.getPiIntegrationResultsForNode(nodeNr);

        for (int i = 0; i < steps; i++) {
            double t = parentTime + (i + 0.5) * stepSize;  // Midpoint of interval
            model.setInterpolatedTime(t);
            double[] state = model.getInterpolatedState();  // [π0, π1, ge0, ge1]

            totalTimeInType[0] += state[0] * stepSize;
            totalTimeInType[1] += state[1] * stepSize;
            }
        }
        // Tree type statistics from BDMM-Prime stochastic mapping method (mapper_skyline.xml)
        double typedTreeLength_0 = 1.207;
        double typedTreeLength_1 = 0.793;
        double tolerance = 5e-3;

        assertEquals("Type 0 edge length mismatch", typedTreeLength_0, totalTimeInType[0], tolerance);
        assertEquals("Type 1 edge length mismatch", typedTreeLength_1, totalTimeInType[1], tolerance);
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

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density
        );

        double tolerance = 1e-3;
        for (int nodeNr = 0; nodeNr <= 1; nodeNr++) {
            Node node = tree.getNode(nodeNr);

            Executor pool = ForkJoinPool.commonPool();
            double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

            MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                    new MultiTypeHiddenEventsIntegrator(
                            parameterization, tree, p0geResults,
                            1e-6, 1e-6,
                            false, true, pool, weightOfNodeSubTree, 0.1
                    );

            multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getDoubleValues(),
                    parameterization, 0.0);

            double multiTypeResult = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr)[0];
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

        // Integrate p0/ge system
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

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density
        );


        double tolerance = 1e-3;
        for (int nodeNr = 0; nodeNr <= 1; nodeNr++) {
            Node node = tree.getNode(nodeNr);

            Executor pool = ForkJoinPool.commonPool();
            double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

            MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                    new MultiTypeHiddenEventsIntegrator(
                            parameterization, tree, p0geResults,
                            1e-6, 1e-6,
                            false, true, pool, weightOfNodeSubTree,
                            0.1
                    );

            multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getDoubleValues(),
                    parameterization, 0.0);

            double multiTypeResult = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr)[0];
            double singleTypeResult = bsp.getExpNrHiddenEventsForBranch(node);

            System.out.printf("Node %d: multi-type = %.10f, single-type = %.10f%n",
                    nodeNr, multiTypeResult, singleTypeResult);

            assertEquals("Mismatch at node " + nodeNr + "single-type result = " + singleTypeResult
                            + " does not match multi-type result = " +  multiTypeResult,
                    singleTypeResult, multiTypeResult, tolerance);

        }
    }

}
