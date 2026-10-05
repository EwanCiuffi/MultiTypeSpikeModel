package test.multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.tree.*;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.SimplexParam;
import test.multitypespike.Params;
import multitypespike.distribution.BranchSpikePrior;
import multitypespike.distribution.MultiTypeHiddenEventsIntegrator;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import bdmmprime.mapping.TypeMappedTree;
import beast.base.util.Randomizer;
import org.junit.jupiter.api.Test;

import java.util.concurrent.Executor;
import java.util.concurrent.ForkJoinPool;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class MultiTypeHiddenEventsTest {

    // Test multi-type hidden events expectation calculation for non-zero birth rate among demes
    @Test
    public void multiTypeEventsBirthAmongDemesTest1() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("1.0");
        SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        Params.real("3.0 3.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.5 0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        Params.real("1.5 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        Params.real("0.0 0.0"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "rhoSampling", new TimedParameter(Params.asVector(origin),
                        Params.real("0.2 0.0"), 2));

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
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
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

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getValues(), parameterization,0.0, 1.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 3.39;
        double sim_hiddenEvents_1 = 0.0;
        double tolerance = 1e-2;

        assertEquals(sim_hiddenEvents_0, hiddenEvents[0], tolerance, "Type 0 hidden events mismatch");
        assertEquals(sim_hiddenEvents_1, hiddenEvents[1], tolerance, "Type 1 hidden events mismatch");
    }


    // Test multi-type hidden events expectation calculation for non-zero birth rate among demes
    @Test
    public void multiTypeEventsBirthAmongDemesTest2() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("1.0");
        SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        Params.real("3.0 3.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.5 0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        Params.real("1.5 1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        Params.real("0.2 0.4"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "rhoSampling", new TimedParameter(Params.asVector(origin),
                        Params.real("0.2 0.0"), 2));

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
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
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

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getValues(), parameterization,0.0, 1.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 2.02553;
        double sim_hiddenEvents_1 = 1.4639;
        double tolerance = 1e-2;

        assertEquals(sim_hiddenEvents_0, hiddenEvents[0], tolerance, "Type 0 hidden events mismatch");
        assertEquals(sim_hiddenEvents_1, hiddenEvents[1], tolerance, "Type 1 hidden events mismatch");
    }


    // Test multi-type hidden events expectation calculation for the zero migration case
    @Test
    public void multiTypeEventsNoMigrationTest() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("2.0");
        SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        Params.real("3.0 3.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.5 0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        Params.real("0.0 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        Params.real("0.0 0.0"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "rhoSampling", new TimedParameter(Params.asVector(origin),
                        Params.real("0.2 0.0"), 2));

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
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
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

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getValues(), parameterization,1.0, 2.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 3.3935;
        double sim_hiddenEvents_1 = 0.0;
        double tolerance = 1e-2;

        assertEquals(sim_hiddenEvents_0, hiddenEvents[0], tolerance, "Type 0 hidden events mismatch");
        assertEquals(sim_hiddenEvents_1, hiddenEvents[1], tolerance, "Type 1 hidden events mismatch");
    }


    // Test multi-type hidden events expectation calculation for equal migration rates between types
    @Test
    public void multiTypeEventsTest1() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("1.0");
        SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        Params.real("3.0 3.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.5 0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        Params.real("0.0 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        Params.real("0.4 0.4"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "rhoSampling", new TimedParameter(Params.asVector(origin),
                        Params.real("0.2 0.0"), 2));

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
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
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

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getValues(), parameterization,0.0, 1.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 2.56971;
        double sim_hiddenEvents_1 = 1.67156;
        double tolerance = 1e-2;

        assertEquals(sim_hiddenEvents_0, hiddenEvents[0], tolerance, "Type 0 hidden events mismatch");
        assertEquals(sim_hiddenEvents_1, hiddenEvents[1], tolerance, "Type 1 hidden events mismatch");
    }


    // Test multi-type hidden events expectation calculation for unequal migration rates between types
    @Test
    public void multiTypeEventsTest2() {

        String newick = "t1[&state=0]:1.0;";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("1.0");
        SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        Params.real("3.0 3.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.5 0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        Params.real("0.0 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        Params.real("0.9 0.23"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "rhoSampling", new TimedParameter(Params.asVector(origin),
                        Params.real("0.2 0.0"), 2));

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
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
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

        multitypeHiddenEvents.integrateSingleLineage(startTypePriorProbs.getValues(), parameterization,0.0, 1.0);

        double[] hiddenEvents = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events  " + hiddenEvents[1]);

        // Hidden events computed from R simulations (hiddenEventsSim.R)
        double sim_hiddenEvents_0 = 2.78364;
        double sim_hiddenEvents_1 = 1.80522;
        double tolerance = 1e-2;

        assertEquals(sim_hiddenEvents_0, hiddenEvents[0], tolerance, "Type 0 hidden events mismatch");
        assertEquals(sim_hiddenEvents_1, hiddenEvents[1], tolerance, "Type 1 hidden events mismatch");
    }


    @Test
    public void piIntegrationTest() {

        String newick = "(t5[&type=0]:5.7,((t1[&type=0]:1,t2[&type=0]:2):1,"
                + "(t3[&type=1]:3,t4[&type=1]:4):0.5):1.3):0.0;";
        Tree tree = new TreeParser(newick, false, false, true, 0);

        RealScalarParam<NonNegativeReal> origin = Params.scalar("6.0");
        SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, Params.real("1.2 1.2"), 2),
                "deathRate", new SkylineVectorParameter(null, Params.real("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(null, Params.real("0.1 0.1"), 2),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.1"), 2),
                "removalProb", new SkylineVectorParameter(null, Params.real("1.0"), 2)
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

        multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getValues(),
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

            assertEquals(expected[i][0], pi[0], tolerance, "π₀ mismatch at node " + nodeNr);
            assertEquals(expected[i][1], pi[1], tolerance, "π₁ mismatch at node " + nodeNr);
        }
    }


    @Test
    public void piIntegrationTest2() {

        String newick = "(t1[&state=0] : 1.0, t2[&state=1] : 1.0);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("2.0");
        SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        Params.real("3.0 3.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.5 0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        Params.real("0.0 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        Params.real("0.2 0.3"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "rhoSampling", new TimedParameter(Params.asVector(origin),
                        Params.real("0.2"), 2));

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

            multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getValues(),
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
        // Tree type statistics from BDMM-Prime stochastic mapping (TypeMappedTree, 250,000 mappings,
        // origin 2.0 as above, i.e. including the stem)
        double typedTreeLength_0 = 0.9705;
        double typedTreeLength_1 = 1.0295;
        double tolerance = 5e-3;

        assertEquals(typedTreeLength_0, totalTimeInType[0], tolerance, "Type 0 edge length mismatch");
        assertEquals(typedTreeLength_1, totalTimeInType[1], tolerance, "Type 1 edge length mismatch");
    }


    @Test
    public void piIntegrationBirthAmongDemesTest() {

        String newick = "(t1[&state=0] : 1.0, t2[&state=1] : 1.0);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("2.0");
        SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        Params.real("3.0 3.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.5 0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        Params.real("1.5 2.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        Params.real("0.2 0.3"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        Params.real("0.0"), 2),
                "rhoSampling", new TimedParameter(Params.asVector(origin),
                        Params.real("0.2"), 2));

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

            multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getValues(),
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
        // Tree type statistics from BDMM-Prime stochastic mapping (TypeMappedTree, 250,000 mappings,
        // origin 2.0 as above, i.e. including the stem)
        double typedTreeLength_0 = 1.1674;
        double typedTreeLength_1 = 0.8326;
        double tolerance = 5e-3;

        System.out.println("Type 0 edge length = " +  totalTimeInType[0]);
        System.out.println("Type 1 edge length = " +  totalTimeInType[1]);

        assertEquals(typedTreeLength_0, totalTimeInType[0], tolerance, "Type 0 edge length mismatch");
        assertEquals(typedTreeLength_1, totalTimeInType[1], tolerance, "Type 1 edge length mismatch");
    }


    @Test
    public void piIntegrationSkylineTest() {

    String newick = "(t1[&state=0] : 1.0, t2[&state=1] : 1.0);";
    Tree tree = new TreeParser(newick, false, false, true, 0);
    RealScalarParam<NonNegativeReal> origin = Params.scalar("1.0001");
    SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

    Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
            "birthRate", new SkylineVectorParameter(
                        Params.real("0.33 0.66"),
                        Params.real("3.0 2.4 0.5 1.5 1.0 2.0"), 2),
            "deathRate", new SkylineVectorParameter(
                        Params.real("0.1 0.2"),
                                Params.real("1.0 0.7 0.1 0.2 0.4 0.7"), 2),
            "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                                Params.real("0.0 0.0"), 2),
            "migrationRate", new SkylineMatrixParameter(
                        Params.real("0.7 0.9"),
                                Params.real("0.2 0.3 1.0 0.7 0.6 0.2"), 2),
            "samplingRate", new SkylineVectorParameter(
                        null,
                                Params.real("0.0"), 2),
            "removalProb", new SkylineVectorParameter(
                        null,
                                Params.real("0.0"), 2),
            "rhoSampling", new TimedParameter(Params.asVector(origin),
                        Params.real("0.2"), 2));

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

        multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getValues(),
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

        assertEquals(typedTreeLength_0, totalTimeInType[0], tolerance, "Type 0 edge length mismatch");
        assertEquals(typedTreeLength_1, totalTimeInType[1], tolerance, "Type 1 edge length mismatch");
    }

    @Test
    public void hiddenEventsODEAnalyticalSingleTypeTest() {

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
                "birthRateAmongDemes", new SkylineMatrixParameter(null, Params.real("0.0"), 1),
                "migrationRate", new SkylineMatrixParameter(null, Params.real("0"), 1),
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
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
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

            multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getValues(),
                    parameterization, 0.0);

            double multiTypeResult = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr)[0];
            double singleTypeResult = bsp.getExpNrHiddenEventsForBranch(node);

            System.out.printf("Node %d: multi-type = %.10f, single-type = %.10f%n",
                    nodeNr, multiTypeResult, singleTypeResult);

            assertEquals(singleTypeResult, multiTypeResult, tolerance, "Mismatch at node " + nodeNr + "single-type result = " + singleTypeResult
                            + " does not match multi-type result = " +  multiTypeResult);

        }
    }


    @Test
    public void hiddenEventsODESingleTypeSkylineTest() {

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
                "birthRateAmongDemes", new SkylineMatrixParameter(null, Params.real("0.0"), 1),
                "migrationRate", new SkylineMatrixParameter(null, Params.real("0.0"), 1),
                "samplingRate", new SkylineVectorParameter(Params.real("2.0"), Params.real("0.5 1.8"), 1),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0"), 1)
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
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1 0.2 0.7 0.1"),
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

            multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbs.getValues(),
                    parameterization, 0.0);

            double multiTypeResult = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr)[0];
            double singleTypeResult = bsp.getExpNrHiddenEventsForBranch(node);

            System.out.printf("Node %d: multi-type = %.10f, single-type = %.10f%n",
                    nodeNr, multiTypeResult, singleTypeResult);

            assertEquals(singleTypeResult, multiTypeResult, tolerance, "Mismatch at node " + nodeNr + "single-type result = " + singleTypeResult
                            + " does not match multi-type result = " +  multiTypeResult);

        }
    }


    @Test
    public void multiTypeSingleTypeEquivalenceTest() {

        String newick = "(t1[&state=1]:1.5, t2[&state=1]:0.5);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealScalarParam<NonNegativeReal> origin = Params.scalar("2.5");

        SimplexParam startTypePriorProbsSingle = Params.simplex("1.0");
        SimplexParam startTypePriorProbsMulti = Params.simplex("0.5 0.5");

        // Single-type parameterization
        Parameterization paramSingle = new CanonicalParameterization();
        paramSingle.initByName(
                "typeSet", new TypeSet(1),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, Params.real("2.0"), 1),
                "deathRate", new SkylineVectorParameter(null, Params.real("1.0"), 1),
                "birthRateAmongDemes", new SkylineMatrixParameter(null, Params.real("0.0"), 1),
                "migrationRate", new SkylineMatrixParameter(null, Params.real("0.0"), 1),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.5"), 1),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0"), 1)
        );

        BirthDeathMigrationDistribution densitySingle = new BirthDeathMigrationDistribution();
        densitySingle.initByName(
                "parameterization", paramSingle,
                "startTypePriorProbs", startTypePriorProbsSingle,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true
        );
        densitySingle.calculateLogP();

        BranchSpikePrior bspSingle = new BranchSpikePrior();
        bspSingle.initByName(
                "parameterization", paramSingle,
                "tree", tree,
                "spikeShape", Params.positive("1.0"),
                "spikes", Params.real("1.0 0.5 0.1"),
                "startTypePriorProbs", startTypePriorProbsSingle,
                "bdmDistr", densitySingle
        );

        // Multi-type parameterization
        Parameterization paramMulti = new CanonicalParameterization();
        paramMulti.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                // All parameters except migration rates are symmetric across types
                "birthRate", new SkylineVectorParameter(null, Params.real("2.0 2.0"), 2),
                "deathRate", new SkylineVectorParameter(null, Params.real("1.0 1.0"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(null, Params.real("0.0 0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(null, Params.real("0.8 1.2"), 2),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.5 0.5"), 2),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0 0.0"), 2)
        );

        BirthDeathMigrationDistribution densityMulti = new BirthDeathMigrationDistribution();
        densityMulti.initByName(
                "parameterization", paramMulti,
                "startTypePriorProbs", startTypePriorProbsMulti,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true
        );
        densityMulti.calculateLogP();

        ContinuousOutputModel[] p0geResultsMulti = densityMulti.getIntegrationResults();

        Executor pool = ForkJoinPool.commonPool();
        double[] weightOfNodeSubTree = new double[tree.getLeafNodeCount() * 2];

        MultiTypeHiddenEventsIntegrator multitypeHiddenEvents =
                new MultiTypeHiddenEventsIntegrator(
                        paramMulti, tree, p0geResultsMulti,
                        1e-6, 1e-6,
                        false, true, pool, weightOfNodeSubTree, 0.1
                );

        multitypeHiddenEvents.integrateHiddenEvents(startTypePriorProbsMulti.getValues(), paramMulti, 0.0);

        // Equivalence comparison
        double tolerance = 1e-3;
        for (int nodeNr = 0; nodeNr <= 1; nodeNr++) {
            Node node = tree.getNode(nodeNr);

            // Get multi-type expectations and sum them across the 2 types
            double[] multiTypeHiddenEventsArr = multitypeHiddenEvents.getExpNrHiddenEventsForNode(nodeNr);
            double multiTypeResultTotal = multiTypeHiddenEventsArr[0] + multiTypeHiddenEventsArr[1];

            // Calculate single-type expectation directly
            double singleTypeResult = bspSingle.getExpNrHiddenEventsForBranch(node);

            System.out.printf("Node %d: multi-type total = %.10f, single-type = %.10f%n",
                    nodeNr, multiTypeResultTotal, singleTypeResult);

            assertEquals(singleTypeResult, multiTypeResultTotal, tolerance, "Mismatch at node " + nodeNr + ": single-type result = " + singleTypeResult
                            + " does not match multi-type total result = " +  multiTypeResultTotal);
        }
    }

    // π at the root must include propagation along the stem from the origin, where
    // startTypePriorProbs applies. Reference: root-type frequencies from BDMM-Prime stochastic mapping.
    @Test
    public void piAtRootStemPropagationTest() {

        String newick = "(t5[&type=0]:5.7,((t1[&type=1]:1,t2[&type=1]:2):1,"
                + "(t3[&type=1]:3,t4[&type=1]:4):0.5):1.3):0.0;";
        Tree tree = new TreeParser(newick, false, false, true, 0);

        RealScalarParam<NonNegativeReal> origin = Params.scalar("8.0");
        SimplexParam startTypePriorProbs = Params.simplex("1.0 0.0");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, Params.real("1.2 1.2"), 2),
                "deathRate", new SkylineVectorParameter(null, Params.real("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(null, Params.real("0.3 0.3"), 2),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.1"), 2),
                "removalProb", new SkylineVectorParameter(null, Params.real("1.0"), 2)
        );

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

        MultiTypeHiddenEventsIntegrator integrator = new MultiTypeHiddenEventsIntegrator(
                parameterization, tree, density.getIntegrationResults(),
                1e-6, 1e-6, false, false, null, new double[tree.getNodeCount()], 0.1);
        integrator.integrateHiddenEvents(startTypePriorProbs.getValues(), parameterization, 0.0);
        double[] piRoot = integrator.getPiAtNode(tree.getRoot().getNr());

        Randomizer.setSeed(42);
        TypeMappedTree mappedTree = new TypeMappedTree();
        mappedTree.initByName(
                "bdmmDistrib", density,
                "startTypePriorProbs", startTypePriorProbs,
                "typeLabel", "type",
                "untypedTree", tree,
                "remapOnLog", true);

        int nMappings = 20000;
        int nRootType1 = 0;
        for (int sample = 1; sample <= nMappings; sample++) {
            mappedTree.remapForLog(sample);
            Node node = mappedTree.getRoot();
            while (node.getChildCount() == 1)
                node = node.getChild(0);
            if ((int) node.getMetaData("type") == 1)
                nRootType1 += 1;
        }
        double expectedPi1 = (double) nRootType1 / nMappings;

        System.out.printf("Root -> π₁=%.4f (stochastic mapping %.4f)%n", piRoot[1], expectedPi1);

        assertEquals(expectedPi1, piRoot[1], 0.015, "π₁ mismatch at root");
        assertEquals(1.0, piRoot[0] + piRoot[1], 1e-8, "π at root not normalised");
    }

    // Two-tip tree used by the tests below; π trajectories are stored.
    private MultiTypeHiddenEventsIntegrator integrateTwoTipTree(Tree tree, String originValue,
                                                               String birthAmongDemes) {
        RealScalarParam<NonNegativeReal> origin = Params.scalar(originValue);
        SimplexParam startTypePriorProbs = Params.simplex("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, Params.real("3.0 3.0"), 2),
                "deathRate", new SkylineVectorParameter(null, Params.real("0.5 0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(null, Params.real(birthAmongDemes), 2),
                "migrationRate", new SkylineMatrixParameter(null, Params.real("0.2 0.3"), 2),
                "samplingRate", new SkylineVectorParameter(null, Params.real("0.0"), 2),
                "removalProb", new SkylineVectorParameter(null, Params.real("0.0"), 2),
                "rhoSampling", new TimedParameter(Params.asVector(origin), Params.real("0.2"), 2));

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

        MultiTypeHiddenEventsIntegrator integrator = new MultiTypeHiddenEventsIntegrator(
                parameterization, tree, density.getIntegrationResults(),
                1e-6, 1e-6, true, false, null, new double[tree.getNodeCount()], 0.1);
        integrator.integrateHiddenEvents(startTypePriorProbs.getValues(), parameterization, 0.0);
        return integrator;
    }

    // With the root at the origin (no stem), startTypePriorProbs applies at the root node,
    // which must still be conditioned on the data below it.
    @Test
    public void piIntegrationNoStemTest() {

        Tree tree = new TreeParser("(t1[&state=0] : 1.0, t2[&state=1] : 1.0);", false, false, true, 0);
        MultiTypeHiddenEventsIntegrator integrator = integrateTwoTipTree(tree, "1.0", "0.0 0.0");

        double[] totalTimeInType = new double[2];
        int steps = 400;
        double stepSize = 1.0 / steps;
        for (int nodeNr = 0; nodeNr < 2; nodeNr++) {
            ContinuousOutputModel model = integrator.getPiIntegrationResultsForNode(nodeNr);
            for (int i = 0; i < steps; i++) {
                model.setInterpolatedTime((i + 0.5) * stepSize);
                double[] state = model.getInterpolatedState();
                totalTimeInType[0] += state[0] * stepSize;
                totalTimeInType[1] += state[1] * stepSize;
            }
        }

        // BDMM-Prime stochastic mapping (TypeMappedTree, 250,000 mappings)
        assertEquals(0.9339, totalTimeInType[0], 5e-3, "Type 0 edge length mismatch");
        assertEquals(1.0661, totalTimeInType[1], 5e-3, "Type 1 edge length mismatch");
    }

    // With birth among demes, the two daughter lineages can start in different types.
    @Test
    public void piDaughterStartBirthAmongDemesTest() {

        Tree tree = new TreeParser("(t1[&state=0] : 1.0, t2[&state=1] : 1.0);", false, false, true, 0);
        MultiTypeHiddenEventsIntegrator integrator = integrateTwoTipTree(tree, "2.0", "1.5 2.0");

        // BDMM-Prime stochastic mapping (TypeMappedTree, 250,000 mappings): P(type 1) at the
        // start of each daughter edge; the parent's own P(type 1) at the root is 0.3662
        double[] expected = {0.3312, 0.4338};

        for (int nodeNr = 0; nodeNr < 2; nodeNr++) {
            ContinuousOutputModel model = integrator.getPiIntegrationResultsForNode(nodeNr);
            model.setInterpolatedTime(1.0);
            double pi1 = model.getInterpolatedState()[1];
            assertEquals(expected[nodeNr], pi1, 5e-3, "Daughter start π₁ mismatch for node " + nodeNr);
        }
        assertEquals(0.3662, integrator.getPiAtNode(2)[1], 5e-3, "π₁ mismatch at root");
    }

}
