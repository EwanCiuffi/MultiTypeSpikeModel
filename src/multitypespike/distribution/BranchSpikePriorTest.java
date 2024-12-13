package test.multitypespike.distribution;

import bdmmprime.parameterization.*;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.*;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.coalescent.RandomTree;
import beast.base.inference.parameter.RealParameter;
import gammaspike.distribution.StumpedTreePrior;
import gammaspike.tree.Stubs;
import mutitypespike.distribution.BranchSpikePrior;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertEquals;

public class BranchSpikePriorTest {

    /**
     * Single-type test for expected number of hidden events on a branch with no rate shifts
     * Compares with the output of GammaSpike Model - 11/12/2024
     */
    @Test
    public void noRateShiftsSingleTypeCaseTest() {

        String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";
        TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
        Tree myTree = treeParser;

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
        bsp.initByName("parameterization", parameterization, "tree", myTree, "gammaShape", "1.0", "spikes", "1.0 0.5 0.1");
        Node node = myTree.getNode(5);

        /*
        Stubs stub = new Stubs();
        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();
        StumpedTreePrior stp = new StumpedTreePrior();
        stp.initByName("lambda", "0.75", "r0", "2.5", "samplingProportion", "0.25", "tree", myTree);
        stub.initByName("tree", myTree, "prior", stp);
        gamma_bsp.initByName("spikes", "1.0 0.5 0.1", "shape", "1.0", "stubs", stub, "tree", myTree);

        stp.getMeanStubNumber(node.getHeight(), node.getParent().getHeight())  = 0.186082405828208.
        gamma_bsp.calculateLogP() = -1.9352885377590174.
        */

        assertEquals(0.186082405828208, bsp.getExpNrHiddenEventsForBranch(node), 1e-5);
        assertEquals(-1.9352885377590174, bsp.calculateLogP(), 1e-5);

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







    public static void main(String[] args) {

        String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";
        TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
        Tree myTree = treeParser;

        RealParameter originParam = new RealParameter("2.0");

        Parameterization parameterization = new CanonicalParameterization();
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
//                        new RealParameter("0.13"), 1),
//                "rhoSampling", new TimedParameter(
//                        originParam,
//                        new RealParameter("1.0"))
//        );

//        Node node = myTree.getNode(5);
//        BranchSpikePrior bsp = new BranchSpikePrior();
//        bsp.initByName("parameterization", parameterization, "tree", myTree, "gammaShape", "1.0", "spikes", "1.0 0.5 0.1");
//        System.out.println(bsp.getExpNrHiddenEventsForInterval(0,0.5));
//        System.out.println(bsp.getExpNrHiddenEventsForBranch(node));

//        System.out.println(parameterization.getNodeTime(node,0));
//        System.out.println(parameterization.getNodeTime(node.getParent(),0));
//
//        System.out.println(node.getParent().getHeight());
//        System.out.println(node.getHeight());


/*
        Stubs stub = new Stubs();
        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();
        StumpedTreePrior stp = new StumpedTreePrior();
        stp.initByName("lambda", "0.75", "r0", "2.5", "samplingProportion", "0.25", "tree", myTree);
        stub.initByName("tree", myTree, "prior", stp);
        gamma_bsp.initByName("spikes", "1.0 0.5 0.1", "shape", "1.0", "stubs", stub, "tree", myTree);
//        System.out.println(stp.getMeanStubNumber(node.getHeight(), node.getParent().getHeight()));
//        System.out.println(gamma_bsp.calculateLogP());

*/
}


}
