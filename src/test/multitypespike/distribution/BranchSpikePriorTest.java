package test.multitypespike.distribution;

import bdmmprime.parameterization.*;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import gammaspike.distribution.StumpedTreePrior;
import gammaspike.tree.Stubs;
import mutitypespike.distribution.BranchSpikePrior;
import org.junit.Test;
import static junit.framework.Assert.assertEquals;

public class BranchSpikePriorTest {


    /**
     * Single-type test for expected number of hidden events on a branch with no rate shifts
     * Compares with the output of GammaSpike Model
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

        assertEquals(0.186082405828208, bsp.getExpNrHiddenEventsForBranch(node), 1e-5);
        assertEquals(-1.9352885377590174, bsp.calculateLogP(), 1e-5);

    }





    public static void main(String[] args) {

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
                        new RealParameter("0.13"), 1),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("1.0"))
        );


        Node node = myTree.getNode(5);
//        BranchSpikePrior bsp = new BranchSpikePrior();
//        bsp.initByName("parameterization", parameterization, "tree", myTree, "gammaShape", "1.0", "spikes", "1.0 0.5 0.1");
//        System.out.println(bsp.getExpNrHiddenEventsForInterval(0,0.5));
//        System.out.println(bsp.getExpNrHiddenEventsForBranch(node));

//        System.out.println(parameterization.getNodeTime(node,0));
//        System.out.println(parameterization.getNodeTime(node.getParent(),0));
//
//        System.out.println(node.getParent().getHeight());
//        System.out.println(node.getHeight());


        Stubs stub = new Stubs();

        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();

//        System.out.println(stub.sampleNStubsOnBranch(5,6));
//        System.out.println(gamma_bsp.calculateLogP());

        StumpedTreePrior stp = new StumpedTreePrior();
        stp.initByName("lambda", "0.75", "r0", "2.5", "samplingProportion", "0.25", "tree", myTree);


        stub.initByName("tree", myTree, "prior", stp);
        gamma_bsp.initByName("spikes", "1.0 0.5 0.1", "shape", "1.0", "stubs", stub, "tree", myTree);

//        System.out.println(stp.getMeanStubNumber(node.getHeight(), node.getParent().getHeight()) );

        System.out.println(gamma_bsp.calculateLogP());

    }



}
