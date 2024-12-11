package test.bdmspike.distribution;

import bdmmprime.parameterization.*;
//import bdmspike.distribution.BranchSpikePrior;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import gammaspike.distribution.BranchSpikePrior;
import gammaspike.distribution.StumpedTreePrior;
import gammaspike.tree.Stubs ;

public class BranchSpikePriorTest {

//    public void noRateChangesCaseTest() {
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

//        System.out.println(node.getParent().getHeight());
//        System.out.println(node.getHeight());


        Stubs stub = new Stubs();

        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();

//        System.out.println(stub.sampleNStubsOnBranch(5,6));
//        System.out.println(gamma_bsp.calculateLogP());
//        assertEquals(, -1.5999999999999988);
        StumpedTreePrior stp = new StumpedTreePrior();
        stp.initByName("lambda", "0.75", "r0", "2.5", "samplingProportion", "0.25", "tree", myTree);
//        System.out.println(stp.getMu());

        System.out.println(stp.getPsi());
//
//        stub.initByName("tree", myTree, "prior", stp);
//        gamma_bsp.initByName("spikes", "1.0 0.5 0.1", "shape", "1.0", "stubs", stub, "tree", myTree);

        System.out.println(stp.getMeanStubNumber(node.getHeight(), node.getParent().getHeight()) );

    }

}
