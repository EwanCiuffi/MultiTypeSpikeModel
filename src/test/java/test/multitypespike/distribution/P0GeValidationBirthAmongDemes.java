package test.multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

public class P0GeValidationBirthAmongDemes {

    public static void main(String[] args) {

        String newick = "(t1[&state=0] : 1.5, t2[&state=1] : 0.5);";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealParameter origin = new RealParameter("2.5");
        RealParameter startTypePriorProbs = new RealParameter("0.5 0.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1.5 1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.1 0.2"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0"), 2));

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

        int nodeNr = 0;
        Node node = tree.getNode(nodeNr);
        double nodeTime = parameterization.getNodeTime(node, 0);
        double end =  node.isRoot() ? parameterization.getNodeTime(node, 0): parameterization.getNodeTime(node.getParent(), 0);

//        System.out.println("nodeTime = "+ nodeTime);
//        System.out.println("endTime = "+ end);

        int steps = 50;
        double stepSize = (nodeTime - end) / (steps - 1);

        ContinuousOutputModel model = p0geResults[nodeNr];

        System.out.println("time\tp0[0]\tp0[1]\tge[0]\tge[1]");
        for (int i = 0; i < steps; i++) {
            double t = nodeTime - i * stepSize;
            model.setInterpolatedTime(t);
            double[] state = model.getInterpolatedState();  // [p0[0], p0[1], ge[0], ge[1]]

            System.out.printf("%.6f\t%.10f\t%.10f\t%.10f\t%.10f%n",
                    t, state[0], state[1], state[2], state[3]);
        }

        try (PrintWriter writer = new PrintWriter(new File("/Users/ewan/code/beast_and_friends/MultiTypeSpikeModel/src/test/multitypespike/distribution/p0ge_output_birthAmongDemes.csv"))) {
            writer.println("time\tp0_0\tp0_1\tge_0\tge_1");
            for (int i = 0; i < steps; i++) {
                double t = nodeTime - i * stepSize;
                model.setInterpolatedTime(t);
                double[] state = model.getInterpolatedState();

                writer.printf("%.6f\t%.10f\t%.10f\t%.10f\t%.10f%n",
                        t, state[0], state[1], state[2], state[3]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}




