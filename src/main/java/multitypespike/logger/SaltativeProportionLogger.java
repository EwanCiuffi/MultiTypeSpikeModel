package multitypespike.logger;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import multitypespike.clockmodel.PunctuatedClockModel;

import java.io.PrintStream;


@Description("Logs the proportion of total evolutionary distance across the tree that is " +
        "attributable to punctuated evolution, both in total and broken down per spike type.")
public class SaltativeProportionLogger extends BEASTObject implements Loggable {

    final public Input<Tree> treeInput = new Input<>("tree", "tree", Input.Validate.REQUIRED);
    final public Input<PunctuatedClockModel> clockInput = new Input<>("clock", "punctuated clock model", Input.Validate.REQUIRED);

    private int nTypes;

    @Override
    public void initAndValidate() {
        nTypes = clockInput.get().nTypes;
    }

    @Override
    public void init(PrintStream out) {
        // Log the total proportion
        out.print(this.clockInput.get().getID() + ".ProportionOfSaltation\t");

        // Dynamically append type-specific columns for multi-type analyses
        if (nTypes > 1) {
            for (int i = 0; i < nTypes; i++) {
                out.print(this.clockInput.get().getID() + ".ProportionOfSaltation_type" + i + "\t");
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        double totalChange = 0;
        double abruptChange = 0;
        double[] abruptChangePerType = new double[nTypes];

        for (Node node : treeInput.get().getNodesAsArray()) {
            if (node.getLength() <= 0 || node.isDirectAncestor() || node.isRoot()) {
                continue;
            }

            // Total distance on branch
            double branchTime = node.getLength();
            double branchRate = clockInput.get().getRateForBranch(node);
            double branchDistance = branchTime * branchRate;

            totalChange += branchDistance;

            // Total abrupt change on this branch
            double abruptChangeBranch = clockInput.get().getSpikeSize(node);
            abruptChange += abruptChangeBranch;

            // Type-specific abrupt changes
            if (nTypes > 1) {
                for (int i = 0; i < nTypes; i++) {
                    abruptChangePerType[i] += clockInput.get().getSpikeSize(node.getNr(), i);
                }
            }
        }

        // Prevent division by zero if there's no tree length/rates
        if (totalChange == 0) {
            out.print("0.0\t");
            if (nTypes > 1) {
                for (int i = 0; i < nTypes; i++) {
                    out.print("0.0\t");
                }
            }
            return;
        }

        // Print total proportion
        out.print((abruptChange / totalChange) + "\t");

        // Print type-specific proportions
        if (nTypes > 1) {
            for (int i = 0; i < nTypes; i++) {
                out.print((abruptChangePerType[i] / totalChange) + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) {
    }
}