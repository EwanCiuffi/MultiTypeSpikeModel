package multitypespike.logger;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import multitypespike.distribution.BranchSpikePrior;

import java.io.PrintStream;

@Description("Logs node type probabilities, either all types or just the maximum.")
public class NodeTypeProbabilityLogger extends CalculationNode implements Function, Loggable {

    final public Input<BranchSpikePrior> branchSpikePriorInput =
            new Input<>("branchSpikePrior", "Branch spike prior", Validate.REQUIRED);

    final public Input<Boolean> logMaxInput =
            new Input<>("logMax", "If true, only log the max type and its probability. If false, log all types.", false);

    private int nTypes, nodeCount;
    private boolean logMax;

    @Override
    public void initAndValidate() {
        nTypes = branchSpikePriorInput.get().nTypes;
        nodeCount = branchSpikePriorInput.get().nodeCount;
        logMax = logMaxInput.get();

        if (nTypes <= 1) {
            throw new RuntimeException("BranchTypeProbabilityLogger is intended only for multi-type analyses (nTypes > 1).");
        }
    }

    @Override
    public int getDimension() {
        if (logMax) {
            return 2 * nodeCount; // Only maxType and maxProb per node
        } else {
            return nTypes * nodeCount; // All types per node
        }
    }

    @Override
    public double getArrayValue(int dim) {
        if (logMax) {
            int nodeNr = dim / 2;
            int index = dim % 2;

            // Find the maximum for this node
            int maxType = 0;
            double maxProb = branchSpikePriorInput.get().getPiVals(nodeNr, 0);
            for (int t = 1; t < nTypes; t++) {
                double prob = branchSpikePriorInput.get().getPiVals(nodeNr, t);
                if (prob > maxProb) {
                    maxProb = prob;
                    maxType = t;
                }
            }
            // Index 0 is maxType, Index 1 is maxProb
            return (index == 0) ? maxType : maxProb;
        } else {
            int nodeNr = dim / nTypes;
            int type = dim % nTypes;
            return branchSpikePriorInput.get().getPiVals(nodeNr, type);
        }
    }


    @Override
    public void init(PrintStream out) {
        String id = this.getID();
        if (id == null || id.isEmpty()) id = "branchTypeProb";

        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            if (logMax) {
                out.print(id + ".node" + nodeNr + ".maxType\t");
                out.print(id + ".node" + nodeNr + ".maxProb\t");
            } else {
                for (int t = 0; t < nTypes; t++) {
                    out.print(id + ".node" + nodeNr + ".type" + t + "\t");
                }
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (int i = 0; i < this.getDimension(); i++) {
            out.print(this.getArrayValue(i) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {}
}