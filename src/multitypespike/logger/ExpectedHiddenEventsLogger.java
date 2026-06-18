package multitypespike.logger;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import multitypespike.distribution.BranchSpikePrior;

import java.io.PrintStream;


@Description("Logs the expected number of hidden speciation events per branch")
public class ExpectedHiddenEventsLogger extends CalculationNode implements Function, Loggable {
    final public Input<BranchSpikePrior> branchSpikePriorInput =
            new Input<>("branchSpikePrior", "Branch spike prior", Validate.REQUIRED);
    final public Input<Boolean> logPerTypeInput = new Input<>(
            "logPerType",
            "If true, log expected hidden events of each type separately for multi-type models; " +
                    "if false, log totals per node (sum across types).",
            false  // default: sum across types
    );

    private boolean logPerType;
    private int nTypes, nodeCount;

    @Override
    public void initAndValidate() {
        logPerType = logPerTypeInput.get();
        nTypes = branchSpikePriorInput.get().nTypes;
        nodeCount = branchSpikePriorInput.get().nodeCount;

        if (nTypes == 1 && logPerType) throw new RuntimeException("logPerType cannot be true for single-type models.");
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            out.print(this.getArrayValue(i) + "\t");
        }
    }


    @Override
    public int getDimension() {
        if(logPerType) return nTypes * nodeCount;
        else return nodeCount;
    }

    @Override
    public double getArrayValue(int dim) {
        if(nTypes == 1 || logPerType) {
            return branchSpikePriorInput.get().getExpectedHiddenEvents(dim);
        } else {
            // Sum across all types for each node
            double sum = 0.0;
            for (int type = 0; type < nTypes; type ++) {
                sum += branchSpikePriorInput.get().getExpectedHiddenEvents(dim, type);
            }
            return sum;
        }
    }

    @Override
    public void init(PrintStream out) {
        String id = this.getID();
        if (id == null || id.isEmpty()) id = "expectedHiddenEvents";

        if (logPerType) {
            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                for (int t = 0; t < nTypes; t++) {
                    out.print(id + ".node" + nodeNr + ".type" + t + "\t");
                }
            }
        } else {
            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                out.print(id + ".node" + nodeNr + "\t");
            }
        }
    }

    public int getNTypes() {
        return nTypes;
    }

    @Override
    public void close(PrintStream out) {
    }

}