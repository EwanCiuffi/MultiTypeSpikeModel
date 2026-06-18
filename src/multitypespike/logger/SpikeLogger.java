package multitypespike.logger;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import multitypespike.clockmodel.PunctuatedClockModel;


@Description("Logs spike values per branch scaled by the spike mean parameter")
public class SpikeLogger extends CalculationNode implements Function, Loggable {
    final public Input<PunctuatedClockModel> clockModelInput =
            new Input<>("clock", "Punctuated clock model input", Validate.REQUIRED);
    final public Input<Boolean> logPerTypeInput = new Input<>(
            "logPerType","If true, log spikes of each type separately for multi-type models; " +
            "if false, log totals per node (sum across types).",false); // default: sum across types

    private int nTypes, nodeCount;
    private boolean logPerType;

    @Override
    public void initAndValidate() {
        nTypes = clockModelInput.get().nTypes;
        nodeCount = clockModelInput.get().nodeCount;
        logPerType = logPerTypeInput.get();
    }

    @Override
    public int getDimension() {
        if(logPerType) return nTypes * nodeCount;
        else return nodeCount;
    }


    @Override
    public double getArrayValue(int dim) {
        if(nTypes == 1) {
            return clockModelInput.get().getSpikeSize(dim);
        } else if(!logPerType) {
            return clockModelInput.get().getSpikeSize(dim);
        } else {
            int type = dim % nTypes;
            int nodeNr = dim / nTypes;
            return clockModelInput.get().getSpikeSize(nodeNr, type);
        }
    }


    /**
     * Loggable interface implementation follows
     */
    @Override
    public void init(PrintStream out) {
        String id = this.getID();
        if (id == null || id.isEmpty()) id = "scaledSpike";

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

    @Override
    public void log(long sampleNr, PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            out.print(this.getArrayValue(i) + "\t");
        }
    }

    public int getNTypes() {
        return nTypes;
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

}
