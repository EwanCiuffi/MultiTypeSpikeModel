package multitypespike.clockmodel;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;

// Based on <GammaSpikeModel>  Copyright (C) <2025>  <Jordan Douglas>

@Description("Scales each spike by the spikeMean parameter")
public class SpikeSize extends CalculationNode implements Function, Loggable {
    final public Input<RealParameter> spikesInput = new Input<>("spikes", "argument to be summed", Validate.REQUIRED);
    final public Input<RealParameter> spikeMeanInput = new Input<>("spikeMean", "spike mean parameter", Validate.REQUIRED);


    @Override
    public void initAndValidate() {

    }

    @Override
    public int getDimension() {
        return spikesInput.get().getDimension();
    }

    @Override
    public double getArrayValue() {
        return getArrayValue(0);
    }



    @Override
    public double getArrayValue(int dim) {
        return spikesInput.get().getValue(dim) * spikeMeanInput.get().getValue();
    }



    /**
     * Loggable interface implementation follows
     */
    @Override
    public void init(PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            String id = this.getID();
            if (id == null || id.equals("")) id = "weightedSpike";
            out.print(id + "." + i + "\t");
        }
    }

    @Override
    public void log(long sampleNr, PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            out.print(this.getArrayValue(i) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

} // class Sum
