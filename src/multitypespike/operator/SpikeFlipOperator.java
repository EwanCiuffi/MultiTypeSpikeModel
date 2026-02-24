package multitypespike.operator;

import java.util.ArrayList;
import java.util.List;

import bdmmprime.parameterization.Parameterization;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;



@Description("Flips a spike from zero to non zero, or vice versa")
public class SpikeFlipOperator extends Operator {

    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch", Input.Validate.REQUIRED);
    final public Input<Double> spikeMeanInput = new Input<>("spikeMean", "mean of the exponential distribution for proposing spikes", 0.005);
    final public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-prime parameterization object (see BDMM-prime package for available parameterizations)",
            Input.Validate.REQUIRED);
    final public Input<Boolean> flipAcrossTypesInput = new Input<>("flipAcrossTypes", "if true, flip all " +
            "spikes for a node across types; if false, flip one spike (default false)", false);

    private int nTypes;
    private int nodeCount;
    private boolean flipAcrossTypes;

    @Override
    public void initAndValidate() {
        nTypes = parameterizationInput.get().getNTypes();
        nodeCount = spikesInput.get().getDimension()/nTypes;
        flipAcrossTypes = flipAcrossTypesInput.get();
    }


    @Override
    public double proposal() {


        double spikeMean = spikeMeanInput.get();
        double lambda = 1.0/spikeMean;
        RealParameter spikes = spikesInput.get();

        if (!flipAcrossTypes) {

            // Sample an index
            final int index = Randomizer.nextInt(spikes.getDimension());
            double sOld = spikes.getValue(index);
            double sNew = 0;
            double pFwd, pRev;


            // Add a spike
            if (sOld == 0) {


                sNew = Randomizer.nextExponential(lambda);
                pRev = 0;
                pFwd = Math.log(lambda) - lambda * sNew; // Exponential density

            }

            // Delete the spike
            else {

                sNew = 0;
                pFwd = 0;
                pRev = Math.log(lambda) - lambda * sOld; // Exponential density

            }

            spikes.setValue(index, sNew);

            return pRev - pFwd;
        }

        // --- MULTI-TYPE ---
        int nodeIndex = Randomizer.nextInt(nodeCount);
        int start = nodeIndex * nTypes; // Determine if we are adding or deleting

        boolean adding = true;
        for (int i = 0; i < nTypes; i++) {
            if (spikes.getValue(start + i) != 0) {
                adding = false;
                break;
            }
        }
        double logHR = 0.0;


        if (adding) { // Add spikes for all types

             for (int i = 0; i < nTypes; i++) {


                 double sNew = Randomizer.nextExponential(lambda);
                 spikes.setValue(start + i, sNew);
                 double pFwd = Math.log(lambda) - lambda * sNew;
                 double pRev = 0;
                 logHR += (pRev - pFwd); }

            } else {
            // Delete spikes for all types


             for (int i = 0; i < nTypes; i++) {


                 double sOld = spikes.getValue(start + i);
                 spikes.setValue(start + i, 0.0);
                 double pFwd = 0;
                 double pRev = Math.log(lambda) - lambda * sOld;
                 logHR += (pRev - pFwd);

             }

        }

        return logHR;

    }


    @Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(spikesInput.get());
        return list;
    }


}
