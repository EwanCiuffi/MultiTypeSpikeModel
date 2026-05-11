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

    final public Input<RealParameter> spikesInput = new Input<>("spikes",
            "spikes associated with each branch", Input.Validate.REQUIRED);
    final public Input<Double> spikeMeanInput = new Input<>("spikeMean",
            "mean of the exponential distribution for proposing spikes", 0.005);
    final public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-prime parameterization object (see BDMM-prime package for available parameterizations)",
            Input.Validate.REQUIRED);
    final public Input<Boolean> flipAcrossTypesInput = new Input<>("flipAcrossTypes",
            "if true, flip all spikes for a node across types; if false, flip one spike (default false)", false);

    private int nTypes;
    private int nodeCount;
    private boolean flipAcrossTypes;

    @Override
    public void initAndValidate() {
        nTypes = parameterizationInput.get().getNTypes();
        nodeCount = spikesInput.get().getDimension() / nTypes;
        flipAcrossTypes = flipAcrossTypesInput.get();
    }

    /**
     * Counts how many of the nTypes spikes for a given node are exactly zero.
     *
     * @param spikes    the spike parameter vector
     * @param nodeStart index of the first spike for this node (= nodeIndex * nTypes)
     * @return number of zero-valued spikes in [nodeStart, nodeStart + nTypes)
     */
    private int countZeroSpikes(RealParameter spikes, int nodeStart) {
        int nZero = 0;
        for (int i = 0; i < nTypes; i++) {
            if (spikes.getValue(nodeStart + i) == 0.0) nZero++;
        }
        return nZero;
    }

    @Override
    public double proposal() {

        final double spikeMean = spikeMeanInput.get();
        final double lambda = 1.0 / spikeMean;
        final RealParameter spikes = spikesInput.get();

        // ---- SINGLE-TYPE: flip one spike chosen uniformly ----
        if (!flipAcrossTypes) {

            final int index = Randomizer.nextInt(spikes.getDimension());
            final double sOld = spikes.getValue(index);

            if (sOld == 0.0) {
                // Birth move: 0 -> Exp(lambda)
                // The index-selection probability 1/D cancels in the HR, so we
                // only need the density ratio.
                // log HR = log p(rev) - log p(fwd)
                //        = log(1) - [log(lambda) - lambda * sNew]
                //        = -(log(lambda) - lambda * sNew)
                final double sNew = Randomizer.nextExponential(lambda);
                spikes.setValue(index, sNew);
                final double logPFwd = Math.log(lambda) - lambda * sNew;
                return -logPFwd;  // pRev = log(1) = 0

            } else {
                // Death move: sOld -> 0
                // log HR = [log(lambda) - lambda * sOld] - log(1)
                //        = log(lambda) - lambda * sOld
                final double logPRev = Math.log(lambda) - lambda * sOld;
                spikes.setValue(index, 0.0);
                return logPRev;   // pFwd = log(1) = 0
            }
        }

        // ---- MULTI-TYPE: flip all spikes for one node simultaneously ----
        //
        // This operator only has well-defined forward and reverse paths when a
        // node is either all-zero (birth) or all-nonzero (death).  A mixed state
        // (e.g. [0.0, 0.5]) can arise at runtime when the single-flip operator
        // runs alongside this one.  If we were to proceed on a mixed node we
        // would overwrite spikes without a valid reverse path, breaking detailed
        // balance.  The correct response is to reject the move immediately.
        //
        // Valid moves and their log Hastings ratios:
        //
        //   Birth (all-zero -> all Exp(lambda)):
        //     log HR = sum_i [ 0 - (log(lambda) - lambda * sNew_i) ]
        //
        //   Death (all-nonzero -> all-zero):
        //     log HR = sum_i [ (log(lambda) - lambda * sOld_i) - 0 ]
        //
        // The node-selection probability 1/nodeCount cancels in both directions.

        final int nodeIndex = Randomizer.nextInt(nodeCount);
        final int start = nodeIndex * nTypes;

        final int nZero = countZeroSpikes(spikes, start);

        // Reject any mixed state: no valid forward/reverse path exists.
        if (nZero != 0 && nZero != nTypes) {
            return Double.NEGATIVE_INFINITY;
        }

        final boolean allZero = (nZero == nTypes);
        double logHR = 0.0;

        if (allZero) {
            // Birth: draw new spike values from Exp(lambda)
            for (int i = 0; i < nTypes; i++) {
                final double sNew = Randomizer.nextExponential(lambda);
                spikes.setValue(start + i, sNew);
                final double logPFwd = Math.log(lambda) - lambda * sNew;
                logHR -= logPFwd;  // pRev term is 0
            }
        } else {
            // Death: set all spikes to zero
            for (int i = 0; i < nTypes; i++) {
                final double sOld = spikes.getValue(start + i);
                spikes.setValue(start + i, 0.0);
                final double logPRev = Math.log(lambda) - lambda * sOld;
                logHR += logPRev;  // pFwd term is 0
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