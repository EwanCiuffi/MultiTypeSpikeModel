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
import org.apache.commons.math.distribution.GammaDistributionImpl;

@Description("Flips a spike from zero to non zero, or vice versa")
public class SpikeFlipOperator extends Operator {

    final public Input<RealParameter> spikesInput = new Input<>("spikes",
            "spikes associated with each branch", Input.Validate.REQUIRED);

    final public Input<RealParameter> spikeShapeInput = new Input<>("spikeShape", "shape parameter for the " +
            "gamma distribution of the spikes.", Input.Validate.REQUIRED);

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


    @Override
    public double proposal() {

        final RealParameter spikes = spikesInput.get();
        final double spikeShape = spikeShapeInput.get().getValue();

        // ---- SINGLE-TYPE: flip one spike chosen uniformly ----
        if (!flipAcrossTypes) {

            final int index = Randomizer.nextInt(spikes.getDimension());
            final double sOld = spikes.getValue(index);

            if (sOld == 0.0) {
                // Birth move: 0 -> Gamma(spikeShape, 1/spikeShape)
                // The index-selection probability 1/D cancels in the HR, so we
                // only need the density ratio.
                // log HR = log p(rev) - log p(fwd)

                // beta = spikeShape instead of 1/spikeShape due to different parameterisation of the Gamma distribution
                double sNew = Randomizer.nextGamma(spikeShape, spikeShape);

                // Ensure we don't propose exactly 0 due to precision
                if (sNew < 1e-9) sNew = 1e-9;

                spikes.setValue(index, sNew);

                // logHR = log(pRev) - log(pFwd) = 0 - logDensity(sNew)
                GammaDistributionImpl gamma = new GammaDistributionImpl(spikeShape, 1.0 / spikeShape);
                return -gamma.logDensity(sNew);

            } else {
                // Death: sOld -> 0

                spikes.setValue(index, 0.0);

                // logHR = log(pRev) - log(pFwd) = logDensity(sOld) - 0
                GammaDistributionImpl gamma = new GammaDistributionImpl(spikeShape, 1.0 / spikeShape);
                return gamma.logDensity(sOld);
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
            GammaDistributionImpl gamma = new GammaDistributionImpl(spikeShape, 1.0 / spikeShape);

            // Birth: draw new spike values from Gamma(spikeShape, 1/spikeShape)
            for (int i = 0; i < nTypes; i++) {

                // beta = spikeShape instead of 1/spikeShape due to different parameterisation of the Gamma distribution
                double sNew = Randomizer.nextGamma(spikeShape, spikeShape);

                if (sNew < 1e-9) sNew = 1e-9;
                spikes.setValue(start + i, sNew);
                logHR -= gamma.logDensity(sNew);
            }
        } else {
            GammaDistributionImpl gamma = new GammaDistributionImpl(spikeShape, 1.0 / spikeShape);

            // Death: set all spikes to zero
            for (int i = 0; i < nTypes; i++) {
                double sOld = spikes.getValue(start + i);
                spikes.setValue(start + i, 0.0);
                logHR += gamma.logDensity(sOld);
            }
        }

        return logHR;
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
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(spikesInput.get());
        return list;
    }
}