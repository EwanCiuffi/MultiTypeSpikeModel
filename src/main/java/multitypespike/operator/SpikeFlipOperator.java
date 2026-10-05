package multitypespike.operator;

import java.util.ArrayList;
import java.util.List;

import bdmmprime.parameterization.Parameterization;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealVector;
import beast.base.util.Randomizer;
import multitypespike.distribution.BranchSpikePrior;

@Description("Flips a spike from zero to non zero, or vice versa")
public class SpikeFlipOperator extends Operator {

    final public Input<RealVectorParam<? extends NonNegativeReal>> spikesInput = new Input<>("spikes",
            "spikes associated with each branch", Input.Validate.REQUIRED);

    final public Input<RealVector<? extends PositiveReal>> spikeShapeInput = new Input<>("spikeShape", "shape parameter for the " +
            "gamma distribution of the spikes.", Input.Validate.REQUIRED);

    final public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-Prime parameterization, giving the number of types.");

    final public Input<BranchSpikePrior> branchSpikePriorInput = new Input<>("branchSpikePrior",
            "Spike prior, giving the number of types when no parameterization is provided.",
            Input.Validate.XOR, parameterizationInput);

    final public Input<Boolean> flipAcrossTypesInput = new Input<>("flipAcrossTypes",
            "if true, flip all spikes for a node across types; if false, flip one spike (default false)", false);

    private int nTypes;
    private int nodeCount;
    private boolean flipAcrossTypes;


    @Override
    public void initAndValidate() {
        flipAcrossTypes = flipAcrossTypesInput.get();
    }

    // The number of types is resolved on first use, once the spike prior has been set up
    private void resolveTypes() {
        nTypes = parameterizationInput.get() != null ? parameterizationInput.get().getNTypes()
                : branchSpikePriorInput.get().nTypes;
        nodeCount = spikesInput.get().size() / nTypes;

        int spikeShapeDim = spikeShapeInput.get().size();
        if (nTypes == 1 && spikeShapeDim > 1) {
            throw new IllegalArgumentException("Single-type model requires exactly one spikeShape parameter.");
        }
        if (nTypes > 1 && spikeShapeDim != 1 && spikeShapeDim != nTypes) {
            throw new IllegalArgumentException("For multi-type models, 'spikeShape' must have dimension 1 (shared) or nTypes (" + nTypes + ").");
        }
    }


    @Override
    public double proposal() {
        if (nTypes == 0) resolveTypes();

        final RealVectorParam<? extends NonNegativeReal> spikes = spikesInput.get();

        // ---- SINGLE-TYPE: flip one spike chosen uniformly ----
        if (!flipAcrossTypes) {

            final int index = Randomizer.nextInt(spikes.size());
            final int type = index % nTypes;
            final double spikeShape = getSpikeShape(type);
            final double sOld = spikes.get(index);

            if (sOld == 0.0) {
                // Birth move: 0 -> Gamma(spikeShape, 1/spikeShape)
                // The index-selection probability 1/D cancels in the HR, so we
                // only need the density ratio.
                // log HR = log p(rev) - log p(fwd)

                // beta = spikeShape instead of 1/spikeShape due to different parameterisation of the Gamma distribution
                double sNew = Randomizer.nextGamma(spikeShape, spikeShape);

                // Ensure we don't propose exactly 0 due to precision
                if (sNew < 1e-9) sNew = 1e-9;

                spikes.set(index, sNew);

                // logHR = log(pRev) - log(pFwd) = 0 - logDensity(sNew)
                return -BranchSpikePrior.logSpikeDensity(sNew, 1, spikeShape);

            } else {
                // Death: sOld -> 0

                spikes.set(index, 0.0);

                // logHR = log(pRev) - log(pFwd) = logDensity(sOld) - 0
                return BranchSpikePrior.logSpikeDensity(sOld, 1, spikeShape);
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
            // Birth: draw new spike values from Gamma(spikeShape, 1/spikeShape)
            for (int i = 0; i < nTypes; i++) {
                double spikeShape = getSpikeShape(i);

                // beta = spikeShape instead of 1/spikeShape due to different parameterisation of the Gamma distribution
                double sNew = Randomizer.nextGamma(spikeShape, spikeShape);

                if (sNew < 1e-9) sNew = 1e-9;
                spikes.set(start + i, sNew);
                logHR -= BranchSpikePrior.logSpikeDensity(sNew, 1, spikeShape);
            }
        } else {
            // Death: set all spikes to zero
            for (int i = 0; i < nTypes; i++) {
                double spikeShape = getSpikeShape(i);
                double sOld = spikes.get(start + i);
                spikes.set(start + i, 0.0);
                logHR += BranchSpikePrior.logSpikeDensity(sOld, 1, spikeShape);
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
    private int countZeroSpikes(RealVectorParam<? extends NonNegativeReal> spikes, int nodeStart) {
        int nZero = 0;
        for (int i = 0; i < nTypes; i++) {
            if (spikes.get(nodeStart + i) == 0.0) nZero++;
        }
        return nZero;
    }

    @Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(spikesInput.get());
        return list;
    }

    public double getSpikeShape(int type) {
        int spikeShapeDim = spikeShapeInput.get().size();
        if (nTypes == 1 || spikeShapeDim == 1) return spikeShapeInput.get().get(0);
        else return spikeShapeInput.get().get(type);
    }
}