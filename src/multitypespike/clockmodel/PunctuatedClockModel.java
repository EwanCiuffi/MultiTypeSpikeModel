package multitypespike.clockmodel;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;



public class PunctuatedClockModel extends BranchRateModel.Base {
    final public Input<Tree> treeInput = new Input<>("tree", "tree input", Input.Validate.REQUIRED);

    final public Input<Function> spikeMeanInput = new Input<>("spikeMean", "mean parameter for each spike", Input.Validate.REQUIRED);

    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree", Input.Validate.REQUIRED);

    final public Input<RealParameter> ratesInput = new Input<>("rates", "Per-branch rate parameters. If nonCentered=false (default), these are the direct lognormal multipliers. " +
            "If true, these are standard N(0,1) values transformed internally.", Input.Validate.OPTIONAL);

    final public Input<BooleanParameter> relaxedInput = new Input<>("relaxed", "if false then use strict clock", Input.Validate.OPTIONAL);

    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "if false then no spikes are inferred", Input.Validate.OPTIONAL);

    final public Input<Boolean> noSpikeOnDatedTipsInput = new Input<>("noSpikeOnDatedTips", "Set to true if dated tips should have a spike of 0", false);

    final public Input<Boolean> nonCenteredInput = new Input<>("nonCentered", "If true, uses non-centered parameterisation where 'rates' are treated as N(0,1) " +
            "and transformed internally to maintain a real-space mean of 1. If false (default), 'rates' are direct multipliers.", false);

    final public Input<RealParameter> rateSDInput = new Input<>("rateSD", "standard deviation of the relaxed-clock lognormal rate distribution. " +
            "Only required when 'nonCentered' is true.", Input.Validate.OPTIONAL);

    public int nTypes, nodeCount;
    int spikeMeanDim, indicatorDim;


    @Override
    public void initAndValidate() {

        if (relaxedInput.get() != null && relaxedInput.get().getValue()) {
            if (ratesInput.get() == null) {
                throw new IllegalArgumentException("If 'relaxed' is true, then the rates input must be provided.");
            }
        }

        if (relaxedInput.get() != null && !relaxedInput.get().getValue()) {
            if (meanRateInput.get() == null) {
                throw new IllegalArgumentException("If 'relaxed' is false, then the clock.rate input must be provided.");
            }
        }

        if (ratesInput.get() != null) {
            ratesInput.get().setDimension(treeInput.get().getNodeCount());
        }

        if (nonCenteredInput.get() && rateSDInput.get() == null) {
            throw new IllegalArgumentException("If 'nonCentered' is true, the 'rateSD' input must be provided.");
        }

        nodeCount = treeInput.get().getNodeCount();
        nTypes = spikesInput.get().getDimension() / nodeCount;

        // Spike mean dimension checks
        spikeMeanDim = spikeMeanInput.get().getDimension();
        if (nTypes == 1 && spikeMeanDim > 1) {
            throw new IllegalArgumentException("Single-type model requires exactly one spikeMean parameter.");
        }
        if (nTypes > 1 && spikeMeanDim != 1 && spikeMeanDim != nTypes) {
            throw new IllegalArgumentException("For multi-type models, 'spikeMean' must have dimension 1 (shared) or nTypes (" + nTypes + ").");
        }

        // Indicator dimension checks
        if (indicatorInput.get() != null) {
            indicatorDim = indicatorInput.get().getDimension();
            if (nTypes == 1 && indicatorDim > 1) {
                throw new IllegalArgumentException("Single-type model requires at most one indicator parameter.");
            }
            if (nTypes > 1 && indicatorDim != 1 && indicatorDim != nTypes) {
                throw new IllegalArgumentException("For multi-type models, 'indicator' must have dimension 1 (shared) or nTypes (" + nTypes + ").");
            }
        }
    }

    /**
     * Get the size of a spike (this will be zero if the node is the root or a sampled ancestor)
     * @param dim
     * @return
     */
    public double getSpikeSize(int dim) {
        Node node = treeInput.get().getNode(dim);
        return getSpikeSize(node);
    }

    /**
     * Get the size of a spike of a particular type
     * @param dim
     * @return
     */
    public double getSpikeSize(int dim, int type) {
        Node node = treeInput.get().getNode(dim);

        // Spike indicator switch
        if (!getIndicator(type)) {
            return 0;
        }

        // Suppress spikes on dated tips if requested
        if (noSpikeOnDatedTipsInput.get() && node.isLeaf() && node.getHeight() > 0) return 0;

        if (node.isRoot() || node.isDirectAncestor()) return 0;

        double spikeMean = getSpikeMean(type);
        double branchSpike = spikesInput.get().getValue(node.getNr() * nTypes + type);

        return branchSpike * spikeMean;
    }


    /**
     * Get the size of a spike (this will be zero if the node is the root or a sampled ancestor)
     * @param node
     * @return
     */
    public double getSpikeSize(Node node) {

        // Suppress spikes on dated tips if requested
        if (noSpikeOnDatedTipsInput.get() && node.isLeaf() && node.getHeight() > 0) return 0;

        if (node.isRoot() || node.isDirectAncestor()) return 0;

        // Compute spike size
        double spikeSum = 0;
        for (int i = 0; i < nTypes; i++) {
            // Only add to the sum if this specific type indicator is on
            if (getIndicator(i)) {
                double spikeMean = getSpikeMean(i);
                double branchSpike = spikesInput.get().getValue(node.getNr() * nTypes + i);
                spikeSum += branchSpike * spikeMean;
            }
        }
        return spikeSum;
    }


    private double getSpikeMean(int type) {
        if (nTypes == 1 || spikeMeanDim == 1) return spikeMeanInput.get().getArrayValue(0);
        else return spikeMeanInput.get().getArrayValue(type);
    }

    private boolean getIndicator(int type) {
        if (indicatorInput.get() == null) return true; // Default to 1 if no indicator input is provided
        if (indicatorDim == 1) return indicatorInput.get().getValue(0);
        return indicatorInput.get().getValue(type);
    }

    @Override
    public double getRateForBranch(Node node) {
        double baseRate = meanRateInput.get().getArrayValue();
        if (node.getLength() <= 0 || node.isRoot() || node.isDirectAncestor()) return baseRate;

        double spikeSize = getSpikeSize(node);

        double effectiveRelaxedRate = getRawRelaxedRate(node);

        // Effective rate takes into account spike and base rate
        double branchDistance = node.getLength() * effectiveRelaxedRate + spikeSize;
        return branchDistance / node.getLength();
    }


    /**
     * Returns the raw relaxed branch rate for a node.
     * This is the value that getRateForBranch would use for the multiplicative
     * contribution to branch distance (i.e. baseRate * relaxed_rate_multiplier).
     */
    private double getRawRelaxedRate(Node node) {
        double baseRate = meanRateInput.get().getArrayValue();
        if (ratesInput.get() == null) return baseRate;
        if (relaxedInput.get() == null || relaxedInput.get().getValue()) {
            return getRateMultiplier(node) * baseRate;
        }
        return baseRate;
    }

    /**
     * Returns the per-branch rate multiplier for a node.
     * Centered (default): Returns the 'rates' parameter directly.
     * Non-centered: Reconstructs the multiplier from z ~ N(0,1) as
     * exp(rateSD * z - rateSD^2 / 2) to maintain a mean of 1 while decoupling
     * the parameters.
     */
    private double getRateMultiplier(Node node) {
        double z = ratesInput.get().getValue(node.getNr());
        if (!nonCenteredInput.get()) {
            return z;
        }
        double sigma = rateSDInput.get().getArrayValue();
        return Math.exp(sigma * z - 0.5 * sigma * sigma);
    }

    // BEAST2 state management

    @Override
    protected boolean requiresRecalculation() {
        if (InputUtil.isDirty(spikesInput) || InputUtil.isDirty(spikeMeanInput) ||
                InputUtil.isDirty(ratesInput) || InputUtil.isDirty(meanRateInput)) {
            return true;
        }
        if (relaxedInput.get() != null && InputUtil.isDirty(relaxedInput)) {
            return true;
        }
        if (indicatorInput.get() != null && InputUtil.isDirty(indicatorInput)) {
            return true;
        }
        if (nonCenteredInput.get() && rateSDInput.get() != null && InputUtil.isDirty(rateSDInput)) {
            return true;
        }
        return false;
    }

    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }


}




