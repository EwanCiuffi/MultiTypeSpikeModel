package multitypespike.clockmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.BoolVector;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;


@Description("Clock model that combines continuous branch rate variation with punctuated spikes of evolution at speciation events")
public class PunctuatedClockModel extends Base {
    final public Input<Tree> treeInput = new Input<>("tree", "tree input", Input.Validate.REQUIRED);

    final public Input<RealVector<? extends NonNegativeReal>> spikeMeanInput = new Input<>("spikeMean", "mean parameter for each spike", Input.Validate.REQUIRED);

    final public Input<RealVector<? extends NonNegativeReal>> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree", Input.Validate.REQUIRED);

    final public Input<RealVectorParam<? extends Real>> ratesInput = new Input<>("rates", "Per-branch rate parameters. If nonCentered=false (default), these are the direct lognormal multipliers. " +
            "If true, these are standard N(0,1) values transformed internally.", Input.Validate.OPTIONAL);

    final public Input<Boolean> relaxedInput = new Input<>("relaxed", "if false then use strict clock (default true)", Input.Validate.OPTIONAL);

    final public Input<BoolVector> indicatorInput = new Input<>("indicator", "if false then no spikes are inferred", Input.Validate.OPTIONAL);

    final public Input<Boolean> noSpikeOnDatedTipsInput = new Input<>("noSpikeOnDatedTips", "Set to true if dated tips should have a spike of 0", false);

    final public Input<Boolean> nonCenteredInput = new Input<>("nonCentered", "If true, uses non-centered parameterisation where relaxed rates are treated as N(0,1) " +
            "and transformed internally to maintain a real-space mean of 1. If false (default), relaxed rates are direct multipliers.", false);

    final public Input<RealScalar<? extends PositiveReal>> rateSDInput = new Input<>("rateSD", "standard deviation of the relaxed-clock lognormal rate distribution. " +
            "Only required when 'nonCentered' is true.", Input.Validate.OPTIONAL);

    public int nTypes, nodeCount;
    int spikeMeanDim, indicatorDim;

    // Per-node spike sums and relaxed rates, rebuilt lazily after any clock input changes
    private double[] spikeSums, relaxedRates;
    private volatile boolean cacheValid = false;


    @Override
    public void initAndValidate() {

        if (relaxedInput.get() != null && relaxedInput.get()) {
            if (ratesInput.get() == null) {
                throw new IllegalArgumentException("If 'relaxed' is true, then the rates input must be provided.");
            }
        }

        if (relaxedInput.get() != null && !relaxedInput.get()) {
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
        nTypes = spikesInput.get().size() / nodeCount;

        // Spike mean dimension checks
        spikeMeanDim = spikeMeanInput.get().size();
        if (nTypes == 1 && spikeMeanDim > 1) {
            throw new IllegalArgumentException("Single-type model requires exactly one spikeMean parameter.");
        }
        if (nTypes > 1 && spikeMeanDim != 1 && spikeMeanDim != nTypes) {
            throw new IllegalArgumentException("For multi-type models, 'spikeMean' must have dimension 1 (shared) or nTypes (" + nTypes + ").");
        }

        // Indicator dimension checks
        if (indicatorInput.get() != null) {
            indicatorDim = indicatorInput.get().size();
            if (nTypes == 1 && indicatorDim > 1) {
                throw new IllegalArgumentException("Single-type model requires at most one indicator parameter.");
            }
            if (nTypes > 1 && indicatorDim != 1 && indicatorDim != nTypes) {
                throw new IllegalArgumentException("For multi-type models, 'indicator' must have dimension 1 (shared) or nTypes (" + nTypes + ").");
            }
        }

        spikeSums = new double[nodeCount];
        relaxedRates = new double[nodeCount];
        cacheValid = false;
    }

    private void updateCache() {
        synchronized (this) {
            if (cacheValid) return;
            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                double spikeSum = 0;
                for (int i = 0; i < nTypes; i++) {
                    if (getIndicator(i))
                        spikeSum += spikesInput.get().get(nodeNr * nTypes + i) * getSpikeMean(i);
                }
                spikeSums[nodeNr] = spikeSum;
                relaxedRates[nodeNr] = getRawRelaxedRate(nodeNr);
            }
            cacheValid = true;
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
        double branchSpike = spikesInput.get().get(node.getNr() * nTypes + type);

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

        if (!cacheValid) updateCache();
        return spikeSums[node.getNr()];
    }


    private double getSpikeMean(int type) {
        if (nTypes == 1 || spikeMeanDim == 1) return spikeMeanInput.get().get(0);
        else return spikeMeanInput.get().get(type);
    }

    private boolean getIndicator(int type) {
        if (indicatorInput.get() == null) return true; // Default to 1 if no indicator input is provided
        if (indicatorDim == 1) return indicatorInput.get().get(0);
        return indicatorInput.get().get(type);
    }

    @Override
    public double getRateForBranch(Node node) {
        double baseRate = meanRateInput.get().get();
        if (node.getLength() <= 0 || node.isRoot() || node.isDirectAncestor()) return baseRate;

        double spikeSize = getSpikeSize(node);

        if (!cacheValid) updateCache();
        double effectiveRelaxedRate = relaxedRates[node.getNr()];

        // Effective rate takes into account spike and base rate
        double branchDistance = node.getLength() * effectiveRelaxedRate + spikeSize;
        return branchDistance / node.getLength();
    }


    /**
     * Returns the raw relaxed branch rate for a node.
     * This is the value that getRateForBranch would use for the multiplicative
     * contribution to branch distance (i.e. baseRate * relaxed_rate_multiplier).
     */
    private double getRawRelaxedRate(int nodeNr) {
        double baseRate = meanRateInput.get().get();
        if (ratesInput.get() == null) return baseRate;
        if (relaxedInput.get() == null || relaxedInput.get()) {
            return getRateMultiplier(nodeNr) * baseRate;
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
    private double getRateMultiplier(int nodeNr) {
        double z = ratesInput.get().get(nodeNr);
        if (!nonCenteredInput.get()) {
            return z;
        }
        double sigma = rateSDInput.get().get();
        return Math.exp(sigma * z - 0.5 * sigma * sigma);
    }

    // BEAST2 state management

    @Override
    protected boolean requiresRecalculation() {
        boolean dirty = InputUtil.isDirty(spikesInput) || InputUtil.isDirty(spikeMeanInput) ||
                InputUtil.isDirty(ratesInput) || InputUtil.isDirty(meanRateInput) ||
                (indicatorInput.get() != null && InputUtil.isDirty(indicatorInput)) ||
                (nonCenteredInput.get() && rateSDInput.get() != null && InputUtil.isDirty(rateSDInput));
        if (dirty) cacheValid = false;
        return dirty;
    }

    @Override
    public void restore() {
        cacheValid = false;
        super.restore();
    }


}




