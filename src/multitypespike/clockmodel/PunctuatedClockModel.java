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
    final public Input<RealParameter> ratesInput = new Input<>("rates", "rates associated with nodes in the tree for sampling of individual rates among branches", Input.Validate.OPTIONAL);
    final public Input<BooleanParameter> relaxedInput = new Input<>("relaxed", "if false then use strict clock", Input.Validate.OPTIONAL);
    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "if false then no spikes are inferred", Input.Validate.OPTIONAL);
    final public Input<Boolean> noSpikeOnDatedTipsInput = new Input<>("noSpikeOnDatedTips", "Set to true if dated tips should have a spike of 0", false);
    final public Input<Boolean> projectRelaxedRatesInput = new Input<>("projectRelaxedRates", "(Experimental feature!) If true, project relaxed rates orthogonally to spikes so spikes explain variance first.", false);

    public int nTypes, nodeCount;
    int spikeMeanDim;

    private double[] projectedRates;
    private boolean projectionDirty = true;

    // Threshold for treating a spike as numerically zero
    private static final double spikeZeroTol = 1e-9;

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

        if (projectRelaxedRatesInput.get() && ratesInput.get() == null) {
            throw new IllegalArgumentException(
                    "'projectRelaxedRates' requires a 'rates' input to be provided.");
        }

        if (ratesInput.get() != null) {
            ratesInput.get().setDimension(treeInput.get().getNodeCount());
        }
        nodeCount = treeInput.get().getNodeCount();
        nTypes = spikesInput.get().getDimension() / nodeCount;

        projectedRates = new double[nodeCount];

        // Spike mean dimension checks
        spikeMeanDim = spikeMeanInput.get().getDimension();
        if (nTypes == 1 && spikeMeanDim > 1) {
            throw new IllegalArgumentException("Single-type model requires exactly one spikeMean parameter.");
        }
        if (nTypes > 1 && spikeMeanDim != 1 && spikeMeanDim != nTypes) {
            throw new IllegalArgumentException("For multi-type models, 'spikeMean' must have dimension 1 (shared) or nTypes (" + nTypes + ").");
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
        if (indicatorInput.get() != null && !indicatorInput.get().getValue()) {
            return 0;
        }

        // Suppress spikes on dated tips if requested
        if (noSpikeOnDatedTipsInput.get()) {
            if (node.isLeaf() && node.getHeight() > 0) {
                return 0;
            }
        }

        if (node.isRoot() || node.isDirectAncestor()) return 0;

        double spikeMean = getSpikeMean(type);
        double branchSpike = spikesInput.get().getValue(node.getNr() * nTypes + type);

        if (branchSpike < spikeZeroTol) return 0;
        else return branchSpike * spikeMean;
    }


    // -------------------------------------------------------------------------
    //  Spike size
    // -------------------------------------------------------------------------

    /**
     * Get the size of a spike (this will be zero if the node is the root or a sampled ancestor)
     * @param node
     * @return
     */
    public double getSpikeSize(Node node) {

        // Spike indicator switch
        if (indicatorInput.get() != null && !indicatorInput.get().getValue()) {
            return 0;
        }

        // Suppress spikes on dated tips if requested
        if (noSpikeOnDatedTipsInput.get()) {
            if (node.isLeaf() && node.getHeight() > 0) {
                return 0;
            }
        }

        if (node.isRoot() || node.isDirectAncestor()) return 0;

        // Compute spike size
        double spikeSum = 0;
        for (int i = 0; i < nTypes; i++) {
            double spikeMean = getSpikeMean(i);
            double branchSpike = spikesInput.get().getValue(node.getNr() * nTypes + i);
            if (branchSpike >= spikeZeroTol) {
                spikeSum += branchSpike * spikeMean;
            }
        }
        return spikeSum;
    }


    private double getSpikeMean(int type) {
        if (nTypes == 1 || spikeMeanDim == 1) return spikeMeanInput.get().getArrayValue(0);
        else return spikeMeanInput.get().getArrayValue(type);
    }


    // -------------------------------------------------------------------------
    // Rate computation
    // -------------------------------------------------------------------------


    @Override
    public double getRateForBranch(Node node) {
        double baseRate = meanRateInput.get().getArrayValue();
        if (node.getLength() <= 0 || node.isRoot() || node.isDirectAncestor()) return baseRate;

        double spikeSize = getSpikeSize(node);

        double effectiveRelaxedRate;
        if (projectRelaxedRatesInput.get()) {
            // Only recompute the projection when inputs have changed
            if (projectionDirty) recomputeProjection();
            effectiveRelaxedRate = projectedRates[node.getNr()];
        } else {
            effectiveRelaxedRate = getRawRelaxedRate(node);
        }

        // Effective rate takes into account spike and base rate
        double branchDistance = node.getLength() * effectiveRelaxedRate + spikeSize;
        return branchDistance / node.getLength();
    }


    /**
     * Returns the raw relaxed branch rate for a node, before any projection.
     * This is the value that getRateForBranch would use for the multiplicative
     * contribution to branch distance (i.e. baseRate * relaxed_rate_multiplier).
     */
    private double getRawRelaxedRate(Node node) {
        double baseRate = meanRateInput.get().getArrayValue();
        if (ratesInput.get() == null) return baseRate;
        if (relaxedInput.get() == null || relaxedInput.get().getValue()) {
            return ratesInput.get().getValue(node.getNr()) * baseRate;
        }
        return baseRate;
    }


    // -------------------------------------------------------------------------
    // Projection
    // -------------------------------------------------------------------------

    private void recomputeProjection() {

        double baseRate = meanRateInput.get().getArrayValue();

        double dotDS = 0.0;  // d · s  where d[i] = r[i] - baseRate
        double dotSS = 0.0;  // s · s

        // Pass 1: accumulate dot products over all active branches
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);
            if (node.isRoot() || node.isDirectAncestor()) {
                projectedRates[nodeNr] = Double.NaN;
                continue;
            }
            double s_i = getSpikeSize(node);
            double d_i = getRawRelaxedRate(node) - baseRate;  // mean-zero deviation
            dotDS += d_i * s_i;
            dotSS += s_i * s_i;
        }

        // When no branch has a spike, the spike vector is zero and projection is
        // undefined — fall back to raw rates (coeff = 0 achieves this).
        double coeff = (dotSS > spikeZeroTol) ? dotDS / dotSS : 0.0;

        // Pass 2: apply projection to each active branch
        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            Node node = treeInput.get().getNode(nodeNr);
            if (node.isRoot() || node.isDirectAncestor()) continue;

            double s_i = getSpikeSize(node);
            double r_i = getRawRelaxedRate(node);

            // r_perp[i] = (r[i] - baseRate) - coeff * s[i] + baseRate
            //           = r[i] - coeff * s[i]
            // For s[i] == 0: r_perp[i] = r[i]  (raw rate, as specified)
            projectedRates[nodeNr] = r_i - coeff * s_i;
        }

        projectionDirty = false;
    }


    // -------------------------------------------------------------------------
    // BEAST2 state management
    // -------------------------------------------------------------------------

    @Override
    protected boolean requiresRecalculation() {
        projectionDirty = true;
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
        return false;
    }

    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        projectionDirty = true;
        super.restore();
    }


}




