package multitypespike.clockmodel;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;


/**
 * Based on <GammaSpikeModel>  Copyright (C) <2025>  <Jordan Douglas>
 */

public class PunctuatedClockModel extends BranchRateModel.Base {
    final public Input<Tree> treeInput = new Input<>("tree", "tree input", Input.Validate.REQUIRED);
    final public Input<Function> spikeMeanInput = new Input<>("spikeMean", "mean parameter for each spike", Input.Validate.REQUIRED);
    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree", Input.Validate.REQUIRED);
    final public Input<RealParameter> ratesInput = new Input<>("rates", "rates associated with nodes in the tree for sampling of individual rates among branches", Input.Validate.OPTIONAL);
    final public Input<BooleanParameter> relaxedInput = new Input<>("relaxed", "if false then use strict clock", Input.Validate.OPTIONAL);
    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "if false then no spikes are inferred", Input.Validate.OPTIONAL);
    final public Input<Boolean> noSpikeOnDatedTipsInput = new Input<>("noSpikeOnDatedTips", "Set to true if dated tips should have a spike of 0", false);

    public int nTypes, nodeCount;
    int spikeMeanDim;
    private static final double SPIKE_ZERO_TOL = 1e-8;

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

        if (branchSpike < SPIKE_ZERO_TOL) return 0;
        else return branchSpike * spikeMean;
    }


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
            spikeSum += branchSpike * spikeMean;
        }
        return spikeSum;
    }


    private double getSpikeMean(int type) {
        if (nTypes == 1 || spikeMeanDim == 1) return spikeMeanInput.get().getArrayValue(0);
        else return spikeMeanInput.get().getArrayValue(type);
    }


    @Override
    public double getRateForBranch(Node node) {

        // Root and sampled ancestors have average rate
        double baseRate = meanRateInput.get().getArrayValue();
        if (node.getLength() <= 0 || node.isRoot() || node.isDirectAncestor()) return baseRate;


        double relaxedBranchRate = getRelaxedBranchRate(node);
        double spikeSize = getSpikeSize(node);
        double branchDistance = node.getLength() * baseRate * relaxedBranchRate + spikeSize;


        // Effective rate takes into account spike and base rate
        return branchDistance / node.getLength();
    }


    public double getRelaxedBranchRate(Node node) {

        if (ratesInput.get() == null) return 1;

        if (relaxedInput.get() == null || relaxedInput.get().getValue()) {
            return ratesInput.get().getValue(node.getNr());
        }

        return 1;
    }

    @Override
    protected boolean requiresRecalculation() {

        if (InputUtil.isDirty(spikesInput) || InputUtil.isDirty(meanRateInput) ||
            InputUtil.isDirty(spikeMeanInput) || InputUtil.isDirty(ratesInput)) {
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

}




