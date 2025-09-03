package multitypespike.clockmodel;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;

// Based on <GammaSpikeModel>  Copyright (C) <2025>  <Jordan Douglas>



public class PunctuatedClockModel extends BranchRateModel.Base {
    final public Input<Tree> treeInput = new Input<>("tree", "tree input", Input.Validate.REQUIRED);
    final public Input<Function> spikeMeanInput = new Input<>("spikeMean", "mean parameter for each spike", Input.Validate.REQUIRED);
    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree", Input.Validate.REQUIRED);
    final public Input<RealParameter> ratesInput = new Input<>("rates", "rates associated with nodes in the tree for sampling of individual rates among branches", Input.Validate.OPTIONAL);
    final public Input<BooleanParameter> relaxedInput = new Input<>("relaxed", "if false then use strict clock", Input.Validate.OPTIONAL);
    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "if false then no spikes are inferred", Input.Validate.OPTIONAL);
    final public Input<Boolean> noSpikeOnDatedTipsInput = new Input<>("noSpikeOnDatedTips", "Set to true if dated tips should have a spike of 0", false);



    @Override
    public void initAndValidate() {
        ratesInput.get().setDimension(treeInput.get().getNodeCount());
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
        double spikeMean = spikeMeanInput.get().getDoubleValues()[0];
        return spikesInput.get().getValue(node.getNr()) * spikeMean;
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
        double effectiveRate = branchDistance / node.getLength();

        return effectiveRate;
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

    public int getSpikeDimension() {
        return spikesInput.get().getDimension();
    }

}



