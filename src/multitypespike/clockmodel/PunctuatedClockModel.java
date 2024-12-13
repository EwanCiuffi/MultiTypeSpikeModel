package multitypespike.clockmodel;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;

// Based on <GammaSpikeModel>  Copyright (C) <2024>  <Jordan Douglas>



public class PunctuatedClockModel extends BranchRateModel.Base {
    final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    final public Input<Function> spikeMeanInput = new Input<>("spikeMean", "mean spike size.", Input.Validate.REQUIRED);
    final public Input<RealParameter> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED);
    final public Input<RealParameter> ratesInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.OPTIONAL);

    final public Input<BooleanParameter> relaxedInput = new Input<>("relaxed", "if false then use strict clock", Input.Validate.OPTIONAL);
    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "burst size is 0 of this is false", Input.Validate.OPTIONAL);

    int nRates;
    double[] ratesArray;
    public void initAndValidate() {
        this.nRates = treeInput.get().getNodeCount();
        this.ratesArray = new double[this.nRates];
    }



    public double getSpikeSize(Node node) {

        if (indicatorInput.get() != null && !indicatorInput.get().getValue()) {
            return 0;
        }

        double spikeMean = spikeMeanInput.get().getDoubleValues()[0];

        return spikesInput.get().getValue(node.getNr()) * spikeMean;

    }

    public double getBranchRate(Node node) {

        if (ratesInput.get() == null) return 1;

        if (relaxedInput.get() == null || relaxedInput.get().getValue()) {
            return ratesInput.get().getValue(node.getNr());
        }

        return 1;


    }
    @Override
    public double getRateForBranch(Node node) {


        // Root has average rate
        double baseRate = meanRateInput.get().getArrayValue();
        if (node.getLength() <= 0 || node.isDirectAncestor() || node.isRoot()) return baseRate;


        double burstRate = getSpikeSize(node);
        double branchRate = getBranchRate(node);
        double branchDistance = node.getLength() * baseRate * branchRate + burstRate;


        // Effective rate takes into account burst and base rate
        double effectiveRate = branchDistance / node.getLength();


        return effectiveRate;
    }

    @Override
    protected boolean requiresRecalculation() {


        if (InputUtil.isDirty(spikesInput) || InputUtil.isDirty(meanRateInput) || InputUtil.isDirty(spikeMeanInput)) {
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


    public double[] getRatesArray() {

        for (int i = 0; i < nRates; i++) {
            Node node = this.treeInput.get().getNode(i);
            ratesArray[i] = this.getRateForBranch(node);
            //Log.warning("rate " + i + " = " + ratesArray[i]);
        }
        return ratesArray;
    }

}
