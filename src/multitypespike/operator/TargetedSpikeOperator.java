package multitypespike.operator;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import multitypespike.distribution.BranchSpikePrior;
import org.apache.commons.math.distribution.GammaDistributionImpl;

@Description("Operator that draws and proposes new spike values on a branch directly from the prior.")
public class TargetedSpikeOperator extends Operator {

    public final Input<RealParameter> spikesInput = new Input<>("spikes",
            "spikes associated with each branch", Input.Validate.REQUIRED);

    public final Input<BranchSpikePrior> branchSpikePriorInput = new Input<>("branchSpikePrior",
            "The prior distribution that calculates the expected number of hidden events.", Input.Validate.REQUIRED);

    public final Input<Tree> treeInput = new Input<>("tree", "tree input.", Input.Validate.REQUIRED);

    private int nTypes;
    final double MAX_CUM_SUM = 0.999;

    @Override
    public void initAndValidate() {
        nTypes = branchSpikePriorInput.get().nTypes;
    }

    @Override
    public double proposal() {
        final RealParameter spikes = spikesInput.get();
        final BranchSpikePrior prior = branchSpikePriorInput.get();

        // Choose a random node
        int nodeNr = Randomizer.nextInt(treeInput.get().getNodeCount());
        Node node = treeInput.get().getNode(nodeNr);

        // Extract old values and fetch shapes
        double[] sOld = new double[nTypes];
        double[] sNew = new double[nTypes];
        double[] shapes = new double[nTypes];
        for (int i = 0; i < nTypes; i++) {
            sOld[i] = spikes.getValue(nodeNr * nTypes + i);
            shapes[i] = prior.getSpikeShape(i);
        }

        // Handle root/direct ancestor nodes
        if (node.isRoot() || node.isDirectAncestor()) {
            double logQOld = 0.0;
            double logQNew = 0.0;
            for (int i = 0; i < nTypes; i++) {
                // Draw directly from Exponential(1.0) / Gamma(1,1)
                sNew[i] = Randomizer.nextExponential(1.0);
                if (sNew[i] == 0.0) sNew[i] = 1e-15;
                spikes.setValue(nodeNr * nTypes + i, sNew[i]);

                // Log density of Exp(1.0) is -x
                logQOld += -sOld[i];
                logQNew += -sNew[i];
            }
            return logQOld - logQNew;
        }

        // Block sample from the exact joint prior
        double[] expHidden = new double[nTypes];
        double[] piVals = new double[nTypes];
        boolean hasFakeParent = node.getParent().isFake();

        for (int i = 0; i < nTypes; i++) {
            expHidden[i] = prior.getExpectedHiddenEvents(nodeNr, i);
            if (nTypes == 1) {
                piVals[i] = 1.0;
            } else {
                piVals[i] = prior.getPiVals(nodeNr, i);
            }
        }

        // Joint categorical sample the type of the observed speciation event
        int obsEventType = -1;
        if (!hasFakeParent) {
            double u = Randomizer.nextDouble();
            double cumsum = 0.0;
            for (int i = 0; i < nTypes; i++) {
                cumsum += Math.max(piVals[i], 0.0);
                if (u <= cumsum) {
                    obsEventType = i;
                    break;
                }
            }
            if (obsEventType == -1) obsEventType = nTypes - 1;
        }

        // Draw new continuous values for all types
        for (int i = 0; i < nTypes; i++) {
            int obsEvent = (i == obsEventType) ? 1 : 0;
            int k = (int) Randomizer.nextPoisson(expHidden[i]);
            int nSpikes = k + obsEvent;

            if (nSpikes == 0) {
                sNew[i] = 0.0;
            } else {
                sNew[i] = Randomizer.nextGamma(shapes[i] * nSpikes, shapes[i]);
                if (sNew[i] == 0.0) sNew[i] = 1e-15;
            }
            spikes.setValue(nodeNr * nTypes + i, sNew[i]);
        }

        // Calculate Hastings Ratio matching the joint node prior
        // Since we drew from the prior, logQ is just the prior density
        double logQOld = calculateNodeLogPrior(sOld, expHidden, piVals, shapes, hasFakeParent);
        double logQNew = calculateNodeLogPrior(sNew, expHidden, piVals, shapes, hasFakeParent);

        return logQOld - logQNew;
    }

    /**
     * Replicates the exact joint prior calculation for a single node.
     * Because the proposal simulates this exactly, the returned Hastings ratio perfectly
     * cancels the prior ratio in the posterior calculation.
     */
    private double calculateNodeLogPrior(double[] s, double[] expHidden, double[] piVals, double[] shapes, boolean hasFakeParent) {
        if (nTypes == 1) return calculateSingleTypeNodeLogPrior(s[0], expHidden[0], shapes[0], hasFakeParent);
        else return calculateMultiTypeNodeLogPrior(s, expHidden, piVals, shapes, hasFakeParent);
    }

    private double calculateSingleTypeNodeLogPrior(double branchSpike, double expNrHiddenEvents, double spikeShape, boolean hasFakeParent) {
        double logP = 0.0;
        boolean isZeroSpike = (branchSpike == 0.0);
        GammaDistributionImpl gamma = new GammaDistributionImpl(1.0, 1.0);

        if (expNrHiddenEvents > 0) {
            double branchP = 0.0;
            double cumsum = 0.0;
            int k = 0;
            double logFactorialK = 0.0;

            while (cumsum < MAX_CUM_SUM || k < 5) {
                double logpk = -expNrHiddenEvents + k * Math.log(expNrHiddenEvents) - logFactorialK;
                cumsum += Math.exp(logpk);

                int nSpikes = hasFakeParent ? k : k + 1;

                if (nSpikes == 0) {
                    if (isZeroSpike) branchP += Math.exp(logpk);
                } else {
                    gamma.setAlpha(spikeShape * nSpikes);
                    gamma.setBeta(1.0 / spikeShape);
                    double gammaLogP = gamma.logDensity(branchSpike);
                    if (!isZeroSpike && Double.isFinite(gammaLogP)) {
                        branchP += Math.exp(logpk + gammaLogP);
                    }
                }
                k++;
                logFactorialK += Math.log(k);
            }
            logP = (branchP > 0) ? Math.log(branchP) : Double.NEGATIVE_INFINITY;
        } else {
            if (!hasFakeParent) {
                if (isZeroSpike) {
                    logP = Double.NEGATIVE_INFINITY;
                } else {
                    gamma.setAlpha(spikeShape);
                    gamma.setBeta(1.0 / spikeShape);
                    logP = gamma.logDensity(branchSpike);
                }
            } else if (!isZeroSpike) {
                logP = Double.NEGATIVE_INFINITY;
            }
        }
        return logP;
    }

    private double calculateMultiTypeNodeLogPrior(double[] s, double[] expHidden, double[] piVals, double[] shapes, boolean hasFakeParent) {
        double logP = 0.0;
        double[] logP0 = new double[nTypes];
        double[] logP1 = new double[nTypes];
        GammaDistributionImpl gamma = new GammaDistributionImpl(1.0, 1.0);

        for (int i = 0; i < nTypes; i++) {
            double branchSpike = s[i];
            boolean isZeroSpike = (branchSpike == 0.0);
            double expNrHiddenEvents = expHidden[i];
            double spikeShape = shapes[i];

            double prob0 = 0.0;
            double prob1 = 0.0;

            if (expNrHiddenEvents > 0) {
                double cumsum = 0.0;
                int k = 0;
                double logFactorialK = 0.0;

                while (cumsum < MAX_CUM_SUM || k < 5) {
                    double logpk = -expNrHiddenEvents + k * Math.log(expNrHiddenEvents) - logFactorialK;
                    double pk = Math.exp(logpk);
                    cumsum += pk;

                    if (k == 0) {
                        if (isZeroSpike) prob0 += pk;
                    } else if (!isZeroSpike) {
                        gamma.setAlpha(spikeShape * k);
                        gamma.setBeta(1.0 / spikeShape);
                        double gammaLogP = gamma.logDensity(branchSpike);
                        if (Double.isFinite(gammaLogP)) prob0 += pk * Math.exp(gammaLogP);
                    }

                    int nSpikes = k + 1;
                    if (!isZeroSpike) {
                        gamma.setAlpha(spikeShape * nSpikes);
                        gamma.setBeta(1.0 / spikeShape);
                        double gammaLogP = gamma.logDensity(branchSpike);
                        if (Double.isFinite(gammaLogP)) prob1 += pk * Math.exp(gammaLogP);
                    }
                    k++;
                    logFactorialK += Math.log(k);
                }
            } else {
                if (isZeroSpike) prob0 += 1.0;
                if (!isZeroSpike) {
                    gamma.setAlpha(spikeShape);
                    gamma.setBeta(1.0 / spikeShape);
                    double gammaLogP = gamma.logDensity(branchSpike);
                    if (Double.isFinite(gammaLogP)) prob1 += Math.exp(gammaLogP);
                }
            }
            logP0[i] = (prob0 > 0) ? Math.log(prob0) : Double.NEGATIVE_INFINITY;
            logP1[i] = (prob1 > 0) ? Math.log(prob1) : Double.NEGATIVE_INFINITY;
        }

        if (hasFakeParent) {
            for (int i = 0; i < nTypes; i++) logP += logP0[i];
        } else {
            double maxLogTerm = Double.NEGATIVE_INFINITY;
            double[] logTerms = new double[nTypes];

            for (int i = 0; i < nTypes; i++) {
                double pi = piVals[i];
                if (pi > 0) {
                    double term = Math.log(pi) + logP1[i];
                    for (int j = 0; j < nTypes; j++) {
                        if (j != i) term += logP0[j];
                    }
                    logTerms[i] = term;
                    if (term > maxLogTerm) maxLogTerm = term;
                } else {
                    logTerms[i] = Double.NEGATIVE_INFINITY;
                }
            }

            if (maxLogTerm == Double.NEGATIVE_INFINITY) {
                logP += Double.NEGATIVE_INFINITY;
            } else {
                double sumExp = 0.0;
                for (int i = 0; i < nTypes; i++) {
                    if (logTerms[i] > Double.NEGATIVE_INFINITY) {
                        sumExp += Math.exp(logTerms[i] - maxLogTerm);
                    }
                }
                logP += maxLogTerm + Math.log(sumExp);
            }
        }
        return logP;
    }

    @Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(spikesInput.get());
        return list;
    }
}