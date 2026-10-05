package multitypespike.logger;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealVector;
import beast.base.util.Randomizer;
import multitypespike.distribution.BranchSpikePrior;

import java.io.PrintStream;
import java.util.Arrays;


@Description("Logs the number of hidden speciation events per branch: either a stochastic sample " +
        "from the posterior conditional distribution (default), or the analytic expected value (logExpectedValue=true)")
public class HiddenEventsLogger extends CalculationNode implements RealVector<Real>, Loggable {
    final public Input<BranchSpikePrior> branchSpikePriorInput =
            new Input<>("branchSpikePrior", "Branch spike prior", Input.Validate.REQUIRED);
    final public Input<Boolean> logPerTypeInput = new Input<>(
            "logPerType","If true, log hidden events of each type separately for multi-type models; " +
                    "if false, log totals per node (sum across types).",false); // default: sum across types
    final public Input<Boolean> logExpectedValueInput = new Input<>(
            "logExpectedValue", "If true, log the analytic expected number of hidden events per branch " +
                    "instead of a stochastic sample.", false);


    protected BranchSpikePrior bsp;
    protected int nTypes, nodeCount;
    protected boolean logPerType;
    protected boolean logExpectedValue;

    // Hidden events sampled for the current log line, per node and type; a new joint sample is
    // drawn once every dimension has been read (log line or tree metadata)
    private double[] sampledEvents;
    private boolean[] dimRead;

    @Override
    public void initAndValidate() {
        bsp = branchSpikePriorInput.get();
        nTypes = bsp.nTypes;
        nodeCount = bsp.nodeCount;
        logPerType = logPerTypeInput.get();
        logExpectedValue = logExpectedValueInput.get();

        if (nTypes == 1 && logPerType) throw new RuntimeException("logPerType cannot be true for single-type models.");
    }

    @Override
    public void log(long sample, PrintStream out) {
        if (!logExpectedValue) sampleHiddenEvents();
        for (int i = 0; i < size(); i ++) {
            out.print(get(i) + "\t");
        }
    }

    @Override
    public int size() {
        if(logPerType) return nTypes * nodeCount;
        else return nodeCount;
    }

    @Override
    public double get(int dim) {
        if (logExpectedValue) {
            if (nTypes == 1 || logPerType) {
                return bsp.getExpectedHiddenEvents(dim);
            } else {
                // Sum across all types for each node
                double sum = 0.0;
                for (int type = 0; type < nTypes; type++) {
                    sum += bsp.getExpectedHiddenEvents(dim, type);
                }
                return sum;
            }
        }

        if (sampledEvents == null || dimRead[dim]) sampleHiddenEvents();
        dimRead[dim] = true;
        if (nTypes == 1 || logPerType) return sampledEvents[dim];

        // Sum across all types for each node
        double sum = 0.0;
        for (int type = 0; type < nTypes; type++) {
            sum += sampledEvents[dim * nTypes + type];
        }
        return sum;
    }

    /**
     * Samples the hidden events on every branch from their conditional distribution given the spikes:
     * first the type of the observed speciation event, then the hidden events of each type.
     */
    private void sampleHiddenEvents() {
        if (sampledEvents == null) {
            sampledEvents = new double[nodeCount * nTypes];
            dimRead = new boolean[size()];
        }
        Arrays.fill(dimRead, false);

        double[] logP0 = new double[nTypes];
        double[] logP1 = new double[nTypes];
        double[] logProbs = new double[2];
        double[] logWeights = new double[nTypes];

        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            Node node = bsp.treeInput.get().getNode(nodeNr);
            if (node.isRoot() || node.isDirectAncestor()) {
                for (int type = 0; type < nTypes; type++) sampledEvents[nodeNr * nTypes + type] = 0;
                continue;
            }

            for (int type = 0; type < nTypes; type++) {
                BranchSpikePrior.branchSpikeLogProbs(getSpike(nodeNr, type),
                        bsp.getExpectedHiddenEvents(nodeNr, type), bsp.getSpikeShape(type), logProbs);
                logP0[type] = logProbs[0];
                logP1[type] = logProbs[1];
            }

            int obsType = -1;
            if (!node.getParent().isFake()) {
                double maxLogWeight = Double.NEGATIVE_INFINITY;
                for (int i = 0; i < nTypes; i++) {
                    double pi = nTypes == 1 ? 1.0 : bsp.getPiVals(nodeNr, i);
                    logWeights[i] = pi > 0 ? Math.log(pi) + logP1[i] : Double.NEGATIVE_INFINITY;
                    for (int j = 0; j < nTypes; j++) {
                        if (j != i) logWeights[i] += logP0[j];
                    }
                    maxLogWeight = Math.max(maxLogWeight, logWeights[i]);
                }
                double[] weights = new double[nTypes];
                for (int i = 0; i < nTypes; i++) weights[i] = Math.exp(logWeights[i] - maxLogWeight);
                obsType = Randomizer.randomChoicePDF(weights);
            }

            for (int type = 0; type < nTypes; type++) {
                int nObs = type == obsType ? 1 : 0;
                sampledEvents[nodeNr * nTypes + type] = sampleHiddenEventCount(getSpike(nodeNr, type),
                        bsp.getExpectedHiddenEvents(nodeNr, type), bsp.getSpikeShape(type), nObs,
                        nObs == 1 ? logP1[type] : logP0[type]);
            }
        }
    }

    private double getSpike(int nodeNr, int type) {
        return bsp.spikesInput.get().get(nodeNr * nTypes + type);
    }

    /**
     * Inverse-CDF draw of the number of hidden events k given the spike, where the spike is made up of
     * k + nObs increments; logTotal is the normalising constant from BranchSpikePrior.branchSpikeLogProbs.
     */
    public static int sampleHiddenEventCount(double spike, double expNrHiddenEvents, double spikeShape,
                                              int nObs, double logTotal) {
        if (!(expNrHiddenEvents > 0.0) || logTotal == Double.NEGATIVE_INFINITY) return 0;

        double logU = Math.log(Randomizer.nextDouble()) + logTotal;
        double logMu = Math.log(expNrHiddenEvents);
        double logPk = -expNrHiddenEvents;
        double logCumSum = Double.NEGATIVE_INFINITY;

        int k = 0;
        for (; k < 100000; k++) {
            int nIncrements = k + nObs;
            double logDensity;
            if (spike == 0.0) logDensity = nIncrements == 0 ? 0.0 : Double.NEGATIVE_INFINITY;
            else logDensity = nIncrements == 0 ? Double.NEGATIVE_INFINITY
                    : BranchSpikePrior.logSpikeDensity(spike, nIncrements, spikeShape);

            logCumSum = BranchSpikePrior.logAdd(logCumSum, logPk + logDensity);
            if (logCumSum >= logU) break;

            logPk += logMu - Math.log(k + 1);
        }
        return k;
    }

    @Override
    public Real getDomain() {
        return Real.INSTANCE;
    }

    @Override
    public List<Double> getElements() {
        List<Double> elements = new ArrayList<>(size());
        for (int i = 0; i < size(); i++) elements.add(get(i));
        return elements;
    }

    @Override
    public void init(PrintStream out) {
        String id = this.getID();
        if (id == null || id.isEmpty()) id = logExpectedValue ? "expectedHiddenEvents" : "nHiddenEvents";

        if (logPerType) {
            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                for (int t = 0; t < nTypes; t++) {
                    out.print(id + ".node" + nodeNr + ".type" + t + "\t");
                }
            }
        } else {
            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                out.print(id + ".node" + nodeNr + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) {
    }

}
