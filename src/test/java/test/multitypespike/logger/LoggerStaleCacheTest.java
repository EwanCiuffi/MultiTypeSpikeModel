package test.multitypespike.logger;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import multitypespike.clockmodel.PunctuatedClockModel;
import multitypespike.logger.SaltativeProportionLogger;
import multitypespike.logger.SpikeLogger;
import org.junit.jupiter.api.Test;
import test.multitypespike.Params;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * When the clock lies outside the posterior (sampleFromPrior="true"), it is never asked to recalculate, so its
 * cache is not invalidated when the spikes or rates change. The loggers must still report the current state.
 */
public class LoggerStaleCacheTest {

    static final double[] SPIKE_MEANS = {0.1, 0.3};

    @Test
    public void loggersFollowParametersTheClockWasNotToldAbout() {
        TreeParser tree = new TreeParser();
        tree.initByName("newick", "((A:1.0,B:2.0):1.0,C:3.0);", "IsLabelledNewick", true, "adjustTipHeights", false);
        int nodeCount = tree.getNodeCount();

        double[] initialSpikes = new double[2 * nodeCount];
        java.util.Arrays.fill(initialSpikes, 1.0);
        RealVectorParam<NonNegativeReal> spikes = new RealVectorParam<>(initialSpikes, NonNegativeReal.INSTANCE);
        RealVectorParam<Real> rates = new RealVectorParam<>(new double[]{1.0}, Real.INSTANCE);
        PunctuatedClockModel clock = new PunctuatedClockModel();
        clock.initByName("tree", tree, "spikes", spikes, "spikeMean", Params.real("0.1 0.3"),
                "rates", rates, "clock.rate", new RealScalarParam<>(0.5, PositiveReal.INSTANCE));

        SaltativeProportionLogger proportion = new SaltativeProportionLogger();
        proportion.initByName("tree", tree, "clock", clock);
        SpikeLogger spikeLogger = new SpikeLogger();
        spikeLogger.initByName("clock", clock);

        logged(proportion);                         // fills the clock's cache
        for (int i = 0; i < spikes.size(); i++) spikes.set(i, 0.5 + i);
        for (int i = 0; i < rates.size(); i++) rates.set(i, 0.2 * (i + 1));

        double total = 0, abrupt = 0;
        double[] abruptPerType = new double[2];
        for (Node node : tree.getNodesAsArray()) {
            if (node.isRoot()) {
                assertEquals(0.0, spikeLogger.get(node.getNr()), 0.0);
                continue;
            }
            double spikeTotal = 0;
            for (int type = 0; type < 2; type++) {
                double spike = spikes.get(2 * node.getNr() + type) * SPIKE_MEANS[type];
                abruptPerType[type] += spike;
                spikeTotal += spike;
            }
            assertEquals(spikeTotal, spikeLogger.get(node.getNr()), 1e-12);
            abrupt += spikeTotal;
            total += node.getLength() * 0.5 * rates.get(node.getNr()) + spikeTotal;
        }

        String[] values = logged(proportion);
        assertEquals(abrupt / total, Double.parseDouble(values[0]), 1e-12);
        assertEquals(abruptPerType[0] / total, Double.parseDouble(values[1]), 1e-12);
        assertEquals(abruptPerType[1] / total, Double.parseDouble(values[2]), 1e-12);
    }

    static String[] logged(SaltativeProportionLogger logger) {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        logger.log(0, new PrintStream(bytes));
        return bytes.toString().trim().split("\t");
    }
}
