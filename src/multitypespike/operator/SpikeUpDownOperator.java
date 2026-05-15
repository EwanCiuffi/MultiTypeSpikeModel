package multitypespike.operator;

import bdmmprime.parameterization.Parameterization;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.kernel.BactrianScaleOperator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Scales up non-zero spikes on a branch and the corresponding branch rate down or vice versa." +
        "For multi-type analyses spikes are scaled uniformly across all types")
public class SpikeUpDownOperator extends BactrianScaleOperator {

    final public Input<RealParameter> spikesInput = new Input<>("spikes",
            "spikes associated with each branch", Input.Validate.REQUIRED);
    final public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-prime parameterization object (see BDMM-prime package for available parameterizations)",
            Input.Validate.REQUIRED);
    final public Input<Boolean> reverseInput = new Input<>("reverseDirection",
            "If true, moves spikes up and corresponding rates down", false);

    private int nTypes;
    private boolean reverse;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        reverse = reverseInput.get();
        nTypes = parameterizationInput.get().getNTypes();

        // Validate index layout: spikes.dimension must equal rates.dimension * nTypes
        // so that spike index = rateIndex * nTypes + typeIndex is well-defined.
        final RealParameter rates = parameterInput.get();
        final RealParameter spikes = spikesInput.get();
        if (spikes.getDimension() != rates.getDimension() * nTypes) {
            throw new IllegalArgumentException(
                    "SpikeUpDownOperator: spikes.dimension (" + spikes.getDimension() +
                            ") must equal rates.dimension (" + rates.getDimension() +
                            ") * nTypes (" + nTypes + ").");
        }
    }

    @Override
    public double proposal() {

        final RealParameter rates = parameterInput.get();
        final RealParameter spikes = spikesInput.get();

        // Sample a rate index uniformly.
        final int index = Randomizer.nextInt(rates.getDimension());
        final double r = rates.getValue(index);

        final double scale = getScaler(index, r);

        // ---- Propose new rate ----
        final double rNew;
        if (!reverse) {
            rNew = r * scale;       // rate goes up
        } else {
            rNew = r / scale;       // rate goes down
        }

        if (rNew < rates.getLower() || rNew > rates.getUpper()) {
            return Double.NEGATIVE_INFINITY;
        }
        rates.setValue(index, rNew);

        // ---- Propose new spikes ----
        // Spikes are laid out as: spikeIndex = rateIndex * nTypes + typeIndex.
        // Only non-zero spikes are scaled; zero spikes (switched off by
        // SpikeFlipOperator) remain at zero and do not contribute to the Jacobian.
        //
        // Jacobian for the full transformation:
        //   r  -> r * scale          contributes +log(scale)
        //   s_i -> s_i / scale       contributes -log(scale)  (per non-zero spike)
        //   s_i = 0 -> 0             contributes nothing
        //
        // log|J| = log(scale) - nonZeroCount * log(scale)
        //        = (1 - nonZeroCount) * log(scale)          [reverse=false]
        //
        // For reverse=true (r/scale, s*scale):
        //   log|J| = -log(scale) + nonZeroCount * log(scale)
        //          = (nonZeroCount - 1) * log(scale)

        int nonZeroCount = 0;

        for (int i = 0; i < nTypes; i++) {
            final int spikeIndex = index * nTypes + i;
            final double s = spikes.getValue(spikeIndex);

            if (s == 0.0) {
                // Spike is off; leave it at zero and skip.
                continue;
            }

            nonZeroCount++;

            final double sNew;
            if (!reverse) {
                sNew = s / scale;   // spike goes down
            } else {
                sNew = s * scale;   // spike goes up
            }

            if (sNew < spikes.getLower() || sNew > spikes.getUpper()) {
                return Double.NEGATIVE_INFINITY;
            }
            spikes.setValue(spikeIndex, sNew);
        }

        // ---- Log Hastings ratio = log Jacobian ----
        if (!reverse) {
            return (1 - nonZeroCount) * Math.log(scale);
        } else {
            return (nonZeroCount - 1) * Math.log(scale);
        }
    }
}