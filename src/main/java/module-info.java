open module multitypespike {
    requires beast.pkgmgmt;
    requires beast.base;
    requires static beast.fx;
    requires static javafx.controls;
    requires sampled.ancestors;
    requires bdmmprime;
    requires commons.math3;

    exports multitypespike.clockmodel;
    exports multitypespike.distribution;
    exports multitypespike.logger;
    exports multitypespike.operator;

    provides beast.base.core.BEASTInterface with
        multitypespike.clockmodel.PunctuatedClockModel,
        multitypespike.distribution.BranchSpikePrior,
        multitypespike.logger.SpikeLogger,
        multitypespike.logger.HiddenEventsLogger,
        multitypespike.logger.SaltativeProportionLogger,
        multitypespike.logger.NodeTypeProbabilityLogger,
        multitypespike.operator.SpikeFlipOperator,
        multitypespike.operator.SpikeUpDownOperator,
        multitypespike.operator.TargetedSpikeOperator;
}
