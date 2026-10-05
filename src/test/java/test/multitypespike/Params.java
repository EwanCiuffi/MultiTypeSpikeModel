package test.multitypespike;

import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.spec.type.RealScalar;

/** Builds typed BEAST 2.8 parameters from space-separated value strings. */
public class Params {

    private static double[] parse(String values) {
        String[] parts = values.trim().split("\\s+");
        double[] result = new double[parts.length];
        for (int i = 0; i < parts.length; i++) result[i] = Double.parseDouble(parts[i]);
        return result;
    }

    public static RealVectorParam<NonNegativeReal> real(String values) {
        return new RealVectorParam<>(parse(values), NonNegativeReal.INSTANCE);
    }

    public static RealVectorParam<PositiveReal> positive(String values) {
        return new RealVectorParam<>(parse(values), PositiveReal.INSTANCE);
    }

    public static RealScalarParam<NonNegativeReal> scalar(String value) {
        return new RealScalarParam<>(Double.parseDouble(value.trim()), NonNegativeReal.INSTANCE);
    }

    public static SimplexParam simplex(String values) {
        return new SimplexParam(parse(values));
    }

    /** One-element vector holding the value of a scalar (e.g. the process length as a change time). */
    public static RealVectorParam<NonNegativeReal> asVector(RealScalar<? extends NonNegativeReal> scalar) {
        return new RealVectorParam<>(new double[]{scalar.get()}, NonNegativeReal.INSTANCE);
    }
}
