package multitypespike.distribution;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import beast.base.evolution.tree.Node;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;


public class MultiTypeHiddenEventsODEs implements FirstOrderDifferentialEquations {

    public double[] expNrHiddenEvents;

    public double[][] b;
    public double[][][] b_ij;
    protected int interval;

    private final int nTypes;
    private final int nodeNr;
    public double[] intervalEndTimes;

    private final ContinuousOutputModel[] p0geComArray;
    private final ContinuousOutputModel piCom;

    protected double integrationMinStep, integrationMaxStep;

    protected FirstOrderIntegrator integrator;

    public MultiTypeHiddenEventsODEs(int nodeNr, Parameterization parameterization, ContinuousOutputModel[] p0geComArray,
                                     ContinuousOutputModel piCom,  double absoluteTolerance, double relativeTolerance) {

        this.nodeNr = nodeNr;
        this.b = parameterization.getBirthRates();
        this.b_ij = parameterization.getCrossBirthRates();
        this.nTypes = parameterization.getNTypes();
        this.intervalEndTimes = parameterization.getIntervalEndTimes();
        this.p0geComArray = p0geComArray;
        this.piCom = piCom;
        integrationMinStep = parameterization.getTotalProcessLength() * 1e-100;
        integrationMaxStep= parameterization.getTotalProcessLength() / 10;

        this.integrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep,
                absoluteTolerance, relativeTolerance);
    }

    @Override
    public int getDimension() {
        return this.nTypes;
    }

    private ContinuousOutputModel getP0GeModel(int nodeNr) {
        return p0geComArray[nodeNr];
    }


    public void setInterval(int interval) {this.interval = interval;}

    public double[] getP0(int nodeNr, double time) {
        ContinuousOutputModel p0geCom = getP0GeModel(nodeNr);
        p0geCom.setInterpolatedTime(time);
        double[] p0ge = p0geCom.getInterpolatedState();

        return p0ge;
    }

    public double[] getPi(double time) {
        piCom.setInterpolatedTime(time);
        double[] p0ge = piCom.getInterpolatedState();

        return p0ge;
    }


    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        double[] piValues = getPi(t);
        double[] p0Values = getP0(nodeNr, t);

        for (int i = 0; i < nTypes; i++) {
            double sum = 2 * b[interval][i] * p0Values[i];

            for (int j = 0; j < nTypes; j++) {
                if (j != i) {
                    sum += b_ij[interval][i][j] * p0Values[j];
                }
            }
            yDot[i] = piValues[i] * sum;
        }
    }

    public void integrate(double[] expNrHiddenEvents, double tStart, double tEnd) {
        integrator.integrate(this, tStart, expNrHiddenEvents, tEnd, expNrHiddenEvents);
    }

    public void integrateForBranch(Node node, Parameterization parameterization, double finalSampleOffset) {
        this.expNrHiddenEvents = new double[nTypes];

        if (node.isRoot() || node.isDirectAncestor()) return;

        int startInterval = parameterization.getNodeIntervalIndex(node.getParent(), finalSampleOffset);
        int endInterval = parameterization.getNodeIntervalIndex(node, finalSampleOffset);
        double t0 = parameterization.getNodeTime(node.getParent(), finalSampleOffset);
        double tEnd = parameterization.getNodeTime(node, finalSampleOffset);

        for (int k = startInterval; k < endInterval; k++) {
            double t1 = intervalEndTimes[k];
            if (Utils.lessThanWithPrecision(t0, t1)) {
                setInterval(k);
                integrate(expNrHiddenEvents, t0, t1);
            }
            t0 = t1;
        }

        if (Utils.lessThanWithPrecision(t0, tEnd)) {
            setInterval(endInterval);
            integrate(expNrHiddenEvents, t0, tEnd);
        }
    }


}
