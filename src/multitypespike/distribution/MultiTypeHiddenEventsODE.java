package multitypespike.distribution;

import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


public class MultiTypeHiddenEventsODE implements FirstOrderDifferentialEquations {

    private final double lambda_i;
    private final double[][] lambda_ij;
    private final int nTypes;
    private final PiState piState;
    private final int nodeNr;
    private final ContinuousOutputModel[] p0geComArray;


    public MultiTypeHiddenEventsODE(int nodeNr, double lambda_i, double[][] lambda_ij,
                                    int nTypes, PiState piState, ContinuousOutputModel[] p0geComArray) {
        this.lambda_i = lambda_i;
        this.lambda_ij = lambda_ij;
        this.nTypes = nTypes;
        this.piState = piState;
        this.nodeNr = nodeNr;
        this.p0geComArray = p0geComArray;
    }


    @Override
    public int getDimension() {
        return this.nTypes;
    }


    private ContinuousOutputModel getP0GeModel(int nodeNr) {
        return p0geComArray[nodeNr];
    }


    public double[] getP0(int nodeNr, double time) {
        ContinuousOutputModel p0geCom = getP0GeModel(nodeNr);
        p0geCom.setInterpolatedTime(time);
        double[] p0ge = p0geCom.getInterpolatedState();

        return p0ge;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        double[] piValues = piState.getPiAtTime(nodeNr, t);
        double[] p0Values = getP0(nodeNr, t);

        for (int i = 0; i < nTypes; i++) {
            double sum = 2 * lambda_i * p0Values[i];

            for (int j = 0; j < nTypes; j++) {
                if (j != i) {
                    sum += lambda_ij[i][j] * p0Values[j];
                }
            }
            yDot[i] = piValues[i] * sum;
        }
    }
}
