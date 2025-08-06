package multitypespike.distribution;

import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class TypedHiddenEventsODE implements FirstOrderDifferentialEquations {

    private final int i;
    private final double lambda_i;
    private final double[][] lambda_ij;
    private final int nTypes;
    private final ContinuousOutputModel p0geCom;
    private final PiState piState;
    private int nodeNr;

    public TypedHiddenEventsODE(int i, int nodeNr, double lambda_i, double[][] lambda_ij, int nTypes,
                                ContinuousOutputModel p0geCom, PiState piState) {
        this.i = i;
        this.lambda_i = lambda_i;
        this.lambda_ij = lambda_ij;
        this.nTypes = nTypes;
        this.p0geCom = p0geCom;
        this.piState = piState;
        this.nodeNr = nodeNr;
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        p0geCom.setInterpolatedTime(t);
        double[] p0 = p0geCom.getInterpolatedState();
        double pi_i = piState.getPiAtTime(nodeNr, t)[i];

        double sum = 2 * pi_i * lambda_i * p0[i];
        for (int j = 0; j < nTypes; j++) {
            if (j != i) {
                sum += pi_i * lambda_ij[i][j] * p0[j];
            }
        }

        yDot[0] = sum;
    }
}
