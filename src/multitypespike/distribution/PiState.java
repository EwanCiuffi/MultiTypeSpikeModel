package multitypespike.distribution;

import org.apache.commons.math3.ode.ContinuousOutputModel;

/**
 * Class containing the values of pi.
 */
public class PiState {

    int dimension;
    protected double[] pi;
    protected ContinuousOutputModel[] integrationResults;

    public PiState(int nTypes) {
        dimension = nTypes;
        pi = new double[nTypes];
    }

//    public ContinuousOutputModel[] getIntegrationResults() {
//        return integrationResults;
//    }

//    public void setIntegrationResults(ContinuousOutputModel[] comArray) {
//        this.integrationResults = comArray;
//    }

//    public double[] getPiAtTime(int nodeNr, double time) {
//        integrationResults[nodeNr].setInterpolatedTime(time);
//        return integrationResults[nodeNr].getInterpolatedState();
//    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();

        for (int type=0; type<dimension; type++) {
            if (type>0)
                sb.append(" ");

            sb.append("pi[").append(type).append("]=").append(pi[type]);
        }

        return sb.toString();
    }
}
