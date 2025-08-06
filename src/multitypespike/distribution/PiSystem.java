package multitypespike.distribution;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

public class PiSystem implements FirstOrderDifferentialEquations {

    private final ContinuousOutputModel p0geCom;

    public double[][] b;
    public double[][][] M, b_ij;

    public double totalProcessLength;

    public int nTypes;
    public int nIntervals;
    public double[] intervalEndTimes;

    protected int interval;

    protected FirstOrderIntegrator piIntegrator;

    protected double integrationMinStep, integrationMaxStep;

    public ContinuousOutputModel[] integrationResults;

    public PiSystem(Parameterization parameterization, Tree tree, ContinuousOutputModel p0geCom, double absoluteTolerance, double relativeTolerance) {

        this.p0geCom = p0geCom;
        this.b = parameterization.getBirthRates();
        this.M = parameterization.getMigRates();
        this.b_ij = parameterization.getCrossBirthRates();

        this.totalProcessLength = parameterization.getTotalProcessLength();

        this.nTypes = parameterization.getNTypes();
        this.nIntervals = parameterization.getTotalIntervalCount();

        this.intervalEndTimes = parameterization.getIntervalEndTimes();

        integrationMinStep = parameterization.getTotalProcessLength() * 1e-100;
        integrationMaxStep= parameterization.getTotalProcessLength() / 10;

        integrationResults = new ContinuousOutputModel[tree.getNodeCount()];

        this.piIntegrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep,
                absoluteTolerance, relativeTolerance);
    }


    public void setInterval(int interval) {this.interval = interval;}

    public int getDimension() {
        return this.nTypes;
    }

    private double[] getP0Ge(double time) {
        //  p0:  (0 .. dim-1)
        //  ge: (dim .. 2*dim-1)

        p0geCom.setInterpolatedTime(time);
        double[] p0ge = p0geCom.getInterpolatedState();

        // Trim away small negative values due to numerical integration errors
        for (int i=0; i<p0ge.length; i++) {
            if (p0ge[i] < 0)
                p0ge[i] = 0.0;
        }
        return p0ge;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        double[] p0ge = getP0Ge(t);

        for (int i = 0; i< nTypes; i++) {

            yDot[i] = 2 * y[i] * (b[interval][i] * p0ge[i] + M[interval][i][i]);

            for (int j = 0; j < nTypes; j++) {
                if (j == i) continue;

                yDot[i] += ((b_ij[interval][i][j] * p0ge[i] + M[interval][i][j]) * (p0ge[nTypes + j] / p0ge[nTypes + i])) * y[i];
                yDot[i] -= ((b_ij[interval][j][i] * p0ge[j] + M[interval][j][i]) * (p0ge[nTypes + i] / p0ge[nTypes + j])) * y[j];

            }
        }
    }


    public void setInitialConditionsForPi(PiState state, double[] startTypePriorProbs) {
        double total = 0.0;
        double[] p0geInit = getP0Ge(0);

        for (int type = 0; type < nTypes; type++) {
            state.pi[type] = p0geInit[type + nTypes] * startTypePriorProbs[type];
            total += state.pi[type];
        }

        for (int type = 0; type < nTypes; type++) {
            state.pi[type] /= total;
        }
    }


    /**
     * Perform the integration of pi with initial conds between tStart and tEnd
     * @param state
     * @param tStart
     * @param tEnd
     * @return result of integration
     */
    public void integrate(PiState state, double tStart, double tEnd) {
        piIntegrator.integrate(this, tStart, state.pi, tEnd, state.pi);
    }


    public void integratePiAlongEdge(Node node, double tStart, PiState state,
                                  Parameterization parameterization, double finalSampleOffset) {

        ContinuousOutputModel com = new ContinuousOutputModel();
        piIntegrator.clearStepHandlers();
        piIntegrator.addStepHandler(com);

        double thisTime = tStart;
        double tEnd = parameterization.getNodeTime(node, finalSampleOffset);
        int thisInterval = parameterization.getIntervalIndex(thisTime);
        int endInterval = parameterization.getNodeIntervalIndex(node, finalSampleOffset);

        while (thisInterval > endInterval) {

            double nextTime = intervalEndTimes[thisInterval-1];

            if (Utils.lessThanWithPrecision(nextTime , thisTime)) {
                setInterval(thisInterval);
                integrate(state, thisTime, nextTime);
            }

            thisTime = nextTime;
            thisInterval -= 1;

        }

        if (Utils.greaterThanWithPrecision(thisTime, tEnd)) {
            setInterval(thisInterval);
            integrate(state, thisTime, tEnd);
        }

        // Store integration results for each edge
        integrationResults[node.getNr()] = com;
        state.setIntegrationResults(integrationResults);
    }


    public void integratePiAtNode(Node node, double parentTime, PiState state,
                                  Parameterization parameterization, double finalSampleOffset) {

        // Get the time of this node (child of parent node)
        double nodeTime = parameterization.getNodeTime(node, finalSampleOffset);

        // Skip integration on origin and sampled ancestor edges
        if (!(node.isRoot() || node.isDirectAncestor())) {
            // Determine intervals and integrate from parent to this node
            integratePiAlongEdge(node, parentTime, state, parameterization, finalSampleOffset);
        }

        // Recurse to children
        for (Node child : node.getChildren()) {
            PiState childState = new PiState(nTypes);
            // Copy current state as starting condition for child
            System.arraycopy(state.pi, 0, childState.pi, 0, nTypes);

            integratePiAtNode(child, nodeTime, childState, parameterization, finalSampleOffset);
        }
    }


    public void integratePi(Tree tree, PiState state, double[] startTypePriorProbs,
                            Parameterization parameterization, double finalSampleOffset) {

        Node root = tree.getRoot();

        // Set initial conditions at root
        setInitialConditionsForPi(state, startTypePriorProbs);

        // Start recursive integration from root
        integratePiAtNode(root, 0, state,
                parameterization, finalSampleOffset);
    }

}