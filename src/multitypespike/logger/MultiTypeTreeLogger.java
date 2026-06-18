package multitypespike.logger;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.TreeWithMetaDataLogger;
import beast.base.inference.StateNode;

import java.io.PrintStream;
import java.util.List;

@Description("Logs tree annotated with multi-type metadata")
public class MultiTypeTreeLogger extends TreeWithMetaDataLogger {

    @Override
    public void log(long sample, PrintStream out) {
        Tree tree = (Tree) treeInput.get().getCurrent();
        List<Function> metadata = parameterInput.get();

        // Resolve StateNodes to their current values
        for (int i = 0; i < metadata.size(); i++) {
            if (metadata.get(i) instanceof StateNode) {
                metadata.set(i, (Function) ((StateNode) metadata.get(i)).getCurrent());
            }
        }

        out.print("tree STATE_" + sample + " = ");
        if (sortTreeInput.get()) {
            tree.getRoot().sort();
        }

        out.print(toNewick(tree.getRoot(), metadata, clockModelInput.get()));
        out.print(";");
    }

    protected String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
        return customToNewick(node, metadataList, branchRateModel);
    }

    private String customToNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
        StringBuilder buf = new StringBuilder();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(customToNewick(node.getLeft(), metadataList, branchRateModel));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(customToNewick(node.getRight(), metadataList, branchRateModel));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }

        StringBuilder metaBuf = new StringBuilder();
        if (!metadataList.isEmpty()) {
            metaBuf.append("[&");
            for (int i = 0; i < metadataList.size(); i++) {
                Function m = metadataList.get(i);
                String id = ((BEASTObject) m).getID();
                int nodeNr = node.getNr();

                if (m instanceof BranchTypeProbabilityLogger) {
                    BranchTypeProbabilityLogger bpl = (BranchTypeProbabilityLogger) m;
                    int nTypes = bpl.getNTypes();

                    // 1. Log the full array of probabilities
                    metaBuf.append(id).append("={");
                    for (int t = 0; t < nTypes; t++) {
                        metaBuf.append(bpl.getProbability(nodeNr, t));
                        if (t < nTypes - 1) metaBuf.append(",");
                    }
                    metaBuf.append("}");

                    // 2. Calculate the maxType and maxProb
                    int maxType = 0;
                    double maxProb = bpl.getProbability(nodeNr, 0);
                    for (int t = 1; t < nTypes; t++) {
                        double p = bpl.getProbability(nodeNr, t);
                        if (p > maxProb) {
                            maxProb = p;
                            maxType = t;
                        }
                    }

                    // 3. Inject maxType and maxProb as completely independent keys for FigTree
                    metaBuf.append(",").append(id).append("_maxType=").append(maxType);
                    metaBuf.append(",").append(id).append("_maxProb=").append(maxProb);

                } else if (m instanceof ExpectedHiddenEventsLogger || m instanceof SpikeLogger) {
                    // Standard multi-type array logging for other custom loggers
                    int nTypes = (m instanceof ExpectedHiddenEventsLogger) ?
                            ((ExpectedHiddenEventsLogger) m).getNTypes() :
                            ((SpikeLogger) m).getNTypes();

                    metaBuf.append(id).append("={");
                    for (int t = 0; t < nTypes; t++) {
                        metaBuf.append(m.getArrayValue(nodeNr * nTypes + t));
                        if (t < nTypes - 1) metaBuf.append(",");
                    }
                    metaBuf.append("}");
                } else {
                    // Fallback to standard single-value BEAST logger
                    metaBuf.append(id).append("=").append(m.getArrayValue(nodeNr));
                }

                if (i < metadataList.size() - 1) metaBuf.append(",");
            }
            metaBuf.append("]");
        }

        buf.append(metaBuf);
        buf.append(":");
        buf.append(node.getLength());
        return buf.toString();
    }
}