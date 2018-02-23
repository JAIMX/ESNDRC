package bap.bapNodeComparators;

import java.util.Comparator;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;

public class NodeBoundbapNodeComparatorForLB implements Comparator<BAPNode>{

    @Override
    public int compare(BAPNode o1, BAPNode o2) {
//        return -Integer.compare(o1.nodeID, o2.nodeID);
        return Double.compare(o1.getBound(), o2.getBound());
    }
}
