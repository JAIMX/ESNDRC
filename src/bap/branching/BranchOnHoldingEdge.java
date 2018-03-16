package bap.branching;

import java.util.Arrays;
import java.util.List;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import bap.branching.branchingDecisions.RoundHoldingEdge;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import model.SNDRC;
import model.SNDRC.Edge;

public class BranchOnHoldingEdge extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem> {

    private int branchTime;
    private double thresholdValue; // fractional in the specific pricing problem
    double branchValue;

    public BranchOnHoldingEdge(SNDRC modelData, List<SNDRCPricingProblem> pricingProblems, double thresholdValue) {
        super(modelData, pricingProblems);
        this.thresholdValue = thresholdValue;
    }

    @Override
    protected boolean canPerformBranching(List<Cycle> solution) {
        // Reset values
        branchValue = 0;

        // calculate the number of holding arcs starting from each time point
        Double[] holdingArcTotalValue = new Double[dataModel.timePeriod];
        for (int t = 0; t < dataModel.timePeriod; t++) {
            holdingArcTotalValue[t] = (double) 0;
        }

        for (Cycle cycle : solution) {
            for (int edgeIndex : cycle.edgeIndexSet) {
                Edge edge = dataModel.edgeSet.get(edgeIndex);
                if (edge.edgeType == 1) { // holding arc
                    holdingArcTotalValue[edge.t1] += cycle.value;
                }
            }
        }

        // check if all values are integers
        boolean isAllInteger = true;
        double bestDifference = 1;

        // Select the time point closest to threshold value
        for (int time = 0; time < dataModel.timePeriod; time++) {
            double value = holdingArcTotalValue[time];

            if (MathProgrammingUtil.isFractional(value)) {
                isAllInteger = false;
                double decimalPart = value - Math.floor(value);
                if (Math.abs(thresholdValue - decimalPart) < bestDifference) {
                    branchTime = time;
                    branchValue = value;
                    bestDifference = Math.abs(thresholdValue - decimalPart);
                }
            }
        }

        return (!isAllInteger);

    }

    protected List<BAPNode<SNDRC, Cycle>> getBranches(BAPNode<SNDRC, Cycle> parentNode) {

        // Branch 1:round down to the nearest integer
        RoundHoldingEdge branchingDecision1 = new RoundHoldingEdge(0, branchTime, branchValue);
        BAPNode<SNDRC, Cycle> node1 = this.createBranch(parentNode, branchingDecision1, parentNode.getSolution(),
                parentNode.getInequalities());

        // Branch 2:round up to the nearest integer
        RoundHoldingEdge branchingDecision2 = new RoundHoldingEdge(1, branchTime, branchValue);
        BAPNode<SNDRC, Cycle> node2 = this.createBranch(parentNode, branchingDecision1, parentNode.getSolution(),
                parentNode.getInequalities());

        return Arrays.asList(node2, node1);
    }

}
