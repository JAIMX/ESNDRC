package bap.branching;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    private double branchValue;
    private SNDRCPricingProblem associatedPricingProblem;

    public BranchOnHoldingEdge(SNDRC modelData, List<SNDRCPricingProblem> pricingProblems, double thresholdValue) {
        super(modelData, pricingProblems);
        this.thresholdValue = thresholdValue;
    }

    @Override
    protected boolean canPerformBranching(List<Cycle> solution) {
        // Reset values
        branchValue = 0;

        // calculate the number of holding arcs starting from each time point
        Map<SNDRCPricingProblem,Double[]> holdingArcTotalValue=new HashMap<>();
        for(SNDRCPricingProblem associatedPricingProblem:pricingProblems){
        	Double[] temp=new Double[dataModel.timePeriod];
        	
            for (int t = 0; t < dataModel.timePeriod; t++) {
                temp[t] = (double) 0;
            }
            holdingArcTotalValue.put(associatedPricingProblem, temp);
        }


        for (Cycle cycle : solution) {
        	SNDRCPricingProblem associatedPricingProblem=cycle.associatedPricingProblem;
        	Double[] temp=holdingArcTotalValue.get(associatedPricingProblem);
        	
            for (int edgeIndex : cycle.edgeIndexSet) {
                Edge edge = dataModel.edgeSet.get(edgeIndex);
                if (edge.edgeType == 1) { // holding arc
                    temp[edge.t1] += cycle.value;
                }
            }
        }

        // check if all values are integers
        boolean isAllInteger = true;
        double bestDifference = 1;

        // Select the time point closest to threshold value
//        for (int time = 0; time < dataModel.timePeriod; time++) {
//            double value = holdingArcTotalValue[time];
//
//            if (MathProgrammingUtil.isFractional(value)) {
//                isAllInteger = false;
//                double decimalPart = value - Math.floor(value);
//                if (Math.abs(thresholdValue - decimalPart) < bestDifference) {
//                    branchTime = time;
//                    branchValue = value;
//                    bestDifference = Math.abs(thresholdValue - decimalPart);
//                }
//            }
//        }
        
        // Select the time point closest to threshold value
        for(SNDRCPricingProblem pricingProblem:pricingProblems){
        	Double[] temp=holdingArcTotalValue.get(pricingProblem);
        	
			for (int time = 0; time < dataModel.timePeriod; time++) {
				double value = temp[time];

				if (MathProgrammingUtil.isFractional(value)) {
					isAllInteger = false;
					double decimalPart = value - Math.floor(value);
					if (Math.abs(thresholdValue - decimalPart) < bestDifference) {
						this.associatedPricingProblem=pricingProblem;
						branchTime = time;
						branchValue = value;
						bestDifference = Math.abs(thresholdValue - decimalPart);
					}
				}
			}
        	
        	
        }

        return (!isAllInteger);

    }

    protected List<BAPNode<SNDRC, Cycle>> getBranches(BAPNode<SNDRC, Cycle> parentNode) {

        // Branch 1:round down to the nearest integer
        RoundHoldingEdge branchingDecision1 = new RoundHoldingEdge(0, branchTime, branchValue,associatedPricingProblem);
        BAPNode<SNDRC, Cycle> node1 = this.createBranch(parentNode, branchingDecision1, parentNode.getSolution(),
                parentNode.getInequalities());

        // Branch 2:round up to the nearest integer
        RoundHoldingEdge branchingDecision2 = new RoundHoldingEdge(1, branchTime, branchValue,associatedPricingProblem);
        BAPNode<SNDRC, Cycle> node2 = this.createBranch(parentNode, branchingDecision2, parentNode.getSolution(),
                parentNode.getInequalities());

        return Arrays.asList(node2, node1);
    }

}
