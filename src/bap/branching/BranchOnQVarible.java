package bap.branching;

import java.util.*;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import bap.branching.branchingDecisions.RoundQ;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import model.SNDRC;

/**
 * Class which creates new branches in the branch-and-price tree. This
 * particular class branches on the variable q.More precisely, the class checks
 * whether there is a fractional q in the solution. We branch by identifying
 * variable q with fractional value closest to a threshold value and round it
 * down and up to the nearest integer.
 * 
 * @author sxx
 *
 */
public class BranchOnQVarible extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem> {

	private int capacityTypeBranching, originNodeBranching; // qso to branch
	private SNDRCPricingProblem pricingProblemForCycle = null;// qso is
	private double thresholdValue;														// fractional in
																// the specific
																// pricing
																// problem

	public BranchOnQVarible(SNDRC modelData, List<SNDRCPricingProblem> pricingProblems,double thresholdValue) {
		super(modelData, pricingProblems);
		this.thresholdValue=thresholdValue;
	}


	/**
	 * Determine on which q varible from the specific capacity type and
	 * originNode we are going to branch.
	 * 
	 * @param solution
	 *            Fractional column generation solution
	 * @return true if a fractional q exists
	 */
	@Override
	protected boolean canPerformBranching(List<Cycle> solution) {
		//Reset values
		pricingProblemForCycle=null;
		capacityTypeBranching=-1;
		originNodeBranching=-1;
		double bestQValue=0;
		
		//For each capacity type and origin node, determine whether there's a fractional q
		Map<SNDRCPricingProblem, Double> qValueMap=new HashMap<>();
		for(SNDRCPricingProblem pricingProblem: pricingProblems){
			qValueMap.put(pricingProblem, (double)dataModel.vehicleLimit[pricingProblem.capacityTypeS][pricingProblem.originNodeO]);
		}
		
		//Aggregate q values
		for(Cycle cycle:solution){
			SNDRCPricingProblem associatedPricingProblem=cycle.associatedPricingProblem;
			Double tempValue=qValueMap.get(associatedPricingProblem);
			qValueMap.put(associatedPricingProblem, tempValue-cycle.value);
		}
		
		//Select the variable qso closest to threshold value
		for(SNDRCPricingProblem pricingProblem:pricingProblems){
			double value=qValueMap.get(pricingProblem);
			if(Math.abs(thresholdValue-value)<Math.abs(thresholdValue-bestQValue)){
				pricingProblemForCycle=pricingProblem;
				capacityTypeBranching=pricingProblem.capacityTypeS;
				originNodeBranching=pricingProblem.originNodeO;
				bestQValue=value;
			}
		}
		
		return MathProgrammingUtil.isFractional(bestQValue);
		
		
	}
	
	
	/**
	 * Create the branches
	 */
	@Override
	protected List<BAPNode<SNDRC,Cycle>> getBranches(BAPNode<SNDRC,Cycle> parentNode){
		//Branch 1:round q down to the nearest integer
		RoundQ branchingDecision1=new RoundQ(0,capacityTypeBranching,originNodeBranching);
		BAPNode<SNDRC,Cycle> node1=this.createBranch(parentNode, branchingDecision1, parentNode.getSolution(), parentNode.getInequalities());
		
		//Branch 2:round q up to the nearest integer
		RoundQ branchingDecision2=new RoundQ(1,capacityTypeBranching,originNodeBranching);
		BAPNode<SNDRC,Cycle> node2=this.createBranch(parentNode, branchingDecision2, parentNode.getSolution(), parentNode.getInequalities());
		
		return Arrays.asList(node2,node1);
	}
}