package bap.branching;

import java.util.*;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import bap.branching.branchingDecisions.RoundQ;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import jdk.nashorn.internal.runtime.regexp.joni.Config;
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

//	private int capacityTypeBranching, originNodeBranching; // qso to branch
	private SNDRCPricingProblem pricingProblemForCycle = null;// qso is
	private double thresholdValue; // fractional in the specific pricing problem
	double bestQValue;

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
		bestQValue=0;
		
		//For each capacity type and origin node, determine whether there's a fractional q
		Map<SNDRCPricingProblem, Double> qValueMap=new HashMap<>();
		for (SNDRCPricingProblem pricingProblem : pricingProblems) {
			qValueMap.put(pricingProblem,
					(double) dataModel.vehicleLimit[pricingProblem.capacityTypeS][pricingProblem.originNodeO]);
		}

		// Aggregate q values
		for (Cycle cycle : solution) {
			SNDRCPricingProblem associatedPricingProblem = cycle.associatedPricingProblem;
			Double tempValue = qValueMap.get(associatedPricingProblem);
			qValueMap.put(associatedPricingProblem, tempValue - cycle.value);
		}
		
		
		boolean isAllInteger=true;
		double bestDifference=1;
		
		//Select the variable qso closest to threshold value
		for(SNDRCPricingProblem pricingProblem:pricingProblems) {
			double value=qValueMap.get(pricingProblem);
			
			if(MathProgrammingUtil.isFractional(value)) {
				isAllInteger=false;
				double decimalQ=value-Math.floor(value);
				if (Math.abs(thresholdValue - decimalQ) < bestDifference) {
					pricingProblemForCycle = pricingProblem;
					bestQValue = value;
					bestDifference=Math.abs(thresholdValue-decimalQ);
				}
			}
		}
		
	
		
		return (!isAllInteger);
		
		
	}
	
	
	/**
	 * Create the branches
	 */
	@Override
	protected List<BAPNode<SNDRC,Cycle>> getBranches(BAPNode<SNDRC,Cycle> parentNode){
//		List<Cycle> initialSolution=new ArrayList<Cycle>();
		
		//Branch 1:round q down to the nearest integer
		RoundQ branchingDecision1=new RoundQ(0,pricingProblemForCycle,bestQValue);
		BAPNode<SNDRC,Cycle> node1=this.createBranch(parentNode, branchingDecision1, parentNode.getSolution(), parentNode.getInequalities());
//		BAPNode<SNDRC,Cycle> node1=this.createBranch(parentNode, branchingDecision1,initialSolution , parentNode.getInequalities());
		
		
//		initialSolution=new ArrayList<Cycle>();
		//Branch 2:round q up to the nearest integer
		RoundQ branchingDecision2=new RoundQ(1,pricingProblemForCycle,bestQValue);
		BAPNode<SNDRC,Cycle> node2=this.createBranch(parentNode, branchingDecision2, parentNode.getSolution(), parentNode.getInequalities());
//		BAPNode<SNDRC,Cycle> node2=this.createBranch(parentNode, branchingDecision2, initialSolution, parentNode.getInequalities());
		
		return Arrays.asList(node2,node1);
	}
}