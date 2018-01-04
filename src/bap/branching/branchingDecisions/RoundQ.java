package bap.branching.branchingDecisions;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import model.SNDRC;

/**
 * Round the branching variable q up and down to the nearest integer
 * @author Helen
 *
 */
public class RoundQ implements BranchingDecision<SNDRC, Cycle> {
	
	public int roundUpOrDown;//0:round down;1:round up
//	public int branchingCapacityType,branchingOriginNode;
	public SNDRCPricingProblem associatedPricingProblem;
	public double qValue;
	
	public RoundQ(int roundUpOrDown,SNDRCPricingProblem associatedPricingProblem,double qValue) {
		this.roundUpOrDown=roundUpOrDown;
//		this.branchingCapacityType=branchingCapacityType;
//		this.branchingOriginNode=branchingOriginNode;
		this.associatedPricingProblem=associatedPricingProblem;
		this.qValue=qValue;
	}
	

	@Override
	public boolean inEqualityIsCompatibleWithBranchingDecision(AbstractInequality inequality) {
		return true;
	}
	
	@Override
    public boolean columnIsCompatibleWithBranchingDecision(Cycle cycle) {
		return true;
	}
	
	@Override
	public String toString() {
		if(roundUpOrDown==0) {
			return "Round variable q(capcityType: "+associatedPricingProblem.capacityTypeS+" originNode: "+associatedPricingProblem.originNodeO+") down";
		}else {
			return "Round variable q(capcityType: "+associatedPricingProblem.capacityTypeS+" originNode: "+associatedPricingProblem.originNodeO+") up";
		}
		
	}
	
}
