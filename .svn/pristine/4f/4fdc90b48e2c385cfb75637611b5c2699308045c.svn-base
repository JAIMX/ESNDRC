package bap.branching.branchingDecisions;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

import cg.Cycle;
import model.SNDRC;

public class RoundServiceEdge implements BranchingDecision<SNDRC, Cycle>{
	
	public int roundUpOrDown;//0:round down;1:round up
	public int branchEdgeIndex;
	public double branchEdgeValue;
	
	public RoundServiceEdge(int roundUpOrDown, int branchEdgeIndex, double branchEdgeValue) {
		this.roundUpOrDown=roundUpOrDown;
		this.branchEdgeIndex=branchEdgeIndex;
		this.branchEdgeValue=branchEdgeValue;
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
			return "Round service edge(edgeIndex: "+branchEdgeIndex+") down";
		}else {
			return "Round service edge(edgeIndex: "+branchEdgeIndex+") up";
		}
		
	}
	

}
