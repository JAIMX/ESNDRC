package cg;

import java.util.Set;

import org.jorlib.frameworks.columnGeneration.colgenMain.AbstractColumn;
import model.SNDRC;

public final class Cycle  extends AbstractColumn<SNDRC, SNDRCPricingProblem>{

	/** Index of edges in a path**/
	public final Set<Integer> edgeIndexSet;
	public final int startTime;
	public final double cost;// parameter in the objective expression
	
	
	public Cycle(SNDRCPricingProblem associatedPricingProblem, boolean isArtificial, String creator, Set<Integer> edgeIndexSet,double cost,int startTime) {
		super(associatedPricingProblem,isArtificial,creator);
		this.edgeIndexSet=edgeIndexSet;
		this.cost=cost;
		this.startTime=startTime;
	}
	
	@Override
	public boolean equals(Object o) {
		if(this==o) {return true;}else {
			if(!(o instanceof Cycle)) {
				return false;
			}
		}
		
		Cycle other=(Cycle)o;
		return this.edgeIndexSet.equals(other.edgeIndexSet)&& (this.isArtificialColumn==other.isArtificialColumn)&&(this.associatedPricingProblem==other.associatedPricingProblem)&&(this.startTime==other.startTime);
		
	}
	
	@Override
	public int hashCode() {
		return edgeIndexSet.hashCode();
	}
	
	@Override
	public String toString() {
		return "artificial: "+isArtificialColumn+" set: "+edgeIndexSet.toString()+" start node= "+associatedPricingProblem.originNodeO+" start time= "+this.startTime;
	}
	
}
