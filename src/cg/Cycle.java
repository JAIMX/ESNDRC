package cg;

import java.util.Set;

import org.jorlib.frameworks.columnGeneration.colgenMain.AbstractColumn;
import model.SNDRC;

public final class Cycle  extends AbstractColumn<SNDRC, SNDRCPricingProblem>{

	/** Index of edges in a path**/
	public final Set<Integer> edgeIndexSet;
	public final int startNode,startTime;
	
	
	public Cycle(SNDRCPricingProblem associatedPricingProblem, boolean isArtificial, String creator, Set<Integer> edgeIndexSet, int startNode,int startTime) {
		super(associatedPricingProblem,isArtificial,creator);
		this.edgeIndexSet=edgeIndexSet;
		this.startNode=startNode;
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
		return this.edgeIndexSet.equals(other.edgeIndexSet)&& this.isArtificialColumn==other.isArtificialColumn;
		
	}
	
	@Override
	public int hashCode() {
		return edgeIndexSet.hashCode();
	}
	
	@Override
	public String toString() {
		return "artificial: "+isArtificialColumn+" set: "+edgeIndexSet.toString()+" start node= "+this.startNode+" start time= "+this.startTime;
	}
	
}
