package cg;

import java.util.HashSet;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.colgenMain.AbstractColumn;

import com.oracle.webservices.internal.api.databinding.DatabindingMode;

import model.SNDRC;

public final class Cycle  extends AbstractColumn<SNDRC, SNDRCPricingProblem>{

	/** Index of edges in a path**/
	public final Set<Integer> edgeIndexSet;
	public final int startTime;
	public final double cost;// parameter in the objective expression
	public final int ifForResourceBoundConstraints;//0: no; 1:yes; 2:for holding branch constraints
	public final int[] pattern; //record the using numbers of each service
	
	
	public Cycle(SNDRCPricingProblem associatedPricingProblem, boolean isArtificial, String creator, Set<Integer> edgeIndexSet,double cost,int startTime,int ifForResourceBoundConstraints,int[] pattern) {
		super(associatedPricingProblem,isArtificial,creator);
		this.edgeIndexSet=edgeIndexSet;
		this.cost=cost;
		this.startTime=startTime;
		this.ifForResourceBoundConstraints=ifForResourceBoundConstraints;
		
		this.pattern=new int[pattern.length];
		System.arraycopy(pattern,0, this.pattern, 0, pattern.length);
	}
	

	
	@Override
	public boolean equals(Object o) {
		if(this==o) {return true;}else {
			if(!(o instanceof Cycle)) {
				return false;
			}
		}
		
		Cycle other=(Cycle)o;
		return this.edgeIndexSet.equals(other.edgeIndexSet)&& (this.isArtificialColumn==other.isArtificialColumn)&&(this.associatedPricingProblem==other.associatedPricingProblem)&&(this.ifForResourceBoundConstraints==other.ifForResourceBoundConstraints);
		
	}
	
	@Override
	public int hashCode() {
		return edgeIndexSet.hashCode();
	}
	
	@Override
	public String toString() {
		return "artificial: "+isArtificialColumn+" set: "+edgeIndexSet.toString()+" start node= "+associatedPricingProblem.originNodeO+" start time= "+this.startTime+" capacity type= "+associatedPricingProblem.capacityTypeS;
	}
	
//	public Cycle copy(Cycle cycle) {
//		Set<Integer> edgeIndexSet=new HashSet<>();
//		
//		for(int edgeIndex:cycle.edgeIndexSet) {
//			edgeIndexSet.add(edgeIndex);
//		}
//		
//		Cycle newCycle=new Cycle(cycle.associatedPricingProblem, cycle.isArtificialColumn, cycle.creator, edgeIndexSet, cycle.cost, cycle.startTime, cycle.ifForResourceBoundConstraints)
//	}
	
}
