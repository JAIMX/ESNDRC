package cg.master.cuts;

import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractCutGenerator;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

public final class StrongInequality extends AbstractInequality {

	public final int edgeIndex;
	public final int commodity;
	
	public StrongInequality(AbstractCutGenerator maintainingGenerator, int edgeIndex,int commodity) {
		super(maintainingGenerator);
		this.edgeIndex=edgeIndex;
		this.commodity=commodity;
	}
	
	@Override
	public boolean equals(Object o) {
		if(this==o)
			return true;
		else if(!(o instanceof StrongInequality))
			return false;
		StrongInequality other=(StrongInequality)o;
		return (this.edgeIndex==other.edgeIndex&&this.commodity==other.commodity);
	}
	
	@Override
	public int hashCode() {
		return this.edgeIndex*19+this.commodity;
	}
	
	@Override 
	public String toString() {
		return "edgeIndex= "+edgeIndex+" commodity= "+commodity;
	}
}
