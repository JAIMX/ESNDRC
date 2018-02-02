package cg;

import java.util.HashSet;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblem;
import model.SNDRC;

public final class SNDRCPricingProblem extends AbstractPricingProblem<SNDRC> {

    /**
     * Create a new Pricing Problem
     *
     * @param dataModel Data model
     * @param name      Name of the pricing problem
     */
	
	public final int capacityTypeS;
	public final int originNodeO;
	public Set<Cycle> fixCycleSet;
	
	public SNDRCPricingProblem(SNDRC modelData,String name,int capacityTypeS,int originNodeO){
		super(modelData,name);
		this.capacityTypeS=capacityTypeS;
		this.originNodeO=originNodeO;
		fixCycleSet=new HashSet<>();
	}
	
}
