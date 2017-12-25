package cg;

import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblem;
import model.SNDRC;

public class SNDRCPricingProblem extends AbstractPricingProblem<SNDRC> {

    /**
     * Create a new Pricing Problem
     *
     * @param dataModel Data model
     * @param name      Name of the pricing problem
     */
	public SNDRCPricingProblem(SNDRC dataModel,String name){
		super(dataModel,name);
	}
	
}
