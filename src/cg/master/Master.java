package cg.master;

import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import ilog.concert.*;
import ilog.cplex.*;
import model.SNDRC;

public final class Master extends AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> {

	private IloObjective obj;
	private IloRange[] flowBalanceConstraints,weakForcingConstraints,resourceBoundConstraints;
	
	public Master(SNDRC dataModel, SNDRCPricingProblem pricingProblem){
		super(dataModel,pricingProblem,OptimizationSense.MINIMIZE);
		System.out.println("Master constructor. Columns: "+masterData.getNrColumns());
	}
	
	/**
     * Builds the master model
     * @return Returns a MasterData object which is a data container for information coming from the master problem
     */
	@Override
	protected SNDRCMasterData buildModel(){
		IloCplex cplex=null;
		
		cplex=new IloCplex();
		cplex.setOut(null);
		cplex.setParam(IloCplex.)
		
		
		
	}
	
}
