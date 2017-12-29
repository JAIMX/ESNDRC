package cg.master;

import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import ilog.concert.*;
import ilog.cplex.*;
import model.SNDRC;
import model.SNDRC.Edge;

public final class Master extends AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> {

	private IloObjective obj;
	private IloNumVar[][] x;
	private IloRange[][] flowBalanceConstraints;
	private IloRange[] weakForcingConstraints;
	private IloRange[][] resourceBoundConstraints;
	
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
		
		
		try {
			cplex=new IloCplex();
			
			cplex.setOut(null);
			cplex.setParam(IloCplex.IntParam.Threads, config.MAXTHREADS);
			
			//Define the objective
			obj=cplex.addMinimize();
			
			//Define variables x
			x=new IloNumVar[dataModel.numDemand][dataModel.numArc];
			for(int p=0;p<dataModel.numDemand;p++) {
				for(int arc=0;arc<dataModel.numArc;arc++ ) {
					Edge edge=dataModel.edgeSet.get(arc);
					x[p][arc]=cplex.numVar(0,Double.MAX_VALUE, "x"+p+edge.u+"->"+edge.v+"t:"+edge.t1+"->"+edge.t2);
				}
			}
			
			//Define flowBalanceConstraints
			flowBalanceConstraints=new IloRange[dataModel.numDemand][dataModel.abstractNumNode];
			
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for(int p=0;p<dataModel.numDemand;p++) {
				for(int i=0;i<dataModel.abstractNumNode;i++) {
					expr.clear();
					//edges which point from i
					for(int edgeIndex:dataModel.pointToEdgeSet.get(i)) {
						expr.addTerm(1,x[p][edgeIndex]);
					}
					
					//edges which point to i
					for(int edgeIndex:dataModel.pointFromEdgeSet.get(i)) {
						expr.addTerm(-1,x[p][edgeIndex]);
					}
					flowBalanceConstraints[p][i]=cplex.addEq(dataModel.b[p][i], expr);
					
				}
			}
			
			
			// Define weakForcingConstraints
			
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
	}
	
}
