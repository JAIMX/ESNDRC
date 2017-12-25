package cg.master;

import java.util.Map;

import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import model.SNDRC;

/**
 * Container which stores information coming from the master problem. It contains:
 * a reference to the cplex model
 * reference to the pricing problem
 * @author sxx
 *
 */
public class SNDRCMasterData extends MasterData<SNDRC,Cycle,SNDRCPricingProblem,IloNumVar>{
	public final IloCplex cplex;
	
	
	public SNDRCMasterData(IloCplex cplex,Map<SNDRCPricingProblem, OrderedBiMap<Cycle, IloNumVar>> varMap){
		super(varMap);
		this.cplex=cplex;
	}

	
}
