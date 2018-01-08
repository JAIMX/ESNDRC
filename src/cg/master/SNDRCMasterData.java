package cg.master;

import java.util.*;

import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import bap.branching.branchingDecisions.RoundQ;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
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
	public final  List<SNDRCPricingProblem> pricingProblems;
//	public Map<RoundQ,IloRange> qBranchingconstraints;
	public Set<RoundQ> qBranchingSet;
	public Map<SNDRCPricingProblem,IloNumVar> qVaribles;
	
	
	public SNDRCMasterData(IloCplex cplex,List<SNDRCPricingProblem> pricingProblems,Map<SNDRCPricingProblem, OrderedBiMap<Cycle, IloNumVar>> varMap,Map<SNDRCPricingProblem,IloNumVar> qVaribles){
		super(varMap);
		this.cplex=cplex;
		this.pricingProblems=pricingProblems;
//		qBranchingconstraints=new LinkedHashMap<>();
		qBranchingSet=new HashSet<>();
		this.qVaribles=qVaribles;
	}

	
}
