package cg.master;

import java.util.*;

import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import bap.branching.branchingDecisions.RoundQ;
import bap.branching.branchingDecisions.RoundServiceEdge;
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
	
	//branch on q variables	
	public final Set<RoundQ> qBranchingSet;
	public final Map<RoundQ,IloRange> qBranchingconstraints;
	
	//branch on service edge
	public final Set<RoundServiceEdge> serviceEdgeBrachingSet;
	public final Map<RoundServiceEdge,IloRange> ServiceEdgeBranchingConstraints;
	
	
	public SNDRCMasterData(IloCplex cplex,List<SNDRCPricingProblem> pricingProblems,Map<SNDRCPricingProblem, OrderedBiMap<Cycle, IloNumVar>> varMap,Map<RoundQ,IloRange> qBranchingconstraints,Map<RoundServiceEdge,IloRange> ServiceEdgeBranchingConstraints){
		super(varMap);
		this.cplex=cplex;
		this.pricingProblems=pricingProblems;
		
		this.qBranchingSet=new HashSet<>();
		this.qBranchingconstraints=new HashMap<>();
		
		if(qBranchingconstraints!=null) {
			for(RoundQ roundQ:qBranchingconstraints.keySet()) {
				this.qBranchingSet.add(roundQ);
				this.qBranchingconstraints.put(roundQ, qBranchingconstraints.get(roundQ));
			}
		}

		
		this.serviceEdgeBrachingSet=new HashSet<>();
		this.ServiceEdgeBranchingConstraints=new HashMap<>();
		if(ServiceEdgeBranchingConstraints!=null) {
			for(RoundServiceEdge roundServiceEdge:ServiceEdgeBranchingConstraints.keySet()) {
				this.serviceEdgeBrachingSet.add(roundServiceEdge);
				this.ServiceEdgeBranchingConstraints.put(roundServiceEdge, ServiceEdgeBranchingConstraints.get(roundServiceEdge));
			}
		}

	}

	
}
