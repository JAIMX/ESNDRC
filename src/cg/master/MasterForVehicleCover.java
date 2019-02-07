package cg.master;

import java.util.*;

import javax.naming.TimeLimitExceededException;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import bap.branching.branchingDecisions.RoundHoldingEdge;
import bap.branching.branchingDecisions.RoundLocalService;
import bap.branching.branchingDecisions.RoundLocalServiceForAllPricingProblems;
import bap.branching.branchingDecisions.RoundQ;
import bap.branching.branchingDecisions.RoundServiceEdge;
import bap.branching.branchingDecisions.RoundServiceEdgeForAllPricingProblems;
import bap.branching.branchingDecisions.RoundTimeService;
import bap.branching.branchingDecisions.RoundTimeServiceForAllPricingProblems;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import cg.master.cuts.StrongInequality;
import cg.master.cuts.StrongInequalityGenerator;
import ilog.concert.*;
import ilog.cplex.*;
import ilog.cplex.IloCplex.UnknownObjectException;
import model.SNDRC;
import model.SNDRC.Edge;

public final class MasterForVehicleCover extends Master {
	
    public MasterForVehicleCover(SNDRC dataModel, List<SNDRCPricingProblem> pricingProblems,
			CutHandler<SNDRC, SNDRCMasterData> cutHandler, StrongInequalityGenerator cutGen,
			boolean ifContainsCutAtRootNode) {
		super(dataModel, pricingProblems, cutHandler, cutGen, ifContainsCutAtRootNode);
	}
    
    
    @Override
    protected SNDRCMasterData buildModel() {
    	IloCplex cplex = null;

        try {

            cplex = new IloCplex();

            cplex.setOut(null);
            cplex.setParam(IloCplex.IntParam.Threads, config.MAXTHREADS);
            cplex.setParam(IloCplex.Param.Simplex.Tolerances.Markowitz, 0.1);
            List<Map<Integer, IloNumVar>> x; // map:edgeIndex, x variable
            IloNumVar[][] q;

            // Define variables x
            // x = new IloNumVar[dataModel.numDemand][dataModel.numArc];
            x = new ArrayList<>();
            for (int p = 0; p < dataModel.numDemand; p++) {
                Map<Integer, IloNumVar> initialX = new HashMap<Integer, IloNumVar>();
                x.add(initialX);
            }

            // add x variables
            for (int p = 0; p < dataModel.numDemand; p++) {
                for (int edgeIndex : dataModel.edgesForX.get(p)) {
                    Edge edge = dataModel.edgeSet.get(edgeIndex);
                    IloNumVar varX = cplex.numVar(0, dataModel.demandSet.get(p).volume,
                            "x" + p + "," + edge.start + "," + edge.end);
                    x.get(p).put(edgeIndex, varX);
                }
            }

            // Define the objective
            /**
             * Here we assume the cost of edge AT is 0
             */
            IloLinearNumExpr exprObj = cplex.linearNumExpr();
            // for (int p = 0; p < dataModel.numDemand; p++) {
            // for (int edgeIndex = 0; edgeIndex < dataModel.numServiceArc;
            // edgeIndex++) {
            // exprObj.addTerm(dataModel.beta *
            // dataModel.edgeSet.get(edgeIndex).duration,
            // x[p][edgeIndex]);
            // }
            // }

//            for (int p = 0; p < dataModel.numDemand; p++) {
//                Map<Integer, IloNumVar> map = x.get(p);
//                for (int edgeIndex : map.keySet()) {
//                    exprObj.addTerm(dataModel.beta * dataModel.edgeSet.get(edgeIndex).duration, map.get(edgeIndex));
//                }
//            }

            obj = cplex.addMinimize(exprObj);

            // Define flowBalanceConstraints
            flowBalanceConstraints = new IloRange[dataModel.numDemand][dataModel.abstractNumNode];

            IloLinearNumExpr expr = cplex.linearNumExpr();
            for (int p = 0; p < dataModel.numDemand; p++) {
                Map<Integer, IloNumVar> map = x.get(p);

                for (int i = 0; i < dataModel.abstractNumNode; i++) {
                    expr.clear();
                    // edges which point from i
                    for (int edgeIndex : dataModel.pointToEdgeSet.get(i)) {
                        if (map.containsKey(edgeIndex)) {
                            expr.addTerm(1, map.get(edgeIndex));
                        }
                    }

                    // edges which point to i
                    for (int edgeIndex : dataModel.pointFromEdgeSet.get(i)) {
                        if (map.containsKey(edgeIndex)) {
                            expr.addTerm(-1, map.get(edgeIndex));
                        }
                    }
//                    flowBalanceConstraints[p][i] = cplex.addEq(dataModel.b[p][i], expr);
                    flowBalanceConstraints[p][i] = cplex.addLe(0,expr);

                }
            }

            // Define weakForcingConstraints
            weakForcingConstraints = new IloRange[dataModel.numServiceArc];
            for (int arcIndex = 0; arcIndex < dataModel.numServiceArc; arcIndex++) {
                expr.clear();
//                for (int p = 0; p < dataModel.numDemand; p++) {
//                    if (x.get(p).containsKey(arcIndex)) {
//                        expr.addTerm(1, x.get(p).get(arcIndex));
//                    }
//                }

//                weakForcingConstraints[arcIndex] = cplex.addGe(0, expr);
                weakForcingConstraints[arcIndex] = cplex.addGe(-dataModel.flowCover[arcIndex], expr);
            }

            // Define resourceBoundConstraints
            resourceBoundConstraints = new IloRange[dataModel.numOfCapacity][dataModel.numNode];
            q = new IloNumVar[dataModel.numOfCapacity][dataModel.numNode];

            for (int s = 0; s < dataModel.numOfCapacity; s++) {
                for (int o = 0; o < dataModel.numNode; o++) {
                    q[s][o] = cplex.numVar(0, dataModel.vehicleLimit[s][o], "q" + s + "," + o);
                }
            }

            for (int s = 0; s < dataModel.numOfCapacity; s++) {
                for (int o = 0; o < dataModel.numNode; o++) {
                    expr.clear();
                    expr.addTerm(1, q[s][o]);
//                    resourceBoundConstraints[s][o] = cplex.addEq(dataModel.vehicleLimit[s][o], expr);
                    resourceBoundConstraints[s][o] = cplex.addGe(dataModel.vehicleLimit[s][o], expr);
                }
            }
            
            
  //----------------------------------------------modify resource bound constraints------------------------------------------------//
//            resourceBoundConstraints = new IloRange[dataModel.numOfCapacity][dataModel.numNode];
//            
//            for (int s = 0; s < dataModel.numOfCapacity; s++) {
//                for (int o = 0; o < dataModel.numNode; o++) {
//                    expr.clear();
//                    resourceBoundConstraints[s][o] = cplex.addGe(dataModel.vehicleLimit[s][o], expr);
//                }
//            }
            
  //----------------------------------------------modify resource bound constraints------------------------------------------------//

            // q branching constraints
            if (qBranchingSet != null) {
                qBranchingconstraints = new HashMap<>();
                // add q branching constraints to model
                for (RoundQ branch : qBranchingSet) {
                    if (branch.roundUpOrDown == 0) { // round down
                        try {
                            expr = cplex.linearNumExpr();
                            expr.addTerm(1,
                                    q[branch.associatedPricingProblem.capacityTypeS][branch.associatedPricingProblem.originNodeO]);
                            IloRange qBranching = cplex.addGe(Math.floor(branch.qValue), expr,
                                    "round down q[" + branch.associatedPricingProblem.capacityTypeS + "]["
                                            + branch.associatedPricingProblem.originNodeO + "]");
                            qBranchingconstraints.put(branch, qBranching);
                        } catch (IloException e) {
                            // TODO Auto-generated catch block
                            e.printStackTrace();
                        }

                    }

                    if (branch.roundUpOrDown == 1) {// round up
                        try {
                            expr = cplex.linearNumExpr();
                            expr.addTerm(1,
                                    q[branch.associatedPricingProblem.capacityTypeS][branch.associatedPricingProblem.originNodeO]);
                            IloRange qBranching = cplex.addLe(Math.ceil(branch.qValue), expr,
                                    "round up q[" + branch.associatedPricingProblem.capacityTypeS + "]["
                                            + branch.associatedPricingProblem.originNodeO + "]");
                            qBranchingconstraints.put(branch, qBranching);
                        } catch (IloException e) {
                            // TODO Auto-generated catch block
                            e.printStackTrace();
                        }

                    }

                }
            }

            // service edge branching constraints
            if (serviceEdgeBrachingSet != null) {
                for (RoundServiceEdge serviceEdgeBranch : serviceEdgeBrachingSet) {
                    if (serviceEdgeBranch.roundUpOrDown == 0) {// round down
                        expr = cplex.linearNumExpr();
                        IloRange edgeBranching = cplex.addGe(Math.floor(serviceEdgeBranch.branchEdgeValue), expr,
                                "round down service edge: " + serviceEdgeBranch.branchEdgeIndex + "originNode: "
                                        + serviceEdgeBranch.pricingProblem.originNodeO + "capacityType: "
                                        + serviceEdgeBranch.pricingProblem.capacityTypeS);
                        ServiceEdgeBranchingConstraints.put(serviceEdgeBranch, edgeBranching);
                    } else { // round up
                        expr = cplex.linearNumExpr();
                        IloRange edgeBranching = cplex.addLe(Math.ceil(serviceEdgeBranch.branchEdgeValue), expr,
                                "round up service edge: " + serviceEdgeBranch.branchEdgeIndex + "originNode: "
                                        + serviceEdgeBranch.pricingProblem.originNodeO + "capacityType: "
                                        + serviceEdgeBranch.pricingProblem.capacityTypeS);
                        ServiceEdgeBranchingConstraints.put(serviceEdgeBranch, edgeBranching);
                    }
                }
            }

            // service edge branching for all pricing problems constraints
            if (serviceEdge4AllBrachingSet != null) {
                for (RoundServiceEdgeForAllPricingProblems serviceEdgeBranch : serviceEdge4AllBrachingSet) {
                    if (serviceEdgeBranch.roundUpOrDown == 0) {// round down
                        expr = cplex.linearNumExpr();
                        IloRange edgeBranching = cplex.addGe(Math.floor(serviceEdgeBranch.branchEdgeValue), expr,
                                "round down service edge: " + serviceEdgeBranch.branchEdgeIndex);
                        serviceEdge4AllBranchingConstraints.put(serviceEdgeBranch, edgeBranching);
                    } else { // round up
                        expr = cplex.linearNumExpr();
                        IloRange edgeBranching = cplex.addLe(Math.ceil(serviceEdgeBranch.branchEdgeValue), expr,
                                "round up service edge: " + serviceEdgeBranch.branchEdgeIndex);
                        serviceEdge4AllBranchingConstraints.put(serviceEdgeBranch, edgeBranching);
                    }
                }
            }

            // local service branching constraints
            if (localServiceBranchingSet != null) {
                for (RoundLocalService localServiceBranch : localServiceBranchingSet) {
                    if (localServiceBranch.roundUpOrDown == 0) {// round down
                        expr = cplex.linearNumExpr();
                        IloRange serviceBranching = cplex.addGe(Math.floor(localServiceBranch.branchValue), expr,
                                "round down local service: " + localServiceBranch.localServiceIndex);
                        localServiceBranchingConstraints.put(localServiceBranch, serviceBranching);
                    } else { // round up
                        expr = cplex.linearNumExpr();
                        IloRange serviceBranching = cplex.addLe(Math.ceil(localServiceBranch.branchValue), expr,
                                "round up local service: " + localServiceBranch.localServiceIndex);
                        localServiceBranchingConstraints.put(localServiceBranch, serviceBranching);
                    }
                }
            }

            // holding edge branching constraints;
            if (holdingEdgeBranchingSet != null) {
                for (RoundHoldingEdge branch : holdingEdgeBranchingSet) {
                    if (branch.roundUpOrDown == 0) { // round down
                        expr = cplex.linearNumExpr();
                        IloRange holdingBranching = cplex.addGe(Math.floor(branch.branchValue), expr,
                                "round down holding at time: " + branch.branchTime);
                        holdingEdgeBranchingConstraints.put(branch, holdingBranching);
                    } else { // round up
                        expr = cplex.linearNumExpr();
                        IloRange holdingBranching = cplex.addLe(Math.ceil(branch.branchValue), expr,
                                "round up holding at time: " + branch.branchTime);
                        holdingEdgeBranchingConstraints.put(branch, holdingBranching);
                    }
                }
            }

            // local service 4 all branching constraints
            if (localService4AllBranchingSet != null) {
                for (RoundLocalServiceForAllPricingProblems localServiceBranch4All : localService4AllBranchingSet) {
                    if (localServiceBranch4All.roundUpOrDown == 0) {// round
                                                                    // down
                        expr = cplex.linearNumExpr();
                        IloRange serviceBranching = cplex.addGe(Math.floor(localServiceBranch4All.branchValue), expr,
                                "round down local service: " + localServiceBranch4All.localServiceIndex);
                        localService4AllBranchingConstraints.put(localServiceBranch4All, serviceBranching);
                    } else { // round up
                        expr = cplex.linearNumExpr();
                        IloRange serviceBranching = cplex.addLe(Math.ceil(localServiceBranch4All.branchValue), expr,
                                "round up local service: " + localServiceBranch4All.localServiceIndex);
                        localService4AllBranchingConstraints.put(localServiceBranch4All, serviceBranching);
                    }
                }
            }
            
            // time serivice branching constraints
            if(timeServiceBranchingSet!=null){
                for(RoundTimeService timeServiceBranch:timeServiceBranchingSet){
                    if(timeServiceBranch.roundUpOrDown==0){ 
                        // round down
                        expr=cplex.linearNumExpr();
                        IloRange branch=cplex.addGe(Math.floor(timeServiceBranch.branchValue), expr,"round down time service: "+timeServiceBranch.localServiceIndex);
                        timeServiceBranchingConstraints.put(timeServiceBranch, branch);
                    }else{
                        // round up
                        expr=cplex.linearNumExpr();
                        IloRange branch=cplex.addLe(Math.ceil(timeServiceBranch.branchValue), expr,"round up time service: "+timeServiceBranch.localServiceIndex);
                        timeServiceBranchingConstraints.put(timeServiceBranch, branch);
                    }
                }
            }
            
            
            // time serivice 4 all branching constraints
            if(timeService4AllBranchingSet!=null){
                for(RoundTimeServiceForAllPricingProblems timeService4AllBranch:timeService4AllBranchingSet){
                    if(timeService4AllBranch.roundUpOrDown==0){ 
                        // round down
                        expr=cplex.linearNumExpr();
                        IloRange branch=cplex.addGe(Math.floor(timeService4AllBranch.branchValue), expr,"round down time service 4 all: "+timeService4AllBranch.localServiceIndex);
                        timeService4AllBranchingConstraints.put(timeService4AllBranch, branch);
                    }else{
                        // round up
                        expr=cplex.linearNumExpr();
                        IloRange branch=cplex.addLe(Math.ceil(timeService4AllBranch.branchValue), expr,"round up time service: "+timeService4AllBranch.localServiceIndex);
                        timeService4AllBranchingConstraints.put(timeService4AllBranch, branch);
                    }
                }
            }

            Map<SNDRCPricingProblem, OrderedBiMap<Cycle, IloNumVar>> varMap = new LinkedHashMap<>();
            for (SNDRCPricingProblem pricingProblem : pricingProblems)
                varMap.put(pricingProblem, new OrderedBiMap<>());

            long time = System.currentTimeMillis();

            // Create a new data object which will store information from the
            // master. This
            // object automatically be passed to the CutHandler class.
            return new SNDRCMasterData(cplex, pricingProblems, varMap, qBranchingconstraints,
                    ServiceEdgeBranchingConstraints, x, q, serviceEdge4AllBranchingConstraints,
                    localServiceBranchingConstraints, holdingEdgeBranchingConstraints,
                    localService4AllBranchingConstraints,timeServiceBranchingConstraints,timeService4AllBranchingConstraints, dataModel, time);

        } catch (IloException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        return null;
    }
    
    
	@Override
    public Map<SNDRCPricingProblem, Integer> getQVarLimit() {
        return masterData.qVariableLimit;
    }
	
    @Override 
    public List<Map<Integer, Double>> getXValues() throws UnknownObjectException, IloException {
        return masterData.xValues;
    }
    
}