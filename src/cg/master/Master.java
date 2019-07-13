package cg.master;

import java.util.*;

import javax.naming.NoInitialContextException;
import javax.naming.TimeLimitExceededException;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import com.sun.media.sound.ModelDestination;

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
import model.SNDRC.Demand;
import model.SNDRC.Edge;

public class Master extends AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> {

    protected IloObjective obj;
    // private IloNumVar[][] x;
    // private List<Set<Integer>> edgesForX;

    protected  IloRange[][] flowBalanceConstraints;
    protected  IloRange[] weakForcingConstraints;
    protected  IloRange[][] resourceBoundConstraints;
    protected  IloRange[][] storeBoundConstraints;
    protected  IloRange[][] chargeBoundConstraints;
    

    // branch on q variables
    protected  Set<RoundQ> qBranchingSet;
    protected  Map<RoundQ, IloRange> qBranchingconstraints;

    // branch on service edge
    protected  Set<RoundServiceEdge> serviceEdgeBrachingSet;
    protected  Map<RoundServiceEdge, IloRange> ServiceEdgeBranchingConstraints;

    // branch on service edge for all pricing problems
    protected  Set<RoundServiceEdgeForAllPricingProblems> serviceEdge4AllBrachingSet;
    protected  Map<RoundServiceEdgeForAllPricingProblems, IloRange> serviceEdge4AllBranchingConstraints;

    // branch on local service
    protected  Set<RoundLocalService> localServiceBranchingSet;
    protected  Map<RoundLocalService, IloRange> localServiceBranchingConstraints;

    // branch on holding edges
    protected  Set<RoundHoldingEdge> holdingEdgeBranchingSet;
    protected  Map<RoundHoldingEdge, IloRange> holdingEdgeBranchingConstraints;

    // branch on local service 4 all
    protected  Set<RoundLocalServiceForAllPricingProblems> localService4AllBranchingSet;
    protected  Map<RoundLocalServiceForAllPricingProblems, IloRange> localService4AllBranchingConstraints;
    
    //branch on time service
    protected  Set<RoundTimeService> timeServiceBranchingSet;
    protected  Map<RoundTimeService, IloRange> timeServiceBranchingConstraints;
    
    //branch on time service 4 all
    protected  Set<RoundTimeServiceForAllPricingProblems> timeService4AllBranchingSet;
    protected  Map<RoundTimeServiceForAllPricingProblems, IloRange> timeService4AllBranchingConstraints;

    // record the branches leading to a node which needs to be solved by cut
    protected  Set<BranchingDecision> branchSetNodeSolvedByCut;
    protected  StrongInequalityGenerator cutGen;
    protected  boolean ifContainsCut;

    public Master(SNDRC dataModel, List<SNDRCPricingProblem> pricingProblems,
            CutHandler<SNDRC, SNDRCMasterData> cutHandler, StrongInequalityGenerator cutGen,
            boolean ifContainsCutAtRootNode) {
        super(dataModel, pricingProblems, cutHandler, OptimizationSense.MINIMIZE);
        this.qBranchingSet = new HashSet<RoundQ>();
        qBranchingconstraints = new HashMap<>();

        this.serviceEdgeBrachingSet = new HashSet<>();
        this.ServiceEdgeBranchingConstraints = new HashMap<>();

        this.serviceEdge4AllBrachingSet = new HashSet<>();
        this.serviceEdge4AllBranchingConstraints = new HashMap<>();

        this.localServiceBranchingSet = new HashSet<>();
        this.localServiceBranchingConstraints = new HashMap<>();

        this.holdingEdgeBranchingSet = new HashSet<>();
        this.holdingEdgeBranchingConstraints = new HashMap<>();

        this.localService4AllBranchingSet = new HashSet<>();
        this.localService4AllBranchingConstraints = new HashMap<>();
        
        this.timeServiceBranchingSet=new HashSet<>();
        this.timeServiceBranchingConstraints=new HashMap<>();
        
        this.timeService4AllBranchingSet=new HashSet<>();
        this.timeService4AllBranchingConstraints=new HashMap<>();

        this.branchSetNodeSolvedByCut = new HashSet<>();
        this.cutGen = cutGen;

        this.ifContainsCut = ifContainsCutAtRootNode;

        // this.buildModel();

        // System.out.println("Master constructor. Columns: " +
        // masterData.getNrColumns());

    }

    /**
     * Builds the master model
     * 
     * @return Returns a MasterData object which is a data container for
     *         information coming from the master problem
     */
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


            for (int p = 0; p < dataModel.numDemand; p++) {
                Map<Integer, IloNumVar> map = x.get(p);
                for (int edgeIndex : map.keySet()) {
                    exprObj.addTerm(dataModel.beta * dataModel.edgeSet.get(edgeIndex).duration, map.get(edgeIndex));
                }
            }

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
                    flowBalanceConstraints[p][i] = cplex.addEq(dataModel.b[p][i], expr);

                }
            }

            // Define weakForcingConstraints
            weakForcingConstraints = new IloRange[dataModel.numServiceArc];
            for (int arcIndex = 0; arcIndex < dataModel.numServiceArc; arcIndex++) {
                expr.clear();
                for (int p = 0; p < dataModel.numDemand; p++) {
                    if (x.get(p).containsKey(arcIndex)) {
                        expr.addTerm(1, x.get(p).get(arcIndex));
                    }
                }

                weakForcingConstraints[arcIndex] = cplex.addGe(0, expr);
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
            
            //Define storeBoundConstraints
            storeBoundConstraints=new IloRange[dataModel.abstractNumNode][dataModel.timePeriod];
            
            for(int l=0;l<dataModel.abstractNumNode;l++) {
            	for(int t=0;t<dataModel.timePeriod;t++) {
            		expr.clear();
            		int nodeIndex=l*dataModel.timePeriod+t;
            		//search holding arc index
            		int holdingEdgeIndex=-1;
            		for(int i=dataModel.numServiceArc;i<dataModel.numArc;i++) {
            			Edge edge=dataModel.edgeSet.get(i);
            			if(edge.start==nodeIndex) {
            				holdingEdgeIndex=i;
            				break;
            			}
            		}
            		
            		for(int k=0;k<dataModel.numDemand;k++) {
            			Demand demand=dataModel.demandSet.get(k);
            			if(demand.destination!=l) {
            				if(x.get(k).containsKey(holdingEdgeIndex)) {
                				expr.addTerm(1, x.get(k).get(holdingEdgeIndex));
            				}
            			}
            		}
            		
            		storeBoundConstraints[l][t]=cplex.addGe(dataModel.storeLimit[l], expr);
            	}
            }
            
            //Define chargeBoundConstraints
            chargeBoundConstraints=new IloRange[dataModel.abstractNumNode][dataModel.timePeriod];
            
            for(int l=0;l<dataModel.abstractNumNode;l++) {
            	for(int t=0;t<dataModel.timePeriod;t++) {
            		expr.clear();
            		chargeBoundConstraints[l][t]=cplex.addGe(dataModel.chargeLimit[l], expr);
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
            
            // time service branching constraints
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
            
            
            // time service 4 all branching constraints
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

    /**
     * Solve the master problem
     * 
     * @param timeLimit
     *            Future point in time by which the solve procedure must be
     *            completed
     * @return true if the master problem has been solved
     * @throws TimeLimitExceededException
     *             TimeLimitExceededException
     */
    @Override
    protected boolean solveMasterProblem(long timeLimit) {

        try {
            // Set time limit
            double timeRemaining = Math.max(1, (timeLimit - System.currentTimeMillis()) / 1000.0);
            // System.out.println(masterData.cplex.toString());

            masterData.cplex.setParam(IloCplex.DoubleParam.TiLim, timeRemaining); // set
                                                                                  // time
                                                                                  // limit
                                                                                  // in
                                                                                  // seconds

            // Potentially export the model
            if (config.EXPORT_MODEL) {
                masterData.cplex.exportModel(config.EXPORT_MASTER_DIR + "master_" + this.getIterationCount() + ".lp");
                // System.out.println(masterData.cplex.toString());
            }

            // Solve the model
            if (!masterData.cplex.solve() || masterData.cplex.getStatus() != IloCplex.Status.Optimal) {
                if (masterData.cplex.getCplexStatus() == IloCplex.CplexStatus.AbortTimeLim) // Aborted
                                                                                            // due
                                                                                            // to
                                                                                            // time
                                                                                            // limit
                    throw new TimeLimitExceededException();
                else {
                    masterData.cplex.exportModel(config.EXPORT_MASTER_DIR + "check.lp");
                    throw new RuntimeException("Master problem solve failed! Status: " + masterData.cplex.getStatus());
                }

            } else {
                masterData.objectiveValue = masterData.cplex.getObjValue();

                masterData.xValues = new ArrayList<>();
                for (int commodity = 0; commodity < dataModel.numDemand; commodity++) {
                    Map tempMap = new HashMap<>();
                    for (int edgeIndex : masterData.x.get(commodity).keySet()) {
                        tempMap.put(edgeIndex, masterData.cplex.getValue(masterData.x.get(commodity).get(edgeIndex)));
                    }
                    masterData.xValues.add(tempMap);
                }
                
                double[][] temp=new double[dataModel.numOfCapacity][dataModel.numNode];
                for (int s = 0; s < dataModel.numOfCapacity; s++) {
                    for (int o = 0; o < dataModel.numNode; o++) {
                        temp[s][o] = masterData.cplex.getDual(resourceBoundConstraints[s][o]);
                    }
                }
                masterData.resourceBoundConstraintsDual=temp;

                // System.out.println("||-----------------------temp solution
                // out---------------------||");
                // this.printSolution();
                // for (int s = 0; s < dataModel.numOfCapacity; s++) {
                // for (int o = 0; o < dataModel.numNode; o++) {
                // System.out.println("q" + s + "," + o + "=" +
                // masterData.cplex.getValue(masterData.q[s][o]));
                // }
                // }
                //
                // for (int demand = 0; demand < dataModel.numDemand; demand++)
                // {
                // for (int edgeIndex = 0; edgeIndex < dataModel.numArc;
                // edgeIndex++) {
                //
                // if (masterData.x.get(demand).containsKey(edgeIndex)) {
                // if
                // (masterData.cplex.getValue(masterData.x.get(demand).get(edgeIndex))
                // > config.PRECISION) {
                // Edge edge = dataModel.edgeSet.get(edgeIndex);
                // System.out.println("x[" + demand + "]:" + edge.start + "->" +
                // edge.end + "="
                // +
                // masterData.cplex.getValue(masterData.x.get(demand).get(edgeIndex)));
                // }
                // }
                //
                // }
                // }
                // System.out.println();

            }
        } catch (IloException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (TimeLimitExceededException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        return true;
    }

    /**
     * Adds a new column to the master problem
     * 
     * @param column
     *            column to add
     */
    @Override
    public void addColumn(Cycle column) {

        try {
            // Register column with objective
            IloColumn iloColumn = masterData.cplex.column(obj, column.cost);

            // weak forcing constraints
            for (int edgeIndex : column.edgeIndexSet) {
                if (dataModel.edgeSet.get(edgeIndex).edgeType == 0) {
                    iloColumn = iloColumn.and(masterData.cplex.column(weakForcingConstraints[edgeIndex],
                            -dataModel.capacity[column.associatedPricingProblem.capacityTypeS]));
                }
            }

            // resource bound constraints
            if (!column.isArtificialColumn) {
                iloColumn = iloColumn.and(masterData.cplex.column(
                        resourceBoundConstraints[column.associatedPricingProblem.capacityTypeS][column.associatedPricingProblem.originNodeO],
                        1));
            }
            // add artificial cycles for resource bound constraints
            if (column.isArtificialColumn && column.ifForResourceBoundConstraints == 1) {
                IloRange constraint = resourceBoundConstraints[column.associatedPricingProblem.capacityTypeS][column.associatedPricingProblem.originNodeO];
                iloColumn = iloColumn.and(masterData.cplex.column(constraint, 1));
            }
            
            // charge bound constraints
            if(!column.isArtificialColumn) {
            	for(int chargeNodeIndex:column.ifCharge) {
            		int l=chargeNodeIndex/dataModel.timePeriod;
            		int t=chargeNodeIndex%dataModel.timePeriod;
            		iloColumn=iloColumn.and(masterData.cplex.column(chargeBoundConstraints[l][t],1));
            	}
            }
            
            //add artificial cycles(ifForResourceBoundConstraints=1) for charge bound constraints
            if(column.isArtificialColumn&&column.ifForResourceBoundConstraints==1) {
            	for(int t=0;t<dataModel.timePeriod;t++) {
            		iloColumn=iloColumn.and(masterData.cplex.column(chargeBoundConstraints[column.associatedPricingProblem.originNodeO][t],1));
            	}
            }

            // service edge branching constraints
            if (column.isArtificialColumn) {
                for (RoundServiceEdge serviceEdgeBranching : masterData.serviceEdgeBrachingSet) {
                    if (serviceEdgeBranching.roundUpOrDown == 1
                            && column.edgeIndexSet.contains(serviceEdgeBranching.branchEdgeIndex)) {
                        IloRange constraint = masterData.ServiceEdgeBranchingConstraints.get(serviceEdgeBranching);
                        iloColumn = iloColumn.and(masterData.cplex.column(constraint, 1));
                    }
                }
            } else { // not an artificial column
                for (RoundServiceEdge serviceEdgeBranching : masterData.serviceEdgeBrachingSet) {
                    if (column.associatedPricingProblem == serviceEdgeBranching.pricingProblem
                            && column.edgeIndexSet.contains(serviceEdgeBranching.branchEdgeIndex)) {
                        IloRange constraint = masterData.ServiceEdgeBranchingConstraints.get(serviceEdgeBranching);
                        iloColumn = iloColumn.and(masterData.cplex.column(constraint, 1));
                    }
                }
            }

            // service edge for all pricing problems branching constraints
            if (column.isArtificialColumn) {
                for (RoundServiceEdgeForAllPricingProblems serviceEdgeBranching : masterData.serviceEdge4AllBrachingSet) {
                    if (serviceEdgeBranching.roundUpOrDown == 1
                            && column.edgeIndexSet.contains(serviceEdgeBranching.branchEdgeIndex)) {
                        IloRange constraint = masterData.serviceEdge4AllBranchingConstraints.get(serviceEdgeBranching);
                        iloColumn = iloColumn.and(masterData.cplex.column(constraint, 1));
                    }
                }
            } else { // not an artificial column
                for (RoundServiceEdgeForAllPricingProblems serviceEdgeBranching : masterData.serviceEdge4AllBrachingSet) {
                    if (column.edgeIndexSet.contains(serviceEdgeBranching.branchEdgeIndex)) {
                        IloRange constraint = masterData.serviceEdge4AllBranchingConstraints.get(serviceEdgeBranching);
                        iloColumn = iloColumn.and(masterData.cplex.column(constraint, 1));
                    }
                }
            }

            // local service for all cycles branching constraints
            if (column.isArtificialColumn && column.ifForResourceBoundConstraints == 0) {
                // first kind of artificial variables
                for (RoundLocalService localServiceBranching : masterData.localServiceBranchingSet) {
                    if (localServiceBranching.roundUpOrDown == 1) {
                        for (int edgeIndex : column.edgeIndexSet) {
                            Edge edge = dataModel.edgeSet.get(edgeIndex);
                            if (edge.serviceIndex == localServiceBranching.localServiceIndex) {
                                IloRange constraint = masterData.localServiceBranchingConstraints
                                        .get(localServiceBranching);
                                iloColumn = iloColumn.and(masterData.cplex.column(constraint, 1));
                                break;
                            }
                        }
                    }
                }

            } else {
                if (!column.isArtificialColumn) {

                    for (RoundLocalService localServiceBranching : masterData.localServiceBranchingSet) {

                        if (column.associatedPricingProblem == localServiceBranching.associatedPricingProblem) {
                            IloRange constraint = masterData.localServiceBranchingConstraints
                                    .get(localServiceBranching);
                            int count = 0;
                            for (int edgeIndex : column.edgeIndexSet) {
                                Edge edge = dataModel.edgeSet.get(edgeIndex);
                                if (edge.serviceIndex == localServiceBranching.localServiceIndex)
                                    count++;
                            }
                            iloColumn = iloColumn.and(masterData.cplex.column(constraint, count));
                        }

                    }

                }
            }

            // holding edges branching constraints(we use artificial variables
            // of kind 3)
            if (column.isArtificialColumn && column.ifForResourceBoundConstraints == 2) {
                for (RoundHoldingEdge branch : masterData.holdingEdgeBranchingSet) {
                    if (branch.roundUpOrDown == 1) {
                        IloRange constraint = masterData.holdingEdgeBranchingConstraints.get(branch);
                        iloColumn = iloColumn.and(masterData.cplex.column(constraint, 1));
                    }
                }
            }
            if (!column.isArtificialColumn) {
                for (int edgeIndex : column.edgeIndexSet) {
                    Edge edge = dataModel.edgeSet.get(edgeIndex);
                    if (edge.edgeType == 1) {
                        int startTime = edge.t1;
                        for (RoundHoldingEdge branch : masterData.holdingEdgeBranchingSet) {
                            if (branch.branchTime == startTime
                                    && branch.associatedPricingProblem == column.associatedPricingProblem) {
                                IloRange constraint = masterData.holdingEdgeBranchingConstraints.get(branch);
                                iloColumn = iloColumn.and(masterData.cplex.column(constraint, 1));
                            }
                        }
                    }
                }
            }

            // local service for all pricing problems branching constraints
            if (column.isArtificialColumn && column.ifForResourceBoundConstraints == 0) {
                // first kind of artificial variables
                for (RoundLocalServiceForAllPricingProblems localServiceBranching : masterData.localService4AllBranchingSet) {
                    if (localServiceBranching.roundUpOrDown == 1) {
                        for (int edgeIndex : column.edgeIndexSet) {
                            Edge edge = dataModel.edgeSet.get(edgeIndex);
                            if (edge.serviceIndex == localServiceBranching.localServiceIndex) {
                                IloRange constraint = masterData.localService4AllBranchingConstraints
                                        .get(localServiceBranching);
                                iloColumn = iloColumn.and(masterData.cplex.column(constraint, 1));
                                break;
                            }
                        }
                    }
                }

            } else {
                if (!column.isArtificialColumn) {

                    for (RoundLocalServiceForAllPricingProblems localServiceBranching : masterData.localService4AllBranchingSet) {

                        IloRange constraint = masterData.localService4AllBranchingConstraints
                                .get(localServiceBranching);
                        int count = 0;
                        for (int edgeIndex : column.edgeIndexSet) {
                            Edge edge = dataModel.edgeSet.get(edgeIndex);
                            if (edge.serviceIndex == localServiceBranching.localServiceIndex)
                                count++;
                        }
                        iloColumn = iloColumn.and(masterData.cplex.column(constraint, count));
                    }

                }
            }
            
            //time service branching constraints(we use third kind of artificial variable)
            if(column.isArtificialColumn&&column.ifForResourceBoundConstraints==2){
                for(RoundTimeService timeServiceBranching:masterData.timeServiceBranchingSet){
                    if(timeServiceBranching.roundUpOrDown==1){
                        IloRange constraint=masterData.timeServiceBranchingConstraints.get(timeServiceBranching);
                        iloColumn=iloColumn.and(masterData.cplex.column(constraint,1));
                    }
                }
            }else{
                if(!column.isArtificialColumn){
                    
                    for(RoundTimeService timeServiceBranching:masterData.timeServiceBranchingSet){
                        if(column.associatedPricingProblem==timeServiceBranching.associatedPricingProblem){
                            IloRange constraint=masterData.timeServiceBranchingConstraints.get(timeServiceBranching);
                            
                            int count=0;
                            for(int edgeIndex:column.edgeIndexSet){
                                Edge edge=dataModel.edgeSet.get(edgeIndex);
                                
                                if(edge.serviceIndex==timeServiceBranching.localServiceIndex&&timeServiceBranching.timeSet.contains(edge.t1)){
                                    count++;
                                }
                            }
                            
                            iloColumn=iloColumn.and(masterData.cplex.column(constraint,count));
                        }
                    }
                }
            }
            
            
            //time service for all pricing problems branching constraints(we use third kind of artificial variable)
            if(column.isArtificialColumn&&column.ifForResourceBoundConstraints==2){
                for(RoundTimeServiceForAllPricingProblems timeService4AllBranching:masterData.timeService4AllBranchingSet){
                    if(timeService4AllBranching.roundUpOrDown==1){
                        IloRange constraint=masterData.timeService4AllBranchingConstraints.get(timeService4AllBranching);
                        iloColumn=iloColumn.and(masterData.cplex.column(constraint,1));
                    }
                }
            }else{
                if(!column.isArtificialColumn){
                    
                    for(RoundTimeServiceForAllPricingProblems timeService4AllBranching:masterData.timeService4AllBranchingSet){
                            IloRange constraint=masterData.timeService4AllBranchingConstraints.get(timeService4AllBranching);
                            
                            int count=0;
                            for(int edgeIndex:column.edgeIndexSet){
                                Edge edge=dataModel.edgeSet.get(edgeIndex);
                                
                                if(edge.serviceIndex==timeService4AllBranching.localServiceIndex&&timeService4AllBranching.timeSet.contains(edge.t1)){
                                    count++;
                                }
                            }
                            
                            iloColumn=iloColumn.and(masterData.cplex.column(constraint,count));
                        
                    }
                }
            }
            
            

            // for strong cuts , including artificial and non-artificial
            // variables(we use
            // first kind of artificial variables to insert)
            for (StrongInequality inequality : masterData.strongInequalities.keySet()) {
                if (column.edgeIndexSet.contains(inequality.edgeIndex)) {
                    iloColumn = iloColumn.and(masterData.cplex.column(masterData.strongInequalities.get(inequality),
                            -dataModel.demandSet.get(inequality.commodity).volume));
                }
            }

            // Create the variable and store it
            if (column.isArtificialColumn) {
                IloNumVar var = masterData.cplex.numVar(iloColumn, 0, Double.MAX_VALUE,
                        "art_" + column.associatedPricingProblem.capacityTypeS + ","
                                + column.associatedPricingProblem.originNodeO + ","
                                + masterData.getNrColumnsForPricingProblem(column.associatedPricingProblem));
                masterData.cplex.add(var);
                masterData.addColumn(column, var);
            } else {
                IloNumVar var = masterData.cplex.numVar(iloColumn, 0, Double.MAX_VALUE,
                        "z_" + column.associatedPricingProblem.capacityTypeS + ","
                                + column.associatedPricingProblem.originNodeO + ","
                                + masterData.getNrColumnsForPricingProblem(column.associatedPricingProblem));
                masterData.cplex.add(var);
                masterData.addColumn(column, var);
            }

        } catch (IloException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    /**
     * Extracts information from the master problem which is required by the
     * pricing problems, e.g. the reduced costs/dual values
     * Attention: We store the charge related cost on holding arcs of modifiedCosts[]!!!
     * 
     * @param pricingProblem
     *            pricing problem
     */
    @Override
    public void initializePricingProblem(SNDRCPricingProblem pricingProblem) {
        try {
            // double[] modifiedCosts = new double[dataModel.numServiceArc];
            double[] modifiedCosts = new double[dataModel.numArc];

            double modifiedCost = 0;

            // initialize for holding arcs
            for (int edgeIndex = dataModel.numServiceArc; edgeIndex < dataModel.numArc; edgeIndex++) {
                modifiedCosts[edgeIndex] = dataModel.chargeObjPara;
            }
            //charge constraints
            for(int edgeIndex=dataModel.numServiceArc;edgeIndex<dataModel.numArc;edgeIndex++) {
            	Edge edge=dataModel.edgeSet.get(edgeIndex);
            	modifiedCosts[edgeIndex]-=masterData.cplex.getDual(chargeBoundConstraints[edge.u][edge.t1]);
            }

            // week force constraints and vij in objective
            for (int edgeIndex = 0; edgeIndex < dataModel.numServiceArc; edgeIndex++) {
                modifiedCosts[edgeIndex] = masterData.cplex.getDual(weakForcingConstraints[edgeIndex])
                        * dataModel.capacity[pricingProblem.capacityTypeS];
                modifiedCosts[edgeIndex] += dataModel.alpha / (dataModel.speed * dataModel.drivingTimePerDay)
                        * dataModel.edgeSet.get(edgeIndex).duration;
            }

            // service edge branching constraints
            for (RoundServiceEdge serviceEdgeBranch : masterData.serviceEdgeBrachingSet) {
                if (serviceEdgeBranch.pricingProblem == pricingProblem) {
                    IloRange constraint = masterData.ServiceEdgeBranchingConstraints.get(serviceEdgeBranch);
                    double dualValue = masterData.cplex.getDual(constraint);
                    modifiedCosts[serviceEdgeBranch.branchEdgeIndex] -= dualValue;
                }
            }

            // service edge 4 all pricing problems branching constraints
            for (RoundServiceEdgeForAllPricingProblems serviceEdgeBranch : masterData.serviceEdge4AllBrachingSet) {
                IloRange constraint = masterData.serviceEdge4AllBranchingConstraints.get(serviceEdgeBranch);
                double dualValue = masterData.cplex.getDual(constraint);
                modifiedCosts[serviceEdgeBranch.branchEdgeIndex] -= dualValue;
            }

            // local service branching constraints
            for (RoundLocalService localServiceBranch : masterData.localServiceBranchingSet) {

                if (localServiceBranch.associatedPricingProblem == pricingProblem) {
                    IloRange constraint = masterData.localServiceBranchingConstraints.get(localServiceBranch);
                    double dualValue = masterData.cplex.getDual(constraint);

                    for (int serviceEdgeIndex = 0; serviceEdgeIndex < dataModel.numServiceArc; serviceEdgeIndex++) {
                        Edge edge = dataModel.edgeSet.get(serviceEdgeIndex);
                        if (edge.serviceIndex == localServiceBranch.localServiceIndex) {
                            modifiedCosts[serviceEdgeIndex] -= dualValue;
                        }
                    }
                }

            }

            /**
             * Here the index of holding arcs are calculated as :
             * localNode*T+t+numServiceArc
             * 
             * note that this calculation can't be applied to learningUB algorithm
             */
            // holding edges branching constraints
            for (RoundHoldingEdge branch : masterData.holdingEdgeBranchingSet) {

                if (branch.associatedPricingProblem == pricingProblem) {
                    IloRange constraint = masterData.holdingEdgeBranchingConstraints.get(branch);
                    double dualValue = masterData.cplex.getDual(constraint);

                    int time = branch.branchTime;
                    for (int localNode = 0; localNode < dataModel.numNode; localNode++) {
                        int edgeIndex = localNode * dataModel.timePeriod + time + dataModel.numServiceArc;
                        modifiedCosts[edgeIndex] -= dualValue;
                    }
                }

            }

            // local service 4 all branching constraints
            for (RoundLocalServiceForAllPricingProblems localServiceBranch : masterData.localService4AllBranchingSet) {
                IloRange constraint = masterData.localService4AllBranchingConstraints.get(localServiceBranch);
                double dualValue = masterData.cplex.getDual(constraint);

                for (int serviceEdgeIndex = 0; serviceEdgeIndex < dataModel.numServiceArc; serviceEdgeIndex++) {
                    Edge edge = dataModel.edgeSet.get(serviceEdgeIndex);
                    if (edge.serviceIndex == localServiceBranch.localServiceIndex) {
                        modifiedCosts[serviceEdgeIndex] -= dualValue;
                    }
                }

            }
            
            
            // time service branching constraints
            for(RoundTimeService timeServiceBranch:masterData.timeServiceBranchingSet) {
            	
            	if(timeServiceBranch.associatedPricingProblem==pricingProblem) {
            		IloRange constraint=masterData.timeServiceBranchingConstraints.get(timeServiceBranch);
            		double dualValue=masterData.cplex.getDual(constraint);
            		
            		for(int edgeIndex=0;edgeIndex<dataModel.numServiceArc;edgeIndex++) {
            			Edge edge=dataModel.edgeSet.get(edgeIndex);
            			if(edge.serviceIndex==timeServiceBranch.localServiceIndex&&timeServiceBranch.timeSet.contains(edge.t1)) {
            				modifiedCosts[edgeIndex]-=dualValue;
            			}
            		}
            	}
            }
            
            // time service 4 all branching constraints
            for(RoundTimeServiceForAllPricingProblems timeService4AllBranch:masterData.timeService4AllBranchingSet) {
            	
            		IloRange constraint=masterData.timeService4AllBranchingConstraints.get(timeService4AllBranch);
            		double dualValue=masterData.cplex.getDual(constraint);
            		
            		for(int edgeIndex=0;edgeIndex<dataModel.numServiceArc;edgeIndex++) {
            			Edge edge=dataModel.edgeSet.get(edgeIndex);
            			if(edge.serviceIndex==timeService4AllBranch.localServiceIndex&&timeService4AllBranch.timeSet.contains(edge.t1)) {
            				modifiedCosts[edgeIndex]-=dualValue;
            			}
            		}
            	
            }
            
            
            

            // strong cuts
            for (StrongInequality inequality : masterData.strongInequalities.keySet()) {
                IloRange constraint = masterData.strongInequalities.get(inequality);
                double dualValue = masterData.cplex.getDual(constraint);
                modifiedCosts[inequality.edgeIndex] += dataModel.demandSet.get(inequality.commodity).volume * dualValue;
            }

            modifiedCost = dataModel.fixedCost[pricingProblem.originNodeO][pricingProblem.capacityTypeS]
                    - masterData.cplex.getDual(
                            resourceBoundConstraints[pricingProblem.capacityTypeS][pricingProblem.originNodeO]);
            pricingProblem.initPricingProblem(modifiedCosts, modifiedCost);

        } catch (IloException e) {
            e.printStackTrace();
        }
    }

    /**
     * Gets the solution from the master problem
     * 
     * @return Returns all non-zero valued columns from the master problem
     */
    @Override
    public List<Cycle> getSolution() {
        List<Cycle> solution = new ArrayList<>();
        try {
            for (SNDRCPricingProblem pricingProblem : pricingProblems) {
                Cycle[] cycles = masterData.getVarMapForPricingProblem(pricingProblem)
                        .getKeysAsArray(new Cycle[masterData.getNrColumnsForPricingProblem(pricingProblem)]);
                IloNumVar[] vars = masterData.getVarMapForPricingProblem(pricingProblem)
                        .getValuesAsArray(new IloNumVar[masterData.getNrColumnsForPricingProblem(pricingProblem)]);
                

                
                double[] values = masterData.cplex.getValues(vars);
                
  

                // Iterate over each column and add it to the solution if it has
                // a non-zero
                // value
                for (int i = 0; i < cycles.length; i++) {
                    cycles[i].value = values[i];
                    if (values[i] >= config.PRECISION) {
                        solution.add(cycles[i]);
                    }
                }
            }
        } catch (IloException e) {
            e.printStackTrace();
        }
        return solution;
    }

    /**
     * Prints the solution
     */
    @Override
    public void printSolution() {
        List<Cycle> solution = this.getSolution();
        for (Cycle m : solution)
            System.out.println(m + ":" + m.value);
    }

    /**
     * Closes the master problem
     */
    @Override
    public void close() {
        masterData.cplex.end();
    }

    /**
     * Checks whether there are any violated inequalities, thereby invoking the
     * cut handler
     * 
     * @return true if violated inqualities have been found (and added to the
     *         master problem)
     */

    @Override
    public boolean hasNewCuts() {

        
        
        if (masterData.boundRecordForCutControl > 0) {
//            System.out.println("boundRecord=" + masterData.boundRecordForCutControl);
//            System.out.println("obj=" + masterData.objectiveValue);
//            System.out.println();

//            if (masterData.ifBoundHasChangedByCut) {
//                if ((masterData.objectiveValue - masterData.boundRecordForCutControl)
//                        / masterData.boundRecordForCutControl < 0.001) {
//                    return false;
//                }
//            }
//
//            if ((masterData.objectiveValue - masterData.boundRecordForCutControl)
//                    / masterData.boundRecordForCutControl > 0.00001) {
//                masterData.ifBoundHasChangedByCut = true;
//            }
            
            

        }

        long time = System.currentTimeMillis();
        if (time - masterData.masterDataBuildTime > 1800000) { // half an hour
            return false;
        }

        masterData.boundRecordForCutControl = masterData.objectiveValue;
        // For convenience, we will precompute values required by the
        // StrongInequalityGenerator class
        // and store it in the masterData object. For each edge, we need to know
        // how
        // often it is used.
        masterData.edgeValueMap.clear();

        for (Cycle cycle : this.getSolution()) {
            for (int edgeIndex : cycle.edgeIndexSet) {
                if (dataModel.edgeSet.get(edgeIndex).edgeType == 0) {
                    if (masterData.edgeValueMap.containsKey(edgeIndex)) {
                        masterData.edgeValueMap.put(edgeIndex, masterData.edgeValueMap.get(edgeIndex) + cycle.value);
                    } else {
                        masterData.edgeValueMap.put(edgeIndex, cycle.value);
                    }
                }
            }
        }

        return super.hasNewCuts();
    }

    /**
     * Listen to branching decisions
     * 
     * @param bd
     *            Branching decision
     */
    @Override
    public void branchingDecisionPerformed(BranchingDecision bd) {

        if (bd instanceof RoundQ) {
            qBranchingSet.add((RoundQ) bd);
        }

        if (bd instanceof RoundServiceEdge) {
            serviceEdgeBrachingSet.add((RoundServiceEdge) bd);
        }

        if (bd instanceof RoundServiceEdgeForAllPricingProblems) {
            serviceEdge4AllBrachingSet.add((RoundServiceEdgeForAllPricingProblems) bd);
        }

        if (bd instanceof RoundLocalService) {
            localServiceBranchingSet.add((RoundLocalService) bd);
        }

        if (bd instanceof RoundHoldingEdge) {
            holdingEdgeBranchingSet.add((RoundHoldingEdge) bd);
        }

        if (bd instanceof RoundLocalServiceForAllPricingProblems) {
            localService4AllBranchingSet.add((RoundLocalServiceForAllPricingProblems) bd);
        }
        
        if (bd instanceof RoundTimeService) {
            timeServiceBranchingSet.add( (RoundTimeService) bd);
        }
        
        if (bd instanceof RoundTimeServiceForAllPricingProblems) {
            timeService4AllBranchingSet.add( (RoundTimeServiceForAllPricingProblems) bd);
        }

        if (branchSetNodeSolvedByCut.contains(bd)) {
            cutHandler.addCutGenerator(cutGen);
            ifContainsCut = true;
        } else {
            if (ifContainsCut) {
                cutHandler.removeCutGenerator(cutGen);
                ifContainsCut = false;
            }
        }

        // destroy the master and rebuild it
        this.close();
        masterData = this.buildModel();
        cutHandler.setMasterData(masterData);

    }

    /**
     * Undo branching decisions during backtracking in the Branch-and-Price tree
     * 
     * @param bd
     *            Branching decision
     */
    @Override
    public void branchingDecisionReversed(BranchingDecision bd) {

        if (bd instanceof RoundQ) {
            // masterData.cplex.remove(masterData.qBranchingconstraints.get(bd));
            // masterData.qBranching.remove(bd);
            try {
                masterData.cplex.remove(masterData.qBranchingconstraints.get(bd));
                masterData.qBranchingSet.remove(bd);
                masterData.qBranchingconstraints.remove(bd);
                this.qBranchingconstraints.remove(bd);
                this.qBranchingSet.remove(bd);
            } catch (IloException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        if (bd instanceof RoundServiceEdge) {
            try {
                masterData.cplex.remove(masterData.ServiceEdgeBranchingConstraints.get(bd));
                masterData.serviceEdgeBrachingSet.remove(bd);
                masterData.ServiceEdgeBranchingConstraints.remove(bd);
                this.ServiceEdgeBranchingConstraints.remove(bd);
                this.serviceEdgeBrachingSet.remove(bd);

            } catch (IloException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        if (bd instanceof RoundServiceEdgeForAllPricingProblems) {
            try {
                masterData.cplex.remove(masterData.serviceEdge4AllBranchingConstraints.get(bd));
                masterData.serviceEdge4AllBrachingSet.remove(bd);
                masterData.serviceEdge4AllBranchingConstraints.remove(bd);
                this.serviceEdge4AllBranchingConstraints.remove(bd);
                this.serviceEdge4AllBrachingSet.remove(bd);

            } catch (IloException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        if (bd instanceof RoundLocalService) {
            try {
                masterData.cplex.remove(masterData.localServiceBranchingConstraints.get(bd));
                masterData.localServiceBranchingSet.remove(bd);
                masterData.localServiceBranchingConstraints.remove(bd);
                this.localServiceBranchingConstraints.remove(bd);
                this.localServiceBranchingSet.remove(bd);

            } catch (IloException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        if (bd instanceof RoundHoldingEdge) {
            try {
                masterData.cplex.remove(masterData.holdingEdgeBranchingConstraints.get(bd));
                masterData.holdingEdgeBranchingSet.remove(bd);
                masterData.holdingEdgeBranchingConstraints.remove(bd);
                this.holdingEdgeBranchingConstraints.remove(bd);
                this.holdingEdgeBranchingSet.remove(bd);

            } catch (IloException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        if (bd instanceof RoundLocalServiceForAllPricingProblems) {
            try {
                masterData.cplex.remove(masterData.localService4AllBranchingConstraints.get(bd));
                masterData.localService4AllBranchingSet.remove(bd);
                masterData.localService4AllBranchingConstraints.remove(bd);
                this.localService4AllBranchingConstraints.remove(bd);
                this.localService4AllBranchingSet.remove(bd);

            } catch (IloException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        
        if (bd instanceof RoundTimeService) {
            try {
                masterData.cplex.remove(masterData.timeServiceBranchingConstraints.get(bd));
                masterData.timeServiceBranchingSet.remove(bd);
                masterData.timeServiceBranchingConstraints.remove(bd);
                this.timeServiceBranchingConstraints.remove(bd);
                this.timeServiceBranchingSet.remove(bd);

            } catch (IloException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        
        if (bd instanceof RoundTimeServiceForAllPricingProblems) {
            try {
                masterData.cplex.remove(masterData.timeService4AllBranchingConstraints.get(bd));
                masterData.timeService4AllBranchingSet.remove(bd);
                masterData.timeService4AllBranchingConstraints.remove(bd);
                this.timeService4AllBranchingConstraints.remove(bd);
                this.timeService4AllBranchingSet.remove(bd);

            } catch (IloException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        if (branchSetNodeSolvedByCut.contains(bd)) {
            if (ifContainsCut) {
                cutHandler.removeCutGenerator(cutGen);
                ifContainsCut = false;
            }

            branchSetNodeSolvedByCut.remove(bd);
        }

    }

    public void addFixVarConstraint(Cycle cycle) throws IloException {

        IloLinearNumExpr expr = masterData.cplex.linearNumExpr();
        expr.addTerm(masterData.getVar(cycle.associatedPricingProblem, cycle), 1);
        IloRange fixVarConstraint = masterData.cplex.addEq(expr, Math.ceil(cycle.value));

        masterData.fixVarConstraints.put(cycle, fixVarConstraint);
        SNDRCPricingProblem pricingProblem = pricingProblems
                .get(cycle.associatedPricingProblem.capacityTypeS * dataModel.numNode
                        + cycle.associatedPricingProblem.originNodeO);
        pricingProblem.fixCycleSet.add(cycle);

        // System.out.println(masterData.cplex.toString());
        // return fixVarConstraint;

    }

    public void removeFixVarConstraint() throws IloException {
        // for(IloRange constraint:constraints) {
        // masterData.cplex.remove(constraint);
        // }

        for (Cycle cycle : masterData.fixVarConstraints.keySet()) {
            masterData.cplex.remove(masterData.fixVarConstraints.get(cycle));
        }

        for (SNDRCPricingProblem problem : pricingProblems) {
            problem.fixCycleSet = new HashSet<>();
        }

        masterData.fixVarConstraints = new HashMap<>();
    }

    public Map<SNDRCPricingProblem, Integer> getQVarLimit() {
        return masterData.qVariableLimit;
    }

    public boolean CheckFeasibility() throws IloException {
        return masterData.cplex.solve();
    }

    // output the current masterData.cplex.tostring file
    public void Output(int nodeIndex) throws IloException {
        masterData.cplex.exportModel(config.EXPORT_MASTER_DIR + nodeIndex + ".lp");
    }

    public List<Map<Integer, Double>> getXValues() throws UnknownObjectException, IloException {

        return masterData.xValues;
    }

    public void AddBranchDecisionForCut(BranchingDecision<SNDRC, Cycle> bd) {
        branchSetNodeSolvedByCut.add(bd);
    }

    public StrongInequalityGenerator getStrongInequalityGenerator() {
        return cutGen;
    }
    
	/**
	 * To compute a bound on the optimal solution of the relaxed master problem, multiple components
	 * are required, including information from the master problem. This function returns that information.
	 * @return value originating from the master problem which is required to calculate a bound on the optimal objective of the master problem
	 */
    @Override
	public double getBoundComponent(){
		return masterData.objectiveValue;
	}
    
    public double[][] getResourceBoundConstraintsDual(){
        return masterData.resourceBoundConstraintsDual;
    }

}
