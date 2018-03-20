package bap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import org.jorlib.demo.frameworks.columnGeneration.cuttingStockCG.cg.PricingProblem;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.CGListener;
import org.jorlib.frameworks.columnGeneration.colgenMain.ColGen;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import com.sun.media.jfxmedia.events.NewFrameEvent;

import bap.bapNodeComparators.NodeBoundbapNodeComparator;
import bap.bapNodeComparators.NodeBoundbapNodeComparatorForLB;
import bap.branching.branchingDecisions.RoundServiceEdge;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import ilog.concert.IloException;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex.UnknownObjectException;
import model.SNDRC;
import model.SNDRC.Edge;

public class BranchAndPrice<V> extends AbstractBranchAndPrice<SNDRC, Cycle, SNDRCPricingProblem> {

    private double thresholdValue;
    private PriorityQueue<BAPNode<SNDRC, Cycle>> lowBoundQueue;

    private final double probLB;
    private final double c;
    private int nrNonImproForAcce;
    private Map<Cycle, Double> optSolutionValueMap;
    private List<Map<Integer, Double>> optXValues;
    private double[] nodeBoundRecord;
    private int helpOutPut;

    public BranchAndPrice(SNDRC modelData, Master master, List<SNDRCPricingProblem> pricingProblems,
            List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
            List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators,
            double objectiveInitialSolution, double thresholdValue, double probLB, double c) {
        super(modelData, master, pricingProblems, solvers, branchCreators, 0, objectiveInitialSolution);
        // this.warmStart(objectiveInitialSolution,initialSolution);
        this.thresholdValue = thresholdValue;
        lowBoundQueue = new PriorityQueue<>(new NodeBoundbapNodeComparatorForLB());
        this.probLB = probLB;
        this.c = c;
        nrNonImproForAcce = 0;
        optSolutionValueMap = new HashMap<>();
        nodeBoundRecord = new double[5000];
        helpOutPut=0;
    }

    /**
     * Generates an artificial solution. Columns in the artificial solution are
     * of high cost such that they never end up in the final solution if a
     * feasible solution exists, since any feasible solution is assumed to be
     * cheaper than the artificial solution. The artificial solution is used to
     * guarantee that the master problem has a feasible solution.
     *
     * @return artificial solution
     */
    @Override
    protected List<Cycle> generateInitialFeasibleSolution(BAPNode<SNDRC, Cycle> node) {

        List<Cycle> artificalVars = new ArrayList<Cycle>();
        // for weak forcing constraints(ifForResourceBoundConstraints=0)
        for (int edgeIndex = 0; edgeIndex < dataModel.numServiceArc; edgeIndex++) {
            Set<Integer> set = new HashSet<>();
            set.add(edgeIndex);
            Cycle cycle = new Cycle(pricingProblems.get(0), true, "Artificial", set, 100000000, 0, 0);
            artificalVars.add(cycle);
        }

        // for resource bound constraints(ifForResourceBoundConstraints=1)
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            Set<Integer> set = new HashSet<>();
            Cycle cycle = new Cycle(pricingProblem, true, "Artificial", set, 100000000, 0, 1);
            artificalVars.add(cycle);
        }

        // for holding edge branch constraints(ifForResourceBoundConstraints=2)
        Set<Integer> set = new HashSet<>();
        Cycle cycle = new Cycle(pricingProblems.get(0), true, "Artificial", set, 100000000, 0, 2);
        artificalVars.add(cycle);

        return artificalVars;
    }

    /**
     * Checks whether the given node is integer
     * 
     * @param node
     *            Node in the Branch-and-Price tree
     * @return true if the solution is an integer solution
     */
    @Override
    protected boolean isIntegerNode(BAPNode<SNDRC, Cycle> node) {
        List<Cycle> result = node.getSolution();

        boolean out = true;
        for (Cycle cycle : result) {
            if (MathProgrammingUtil.isFractional(cycle.value)) {
                out = false;
                break;
            }
        }

        return out;
    }

    /**
     * Starts running the Branch-and-Price algorithm. Note: In the current
     * version of the code, one should not invoke this function multiple times
     * on the same instance!
     * 
     * @param timeLimit
     *            Future point in time by which the algorithm should finish
     */
    @Override
    public void runBranchAndPrice(long timeLimit) {
        notifier.fireStartBAPEvent(); // Signal start Branch-and-Price process
        this.runtime = System.currentTimeMillis();

        // Check whether an warm start is provided, if not, invoke
        // generateInitialFeasibleSolution
        BAPNode<SNDRC, Cycle> rootNode = queue.peek();
        if (rootNode.getInitialColumns().isEmpty())
            rootNode.addInitialColumns(this.generateInitialFeasibleSolution(rootNode));

        lowBoundQueue.add(rootNode);

        // Start processing nodes until the queue is empty
        while (!queue.isEmpty()) {
            BAPNode<SNDRC, Cycle> bapNode = queue.poll();
            // lowBoundQueue.poll();

            notifier.fireNextNodeEvent(bapNode);
            // lowBoundQueue.poll();
            lowBoundQueue.remove(bapNode);

            // Prune this node if its bound is worse than the best found
            // solution. Since all solutions are integral, we may round up/down,
            // depending on the optimization sense
            if (this.nodeCanBePruned(bapNode)) {
                notifier.firePruneNodeEvent(bapNode, bapNode.getBound());
                nodesProcessed++;
                continue;
            }

            graphManipulator.next(bapNode); // Prepare data structures for the
                                            // next node

            // Generate an initial solution for this node to guarantee that the
            // master problem is feasible
            if (bapNode.nodeID != 0) {
                bapNode.addInitialColumns(this.generateInitialFeasibleSolution(bapNode));
            }

            // Solve the next BAPNode
            try {
                this.solveBAPNode(bapNode, timeLimit);

                // output the model
//                ((Master) master).Output(bapNode.nodeID);
                nodeBoundRecord[bapNode.nodeID] = bapNode.getBound();

//                if (bapNode.nodeID == 0) {
//                    bapNodeSolutionOutput(bapNode);
//                }
                

            } catch (TimeLimitExceededException e) {
                queue.add(bapNode);
                lowBoundQueue.add(bapNode);
                notifier.fireTimeOutEvent(bapNode);
                break;
            }

            // Prune this node if its bound is worse than the best found
            // solution. Since all solutions are integral, we may round up/down,
            // depending on the optimization sense
            if (this.nodeCanBePruned(bapNode)) {
                notifier.firePruneNodeEvent(bapNode, bapNode.getBound());
                nodesProcessed++;
                continue;
            }

            // Check whether the node is infeasible, i.e. whether there are
            // artifical columns in the solution. If so, ignore it and continue
            // with the next node.
            if (this.isInfeasibleNode(bapNode)) {
                notifier.fireNodeIsInfeasibleEvent(bapNode);
                nodesProcessed++;
                continue;
            }

            // If solution is integral, check whether it is better than the
            // current best solution
            if (this.isIntegerNode(bapNode)) {
                int integerObjective = MathProgrammingUtil.doubleToInt(bapNode.getObjective());
                notifier.fireNodeIsIntegerEvent(bapNode, bapNode.getBound(), integerObjective);
                if (optimizationSenseMaster == OptimizationSense.MINIMIZE
                        && integerObjective < this.upperBoundOnObjective) {
                    this.objectiveIncumbentSolution = integerObjective;
                    this.upperBoundOnObjective = integerObjective;
                    this.incumbentSolution = bapNode.getSolution();

                    optSolutionValueMap = new HashMap<>();
                    for (Cycle cycle : incumbentSolution) {
                        optSolutionValueMap.put(cycle, cycle.value);
                    }
                    try {
                        optXValues = ((Master) master).getXValues();
//                        bapNodeSolutionOutput(bapNode);
                    } catch (IloException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                } else if (optimizationSenseMaster == OptimizationSense.MAXIMIZE
                        && integerObjective > this.lowerBoundOnObjective) {
                    this.objectiveIncumbentSolution = integerObjective;
                    this.lowerBoundOnObjective = integerObjective;
                    this.incumbentSolution = bapNode.getSolution();

                    optSolutionValueMap = new HashMap<>();
                    for (Cycle cycle : incumbentSolution) {
                        optSolutionValueMap.put(cycle, cycle.value);
                    }
                    try {
                        optXValues = ((Master) master).getXValues();
                    } catch (IloException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                }
            } else { // We need to branch

                // ----------------------------------------------------------------------------------------------------------------------------------//

                // An acceleration technique for ub
                try {
                    double prob = CalculateProb();
                    double random = Math.random();
                    if (random < prob) {
                        this.AccelerationForUB(bapNode);
                    }

                } catch (IloException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }

//                if(Math.abs(bapNode.getBound()-3636)<0.0001&&helpOutPut==0){
//                    System.out.println("node ID="+bapNode.nodeID+" bound="+bapNode.getBound());
//                    try {
//                        bapNodeSolutionOutput(bapNode);
//                        helpOutPut=1;
//                    } catch (IloException e) {
//                        // TODO Auto-generated catch block
//                        e.printStackTrace();
//                    }
//                    
//                }
//                
//                if(Math.abs(bapNode.getBound()-3661)<0.0001&&helpOutPut==1){
//                    System.out.println("node ID="+bapNode.nodeID+" bound="+bapNode.getBound());
//                    try {
//                        bapNodeSolutionOutput(bapNode);
//                    } catch (IloException e) {
//                        // TODO Auto-generated catch block
//                        e.printStackTrace();
//                    }
//                    helpOutPut=2;
//                }
                

//                if (bapNode.getParentID() >= 0) {
//                    if (Math.abs(bapNode.getBound() - nodeBoundRecord[bapNode.getParentID()]) < 0.000001
//                            && Math.abs(bapNode.getBound() - 2326) < 0.000001){
////                        try {
////                            System.out.println(bapNode.getBranchingDecision().toString());
////                            bapNodeSolutionOutput(bapNode);
////                        } catch (IloException e) {
////                            // TODO Auto-generated catch block
////                            e.printStackTrace();
////                        }
//                        throw new RuntimeException(
//                                "LB doesn't improve!!! " + "node:" + bapNode.nodeID + " " + bapNode.getParentID()
//                                );
//                        
//                    }
//                }

                notifier.fireNodeIsFractionalEvent(bapNode, bapNode.getBound(), bapNode.getObjective());
                List<BAPNode<SNDRC, Cycle>> newBranches = new ArrayList<>();
                for (AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem> bc : branchCreators) {
                    newBranches.addAll(bc.branch(bapNode));
                    if (!newBranches.isEmpty())
                        break;
                }

                if (newBranches.isEmpty())
                    throw new RuntimeException(
                            "BAP encountered fractional solution, but non of the BranchCreators produced any new branches?");
                else {
                    queue.addAll(newBranches);
                    lowBoundQueue.addAll(newBranches);
                    notifier.fireBranchEvent(bapNode, Collections.unmodifiableList(newBranches));
                }
            }

            nodesProcessed++;
        }

        // Update statistics
        if (queue.isEmpty()) { // Problem solved to optimality
            this.isOptimal = true;
            if (optimizationSenseMaster == OptimizationSense.MINIMIZE)
                this.lowerBoundOnObjective = this.objectiveIncumbentSolution;
            else
                this.upperBoundOnObjective = this.objectiveIncumbentSolution;
        } else { // Problem NOT solved to optimality
            this.isOptimal = false;
            if (optimizationSenseMaster == OptimizationSense.MINIMIZE) {
                lowerBoundOnObjective = queue.peek().getBound();
                for (BAPNode bapNode : queue) {
                    lowerBoundOnObjective = Math.min(lowerBoundOnObjective, bapNode.getBound());
                }
            } else {
                upperBoundOnObjective = queue.peek().getBound();
                for (BAPNode bapNode : queue) {
                    upperBoundOnObjective = Math.max(upperBoundOnObjective, bapNode.getBound());
                }
            }
        }
        notifier.fireStopBAPEvent(); // Signal that BAP has been completed
        this.runtime = System.currentTimeMillis() - runtime;
    }

    public Double CalculateProb() {
        double prob = 2 - probLB;
        prob -= (2 - 2 * probLB) / (1 + Math.pow(Math.E, -c * nrNonImproForAcce));
        return prob;
    }

    public void AccelerationForUB(BAPNode<SNDRC, Cycle> bapNode) throws IloException {
        boolean ifFindBetterUB = false;

        List<Cycle> solution = bapNode.getSolution();

        // Set<IloRange> fixVarConstraints=new HashSet<>();

        double objRecord = bapNode.getObjective();
        double boundRecord = bapNode.getBound();
        List<Cycle> solutionRecord = new ArrayList<Cycle>();
        Map<Cycle, Double> mapValue = new HashMap<>();
        for (Cycle cycle : bapNode.getSolution()) {
            solutionRecord.add(cycle);
            mapValue.put(cycle, cycle.value);
        }
        List<AbstractInequality> inequalityRecord = new ArrayList<>();
        for (AbstractInequality inequality : bapNode.getInequalities()) {
            inequalityRecord.add(inequality);
        }

        Map<SNDRCPricingProblem, Integer> cycleLimit = ((Master) master).getQVarLimit();
        Map<SNDRCPricingProblem, Integer> fixUsed = new HashMap<>();
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            fixUsed.put(pricingProblem, 0);
        }

        while (!this.isIntegerNode(bapNode)) {

            boolean ifAllBelowThresholdValue = true;
            boolean ifFindOneToFix = false;
            solution = bapNode.getSolution();

            for (Cycle cycle : solution) {

                if (MathProgrammingUtil.isFractional(cycle.value)) {
                    double decimalValue = cycle.value - (int) cycle.value;
                    if (decimalValue > thresholdValue && cycleLimit.get(cycle.associatedPricingProblem) > fixUsed
                            .get(cycle.associatedPricingProblem)) {
                        fixUsed.put(cycle.associatedPricingProblem, fixUsed.get(cycle.associatedPricingProblem) + 1);
                        ifFindOneToFix = true;
                        ifAllBelowThresholdValue = false;
                        ((Master) master).addFixVarConstraint(cycle);
                    }
                }
            }

            // if all cycles' value are below the threshold value, fix the
            // variable with highest fractional decimal value
            double record = 0;
            Cycle cycleRecord = null;
            if (ifAllBelowThresholdValue) {
                for (Cycle cycle : solution) {
                    if (MathProgrammingUtil.isFractional(cycle.value) && cycleLimit
                            .get(cycle.associatedPricingProblem) > fixUsed.get(cycle.associatedPricingProblem)) {
                        double decimalValue = cycle.value - (int) cycle.value;
                        if (decimalValue > record) {
                            cycleRecord = cycle;
                            record = decimalValue;
                        }
                    }
                }

                if (cycleRecord != null) {
                    fixUsed.put(cycleRecord.associatedPricingProblem,
                            fixUsed.get(cycleRecord.associatedPricingProblem) + 1);
                    ifFindOneToFix = true;
                    ((Master) master).addFixVarConstraint(cycleRecord);
                } else {
                    break;
                }

            }

            // here we should check if the master problem is feasible
            if (((Master) master).CheckFeasibility() == false) {
                break;
            }

            // start to solve master problem
            // bapNode.addInitialColumns(this.generateInitialFeasibleSolution(bapNode));

            ColGen<SNDRC, Cycle, SNDRCPricingProblem> cg = null;
            try {
                List<Cycle> nullList = new ArrayList<>();
                // cg = new ColGen<>(dataModel, master, pricingProblems,
                // solvers, pricingProblemManager, bapNode.getInitialColumns(),
                // objectiveIncumbentSolution, bapNode.getBound()); //Solve the
                // node
                cg = new ColGen<>(dataModel, master, pricingProblems, solvers, pricingProblemManager, nullList,
                        objectiveIncumbentSolution, bapNode.getBound());
                // for(CGListener listener : columnGenerationEventListeners)
                // cg.addCGEventListener(listener);
                try {
                    cg.solve(System.currentTimeMillis() + 360000000);
                } catch (TimeLimitExceededException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            } finally {
                // //Update statistics
                // if(cg != null) {
                // timeSolvingMaster += cg.getMasterSolveTime();
                // timeSolvingPricing += cg.getPricingSolveTime();
                // totalNrIterations += cg.getNumberOfIterations();
                // totalGeneratedColumns += cg.getNrGeneratedColumns();
                // notifier.fireFinishCGEvent(bapNode, cg.getBound(),
                // cg.getObjective(), cg.getNumberOfIterations(),
                // cg.getMasterSolveTime(), cg.getPricingSolveTime(),
                // cg.getNrGeneratedColumns());
                // }
            }
            bapNode.storeSolution(cg.getObjective(), cg.getBound(), cg.getSolution(), cg.getCuts());
            if (this.nodeCanBePruned(bapNode) || this.isInfeasibleNode(bapNode)) {
                break;
            }

            if (this.isIntegerNode(bapNode)) {

                int integerObjective = MathProgrammingUtil.doubleToInt(bapNode.getObjective());
                // notifier.fireNodeIsIntegerEvent(bapNode, bapNode.getBound(),
                // integerObjective);
                if (optimizationSenseMaster == OptimizationSense.MINIMIZE
                        && integerObjective < this.upperBoundOnObjective) {
                    this.objectiveIncumbentSolution = integerObjective;
                    this.upperBoundOnObjective = integerObjective;
                    this.incumbentSolution = new ArrayList<>();
                    for (Cycle cycle : bapNode.getSolution()) {
                        this.incumbentSolution.add(cycle);
                    }
                    optSolutionValueMap = new HashMap<>();
                    for (Cycle cycle : incumbentSolution) {
                        optSolutionValueMap.put(cycle, cycle.value);
                    }
                    try {
                        optXValues = ((Master) master).getXValues();
//                        bapNodeSolutionOutput(bapNode);
                    } catch (IloException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }

                    // deal with nrNonImproForAcce
                    nrNonImproForAcce = 0;
                    ifFindBetterUB = true;

                } else if (optimizationSenseMaster == OptimizationSense.MAXIMIZE
                        && integerObjective > this.lowerBoundOnObjective) {
                    this.objectiveIncumbentSolution = integerObjective;
                    this.lowerBoundOnObjective = integerObjective;
                    this.incumbentSolution = new ArrayList<>();
                    for (Cycle cycle : bapNode.getSolution()) {
                        this.incumbentSolution.add(cycle);
                    }
                    optSolutionValueMap = new HashMap<>();
                    for (Cycle cycle : incumbentSolution) {
                        optSolutionValueMap.put(cycle, cycle.value);
                    }
                    try {
                        optXValues = ((Master) master).getXValues();
                    } catch (IloException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }

                    nrNonImproForAcce = 0;
                    ifFindBetterUB = true;
                }

                break;

            }

        }

        ((Master) master).removeFixVarConstraint();

        // After the acceleration technique and removing the fix variable
        // constraints, we need to solve the master problem again
        for (Cycle cycle : mapValue.keySet()) {
            cycle.value = mapValue.get(cycle);
        }
        bapNode.storeSolution(objRecord, boundRecord, solutionRecord, inequalityRecord);

        if (!ifFindBetterUB)
            nrNonImproForAcce++;

    }

    public PriorityQueue<BAPNode<SNDRC, Cycle>> getLowBoundQueue() {
        return lowBoundQueue;
    }

    public Map<Cycle, Double> GetOptSolutionValueMap() {
        return optSolutionValueMap;
    }

    public List<Map<Integer, Double>> GetOptXValues() {
        return optXValues;
    }
    
    public void bapNodeSolutionOutput(BAPNode<SNDRC, Cycle> bapNode) throws UnknownObjectException, IloException{
        System.out.println("Now the node bound="+bapNode.getBound());
        List<Cycle> rootSolution = bapNode.getSolution();
        for (Cycle cycle : rootSolution) {
            System.out.println(cycle.toString() + ":" + cycle.value);
            for (int edgeIndex : cycle.edgeIndexSet) {
                Edge edge = dataModel.edgeSet.get(edgeIndex);
                if (edge.edgeType == 0) {
                    System.out.print("(" + edge.u + "," + edge.t1 + ")->(" + edge.v + "," + edge.t2 + ") ");
                }
            }
            System.out.println();
        }
        
        List<Map<Integer,Double>> optXValues=((Master) master).getXValues();
        //output x variables
        for(int demand=0;demand<dataModel.numDemand;demand++){
            for(int edgeIndex:optXValues.get(demand).keySet()){
                if(optXValues.get(demand).get(edgeIndex)>0.01){
                    Edge edge=dataModel.edgeSet.get(edgeIndex);
                    System.out.println("x[" + demand + "]:" + edge.start + "->" + edge.end +" "+edge.u+","+edge.t1+"->"+edge.v+","+edge.t2+" = "+optXValues.get(demand).get(edgeIndex));
                }
            }
            System.out.println();
        }
    }

}
