package bap;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.colgenMain.ColGen;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.Configuration;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import com.sun.xml.internal.ws.api.server.ServiceDefinition;

import bap.bapNodeComparators.NodeBoundbapNodeComparator;
import bap.bapNodeComparators.NodeBoundbapNodeComparatorForLB;
import bap.branching.BranchOnLocalService;
import bap.branching.BranchOnLocalServiceForAllPricingProblems;
import bap.branching.BranchOnServiceEdge;
import cg.Cycle;
import cg.ExactPricingProblemSolver;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import cg.master.SNDRCMasterData;
import cg.master.cuts.StrongInequality;
import cg.master.cuts.StrongInequalityGenerator;
import ilog.concert.IloException;
import ilog.cplex.IloCplex.UnknownObjectException;
import logger.BapLoggerB;
import logger.BapLoggerB_M;
import model.SNDRC;
import model.SNDRC.Edge;
import model.SNDRC.Service;

public class BranchAndPriceB_M <V> extends AbstractBranchAndPrice<SNDRC, Cycle, SNDRCPricingProblem>{

    private double thresholdValue;
    private PriorityQueue<BAPNode<SNDRC, Cycle>> lowBoundQueue;

    private final double probLB;
    private final double c;
    private int nrNonImproForAcce;
    private Map<Cycle, Double> optSolutionValueMap;
    private List<Map<Integer, Double>> optXValues;
//    private double[] nodeBoundRecord;
//    private int helpOutPut;
    
    //For learning upper bound
    private int[] cutFrequency;
    private double[] edgeFrequency,accumulatedReducedCost;
    private int nodeFre;
    private double alphaForEdgeFre,leanringCheckPercent;
    private int timeCompress;
    private boolean ifUseLearningUB,ifAccelerationForUB;
    private boolean[] subEdgeRecord;
    private boolean ifOptGetFromSubGraph;
    
    
    
    /**
     * 
     * @param modelData
     * @param master
     * @param pricingProblems
     * @param solvers
     * @param branchCreators
     * @param objectiveInitialSolution
     * @param thresholdValue   for AccelerationForUB()
     * @param probLB  for AccelerationForUB()
     * @param c  for AccelerationForUB()
     * @param nodeFre for learningUB()
     * @param alphaForEdgeFre  parameter for learningUB(), we use it to decide whether a service edge should be included to the subgraph
     * @param leanringCheckPercent  when the new edge subset chosen has this percent of edges different from current subset, we decide to learnUB() 
     */

    public BranchAndPriceB_M(SNDRC modelData, Master master, List<SNDRCPricingProblem> pricingProblems,
            List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
            List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators,
            double objectiveInitialSolution, double thresholdValue, double probLB, double c,int nodeFre,double alphaForEdgeFre,int timeCompress,double leanringCheckPercent,boolean ifUseLearningUB,boolean ifAccelerationForUB) {
        super(modelData, master, pricingProblems, solvers, branchCreators, 0, objectiveInitialSolution);
        // this.warmStart(objectiveInitialSolution,initialSolution);
        this.thresholdValue = thresholdValue;
        lowBoundQueue = new PriorityQueue<>(new NodeBoundbapNodeComparatorForLB());
        this.probLB = probLB;
        this.c = c;
        nrNonImproForAcce = 0;
        optSolutionValueMap = new HashMap<>();
//        nodeBoundRecord = new double[10000];
//        helpOutPut=0;
        
        edgeFrequency=new double[modelData.numServiceArc];
        cutFrequency=new int[modelData.numServiceArc];
        accumulatedReducedCost=new double[modelData.numServiceArc];
        this.nodeFre=nodeFre;
        this.alphaForEdgeFre=alphaForEdgeFre;
        this.timeCompress=timeCompress;
        this.ifUseLearningUB=ifUseLearningUB;
        
        this.subEdgeRecord=new boolean[dataModel.numServiceArc];
        for(int i=0;i<subEdgeRecord.length;i++){
            subEdgeRecord[i]=false;
        }
        this.leanringCheckPercent=leanringCheckPercent;
        
        this.ifAccelerationForUB=ifAccelerationForUB;
        
        this.ifOptGetFromSubGraph=false;
        
        
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

            
            Double parentBound=bapNode.getBound();
            
            // Solve the next BAPNode
            try {
                
                
                this.solveBAPNode(bapNode, timeLimit);

                // output the model
//                ((Master) master).Output(bapNode.nodeID);
//                nodeBoundRecord[bapNode.nodeID] = bapNode.getBound();

                

            } catch (TimeLimitExceededException e) {
                queue.add(bapNode);
                lowBoundQueue.add(bapNode);
                notifier.fireTimeOutEvent(bapNode);
                break;
            }
            
            
            
            if(this.ifUseLearningUB==true){
              //Collect the information for learning upper bound
                if((!this.nodeCanBePruned(bapNode))&&(!this.isInfeasibleNode(bapNode))){  //bapNode is integer or fractional
                    
                    //for edgeFrequency
                    List<Cycle> solution=bapNode.getSolution();
                    for(Cycle cycle:solution){
                        double value=cycle.value;
                        for(int edgeIndex:cycle.edgeIndexSet){
                            if(edgeIndex<dataModel.numServiceArc){
                                edgeFrequency[edgeIndex]+=value;
                            }
                        }
                    }

                    
                    
                    //for cutFrequence
                    List<AbstractInequality> initialCutSet=bapNode.getInitialInequalities();
                    List<AbstractInequality> cutSet=bapNode.getInequalities();
                    for(AbstractInequality cut:cutSet){
                        if(!initialCutSet.contains(cut)){
                            
                            if(cut instanceof StrongInequality) {
                                StrongInequality  strongInequality = (StrongInequality) cut;
                                cutFrequency[strongInequality.edgeIndex]++;
                            }
                            
                            
                        }
                    }
                    
                    
                    //for accumulatedReducedCost
                    
                    //pick the most used pricing problem
                    Map<SNDRCPricingProblem,Double> temp=new HashMap<>();
                    for(Cycle cycle:solution){
                        if(!temp.keySet().contains(cycle.associatedPricingProblem)){
                            temp.put(cycle.associatedPricingProblem, cycle.value);
                        }else{
                            temp.put(cycle.associatedPricingProblem, temp.get(cycle.associatedPricingProblem)+cycle.value);
                        }
                    }
                    
                    SNDRCPricingProblem mostUsedPricingProblem=null;
                    Double record=Double.MIN_VALUE;
                    
                    for(SNDRCPricingProblem pricingProblem:temp.keySet()){
                        Double count=temp.get(pricingProblem);
                        if(count>record){
                            mostUsedPricingProblem=pricingProblem;
                            record=count;
                        }
                    }
                    
                    for(int edgeIndex=0;edgeIndex<dataModel.numService;edgeIndex++){
                        double reducedCost=mostUsedPricingProblem.dualCost;
                        double[] reducedCosts=mostUsedPricingProblem.dualCosts;
                        

                        accumulatedReducedCost[edgeIndex]+=reducedCost+reducedCosts[edgeIndex]-dataModel.fixedCost[mostUsedPricingProblem.originNodeO][mostUsedPricingProblem.capacityTypeS];
                    }
                    
                }
                
                if(this.nodesProcessed % nodeFre==0&&this.nodesProcessed!=0&&ifUseLearningUB){
                    try {
                        LearningUB();
                    } catch (TimeLimitExceededException | IloException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                    
//                    this.ifUseLearningUB=false;
                    
                    edgeFrequency=new double[dataModel.numServiceArc];
                    cutFrequency=new int[dataModel.numServiceArc];
                    accumulatedReducedCost=new double[dataModel.numServiceArc];
                    
                }
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
                    
                    this.ifOptGetFromSubGraph=false;
                    
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

                 //An acceleration technique for ub
                
                if(ifAccelerationForUB){
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
                }



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
                    
                    //if node bound doesn't improve, we record its two children by leading branch, add these branch to master
//                    if(Math.abs(parentBound-bapNode.getBound())<0.00001){
//                        for(BAPNode<SNDRC, Cycle> child:newBranches){
//                            ((Master) master).AddBranchDecisionForCut(child.getBranchingDecision());
//                        }
//                    }
                    
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
                        
//                        break;  //fix only one cycle each time
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
                    this.ifOptGetFromSubGraph=false;

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
        System.out.println();
        
//        List<Map<Integer,Double>> optXValues=((Master) master).getXValues();
//        //output x variables
//        for(int demand=0;demand<dataModel.numDemand;demand++){
//            for(int edgeIndex:optXValues.get(demand).keySet()){
//                if(optXValues.get(demand).get(edgeIndex)>0.01){
//                    Edge edge=dataModel.edgeSet.get(edgeIndex);
//                    System.out.println("x[" + demand + "]:" + edge.start + "->" + edge.end +" "+edge.u+","+edge.t1+"->"+edge.v+","+edge.t2+" = "+optXValues.get(demand).get(edgeIndex));
//                }
//            }
//            System.out.println();
//        }
    }
    
    /**
     * According to the statistic information from edgeFrequency,cutFrequency and accumulatedReducedCost, 
     * we pick up some important edges to build a subgraph. If the problem is infeasible, we adjust some parameters to enlarge the edge set.
     * @throws IloException 
     * @throws TimeLimitExceededException 
     */
    public void LearningUB() throws TimeLimitExceededException, IloException{
        
        Set<Integer> serviceEdgeSet=new TreeSet<>();
//        int startTime=(int) (Math.random()*this.timeCompress);
        
        
        Set<Integer> finalSet=new HashSet<>();
        
        for(int startTime=0;startTime<this.timeCompress;startTime+=100){
            
            serviceEdgeSet=new HashSet<>();
            
            for(int serviceIndex=0;serviceIndex<dataModel.numService;serviceIndex++){
                Service service=dataModel.serviceSet.get(serviceIndex);
                
                
                for(int time=0;time<dataModel.timePeriod;time+=timeCompress){
                   int time0=time;
                   int timeLast;
                   int time1=time0+timeCompress-1;
                   if(time1>=dataModel.timePeriod){
                       time1=dataModel.timePeriod-1;
                   }
                   timeLast=time1-time0+1;
                   
                   time0=(time0+startTime)%dataModel.timePeriod;
                   time1=(time1+startTime)%dataModel.timePeriod;
                   
                   Set<Integer> timeSet=new HashSet<>();
                   for(int t=0;t<timeLast;t++){
                       int currentTime=(time0+t)%dataModel.timePeriod;
                       timeSet.add(currentTime);
                   }
                   
                   // edgeIndex=serviceIndex*timePeriod+t
                   double sum=0;
                   int maxEdgeIndex=-1;
                   double record=Double.MIN_VALUE;
                   Set<Integer> edgeIndexSet=new HashSet<>();
                   
                   for(int t:timeSet){
                       int tempEdgeIndex=serviceIndex*dataModel.timePeriod+t;
                       edgeIndexSet.add(tempEdgeIndex);
                       
                       sum+=edgeFrequency[tempEdgeIndex];
                       if(record<edgeFrequency[tempEdgeIndex]){
                           record=edgeFrequency[tempEdgeIndex];
                           maxEdgeIndex=tempEdgeIndex;
                       }
                   }
                   
                   if(sum>alphaForEdgeFre*nodeFre*timeLast){
                       serviceEdgeSet.add(maxEdgeIndex);
                       
                       
                       
                       edgeIndexSet.remove(maxEdgeIndex);
                       for(int edgeIndex:edgeIndexSet){
                           if((edgeFrequency[maxEdgeIndex]-edgeFrequency[edgeIndex])/edgeFrequency[maxEdgeIndex]<0.1){
                               serviceEdgeSet.add(edgeIndex);
                           }
                       }
                       
                       
                   }
                   
                   
                       
                }
                
            }
            
            
            
            for(int edgeIndex:serviceEdgeSet){
                if(!finalSet.contains(edgeIndex)){
                    finalSet.add(edgeIndex);
                }
            }
            
        }
        
        serviceEdgeSet=finalSet;

        

        
        

        
        
        //if the serviceEdgeSet has a percentage of over leanringCheckPercent different from last subgraph, we will solve the new sub problem
        int count=0;
        for(int edgeIndex=0;edgeIndex<dataModel.numServiceArc;edgeIndex++){
            if(subEdgeRecord[edgeIndex]&&!serviceEdgeSet.contains(edgeIndex)){
                count++;
            }
            if(!subEdgeRecord[edgeIndex]&&serviceEdgeSet.contains(edgeIndex)){
                count++;
            }
            
        }
        
        if(count>=dataModel.numServiceArc*leanringCheckPercent){
            
            for(int edgeIndex=0;edgeIndex<dataModel.numServiceArc;edgeIndex++){
                if(serviceEdgeSet.contains(edgeIndex)){
                    subEdgeRecord[edgeIndex]=true;
                }else{
                    subEdgeRecord[edgeIndex]=false;
                }
            }
            
            
          // output learning information
          System.out.println("Yes");
//          System.out.println(Arrays.toString(edgeFrequency));
//          System.out.println(serviceEdgeSet.toString());
          System.out.println(serviceEdgeSet.size());
          System.out.println();
//          
//          for(int edgeIndex=0;edgeIndex<dataModel.numServiceArc;edgeIndex=edgeIndex+3){
//              System.out.println(edgeIndex+"-"+(edgeIndex+2)+": ");
//              for(int i=0;i<3;i++){
//                  int currentEdgeIndex=edgeIndex+i;
//                  if(currentEdgeIndex<dataModel.numServiceArc){
//                      System.out.print(edgeFrequency[currentEdgeIndex]+" ");
//                  }
//              }
//              System.out.println();
//              
//              for(int i=0;i<3;i++){
//                  int currentEdgeIndex=edgeIndex+i;
//                  if(currentEdgeIndex<dataModel.numServiceArc){
//                      System.out.print(cutFrequency[currentEdgeIndex]+" ");
//                  }
//              }
//              System.out.println();
//              
//              for(int i=0;i<3;i++){
//                  int currentEdgeIndex=edgeIndex+i;
//                  if(currentEdgeIndex<dataModel.numServiceArc){
//                      System.out.print(accumulatedReducedCost[currentEdgeIndex]+" ");
//                  }
//              }
//              
//              System.out.println();
//              System.out.println();
//          }
          
          
            
            // after the built of serviceEdgeSet, we set up a new sub problem and solve it by branch and price
            SNDRC subGraph=new SNDRC(dataModel,serviceEdgeSet);
            
//            subGraph.isFeasibleForX=true;
//            System.out.println(subGraph.isFeasibleForX);
//            subGraph.isFeasibleForX=true;
            if(subGraph.isFeasibleForX){
                
                
//--------------------------------------------------------------------------------------------------------------------------------------                
//                long time0 = System.currentTimeMillis();
//                ColumnGenerationBasedHeuristic solver = new ColumnGenerationBasedHeuristic(subGraph, 0.65, true);
//                solver.Solve();
//                long time1 = System.currentTimeMillis();
//                System.out.println("Total time= " + (time1 - time0));
//--------------------------------------------------------------------------------------------------------------------------------------     
                
                
                this.ifUseLearningUB=false;
                
                //output subEdgeSet
//                Set tempSet=new TreeSet<>();
//                for(int edgeIndex:serviceEdgeSet){
//                    tempSet.add(edgeIndex);
//                }
//                System.out.println(tempSet.toString());
                
                //Create the pricing problems
                List<SNDRCPricingProblem> subPricingProblems=new LinkedList<SNDRCPricingProblem>();
                for(int capacityType=0;capacityType<subGraph.numOfCapacity;capacityType++) {
                    for(int originNode=0;originNode<subGraph.numNode;originNode++) {
                        String name="capacity type: "+capacityType+" origin node: "+originNode;
                        SNDRCPricingProblem subPricingProblem=new SNDRCPricingProblem(subGraph,name,capacityType,originNode);
                        subPricingProblems.add(subPricingProblem);
                    }
                }
                
              //Create a cutHandler
                CutHandler<SNDRC, SNDRCMasterData> subCutHandler=new CutHandler<>();
                StrongInequalityGenerator subCutGen=new StrongInequalityGenerator(subGraph,subPricingProblems,0);
//              subCutHandler.addCutGenerator(subCutGen);
                
              //Create the Master Problem
                Master subMaster=new Master(subGraph,subPricingProblems,subCutHandler,subCutGen,false);
                
               //Define which solvers to use
                List<Class<?extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> subSolvers=Collections.singletonList(ExactPricingProblemSolver.class);
                
                //Define one or more Branch creators
                List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators=Arrays.asList( new BranchOnLocalServiceForAllPricingProblems(subGraph, subPricingProblems, 0.5),new BranchOnLocalService(subGraph, subPricingProblems, 0.5),new BranchOnServiceEdge(subGraph, subPricingProblems, 0.5));
                
              //Create a Branch-and-Price instance
                BranchAndPriceB_M subBap=new BranchAndPriceB_M(subGraph, subMaster, subPricingProblems, subSolvers, branchCreators,this.objectiveIncumbentSolution,0.6,0.3,0.1,20,0.5,5,0.1,false,true);
//              bap.setNodeOrdering(new BFSbapNodeComparator());
                subBap.setNodeOrdering(new NodeBoundbapNodeComparator());
                
                BapLoggerB_M logger=new BapLoggerB_M(subBap, new File("./output/subBAPlogger.log"));
                
//                subBap.runBranchAndPrice(System.currentTimeMillis()+14400000L); // 4 hours
                
                //the running time of learningUB() is 5h-time has ran
                long runtime=System.currentTimeMillis()-this.runtime;
                long timeLeft=18000000-runtime;
//                long timeLeft=20000-runtime;
//                System.out.println("runTime="+runtime);
//                System.out.println("timeLeft="+timeLeft);
                subBap.runBranchAndPrice(System.currentTimeMillis()+timeLeft);
                
                if(subBap.hasSolution()){
                    if(subBap.objectiveIncumbentSolution<this.objectiveIncumbentSolution){
                        
                        this.objectiveIncumbentSolution = subBap.objectiveIncumbentSolution;
                        this.upperBoundOnObjective = subBap.objectiveIncumbentSolution;
                        this.incumbentSolution = subBap.getSolution();

                        optSolutionValueMap = new HashMap<>();
                        optSolutionValueMap=subBap.optSolutionValueMap;
                        
                        optXValues=subBap.optXValues;
                        
                        this.ifOptGetFromSubGraph=true;
                    }
                }
                
                subBap.close();
                subCutHandler.close();
            }
            
            
        }else{
            System.out.println("No");
        }
        
    }
    
 
   public boolean GetIfOptGetFromSubGraph(){
	   return ifOptGetFromSubGraph;
   }
    
}


