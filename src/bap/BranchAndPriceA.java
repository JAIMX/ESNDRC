package bap;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Scanner;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.colgenMain.ColGen;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;

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
import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;
import logger.BapLoggerA;
import model.SNDRC;
import model.SNDRC.Demand;
import model.SNDRC.Edge;
import model.SNDRC.Service;

public class BranchAndPriceA <V> extends AbstractBranchAndPrice<SNDRC, Cycle, SNDRCPricingProblem>{

    private double thresholdValue;
    private PriorityQueue<BAPNode<SNDRC, Cycle>> lowBoundQueue;

    private final double probLB;
    private final double c;
    private int nrNonImproForAcce;
    private Map<Cycle, Double> optSolutionValueMap;
    private List<Map<Integer, Double>> optXValues;
    private List<Map<Integer, Double>> xValuesForRootLP;
//    private double[] nodeBoundRecord;
//    private int helpOutPut;
    
    //For learning upper bound
    private int[] cutFrequency;
    private double[] edgeFrequency,accumulatedReducedCost;
    private int nodeFre;
    private double alphaForEdgeFre;
    private int timeCompress;
    private boolean ifUseLearningUB;
    private final boolean ifOptGetFromSubGraph;
    private int cycleRecordForIntesification;
    private Set<Cycle> cycleSetRecord;
    private int freqForInten0=5;
    private int freqForIntensification=freqForInten0;
    
    
    
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
     */

    public BranchAndPriceA(SNDRC modelData, Master master, List<SNDRCPricingProblem> pricingProblems,
            List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
            List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators,
            double objectiveInitialSolution, double thresholdValue, double probLB, double c,int nodeFre,double alphaForEdgeFre,int timeCompress,boolean ifUseLearningUB) {
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

    	int[] temp=new int[dataModel.numService];
        List<Cycle> artificalVars = new ArrayList<Cycle>();
        // for weak forcing constraints(ifForResourceBoundConstraints=0)
        for (int edgeIndex = 0; edgeIndex < dataModel.numServiceArc; edgeIndex++) {
            Set<Integer> set = new HashSet<>();
            HashSet<Integer> set2=new HashSet<>();
            set.add(edgeIndex);
            Cycle cycle = new Cycle(pricingProblems.get(0), true, "Artificial", set, 100000000, 0, 0,temp,set2);
            artificalVars.add(cycle);
        }

        // for resource bound constraints(ifForResourceBoundConstraints=1)
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            Set<Integer> set = new HashSet<>();
            HashSet<Integer> set2=new HashSet<>();
            Cycle cycle = new Cycle(pricingProblem, true, "Artificial", set, 100000000, 0, 1,temp,set2);
            artificalVars.add(cycle);
        }

        // for holding edge branch constraints(ifForResourceBoundConstraints=2)
        Set<Integer> set = new HashSet<>();
        HashSet<Integer> set2=new HashSet<>();
        Cycle cycle = new Cycle(pricingProblems.get(0), true, "Artificial", set, 100000000, 0, 2,temp,set2);
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
        this.cycleSetRecord=new HashSet<Cycle>();
        Map<Cycle,Double> lpSumColumnCount=new HashMap<>();

        
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
                
                if(bapNode.nodeID==0){
                    System.out.println("root node bound= "+bapNode.getBound());
                    
                    //explore the performance on non-improvement LP solutions
//                    List<Cycle> solution =bapNode.getSolution();
//                    for(Cycle cycle:solution){
//                        System.out.println(cycle);
//                        System.out.println(out(cycle) + ":" + cycle.value);
//                        System.out.println();
//                    }
                    

                    
                    xValuesForRootLP=new ArrayList<>();
                    xValuesForRootLP=((Master) master).getXValues();
                    
                    
                    // output x variables
//                    for (int demand = 0; demand < dataModel.numDemand; demand++) {
//                        for (int edgeIndex : xValuesForRootLP.get(demand).keySet()) {
//                            if (xValuesForRootLP.get(demand).get(edgeIndex) > 0.01) {
//                                Edge edge;
//
//                                if (!ifOptGetFromSubGraph) {
//                                    edge = dataModel.edgeSet.get(edgeIndex);
//                                } else {
//                                    edge = dataModel.subEdgeSet.get(edgeIndex);
//                                }
//
//                                
//                                if(edge.edgeType==0){
//                                    System.out.println("x[" + demand + "]:" + edge.u + "," + edge.t1 + "->" + edge.v + "," + edge.t2
//                                            + "= " + xValuesForRootLP.get(demand).get(edgeIndex) + " " + edge.duration);
//                                }
//
//                            }
//                        }
//                        System.out.println();
//                    }
//                    
//                    System.out.println("End");
                    
                }

                // output the model
//                ((Master) master).Output(bapNode.nodeID);
//                nodeBoundRecord[bapNode.nodeID] = bapNode.getBound();

                

            } catch (IloException | TimeLimitExceededException e) {
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
                
                if(this.nodesProcessed % nodeFre==0&&this.nodesProcessed!=0){
                    LearningUB();
                    
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
            
            //record cycles for intensification
            cycleRecordForIntesification++;
            for(Cycle cycle:bapNode.getSolution()){
            	if(lpSumColumnCount.containsKey(cycle)){
            		lpSumColumnCount.put(cycle, lpSumColumnCount.get(cycle)+cycle.value);
            	}else{
            		lpSumColumnCount.put(cycle, cycle.value);
            	}
            	

            }
            for(SNDRCPricingProblem pri:pricingProblems) {
            	for(Cycle cycle:master.getColumns(pri)) {
                	if(!cycleSetRecord.contains(cycle)) {
                		cycleSetRecord.add(cycle);           		
                	}
            	}
            }
            
            if(cycleRecordForIntesification>=freqForIntensification){
            	
            	System.out.println("Before intensification, we have "+cycleSetRecord.size()+" columns");
            	System.out.println("After pick up:");
            	Set<Cycle> subCycleSet=subCycle(lpSumColumnCount,cycleSetRecord);
            	try {
					Intensification(subCycleSet);
					break;
				} catch (IloException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
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

                 //An acceleration technique for ub
                
                if(!ifUseLearningUB){
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
                    if(Math.abs(parentBound-bapNode.getBound())<0.00001){
                        for(BAPNode<SNDRC, Cycle> child:newBranches){
                            ((Master) master).AddBranchDecisionForCut(child.getBranchingDecision());
                        }
                        
                    }
                    
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
    
    public Set<Cycle> subCycle(Map<Cycle, Double>lpSumColumnCount, Set<Cycle> cycleSet){

    	//note that we store artificial variables
    	Set<Cycle> artificialCycle=new HashSet<>();
    	for(Cycle cycle:cycleSet){
    		if(cycle.isArtificialColumn){
    			artificialCycle.add(cycle);
    		}
    	}
    	cycleSet.removeAll(artificialCycle);
    	
//    	System.out.println(lpSumColumnCount.toString());
//    	System.out.println(cycleSet.size());
//    	System.out.println(artificialCycle.size());
    	
    	//classify by pattern,we transfer int[] to string to recognize
    	Map<String,ArrayList<Cycle>> group=new HashMap<>();
    	for(Cycle cycle:cycleSet){
    		String string=Arrays.toString(cycle.pattern);
    		if(!group.containsKey(string)){
    			ArrayList<Cycle> list=new ArrayList<>();
    			list.add(cycle);
    			group.put(string, list);
    		}else{
    			group.get(string).add(cycle);
    		}
    	}
    	
//    	System.out.println("group");
//    	for(String string:group.keySet()){
//    		System.out.println(string+" "+group.get(string).size());
//    	}
    	
    	//calculate distances in each group and record them in priority queue
    	PriorityQueue<CyclePair> pq=new PriorityQueue<CyclePair>();
    	for(String string:group.keySet()){
    		ArrayList<Cycle> list=group.get(string);
    		for(int i=0;i<list.size()-1;i++){
    			for(int j=i+1;j<list.size();j++){
    				Cycle cycle1=list.get(i);
    				Cycle cycle2=list.get(j);
    				CyclePair pair=new CyclePair(cycle1, cycle2, calDistance(cycle1,cycle2));
    				pq.add(pair);
    			}
    		}
    	}
    	
    	//delete cycles
    	while(cycleSet.size()>1000&&pq.size()>0){
    		CyclePair pair=pq.poll();
    		if(cycleSet.contains(pair.cycle1)&&cycleSet.contains(pair.cycle2)){
    			if(lpSumColumnCount.containsKey(pair.cycle1)&&lpSumColumnCount.containsKey(pair.cycle2)) {
            		if(lpSumColumnCount.get(pair.cycle1)>lpSumColumnCount.get(pair.cycle2)){
            			cycleSet.remove(pair.cycle2);
            		}else{
            			cycleSet.remove(pair.cycle1);
            		}
    			}else {
    				boolean temp=false;
    				if(lpSumColumnCount.containsKey(pair.cycle1)){
    					cycleSet.remove(pair.cycle2);
    					temp=true;
    				}
    				if(lpSumColumnCount.containsKey(pair.cycle2)){
    					cycleSet.remove(pair.cycle1);
    					temp=true;
    				}
    				if(!temp&&pair.cycle1.value>pair.cycle2.value) {
    					cycleSet.remove(pair.cycle2);
    				}else {
    					if(!temp){
        					cycleSet.remove(pair.cycle1);
    					}
    				}
    			}
    		}
    	}
    	
//    	cycleSet.addAll(artificialCycle);
    	return cycleSet;
    	
    }
    
    public double calDistance(Cycle cycle1,Cycle cycle2){
    	//for each service,we record the start time and calculate the difference
    	ArrayList<Integer>[] record1=(ArrayList<Integer>[]) new ArrayList[dataModel.numService];
    	ArrayList<Integer>[] record2=(ArrayList<Integer>[]) new ArrayList[dataModel.numService];
    	for(int edgeIndex:cycle1.edgeIndexSet){
    		Edge edge=dataModel.edgeSet.get(edgeIndex);
    		if(edge.edgeType==0){ //service edge
    			if(record1[edge.serviceIndex]==null){
    				ArrayList<Integer> list=new ArrayList<>();
    				list.add(edge.t1);
    				record1[edge.serviceIndex]=list;
    			}else{
    				record1[edge.serviceIndex].add(edge.t1);
    			}
    		}
    	}
    	
    	for(int edgeIndex:cycle2.edgeIndexSet){
    		Edge edge=dataModel.edgeSet.get(edgeIndex);
    		if(edge.edgeType==0){ //service edge
    			if(record2[edge.serviceIndex]==null){
    				ArrayList<Integer> list=new ArrayList<>();
    				list.add(edge.t1);
    				record2[edge.serviceIndex]=list;
    			}else{
    				record2[edge.serviceIndex].add(edge.t1);
    			}
    		}
    	}
    	
    	double sum=0;
    	int count=0;
    	for(int i=0;i<dataModel.numService;i++){
    		ArrayList<Integer> list1=record1[i];
    		ArrayList<Integer> list2=record2[i];
    		

    		if(list1!=null){
    			Collections.sort(list1);
    			Collections.sort(list2);
        		for(int j=0;j<list1.size();j++){
        			sum+=Math.abs(list1.get(j)-list2.get(j));
        			if(list1.get(j)!=list2.get(j)) count++;
        		}
    		}

    	}
    	
    	return sum/count; 	
    	
    }
    
    private class CyclePair implements Comparable<CyclePair>{
    	private Cycle cycle1,cycle2;
    	private double distance;
    	
    	public CyclePair(Cycle cycle1,Cycle cycle2,double distance){
    		this.cycle1=cycle1;
    		this.cycle2=cycle2;
    		this.distance=distance;
    	}
    	
    	@Override
    	public int compareTo(CyclePair that){
            if(this.distance<that.distance-0.001){
    			return -1;
    		}else{
    			if(this.distance>that.distance+0.001){
    				return 1;
    			}
    		}
            return 0;
    	}
    }
    
    public void Intensification(Set<Cycle> cycleSet) throws IloException{
    	
    	System.out.println("==================Intensification===================");
    	System.out.println("We add "+cycleSet.size()+" columns to cplex.");
		IloCplex cplex=new IloCplex();
		
//		FileOutputStream outputStream=new FileOutputStream(new File("./output/cplexOut.txt"));
//		cplex.setOut(outputStream);
		
		cplex.setParam(IloCplex.IntParam.Threads, 4);
		cplex.setParam(IloCplex.Param.Simplex.Tolerances.Markowitz, 0.1);
		cplex.setParam(IloCplex.DoubleParam.TiLim, 7200); //2 hours
		cplex.setParam(IloCplex.DoubleParam.EpGap,0.005);
		List<Map<Integer,IloNumVar>> x; //map:edgeIndex, x variable
//		Map<Path,IloNumVar> pathVarMap=new HashMap<>();
		
		// Define variables x
		x=new ArrayList<Map<Integer,IloNumVar>>();
		for(int p=0;p<dataModel.numDemand;p++) {
			Map<Integer,IloNumVar> initialX=new HashMap<Integer,IloNumVar>();
			x.add(initialX);
		}
		
		// add x variables
		for(int p=0;p<dataModel.numDemand;p++) {
			for(int edgeIndex:dataModel.edgesForX.get(p)) {
				Edge edge=dataModel.edgeSet.get(edgeIndex);
				IloNumVar varX=cplex.numVar(0, dataModel.demandSet.get(p).volume,"x" +p+","+ edge.start + "," + edge.end );
				x.get(p).put(edgeIndex, varX);
			}
		}
		
		
		// Define the objective
		/**
		 * Here we assume the cost of edge AT is 0
		 */
		IloLinearNumExpr exprObj = cplex.linearNumExpr();
		
		for(int p=0;p<dataModel.numDemand;p++) {
			Map<Integer,IloNumVar> map=x.get(p);
			for(int edgeIndex:map.keySet()) {
				exprObj.addTerm(dataModel.beta*dataModel.edgeSet.get(edgeIndex).duration, map.get(edgeIndex));
			}
		}
		

		IloObjective obj = cplex.addMinimize(exprObj);
		
		
		
		// Define flowBalanceConstraints
		IloRange[][] flowBalanceConstraints = new IloRange[dataModel.numDemand][dataModel.abstractNumNode];

		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int p = 0; p < dataModel.numDemand; p++) {
			Map<Integer,IloNumVar> map=x.get(p);
			
			for (int i = 0; i < dataModel.abstractNumNode; i++) {
				expr.clear();
				// edges which point from i
				for (int edgeIndex : dataModel.pointToEdgeSet.get(i)) {
					if(map.containsKey(edgeIndex)) {
						expr.addTerm(1, map.get(edgeIndex));
					}
				}

				// edges which point to i
				for (int edgeIndex : dataModel.pointFromEdgeSet.get(i)) {
					if(map.containsKey(edgeIndex)) {
						expr.addTerm(-1, map.get(edgeIndex));
					}
				}
				flowBalanceConstraints[p][i] = cplex.addEq(dataModel.b[p][i], expr);

			}
		}
		
		
		// Define weakForcingConstraints
		IloRange[] weakForcingConstraints = new IloRange[dataModel.numServiceArc];
		for (int arcIndex = 0; arcIndex < dataModel.numServiceArc; arcIndex++) {
			expr.clear();
			for (int p = 0; p < dataModel.numDemand; p++) {
				if(x.get(p).containsKey(arcIndex)) {
					expr.addTerm(1, x.get(p).get(arcIndex));
				}
			}

			weakForcingConstraints[arcIndex] = cplex.addGe(0, expr);
		}
		
		
		
		
		// Define resourceBoundConstraints
		IloRange[][] resourceBoundConstraints = new IloRange[dataModel.numOfCapacity][dataModel.numNode];
		for (int s = 0; s < dataModel.numOfCapacity; s++) {
			for (int o = 0; o < dataModel.numNode; o++) {
//				expr.clear();
//				expr.addTerm(1, q[s][o]);
				expr.clear();
//				resourceBoundConstraints[s][o] = cplex.addEq(dataModel.vehicleLimit[s][o], 0);
//				resourceBoundConstraints[s][o]=cplex.range(0,dataModel.vehicleLimit[s][o]);
				resourceBoundConstraints[s][o]=cplex.addRange(0,dataModel.vehicleLimit[s][o]);
			}
		}
		
        //Define storeBoundConstraints
        IloRange[][] storeBoundConstraints=new IloRange[dataModel.numNode][dataModel.timePeriod];
        for(int l=0;l<dataModel.numNode;l++) {
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
        IloRange[][] chargeBoundConstraints=new IloRange[dataModel.numNode][dataModel.timePeriod];
        
        for(int l=0;l<dataModel.numNode;l++) {
        	for(int t=0;t<dataModel.timePeriod;t++) {
        		expr.clear();
        		chargeBoundConstraints[l][t]=cplex.addGe(dataModel.chargeLimit[l], expr);
        	}
        }
		
		
		
		
		// add all columns
        int count=0;
        Map<Cycle,IloNumVar> cycleVarMap=new HashMap<>();
		for(Cycle cycle:cycleSet) {
			// Register column with objective
			IloColumn iloColumn = cplex.column(obj, cycle.cost);

			// weak forcing constraints
			for (int edgeIndex : cycle.edgeIndexSet) {
				if(dataModel.edgeSet.get(edgeIndex).edgeType==0) {
					iloColumn = iloColumn.and(cplex.column(weakForcingConstraints[edgeIndex],
							-dataModel.capacity[cycle.associatedPricingProblem.capacityTypeS]));
				}
			}
			
			
			// resource bound constraints
			iloColumn = iloColumn.and(cplex.column(resourceBoundConstraints[cycle.associatedPricingProblem.capacityTypeS][cycle.associatedPricingProblem.originNodeO],1));
			
			//charge bound constraints
			for(int chargeEdgeIndex:cycle.ifCharge){
        		int chargeNodeIndex=dataModel.edgeSet.get(chargeEdgeIndex).start;
        		int l=chargeNodeIndex/dataModel.timePeriod;
        		int t=chargeNodeIndex%dataModel.timePeriod;
        		iloColumn=iloColumn.and(cplex.column(chargeBoundConstraints[l][t],1));
			}
//			cplex.exportModel("check.lp");
			
			
			// Create the variable and store it
//			IloNumVar var =cplex.intVar(iloColumn, 0, dataModel.vehicleLimit[capacityType][originNode],"z_" + capacityType + ","+originNode+","+count);
			IloNumVar var =cplex.intVar(iloColumn, 0, Integer.MAX_VALUE,"z_" + cycle.associatedPricingProblem.capacityTypeS + ","+cycle.associatedPricingProblem.originNodeO+","+count);
			cycleVarMap.put(cycle, var);
			count++;
			
		}
		
		
		
		cplex.solve();
		
		
        int integerObjective = MathProgrammingUtil.doubleToInt(cplex.getObjValue());
        
        if (integerObjective < this.upperBoundOnObjective) {
            this.objectiveIncumbentSolution = integerObjective;
            this.upperBoundOnObjective = integerObjective;
            this.incumbentSolution = new ArrayList<>();
            for (Cycle cycle : cycleSet) {
            	int value=MathProgrammingUtil.doubleToInt(cplex.getValue(cycleVarMap.get(cycle)));
            	if(value>0){
                    this.incumbentSolution.add(cycle);
            	}
            }
            
            optSolutionValueMap = new HashMap<>();
            for (Cycle cycle : incumbentSolution) {
            	double value=cplex.getValue(cycleVarMap.get(cycle));
                optSolutionValueMap.put(cycle, value);
            }
            
            optXValues=new ArrayList<>();
    		for (int demand = 0; demand < dataModel.numDemand; demand++) {
    			Map<Integer,Double> map=new HashMap<>();
    			for (int edgeIndex: x.get(demand).keySet()) {
    				map.put(edgeIndex, cplex.getValue(x.get(demand).get(edgeIndex)));
    			}
    			optXValues.add(map);
    		}

    		freqForIntensification=freqForInten0;
    		System.out.println("We use intensification finding a better solution: "+cplex.getObjValue());
        }else{
            freqForIntensification=freqForIntensification*2;
        }
        
        //update cycleRecordForIntesification
        cycleRecordForIntesification=0;
        cycleSetRecord.clear();
        cycleSetRecord.addAll(incumbentSolution);
		
    }
    
    public String out(Cycle column) {

        Queue<Edge> path = new PriorityQueue<>();
        Set<Edge> chargeEdgeSet=new HashSet<>();
        for(int edgeIndex:column.edgeIndexSet){
        	if(column.ifCharge.contains(edgeIndex)) chargeEdgeSet.add(dataModel.edgeSet.get(edgeIndex));
        }

        if (!ifOptGetFromSubGraph) {
            for (int edgeIndex : column.edgeIndexSet) {
                path.add(dataModel.edgeSet.get(edgeIndex));
            }
        } else {
            for (int edgeIndex : column.edgeIndexSet) {
                path.add(dataModel.subEdgeSet.get(edgeIndex));
            }
        }

        StringBuilder pathRecord = new StringBuilder();

        Edge edge = null;
        int size = path.size();
        for (int i = 0; i < size; i++) {

            edge = path.poll();
            pathRecord.append('(');
            pathRecord.append(edge.u);
            pathRecord.append(',');
            pathRecord.append(edge.t1);
            pathRecord.append(')');

            pathRecord.append("->");
            if(chargeEdgeSet.contains(edge)){
            	pathRecord.append("charge");
            }

        }

        pathRecord.append('(');
        pathRecord.append(edge.v);
        pathRecord.append(',');
        pathRecord.append(edge.t2);
        pathRecord.append(')');

        return pathRecord.toString();

    }

    public Double CalculateProb() {
        double prob = 2 - probLB;
        prob -= (2 - 2 * probLB) / (1 + Math.pow(Math.E, -c * nrNonImproForAcce));
//        return prob;
        return -1.0;
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
    public List<Map<Integer, Double>> GetxValuesForRootLP() {
        return xValuesForRootLP;
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
     */
    public void LearningUB(){
        
        Set<Integer> serviceEdgeSet=new HashSet<>();
        int startTime=(int) (Math.random()*this.timeCompress);
        
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
               for(int t:timeSet){
                   int tempEdgeIndex=serviceIndex*dataModel.timePeriod+t;
                   sum+=edgeFrequency[tempEdgeIndex];
                   if(record<edgeFrequency[tempEdgeIndex]){
                       record=edgeFrequency[tempEdgeIndex];
                       maxEdgeIndex=tempEdgeIndex;
                   }
               }
               
               if(sum>alphaForEdgeFre*nodeFre*timeLast){
                   serviceEdgeSet.add(maxEdgeIndex);
               }
                   
            }
            
        }
        
        
//        System.out.println(Arrays.toString(edgeFrequency));
//        System.out.println(serviceEdgeSet.size());
        
        
        
        // after the built of serviceEdgeSet, we set up a new sub problem and solve it by branch and price
        SNDRC subGraph=new SNDRC(dataModel,serviceEdgeSet);
        
//        subGraph.isFeasibleForX=true;
//        System.out.println(subGraph.isFeasibleForX);
//        subGraph.isFeasibleForX=true;
        if(subGraph.isFeasibleForX){
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
//          subCutHandler.addCutGenerator(subCutGen);
            
          //Create the Master Problem
            Master subMaster=new Master(subGraph,subPricingProblems,subCutHandler,subCutGen,false);
            
           //Define which solvers to use
            List<Class<?extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> subSolvers=Collections.singletonList(ExactPricingProblemSolver.class);
            
            //Define one or more Branch creators
            List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators=Arrays.asList( new BranchOnLocalServiceForAllPricingProblems(subGraph, subPricingProblems, 0.5),new BranchOnLocalService(subGraph, subPricingProblems, 0.5),new BranchOnServiceEdge(subGraph, subPricingProblems, 0.5));
            
          //Create a Branch-and-Price instance
            BranchAndPriceA subBap=new BranchAndPriceA(subGraph, subMaster, subPricingProblems, subSolvers, branchCreators,this.objectiveIncumbentSolution,0.6,0.3,0.1,20,0.5,3,false);
//          bap.setNodeOrdering(new BFSbapNodeComparator());
            subBap.setNodeOrdering(new NodeBoundbapNodeComparator());
            
            BapLoggerA logger=new BapLoggerA(subBap, new File("./output/subBAPlogger.log"));
            
            subBap.runBranchAndPrice(System.currentTimeMillis()+3600000L); // one hour
            if(subBap.hasSolution()){
                if(subBap.objectiveIncumbentSolution<this.objectiveIncumbentSolution){
                    
                    this.objectiveIncumbentSolution = subBap.objectiveIncumbentSolution;
                    this.upperBoundOnObjective = subBap.objectiveIncumbentSolution;
                    this.incumbentSolution = subBap.getSolution();

                    optSolutionValueMap = new HashMap<>();
                    optSolutionValueMap=subBap.optSolutionValueMap;
                    
                    optXValues=subBap.optXValues;
                    
                    
                }
            }
            
            subBap.close();
            subCutHandler.close();
        }
        

        
    }
    
    
    
    public boolean GetIfOptGetFromSubGraph(){
 	   return ifOptGetFromSubGraph;
    }
    
    
    public void output(String filename,BranchAndPriceA bap) throws FileNotFoundException{
    	PrintWriter out=new PrintWriter(filename);
    	
    	out.println("LP information:");
    	for(int k=0;k<dataModel.numDemand;k++){
    		Map<Integer,Double> map=(Map<Integer, Double>) bap.GetxValuesForRootLP().get(k);
    		out.println(k);
    		for(int edgeIndex:map.keySet()){
    			double value=map.get(edgeIndex);
    			if(value>0.0001&&dataModel.edgeSet.get(edgeIndex).edgeType==0){
    				out.print(edgeIndex+" "+value+" ");
    			}
    		}
    		out.println();
    		
    	}
}
}
