package cg;

import java.util.*;

import org.jorlib.demo.frameworks.columnGeneration.cuttingStockCG.cg.PricingProblem;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;

import cg.master.SNDRCMasterData;
import model.SNDRC;
import model.SNDRC.Edge;

public class ExactPricingProblemSolver extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem> {

    private double[] modifiedCosts;
    private double modifiedCost;
    
    

    /**
     * Creates a new solver instance for a particular pricing problem
     * 
     * @param dataModel
     *            data model
     * @param pricingProblem
     *            pricing problem
     */
    public ExactPricingProblemSolver(SNDRC dataModel, SNDRCPricingProblem pricingProblem) {
        super(dataModel, pricingProblem);
        this.name = "ExactShortestPathSolver";

    }

    /**
     * Main method which solves the pricing problem.
     * 
     * 
     * @return List of columns (cycles) with negative reduced cost.
     * @throws TimeLimitExceededException
     *             TimeLimitExceededException
     */

    @Override
    protected List<Cycle> generateNewColumns() throws TimeLimitExceededException {
        List<Cycle> newRoutes = new ArrayList<>();
        this.objective=Double.MAX_VALUE;

        // explore routes starting from different time
        for (int startTime = 0; startTime < dataModel.timePeriod; startTime++) {

            int originNodeIndex = pricingProblem.originNodeO * dataModel.timePeriod + startTime;
//            System.out.println("originNodeIndex="+originNodeIndex);

            double[][] dpFunction = new double[dataModel.abstractNumNode][dataModel.powerCapacity+1];
            int[][] pathRecord = new int[dataModel.abstractNumNode][dataModel.powerCapacity+1];
            int[][] chargeRecord=new int[dataModel.abstractNumNode][dataModel.powerCapacity+1]; //record the power left at last node
            for (int i = 0; i < dpFunction.length; i++) {
                for(int power=0;power<=dataModel.powerCapacity;power++){
                	dpFunction[i][power]=Double.MAX_VALUE;
                	chargeRecord[i][power]=-1;
                }
            }


            // update for original node
            dpFunction[originNodeIndex][dataModel.powerCapacity]=0;

            int durationLimit = dataModel.timePeriod;
            
//            boolean debug=false;
//            if(originNodeIndex==23&&pricingProblem.capacityTypeS==1) debug=true;

            // update for the following nodes
            for (int currentTime = startTime; currentTime < startTime + dataModel.timePeriod; currentTime++) {

                int time = currentTime % dataModel.timePeriod;

                for (int localNode = 0; localNode < dataModel.numNode; localNode++) {
                    int currentNodeIndex = localNode * dataModel.timePeriod + time;
//                    if(debug) System.out.println("currentNodeIndex="+currentNodeIndex);
                    
                    //check dominated label for currentNode
                    for(int power1=0;power1<dataModel.powerCapacity;power1++){
                    	for(int power2=power1+1;power2<=dataModel.powerCapacity;power2++){
                    		if(dpFunction[currentNodeIndex][power2]<dpFunction[currentNodeIndex][power1]-0.0001){
                    			dpFunction[currentNodeIndex][power1]=Double.MAX_VALUE;
                    		}
                    	}
                    }
                    
                    for(int power=0;power<=dataModel.powerCapacity;power++){
//                    	if(debug) System.out.println("power="+power);
                    	
                        if (dpFunction[currentNodeIndex][power] < Double.MAX_VALUE - 1) {
                        	int holdingEdgeIndex=-1;
                            for (int edgeIndex : dataModel.pointToEdgeSet.get(currentNodeIndex)) {
                                Edge edge = dataModel.edgeSet.get(edgeIndex);
                                
                                /**
                                 * Note that the duration of holding arcs is 0(duration represents the distance)
                                 */
                                if (edge.edgeType == 0) { // service arcs
                                    if (edge.duration < durationLimit|| (edge.duration == durationLimit && edge.end == originNodeIndex)){
                                    	if(edge.duration<=power){
                                    		if(dpFunction[edge.end][power-edge.duration]>dpFunction[currentNodeIndex][power]+modifiedCosts[edgeIndex]){
                                    			dpFunction[edge.end][power-edge.duration]=dpFunction[currentNodeIndex][power]+modifiedCosts[edgeIndex];
                                    			pathRecord[edge.end][power-edge.duration]=edgeIndex;
                                    			chargeRecord[edge.end][power-edge.duration]=-1;
//                                    			if(debug) System.out.println("service arc-We update df "+edge.end+","+(power-edge.duration)+"="+dpFunction[edge.end][power-edge.duration]);
                                    		}
                                    	}
                                    }
                                } else { // holding arcs
                                	holdingEdgeIndex=edgeIndex;
                                    if (durationLimit > 1|| (durationLimit == 1 && localNode == pricingProblem.originNodeO)) {
                                        if (dpFunction[edge.end][power] > dpFunction[currentNodeIndex][power]+0) {
                                        	dpFunction[edge.end][power] = dpFunction[currentNodeIndex][power];
                                            pathRecord[edge.end][power] = edgeIndex;
                                            chargeRecord[edge.end][power]=-1;
//                                            if(debug) System.out.println("holding arc-We update df "+edge.end+","+(power)+"="+dpFunction[edge.end][power]);
                                        }
                                    }
                                }
                            }
                            //charge situation
                            if (durationLimit > 1|| (durationLimit == 1 && localNode == pricingProblem.originNodeO)) {
                            	if(power<dataModel.powerCapacity){
                            		int chargePower=power+dataModel.chargeOnceDistance;
                            		chargePower=Math.min(chargePower, dataModel.powerCapacity);
                            		Edge edge=dataModel.edgeSet.get(holdingEdgeIndex);
                            		if(dpFunction[edge.end][chargePower]>dpFunction[currentNodeIndex][power]+modifiedCosts[holdingEdgeIndex]){
                            			dpFunction[edge.end][chargePower]=dpFunction[currentNodeIndex][power]+modifiedCosts[holdingEdgeIndex];
                            			pathRecord[edge.end][chargePower]=holdingEdgeIndex;
//                            			chargeRecord.put(edge.end, power);
                            			chargeRecord[edge.end][chargePower]=power;
//                            			if(debug)System.out.println("charge arc-We update df "+edge.end+","+(chargePower)+"="+dpFunction[edge.end][chargePower]);
                            		}
                            	}
                            }
                        }
                    }
                }
                durationLimit--;
            }

            double minimalValue=Double.MAX_VALUE;
            int powerRecord=-1;
            for(int power=0;power<=dataModel.powerCapacity;power++){
            	if(dpFunction[originNodeIndex][power]<minimalValue-0.001){
            		minimalValue=dpFunction[originNodeIndex][power];
            		powerRecord=power;
            	}
            }
            
            if (dpFunction[originNodeIndex][powerRecord] < Double.MAX_VALUE - 1) {
                // if(dpFunction[originNodeIndex]+modifiedCost<-config.PRECISION)
                // {
            	boolean checkIfless=false;
            	if(dpFunction[originNodeIndex][powerRecord] + modifiedCost<this.objective-0.0001) checkIfless=true;
            	this.objective=Math.min(this.objective, dpFunction[originNodeIndex][powerRecord] + modifiedCost);
                if (dpFunction[originNodeIndex][powerRecord] + modifiedCost < -0.001) {
                    HashSet<Integer> edgeIndexSet = new HashSet<Integer>();
                    HashSet<Integer> chargeEdgeSet=new HashSet<Integer>();
                    double cost = 0;
                    double totalLength = 0;

                    // create edgeIndexSet for this path
                    int currentNodeIndex = originNodeIndex;
                    int powerLeft=powerRecord;
                    do {
                        // update lastNodeIndex
                    	int record=pathRecord[currentNodeIndex][powerLeft];
                    	Edge edge=dataModel.edgeSet.get(record);
                    	if(edge.edgeType==1&&chargeRecord[currentNodeIndex][powerLeft]>=0){ //charge
                            totalLength += edge.duration;
                    		edgeIndexSet.add(record);
                    		chargeEdgeSet.add(record);
                    		powerLeft=chargeRecord[currentNodeIndex][powerLeft];
                    		currentNodeIndex=edge.start;
                    	}else{
                    		totalLength+=edge.duration;
                    		edgeIndexSet.add(record);
                    		powerLeft+=edge.duration;
                    		currentNodeIndex=edge.start;
                    	}
                    } while (currentNodeIndex != originNodeIndex);

                    cost = dataModel.alpha * totalLength / (dataModel.speed * dataModel.drivingTimePerDay)
                            + dataModel.fixedCost[pricingProblem.originNodeO][pricingProblem.capacityTypeS]+chargeEdgeSet.size()*dataModel.chargeObjPara;
                    

                    // if just the start time is different, we only add one
                    // cycle
                    boolean repeat = false;
                    for (Cycle route : newRoutes) {
                        if (edgeIndexSet.equals(route.edgeIndexSet)) {
                            repeat = true;
                            break;
                        }
                        
                    }
                    
                    //calculate pattern
                    int[] pattern=new int[dataModel.numService];
                    for(int edgeIndex:edgeIndexSet) {
                    	Edge edge=dataModel.edgeSet.get(edgeIndex);
                    	if(edge.edgeType==0) {
                        	pattern[edge.serviceIndex]++;
                    	}
                    }

                    if (!repeat) {
                        Cycle cycle = new Cycle(pricingProblem, false, "exactPricing", edgeIndexSet, cost, startTime,0,pattern,chargeEdgeSet);
                        if (!pricingProblem.fixCycleSet.contains(cycle)){  //&&checkIfless) {
//                        	while(newRoutes.size()>0) newRoutes.remove(newRoutes.size()-1);
                            newRoutes.add(cycle);
                        }
                    }

                }
            }

        }

        return newRoutes;
    }

    /**
     * Update the objective function of the pricing problem with the new pricing
     * information (modified costs). The modified costs are stored in the
     * pricing problem.
     */
    @Override
    protected void setObjective() {
    	//For modifiedCosts, the costs of holding arcs represent the costs of charge activity. The cost of practical holding arcs are still 0.
        modifiedCosts = Arrays.copyOf(pricingProblem.dualCosts, pricingProblem.dualCosts.length);
        modifiedCost = pricingProblem.dualCost;
    }

    /**
     * Close the pricing problem
     */
    @Override
    public void close() {
        // cplex.end();
    }
    
	/**
	 * Returns a bound on the objective of the pricing problem. If the pricing problem is solved to optimality, this function would typically return the objective value.
	 * Alternatively, the objective value of a relaxation of the Pricing Problem may be returned, e.g. the LP relaxation when the Pricing Problem is implemented as a MIP, or the value of a Lagrangian Relaxation.
	 * @return a bound on the objective of the pricing problem)
	 */
	public double getBound(){
	    
		return this.objective;
	}

}
