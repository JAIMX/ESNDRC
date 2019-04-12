package cg;

import java.util.*;

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

            double[] dpFunction = new double[dataModel.abstractNumNode];
            int[] pathRecord = new int[dataModel.abstractNumNode];
            for (int i = 0; i < dpFunction.length; i++) {
                dpFunction[i] = Double.MAX_VALUE;
            }

            // //update for original node
            // //service arcs
            // for(int edgeIndex:dataModel.pointToEdgeSet.get(originNodeIndex)){
            // Edge edge=dataModel.edgeSet.get(edgeIndex);
            // if(edge.edgeType==0) {
            // dpFunction[edge.end]=modifiedCosts[edgeIndex];
            // pathRecord[edge.end]=edgeIndex;
            // }else { // holding arcs
            // dpFunction[edge.end]=0;
            // pathRecord[edge.end]=edgeIndex;
            // }
            // }

            // update for original node
            for (int edgeIndex : dataModel.pointToEdgeSet.get(originNodeIndex)) {
                Edge edge = dataModel.edgeSet.get(edgeIndex);
                dpFunction[edge.end] = modifiedCosts[edgeIndex];
                pathRecord[edge.end] = edgeIndex;
            }

            int durationLimit = dataModel.timePeriod;

            // update for the following nodes
            for (int currentTime = startTime + 1; currentTime < startTime + dataModel.timePeriod; currentTime++) {
                durationLimit--;
                int time = currentTime % dataModel.timePeriod;

                for (int localNode = 0; localNode < dataModel.numNode; localNode++) {
                    int currentNodeIndex = localNode * dataModel.timePeriod + time;
                    if (dpFunction[currentNodeIndex] < Double.MAX_VALUE - 1) {

                        for (int edgeIndex : dataModel.pointToEdgeSet.get(currentNodeIndex)) {
                            Edge edge = dataModel.edgeSet.get(edgeIndex);
                            
                            /**
                             * Note that the duration of holding arcs is 0(duration represents the distance)
                             */
                            
//                            if (edge.edgeType == 0) { // service arcs
//                                if (edge.duration < durationLimit
//                                        || (edge.duration == durationLimit && edge.end == originNodeIndex)) {
//
//                                    if (dpFunction[edge.end] > dpFunction[currentNodeIndex]
//                                            + modifiedCosts[edgeIndex]) {
//                                        dpFunction[edge.end] = dpFunction[currentNodeIndex] + modifiedCosts[edgeIndex];
//                                        pathRecord[edge.end] = edgeIndex;
//                                    }
//
//                                }
//                            } else { // holding arcs
//                                if (durationLimit > 1
//                                        || (durationLimit == 1 && localNode == pricingProblem.originNodeO)) {
//                                    if (dpFunction[edge.end] > dpFunction[currentNodeIndex]) {
//                                        dpFunction[edge.end] = dpFunction[currentNodeIndex];
//                                        pathRecord[edge.end] = edgeIndex;
//                                    }
//                                }
//                            }


                            if (edge.edgeType == 0) { // service arcs
                                if (edge.duration < durationLimit
                                        || (edge.duration == durationLimit && edge.end == originNodeIndex)) {

                                    if (dpFunction[edge.end] > dpFunction[currentNodeIndex]
                                            + modifiedCosts[edgeIndex]) {
                                        dpFunction[edge.end] = dpFunction[currentNodeIndex] + modifiedCosts[edgeIndex];
                                        pathRecord[edge.end] = edgeIndex;
                                    }

                                }
                            } else { // holding arcs
                                if (durationLimit > 1
                                        || (durationLimit == 1 && localNode == pricingProblem.originNodeO)) {
                                    if (dpFunction[edge.end] > dpFunction[currentNodeIndex]+modifiedCosts[edgeIndex]) {
                                        dpFunction[edge.end] = dpFunction[currentNodeIndex]+modifiedCosts[edgeIndex];
                                        pathRecord[edge.end] = edgeIndex;
                                    }
                                }
                            }
                            

                        }


                    }
                }
            }

            if (dpFunction[originNodeIndex] < Double.MAX_VALUE - 1) {
                // if(dpFunction[originNodeIndex]+modifiedCost<-config.PRECISION)
                // {
            	boolean checkIfless=false;
            	if(dpFunction[originNodeIndex] + modifiedCost<this.objective-0.0001) checkIfless=true;
            	this.objective=Math.min(this.objective, dpFunction[originNodeIndex] + modifiedCost);
                if (dpFunction[originNodeIndex] + modifiedCost < -0.1) {
                    Set<Integer> edgeIndexSet = new HashSet<>();
                    double cost = 0;
                    double totalLength = 0;

                    // create edgeIndexSet for this path
                    int currentNodeIndex = originNodeIndex;
                    int lastNodeIndex;

                    do {
                        // update lastNodeIndex
                        Edge edge = dataModel.edgeSet.get(pathRecord[currentNodeIndex]);

                        lastNodeIndex = edge.start;
                        totalLength += edge.duration;
                        edgeIndexSet.add(pathRecord[currentNodeIndex]);

                        // update currentNodeIndex
                        currentNodeIndex = lastNodeIndex;

                    } while (lastNodeIndex != originNodeIndex);

                    cost = dataModel.alpha * totalLength / (dataModel.speed * dataModel.drivingTimePerDay)
                            + dataModel.fixedCost[pricingProblem.originNodeO][pricingProblem.capacityTypeS];

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
                        Cycle cycle = new Cycle(pricingProblem, false, "exactPricing", edgeIndexSet, cost, startTime,0,pattern);
                        if (!pricingProblem.fixCycleSet.contains(cycle)&&checkIfless) {
                        	while(newRoutes.size()>0) newRoutes.remove(newRoutes.size()-1);
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
