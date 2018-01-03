package cg;

import java.util.*;

import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;

import ilog.concert.IloException;
import ilog.concert.IloObjective;
import ilog.cplex.IloCplex;
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
	 * @return List of columns (cycles) with negative reduced cost.
	 * @throws TimeLimitExceededException TimeLimitExceededException
	 */
	
	@Override
	protected List<Cycle> generateNewColumns()throws TimeLimitExceededException {
		List<Cycle> newRoutes=new ArrayList<>();
		
		//explore routes starting from different time
		for(int startTime=0;startTime<dataModel.timePeriod;startTime++){
			
			int originNodeIndex=pricingProblem.originNodeO*dataModel.timePeriod+startTime;
			
			double[] dpFunction=new double[dataModel.abstractNumNode];
			int[] pathRecord=new int[dataModel.abstractNumNode];
			for(int i=0;i<dpFunction.length;i++){
				dpFunction[i]=Double.MAX_VALUE;
			}
			
			//update for original node
			//service arcs
			for(int edgeIndex:dataModel.pointToEdgeSet.get(originNodeIndex)){
				Edge edge=dataModel.edgeSet.get(edgeIndex);
				dpFunction[edge.end]=modifiedCosts[edgeIndex];
				pathRecord[edge.end]=edgeIndex;
			}
			//holding arcs
			dpFunction[(originNodeIndex+1)%dataModel.timePeriod]=0;
			if(startTime==dataModel.timePeriod-1){
				pathRecord[originNodeIndex-dataModel.timePeriod+1]=-1;
			}else pathRecord[originNodeIndex+1]=-1;
			
			
			
			
			
			
			
		}
		
		
		
		return newRoutes;	
	}
	
	
	/**
	 * Update the objective function of the pricing problem with the new pricing information (modified costs).
	 * The modified costs are stored in the pricing problem.
	 */
	@Override
	protected void setObjective() {
		modifiedCosts=Arrays.copyOf(pricingProblem.dualCosts, pricingProblem.dualCosts.length);
		modifiedCost=pricingProblem.dualCost;
	}

	
	
	
}	

