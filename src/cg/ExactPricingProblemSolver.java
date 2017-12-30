package cg;

import java.util.*;

import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;

import ilog.concert.IloException;
import ilog.concert.IloObjective;
import ilog.cplex.IloCplex;
import model.SNDRC;

public class ExactPricingProblemSolver extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem> {

	private double[] modifiedCost;
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
		
		
		
		
		return newRoutes;	
	}
	
	
	/**
	 * Update the objective function of the pricing problem with the new pricing information (modified costs).
	 * The modified costs are stored in the pricing problem.
	 */
	@Override
	protected void setObjective() {
		modifiedCost=Arrays.copyOf(pricingProblem.dualCosts, pricingProblem.dualCosts.length);
	}

	
	
	
}	

