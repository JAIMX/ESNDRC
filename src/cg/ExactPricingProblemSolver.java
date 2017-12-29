package cg;

import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;

import ilog.concert.IloObjective;
import ilog.cplex.IloCplex;
import model.SNDRC;

public class ExactPricingProblemSolver extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem> {

	private IloCplex cplex; // Cplex instance
	private IloObjective obj;// Objective function

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
}
