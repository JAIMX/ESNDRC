package cg;

import java.util.List;

import org.jorlib.frameworks.columnGeneration.colgenMain.ColGen;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemManager;

import cg.master.SNDRCMasterData;
import model.SNDRC;

public class ColGenPlus extends ColGen{

	/**
	 * Create a new column generation instance
	 * @param dataModel data model
	 * @param master master problem
	 * @param pricingProblems pricing problems
	 * @param solvers pricing problem solvers
	 * @param initSolution initial solution
	 * @param cutoffValue cutoff Value. If the master is a minimization problem, the Colgen procedure is terminated if {@code ceil(boundOnMasterObjective) >= cutoffValue}. If the master is a maximization problem, the Colgen procedure is terminated if {@code floor(boundOnMasterObjective) <= cutoffValue}.
	 * @param boundOnMasterObjective Bound on the best attainable objective value from the master problem. Assuming that the master is a minimization problem, the Colgen procedure is terminated if {@code ceil(boundOnMasterObjective) >= cutoffValue}.
	 */
	public ColGenPlus(SNDRC dataModel, 
			AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> master,
			List<SNDRCPricingProblem> pricingProblems,
			List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
			List<Cycle> initSolution,
			int cutoffValue,
		  	double boundOnMasterObjective) {
		super(dataModel,master,pricingProblems,solvers,initSolution,cutoffValue,boundOnMasterObjective);
		
	}
	
	/**
	 * Create a new column generation instance
	 * @param dataModel data model
	 * @param master master problem
	 * @param pricingProblem pricing problem
	 * @param solvers pricing problem solvers
	 * @param initSolution initial solution
	 * @param cutoffValue cutoff Value. If the master is a minimization problem, the Colgen procedure is terminated if {@code ceil(boundOnMasterObjective) >= cutoffValue}. If the master is a maximization problem, the Colgen procedure is terminated if {@code floor(boundOnMasterObjective) <= cutoffValue}.
	 * @param boundOnMasterObjective Bound on the best attainable objective value from the master problem. Assuming that the master is a minimization problem, the Colgen procedure is terminated if {@code ceil(boundOnMasterObjective) >= cutoffValue}.
	 */
	public ColGenPlus(SNDRC dataModel, 
			AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> master,
			SNDRCPricingProblem pricingProblem,
			List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
			List<Cycle> initSolution,
			int cutoffValue,
		  	double boundOnMasterObjective) {
		super(dataModel,master,pricingProblem,solvers,initSolution,cutoffValue,boundOnMasterObjective);
		
	}
	
	/**
	 * Create a new column generation instance
	 * @param dataModel data model
	 * @param master master problem
	 * @param pricingProblem pricing problem
	 * @param solvers pricing problem solvers
	 * @param initSolution initial solution
	 * @param cutoffValue cutoff Value. If the master is a minimization problem, the Colgen procedure is terminated if {@code ceil(boundOnMasterObjective) >= cutoffValue}. If the master is a maximization problem, the Colgen procedure is terminated if {@code floor(boundOnMasterObjective) <= cutoffValue}.
	 * @param boundOnMasterObjective Bound on the best attainable objective value from the master problem. Assuming that the master is a minimization problem, the Colgen procedure is terminated if {@code ceil(boundOnMasterObjective) >= cutoffValue}.
	 */
	public ColGenPlus(SNDRC dataModel, 
			AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> master,
			List<SNDRCPricingProblem> pricingProblems,
			List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
			PricingProblemManager<SNDRC, Cycle, SNDRCPricingProblem> pricingProblemManager,
			List<Cycle> initSolution,
			int cutoffValue,
		  	double boundOnMasterObjective) {
		super(dataModel,master,pricingProblems,solvers,pricingProblemManager,initSolution,cutoffValue,boundOnMasterObjective);
		
	}
	
	
	
	
	@Override
	public void solve(long timeLimit) throws TimeLimitExceededException{
		//set time limit pricing problems
		pricingProblemManager.setTimeLimit(timeLimit);
		colGenSolveTime=System.currentTimeMillis();
		
		boolean foundNewColumns=false; //Identify whether the pricing problem generated new columns
		boolean hasNewCuts; //Identify whether the master problem violates any valid inequalities
		notifier.fireStartCGEvent();
		do{
			nrOfColGenIterations++;
			hasNewCuts=false;
			
			//Solve the master
			this.invokeMaster(timeLimit);

			//We can stop when the optimality gap is closed. We still need to check for violated inequalities though.
			if(Math.abs(objectiveMasterProblem - boundOnMasterObjective)<config.PRECISION){
				//Check whether there are inequalities. Otherwise potentially an infeasible integer solution (e.g. TSP solution with subtours) might be returned.
				if(config.CUTSENABLED){
					long time=System.currentTimeMillis();
					hasNewCuts=master.hasNewCuts();
					masterSolveTime+=(System.currentTimeMillis()-time); //Generating inequalities is considered part of the master problem
					if(hasNewCuts)
						continue;
					else
						break;
				}else
					break;
			}
			
			//Solve the pricing problem and possibly update the bound on the master problem objective
			List<Cycle> newColumns=this.invokePricingProblems(timeLimit); //List containing new columns generated by the pricing problem
			foundNewColumns=!newColumns.isEmpty();

			//Check whether the boundOnMasterObjective exceeds the cutoff value
			if(boundOnMasterExceedsCutoffValue())
				break;
			else if(System.currentTimeMillis() >= timeLimit){ //Check whether we are still within the timeLimit
				notifier.fireTimeLimitExceededEvent();
				throw new TimeLimitExceededException();
			}else if(config.CUTSENABLED && !foundNewColumns){ //Check for inequalities. This can only be done if the master problem hasn't changed (no columns can be added).
				long time=System.currentTimeMillis();
				hasNewCuts=master.hasNewCuts();
				masterSolveTime+=(System.currentTimeMillis()-time); //Generating inequalities is considered part of the master problem
			}
			
		}while(foundNewColumns || hasNewCuts);
		this.boundOnMasterObjective = (optimizationSenseMaster == OptimizationSense.MINIMIZE ? Math.max(this.boundOnMasterObjective, this.objectiveMasterProblem) : Math.min(this.boundOnMasterObjective, this.objectiveMasterProblem)); //When solved to optimality, the bound on the master problem objective equals the objective value.
		colGenSolveTime=System.currentTimeMillis()-colGenSolveTime;
		notifier.fireFinishCGEvent();
	}
	
	
	
	
	
	
	
	
	
}
