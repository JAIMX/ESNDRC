package cg;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jorlib.frameworks.columnGeneration.colgenMain.ColGen;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.pricing.DefaultPricingProblemSolverFactory;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemBundle;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemManager;

import cg.master.SNDRCMasterData;
import model.SNDRC;

public class ColGenConditionalTerminate extends ColGen<SNDRC, Cycle, SNDRCPricingProblem> {
	private long terminateTimeLimit;

    /**
     * Create a new column generation instance
     * 
     * @param dataModel
     *            data model
     * @param master
     *            master problem
     * @param pricingProblems
     *            pricing problems
     * @param solvers
     *            pricing problem solvers
     * @param initSolution
     *            initial solution
     * @param cutoffValue
     *            cutoff Value. If the master is a minimization problem, the
     *            Colgen procedure is terminated if
     *            {@code ceil(boundOnMasterObjective) >= cutoffValue}. If the
     *            master is a maximization problem, the Colgen procedure is
     *            terminated if
     *            {@code floor(boundOnMasterObjective) <= cutoffValue}.
     * @param boundOnMasterObjective
     *            Bound on the best attainable objective value from the master
     *            problem. Assuming that the master is a minimization problem,
     *            the Colgen procedure is terminated if
     *            {@code ceil(boundOnMasterObjective) >= cutoffValue}.
     */
    public ColGenConditionalTerminate(SNDRC dataModel, AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> master,
            List<SNDRCPricingProblem> pricingProblems,
            List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
            List<Cycle> initSolution, int cutoffValue, double boundOnMasterObjective, long terminateTimeLimit) {
        super(dataModel, master, pricingProblems, solvers, initSolution, cutoffValue, boundOnMasterObjective);
        this.terminateTimeLimit=terminateTimeLimit;
    }
    
	/**
	 * Solve the Column Generation problem. First the master problem is solved. Next the pricing problems(s) is (are) solved. To solve the pricing problems, the pricing
	 * solvers are invoked one by one in a hierarchical fashion. First the first solver is invoked to solve the pricing problems. Any new columns generated are immediately returned.
	 * If it fails to find columns, the next solver is invoked and so on. If the pricing problem discovers new columns, they are added to the master problem and the method continues
	 * with the next column generation iteration.<br>
	 * If no new columns are found, the method checks for violated inequalities. If there are violated inequalities, they are added to the master problem and the method continues with the
	 * next column generation iteration.<br>
	 * The solve procedure terminates under any of the following conditions:
	 * <ol>
	 * <li>the solver could not identify new columns</li>
	 * <li>Time limit exceeded</li>
	 * <li>The bound on the best attainable solution to the master problem is worse than the cutoff value. Assuming that the master is a minimization problem, the Colgen procedure is terminated if {@code ceil(boundOnMasterObjective) >= cutoffValue}</li>
	 * <li>The solution to the master problem is provable optimal, i.e the bound on the best attainable solution to the master problem equals the solution of the master problem.</li>
	 * </ol>
	 * @param timeLimit Future point in time (ms) by which the procedure should be finished. Should be defined as: {@code System.currentTimeMilis()+<desired runtime>}
	 * @throws TimeLimitExceededException Exception is thrown when time limit is exceeded
	 */
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
			
			if(System.currentTimeMillis()>=terminateTimeLimit) {
				break;
			}
			
		}while(foundNewColumns || hasNewCuts);
		this.boundOnMasterObjective = (optimizationSenseMaster == OptimizationSense.MINIMIZE ? Math.max(this.boundOnMasterObjective, this.objectiveMasterProblem) : Math.min(this.boundOnMasterObjective, this.objectiveMasterProblem)); //When solved to optimality, the bound on the master problem objective equals the objective value.
		colGenSolveTime=System.currentTimeMillis()-colGenSolveTime;
		notifier.fireFinishCGEvent();
	}
    
    
}
