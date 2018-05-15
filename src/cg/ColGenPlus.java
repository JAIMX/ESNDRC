package cg;

import java.util.List;

import org.jorlib.frameworks.columnGeneration.colgenMain.ColGen;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.pricing.PricingProblemManager;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import cg.master.SNDRCMasterData;
import model.SNDRC;

public class ColGenPlus extends ColGen<SNDRC, Cycle, SNDRCPricingProblem> {
    private double gapTolerance;

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
    public ColGenPlus(SNDRC dataModel, AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> master,
            List<SNDRCPricingProblem> pricingProblems,
            List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
            List<Cycle> initSolution, int cutoffValue, double boundOnMasterObjective, double gapTolerance) {
        super(dataModel, master, pricingProblems, solvers, initSolution, cutoffValue, boundOnMasterObjective);
        this.gapTolerance = gapTolerance;

    }

    /**
     * Create a new column generation instance
     * 
     * @param dataModel
     *            data model
     * @param master
     *            master problem
     * @param pricingProblem
     *            pricing problem
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
    public ColGenPlus(SNDRC dataModel, AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> master,
            SNDRCPricingProblem pricingProblem,
            List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
            List<Cycle> initSolution, int cutoffValue, double boundOnMasterObjective, double gapTolerance) {
        super(dataModel, master, pricingProblem, solvers, initSolution, cutoffValue, boundOnMasterObjective);
        this.gapTolerance = gapTolerance;

    }

    /**
     * Create a new column generation instance
     * 
     * @param dataModel
     *            data model
     * @param master
     *            master problem
     * @param pricingProblem
     *            pricing problem
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
    public ColGenPlus(SNDRC dataModel, AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> master,
            List<SNDRCPricingProblem> pricingProblems,
            List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers,
            PricingProblemManager<SNDRC, Cycle, SNDRCPricingProblem> pricingProblemManager, List<Cycle> initSolution,
            int cutoffValue, double boundOnMasterObjective, double gapTolerance) {
        super(dataModel, master, pricingProblems, solvers, pricingProblemManager, initSolution, cutoffValue,
                boundOnMasterObjective);
        this.gapTolerance = gapTolerance;
    }

    @Override
    public void solve(long timeLimit) throws TimeLimitExceededException {
        // set time limit pricing problems
        pricingProblemManager.setTimeLimit(timeLimit);
        colGenSolveTime = System.currentTimeMillis();

        boolean foundNewColumns = false; // Identify whether the pricing problem
                                         // generated new columns
        boolean hasNewCuts; // Identify whether the master problem violates any
                            // valid inequalities
        notifier.fireStartCGEvent();

        double gap0 = -1;
        double gap1 = -1;
        do {
            nrOfColGenIterations++;
            hasNewCuts = false;

            // Solve the master
            this.invokeMaster(timeLimit);

            // check if the solution is integral
            boolean isIntegral = true;
            boolean isFeasible = true;
            List<Cycle> result = this.getSolution();
            for (Cycle cycle : result) {
                if (MathProgrammingUtil.isFractional(cycle.value)) {
                    isIntegral = false;
                    break;
                }
            }
            for (Cycle cycle : result) {
                if(cycle.isArtificialColumn){
                    isFeasible=false;
                    break;
                }
            }

            // update gap0 and gap1
            gap0 = gap1;
            if (boundOnMasterObjective <= config.PRECISION) {
                gap1 = objectiveMasterProblem - 1;
            } else {
                gap1 = (objectiveMasterProblem - boundOnMasterObjective) / boundOnMasterObjective;
            }

            // after the update of objectiveMasterProblem £¬we check if stop
            // column generation process

            // Case1: boundOnMasterObjective exceeds cutoff value
            if (boundOnMasterExceedsCutoffValue()) {
                break;
            }

            // Case2: objectiveMasterProblem is smaller than cutoff value(here
            // we should note infeasible case,set the objective value of
            // artificial variables big enough compared to cutoff value)
            if (objectiveMasterProblem < cutoffValue - config.PRECISION && gap1 < 0.3 && !isIntegral&&isFeasible) {
                if (config.CUTSENABLED) {
                    long time = System.currentTimeMillis();
                    hasNewCuts = master.hasNewCuts();
                    masterSolveTime += (System.currentTimeMillis() - time); // Generating
                                                                            // inequalities
                                                                            // is
                                                                            // considered
                                                                            // part
                                                                            // of
                                                                            // the
                                                                            // master
                                                                            // problem
                    if (hasNewCuts) {
                        gap0 = -1;
                        gap1 = -1;
                        continue;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }

            // Case3:objectiveMasterProblem is over cutoff value and
            // boundOnMasterObjective below cutoff value
            if (objectiveMasterProblem > cutoffValue + config.PRECISION && !boundOnMasterExceedsCutoffValue()) {
                if (gap0 > 0 && gap0 - gap1 < gapTolerance - config.PRECISION && gap1 < 0.3 && !isIntegral&&isFeasible) {

                    if (config.CUTSENABLED) {
                        long time = System.currentTimeMillis();
                        hasNewCuts = master.hasNewCuts();
                        masterSolveTime += (System.currentTimeMillis() - time); // Generating
                                                                                // inequalities
                                                                                // is
                                                                                // considered
                                                                                // part
                                                                                // of
                                                                                // the
                                                                                // master
                                                                                // problem
                        if (hasNewCuts) {
                            gap0 = -1;
                            gap1 = -1;
                            continue;
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }

                }
            }

            // Case4:solve to optimal
            if (Math.abs(objectiveMasterProblem - boundOnMasterObjective) < config.PRECISION) {
                // Check whether there are inequalities. Otherwise potentially
                // an infeasible integer solution (e.g. TSP solution with
                // subtours) might be returned.
                if (config.CUTSENABLED) {
                    long time = System.currentTimeMillis();
                    hasNewCuts = master.hasNewCuts();
                    masterSolveTime += (System.currentTimeMillis() - time); // Generating
                                                                            // inequalities
                                                                            // is
                                                                            // considered
                                                                            // part
                                                                            // of
                                                                            // the
                                                                            // master
                                                                            // problem
                    if (hasNewCuts)
                        continue;
                    else
                        break;
                } else
                    break;
            }

            // Solve the pricing problem and possibly update the bound on the
            // master problem objective
            List<Cycle> newColumns = this.invokePricingProblems(timeLimit); // List
                                                                            // containing
                                                                            // new
                                                                            // columns
                                                                            // generated
                                                                            // by
                                                                            // the
                                                                            // pricing
                                                                            // problem
            foundNewColumns = !newColumns.isEmpty();

            if (System.currentTimeMillis() >= timeLimit) { // Check whether we
                                                           // are still within
                                                           // the timeLimit
                notifier.fireTimeLimitExceededEvent();
                throw new TimeLimitExceededException();
            }
            // else if(config.CUTSENABLED && !foundNewColumns&&!hasNewCuts) {
            // long time=System.currentTimeMillis();
            // hasNewCuts=master.hasNewCuts();
            // masterSolveTime+=(System.currentTimeMillis()-time);
            // }

        } while (foundNewColumns || hasNewCuts);
        // this.boundOnMasterObjective = (optimizationSenseMaster ==
        // OptimizationSense.MINIMIZE ? Math.max(this.boundOnMasterObjective,
        // this.objectiveMasterProblem) : Math.min(this.boundOnMasterObjective,
        // this.objectiveMasterProblem)); //When solved to optimality, the bound
        // on the master problem objective equals the objective value.
        colGenSolveTime = System.currentTimeMillis() - colGenSolveTime;
        notifier.fireFinishCGEvent();
    }

    @Override
    protected double calculateBoundOnMasterObjective(
            Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>> solver) {
        double[] bounds = pricingProblemManager.getBoundsOnPricingProblems(solver);
        double bound = master.getBoundComponent();
        for (double ele : bounds) {
            bound += ele;
        }

        bound = Math.max(bound, this.boundOnMasterObjective);
        if (bound > master.getBoundComponent() + config.PRECISION) {
            bound = master.getBoundComponent();
        }
        return bound;
    }

}
