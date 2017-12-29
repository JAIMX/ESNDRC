package cg.master;

import java.util.*;

import javax.naming.TimeLimitExceededException;

import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import ilog.concert.*;
import ilog.cplex.*;
import model.SNDRC;
import model.SNDRC.Edge;

public final class Master extends AbstractMaster<SNDRC, Cycle, SNDRCPricingProblem, SNDRCMasterData> {

	private IloObjective obj;
	private IloNumVar[][] x;
	private IloNumVar[][] q;
	private IloRange[][] flowBalanceConstraints;
	private IloRange[] weakForcingConstraints;
	private IloRange[][] resourceBoundConstraints;

	public Master(SNDRC dataModel, SNDRCPricingProblem pricingProblem) {
		super(dataModel, pricingProblem, OptimizationSense.MINIMIZE);
		System.out.println("Master constructor. Columns: " + masterData.getNrColumns());
	}

	/**
	 * Builds the master model
	 * 
	 * @return Returns a MasterData object which is a data container for information
	 *         coming from the master problem
	 */
	@Override
	protected SNDRCMasterData buildModel() {
		IloCplex cplex = null;

		try {
			cplex = new IloCplex();

			cplex.setOut(null);
			cplex.setParam(IloCplex.IntParam.Threads, config.MAXTHREADS);

			// Define variables x
			x = new IloNumVar[dataModel.numDemand][dataModel.numArc];
			for (int p = 0; p < dataModel.numDemand; p++) {
				for (int arc = 0; arc < dataModel.numArc; arc++) {
					Edge edge = dataModel.edgeSet.get(arc);
					x[p][arc] = cplex.numVar(0, Double.MAX_VALUE,
							"x" + p + edge.u + "->" + edge.v + "t:" + edge.t1 + "->" + edge.t2);
				}
			}

			// Define the objective
			/**
			 * Here we assume the cost of edge AT is 0
			 */
			IloLinearNumExpr exprObj = cplex.linearNumExpr();
			for (int p = 0; p < dataModel.numDemand; p++) {
				for (int edgeIndex = 0; edgeIndex < dataModel.numServiceArc; edgeIndex++) {
					exprObj.addTerm(dataModel.beta * dataModel.edgeSet.get(edgeIndex).duration, x[p][edgeIndex]);
				}
			}

			obj = cplex.addMinimize(exprObj);

			// Define flowBalanceConstraints
			flowBalanceConstraints = new IloRange[dataModel.numDemand][dataModel.abstractNumNode];

			IloLinearNumExpr expr = cplex.linearNumExpr();
			for (int p = 0; p < dataModel.numDemand; p++) {
				for (int i = 0; i < dataModel.abstractNumNode; i++) {
					expr.clear();
					// edges which point from i
					for (int edgeIndex : dataModel.pointToEdgeSet.get(i)) {
						expr.addTerm(1, x[p][edgeIndex]);
					}

					// edges which point to i
					for (int edgeIndex : dataModel.pointFromEdgeSet.get(i)) {
						expr.addTerm(-1, x[p][edgeIndex]);
					}
					flowBalanceConstraints[p][i] = cplex.addEq(dataModel.b[p][i], expr);

				}
			}

			// Define weakForcingConstraints
			weakForcingConstraints = new IloRange[dataModel.numServiceArc];
			for (int arcIndex = 0; arcIndex < dataModel.numServiceArc; arcIndex++) {
				expr.clear();
				for (int p = 0; p < dataModel.numDemand; p++) {
					expr.addTerm(1, x[p][arcIndex]);
				}

				weakForcingConstraints[arcIndex] = cplex.addGe(0, expr);
			}

			// Define resourceBoundConstraints
			resourceBoundConstraints = new IloRange[dataModel.numOfCapacity][dataModel.numNode];
			q = new IloNumVar[dataModel.numOfCapacity][dataModel.numNode];

			for (int s = 0; s < dataModel.numOfCapacity; s++) {
				for (int o = 0; o < dataModel.numNode; o++) {
					q[s][o] = cplex.numVar(0, dataModel.vehicleLimit[s][o], "q" + s + "," + o);
				}
			}

			for (int s = 0; s < dataModel.numOfCapacity; s++) {
				for (int o = 0; o < dataModel.numNode; o++) {
					expr.clear();
					expr.addTerm(1, q[s][o]);
					resourceBoundConstraints[s][o] = cplex.addEq(dataModel.vehicleLimit[s][o], expr);
				}
			}

		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Map<SNDRCPricingProblem, OrderedBiMap<Cycle, IloNumVar>> varMap = new LinkedHashMap<>();
		for (SNDRCPricingProblem pricingProblem : pricingProblems)
			varMap.put(pricingProblem, new OrderedBiMap<>());

		// Create a new data object which will store information from the master. This
		// object automatically be passed to the CutHandler class.
		return new SNDRCMasterData(cplex, pricingProblems, varMap);

	}

	/**
	 * Solve the master problem
	 * 
	 * @param timeLimit
	 *            Future point in time by which the solve procedure must be
	 *            completed
	 * @return true if the master problem has been solved
	 * @throws TimeLimitExceededException
	 *             TimeLimitExceededException
	 */
	@Override
	protected boolean solveMasterProblem(long timeLimit){

		try {
			//Set time limit
			double timeRemaining=Math.max(1,(timeLimit-System.currentTimeMillis())/1000.0);
			
			masterData.cplex.setParam(IloCplex.DoubleParam.TiLim, timeRemaining); //set time limit in seconds
			
			//Potentially export the model
			if(config.EXPORT_MODEL) masterData.cplex.exportModel(config.EXPORT_MASTER_DIR+"master_"+this.getIterationCount()+".lp");
			
			//Solve the model
			if(!masterData.cplex.solve() || masterData.cplex.getStatus()!=IloCplex.Status.Optimal){
				if(masterData.cplex.getCplexStatus()==IloCplex.CplexStatus.AbortTimeLim) //Aborted due to time limit
					throw new TimeLimitExceededException();
				else
					throw new RuntimeException("Master problem solve failed! Status: "+masterData.cplex.getStatus());
			}else{
				masterData.objectiveValue=masterData.cplex.getObjValue();
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (TimeLimitExceededException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
		return true;
	}

	/**
	 * Adds a new column to the master problem
	 * @param column column to add
	 */
	@Override
	public void addColumn(Cycle column) {
		
		try {
			//Register column with objective
			IloColumn iloColumn=masterData.cplex.column(obj,column.cost);
			
			//weak forcing constraints
			for(int edgeIndex:column.edgeIndexSet) {
				iloColumn=iloColumn.and(masterData.cplex.column(weakForcingConstraints[edgeIndex],-dataModel.capacity[column.capacityIndex]));
			}
			
			//resource bound constraints
			iloColumn=iloColumn.and(masterData.cplex.column(resourceBoundConstraints[column.capacityIndex][column.startNode],1));
			
			
			//Create the variable and store it
			IloNumVar var=masterData.cplex.numVar(iloColumn, 0, Double.MAX_VALUE, "z_"+masterData.getNrColumnsForPricingProblem(column.associatedPricingProblem));
			masterData.cplex.add(var);
			masterData.addColumn(column,var);
			
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	/**
	 * Extracts information from the master problem which is required by the pricing problems, e.g. the reduced costs/dual values
	 * @param pricingProblem pricing problem
	 */
	@Override
	public void initializePricingProblem(SNDRCPricingProblem pricingProblem) {
		try {
			double[] modifiedCosts=new double[]; //Modified cost for every edge
			
		} catch (IloException e) {
			e.printStackTrace();
		}
	}
	

}
