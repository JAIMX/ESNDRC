package cg.master.cuts;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractCutGenerator;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

import cg.Cycle;
import cg.master.SNDRCMasterData;
import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex.UnknownObjectException;
import model.SNDRC;

public final class StrongInequalityGenerator extends AbstractCutGenerator<SNDRC, SNDRCMasterData> {

	public StrongInequalityGenerator(SNDRC modelData) {
		super(modelData, "strongIneqGenerator");
	}

	/**
	 * Generate inequalities using the data originating from the master problem
	 * 
	 * @return Returns true if a violated inequality has been found
	 */

	@Override
	public List<AbstractInequality> generateInqualities() {
		// Check for violated situations. When found, generate an inequality.
		for (int edgeIndex : masterData.edgeValueMap.keySet()) {
			Double edgeValue = masterData.edgeValueMap.get(edgeIndex);
			if (edgeValue < 1) {
				for (int commodity = 0; commodity < dataModel.numDemand; commodity++) {
					if(masterData.xValues.get(commodity).containsKey(edgeValue)) {
						double xValue = masterData.xValues.get(commodity).get(edgeIndex);
						double demandValue = dataModel.demandSet.get(commodity).volume;

						
						if (xValue > demandValue * edgeValue) {
							StrongInequality inequality = new StrongInequality(this, edgeIndex, commodity);
							this.addCut(inequality);
							return Collections.singletonList(inequality);
						}
					}

				}
			}
		}

		return Collections.emptyList();
	}

	/**
	 * If a violated inequality has been found add it to the master problem.
	 * 
	 * @param strongInequality
	 *            strong inequality
	 */
	private void addCut(StrongInequality strongInequality) {
		if (masterData.strongInequalities.containsKey(strongInequality))
			throw new RuntimeException(
					"Error, duplicate strong cut is being generated! This cut should already exist in the master problem: "
							+ strongInequality);

		// Create the inequality in cplex
		try {
			IloLinearNumExpr expr = masterData.cplex.linearNumExpr();
			expr.addTerm(1, masterData.x.get(strongInequality.commodity).get(strongInequality.edgeIndex));

			// register the columns with this constraints
			for (Cycle cycle : masterData.getColumns()) {
				if (cycle.edgeIndexSet.contains(strongInequality.edgeIndex)) {
					IloNumVar var = masterData.getVar(cycle);
					expr.addTerm(-dataModel.demandSet.get(strongInequality.commodity).volume, var);
				}
			}

			IloRange strongConstraint = masterData.cplex.addGe(0, expr, "strong cut");
			masterData.strongInequalities.put(strongInequality, strongConstraint);

		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/**
	 * Add a strong inequality to the master problem
	 * 
	 * @param cut
	 *            AbstractInequality
	 */
	@Override
	public void addCut(AbstractInequality cut) {
		if (!(cut instanceof StrongInequality))
			throw new IllegalArgumentException("This AbstractCutGenerator can ONLY add StrongInequalities");
		StrongInequality strongInequality = (StrongInequality) cut;
		this.addCut(strongInequality);
	}

	/**
	 * Retuns a list of inequalities that have been generated.
	 * 
	 * @return Retuns a list of inequalities that have been generated.
	 */
	@Override
	public List<AbstractInequality> getCuts() {
		return new ArrayList<>(masterData.strongInequalities.keySet());
	}

	/**
	 * Close the generator
	 */
	@Override
	public void close() {
	} // Nothing to do here

}
