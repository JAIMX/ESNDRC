package cg.master.cuts;

import java.util.Collections;
import java.util.List;

import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractCutGenerator;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

import cg.master.SNDRCMasterData;
import ilog.concert.IloException;
import ilog.cplex.IloCplex.UnknownObjectException;
import model.SNDRC;

public final class StrongInequalityGenerator extends AbstractCutGenerator<SNDRC, SNDRCMasterData>{

	public StrongInequalityGenerator(SNDRC modelData) {
		super(modelData,"strongIneqGenerator");
	}
	
	/**
	 * Generate inequalities using the data originating from the master problem
	 * @return Returns true if a violated inequality has been found
	 */
	
	@Override
	public List<AbstractInequality> generateInqualities() {
		//Check for violated situations. When found, generate an inequality.
		for(int edgeIndex:masterData.edgeValueMap.keySet()) {
			Double edgeValue=masterData.edgeValueMap.get(edgeIndex);
			if(edgeValue<1) {
				for(int commodity=0;commodity<dataModel.numDemand;commodity++) {
					
					
					
					try {
						if(masterData.cplex.getValue(masterData.x.get(commodity).get(edgeIndex))>dataModel.demandSet.get(commodity).volume*edgeValue) {
							StrongInequality inequality=new StrongInequality(this, edgeIndex, commodity);
							this.addCut(inequality);
							return Collections.singletonList(inequality);
						}
					} catch (UnknownObjectException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (IloException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					
					
					
					
				}
			}
		}
		
		return Collections.emptyList();
	}
	
	
	
	/**
	 * If a violated inequality has been found add it to the master problem.
	 * @param strongInequality strong inequality
	 */
	private void addCut(StrongInequality strongInequality){
		if(masterData.strongInequalities.containsKey(strongInequality))
			throw new RuntimeException("Error, duplicate strong cut is being generated! This cut should already exist in the master problem: "+strongInequality);
		
		//Create the inequality in cplex
		
		
		
		
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
}
