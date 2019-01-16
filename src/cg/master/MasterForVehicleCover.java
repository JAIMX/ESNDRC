package cg.master;

import java.util.*;

import javax.naming.TimeLimitExceededException;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.AbstractMaster;
import org.jorlib.frameworks.columnGeneration.master.OptimizationSense;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import bap.branching.branchingDecisions.RoundHoldingEdge;
import bap.branching.branchingDecisions.RoundLocalService;
import bap.branching.branchingDecisions.RoundLocalServiceForAllPricingProblems;
import bap.branching.branchingDecisions.RoundQ;
import bap.branching.branchingDecisions.RoundServiceEdge;
import bap.branching.branchingDecisions.RoundServiceEdgeForAllPricingProblems;
import bap.branching.branchingDecisions.RoundTimeService;
import bap.branching.branchingDecisions.RoundTimeServiceForAllPricingProblems;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import cg.master.cuts.StrongInequality;
import cg.master.cuts.StrongInequalityGenerator;
import ilog.concert.*;
import ilog.cplex.*;
import ilog.cplex.IloCplex.UnknownObjectException;
import model.SNDRC;
import model.SNDRC.Edge;

public final class MasterForVehicleCover extends Master {
	
    public MasterForVehicleCover(SNDRC dataModel, List<SNDRCPricingProblem> pricingProblems,
			CutHandler<SNDRC, SNDRCMasterData> cutHandler, StrongInequalityGenerator cutGen,
			boolean ifContainsCutAtRootNode) {
		super(dataModel, pricingProblems, cutHandler, cutGen, ifContainsCutAtRootNode);
	}
    
	@Override
    public Map<SNDRCPricingProblem, Integer> getQVarLimit() {
        return masterData.qVariableLimit;
    }
	
    @Override 
    public List<Map<Integer, Double>> getXValues() throws UnknownObjectException, IloException {
        return masterData.xValues;
    }
    
}