package bap;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import model.SNDRC;

public class BranchAndPrice extends AbstractBranchAndPrice<SNDRC, Cycle, SNDRCPricingProblem>{

	public BranchAndPrice(SNDRC modelData,Master master,List<SNDRCPricingProblem> pricingProblems,List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle,SNDRCPricingProblem>>> solvers,List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators,double objectiveInitialSolution) {
		super(modelData,master,pricingProblems,solvers,branchCreators,0,objectiveInitialSolution);
//		this.warmStart(objectiveInitialSolution,initialSolution);
	}
	
	
    /**
     * Generates an artificial solution. Columns in the artificial solution are of high cost such that they never end up in the final solution
     * if a feasible solution exists, since any feasible solution is assumed to be cheaper than the artificial solution. The artificial solution is used
     * to guarantee that the master problem has a feasible solution.
     *
     * @return artificial solution
     */
    @Override
    protected List<Cycle> generateInitialFeasibleSolution(BAPNode<SNDRC,Cycle> node) {
    	
    	List<Cycle> artificalVars=new ArrayList<Cycle>();
    	// for weak forcing constraints
    	for(int edgeIndex=0;edgeIndex<dataModel.numServiceArc;edgeIndex++) {
    		Set<Integer> set=new HashSet<>();
    		set.add(edgeIndex);
    		Cycle cycle=new Cycle(pricingProblems.get(0),true,"Artificial",set,100000000,0);
    		artificalVars.add(cycle);
    	}
    	
    	// for service edge branching constraints
    
        return artificalVars;
    }
    
    
    /**
     * Checks whether the given node is integer
     * @param node Node in the Branch-and-Price tree
     * @return true if the solution is an integer solution
     */
    @Override
    protected boolean isIntegerNode(BAPNode<SNDRC,Cycle> node) {
    	List<Cycle> result=node.getSolution();
    	
    	boolean out=true;
    	for(Cycle cycle:result) {
    		if(MathProgrammingUtil.isFractional(cycle.value)) {
    			out=false;
    			break;
    		}
    	}
    	
    	return out;
    }
	
}
