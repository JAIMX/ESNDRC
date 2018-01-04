package bap;

import java.util.List;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import model.SNDRC;

public class BranchAndPrice extends AbstractBranchAndPrice<SNDRC, Cycle, SNDRCPricingProblem>{

	public BranchAndPrice(SNDRC modelData,Master master,List<SNDRCPricingProblem> pricingProblems,List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle,SNDRCPricingProblem>>> solvers,List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators,int objectiveInitialSolution,List<Cycle> initialSolution) {
		super(modelData,master,pricingProblems,solvers,branchCreators,0,objectiveInitialSolution);
		this.warmStart(objectiveInitialSolution,initialSolution);
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
        Matching matching1=new Matching("Artificial", true,	pricingProblems.get(0), incumbentSolution.get(0).edges, incumbentSolution.get(0).succ, objectiveIncumbentSolution);
        Matching matching2=new Matching("Artificial", true,	pricingProblems.get(1), incumbentSolution.get(1).edges, incumbentSolution.get(1).succ, objectiveIncumbentSolution);
        return Arrays.asList(matching1, matching2);
    }
	
}
