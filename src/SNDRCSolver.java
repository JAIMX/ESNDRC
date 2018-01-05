import java.io.File;
import java.util.*;

import org.jorlib.demo.frameworks.columnGeneration.graphColoringBAP.cg.master.Master;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.io.SimpleBAPLogger;
import org.jorlib.frameworks.columnGeneration.io.SimpleDebugger;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;

import bap.BranchAndPrice;
import bap.branching.BranchOnQVarible;
import cg.Cycle;
import cg.ExactPricingProblemSolver;
import cg.SNDRCPricingProblem;
import model.SNDRC;

public class SNDRCSolver {
	public SNDRCSolver(SNDRC dataModel) {
		
		//Create the pricing problems
		List<SNDRCPricingProblem> pricingProblems=new ArrayList<>();
		for(int capacityType=0;capacityType<dataModel.numOfCapacity;capacityType++) {
			for(int originNode=0;originNode<dataModel.numNode;originNode++) {
				String name="capacity type: "+capacityType+" origin node: "+originNode;
				pricingProblems.add(new SNDRCPricingProblem(dataModel,name,capacityType,originNode));
			}
		}
		
		
		//Create the Master Problem
		Master master=new Master(dataModel,pricingProblems);
		
		//Define which solvers to use
		List<Class<?extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers=Collections.singletonList(ExactPricingProblemSolver.class);
		
		//Define one or more Branch creators
		List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators= Collections.singletonList(new BranchOnQVarible(dataModel, pricingProblems,0.5));
		
		//Create a Branch-and-Price instance
		BranchAndPrice bap=new BranchAndPrice(dataModel, master, pricingProblems, solvers, branchCreators,Double.MAX_VALUE);

		//OPTIONAL: Attach a debugger
		SimpleDebugger debugger=new SimpleDebugger(bap, true);

		//OPTIONAL: Attach a logger to the Branch-and-Price procedure.
		SimpleBAPLogger logger=new SimpleBAPLogger(bap, new File("./output/SNDRC.log"));

		//Solve the TSP problem through Branch-and-Price
		bap.runBranchAndPrice(System.currentTimeMillis()+8000000L);
		
		
		//Print solution:
		System.out.println("================ Solution ================");
		System.out.println("BAP terminated with objective : "+bap.getObjective());
		System.out.println("Total Number of iterations: "+bap.getTotalNrIterations());
		System.out.println("Total Number of processed nodes: "+bap.getNumberOfProcessedNodes());
		System.out.println("Total Time spent on master problems: "+bap.getMasterSolveTime()+" Total time spent on pricing problems: "+bap.getPricingSolveTime());
		if(bap.hasSolution()) {
			System.out.println("Solution is optimal: "+bap.isOptimal());
			System.out.println("Columns (only non-zero columns are returned):");
			List<Cycle> solution = bap.getSolution();
			for (Cycle column : solution)
				System.out.println(column);
		}
		
		bap.close();
		
		
	}

}
