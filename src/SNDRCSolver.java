import java.io.File;
import java.io.IOException;
import java.util.*;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.bapNodeComparators.BFSbapNodeComparator;
import org.jorlib.frameworks.columnGeneration.io.SimpleDebugger;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.Configuration;

import bap.BranchAndPrice;
import bap.bapNodeComparators.NodeBoundbapNodeComparator;
import bap.branching.BranchOnLocalService;
import bap.branching.BranchOnQVarible;
import bap.branching.BranchOnServiceEdge;
import bap.branching.BranchOnServiceEdgeForAllPricingProblems;
import cg.Cycle;
import cg.ExactPricingProblemSolver;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import cg.master.SNDRCMasterData;
import cg.master.cuts.StrongInequalityGenerator;
import logger.BapLogger;
import model.SNDRC;
import model.SNDRC.Edge;

public class SNDRCSolver {
	SNDRC dataModel;
	public SNDRCSolver(SNDRC dataModel) {
		this.dataModel=dataModel;
		
		//Create the pricing problems
		List<SNDRCPricingProblem> pricingProblems=new LinkedList<SNDRCPricingProblem>();
		for(int capacityType=0;capacityType<dataModel.numOfCapacity;capacityType++) {
			for(int originNode=0;originNode<dataModel.numNode;originNode++) {
				String name="capacity type: "+capacityType+" origin node: "+originNode;
				SNDRCPricingProblem pricingProblem=new SNDRCPricingProblem(dataModel,name,capacityType,originNode);
				pricingProblems.add(pricingProblem);
			}
		}
		
		
		//Create a cutHandler
		CutHandler<SNDRC, SNDRCMasterData> cutHandler=new CutHandler<>();
		StrongInequalityGenerator cutGen=new StrongInequalityGenerator(dataModel,pricingProblems,1);
//		cutHandler.addCutGenerator(cutGen);
		
		//Create the Master Problem
		Master master=new Master(dataModel,pricingProblems,cutHandler);
		
		//Define which solvers to use
		List<Class<?extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers=Collections.singletonList(ExactPricingProblemSolver.class);
		
		//Define one or more Branch creators
//		List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators=new ArrayList<AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>>();
//		List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators=Arrays.asList(new BranchOnQVarible(dataModel, pricingProblems,0.5),new BranchOnServiceEdgeForAllPricingProblems(dataModel, pricingProblems, 0.5),new BranchOnServiceEdge(dataModel, pricingProblems, 0.5));
		List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators=Arrays.asList(new BranchOnLocalService(dataModel, pricingProblems, 0.5),new BranchOnQVarible(dataModel, pricingProblems,0.5),new BranchOnServiceEdge(dataModel, pricingProblems, 0.5));
//	    List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators=Arrays.asList(new BranchOnQVarible(dataModel, pricingProblems,0.5),new BranchOnServiceEdge(dataModel, pricingProblems, 0.5));
	      
		//Create a Branch-and-Price instance
//		BranchAndPrice bap=new BranchAndPrice(dataModel, master, pricingProblems, solvers, branchCreators,Double.MAX_VALUE,0.5,0.3,0.1);
		BranchAndPrice bap=new BranchAndPrice(dataModel, master, pricingProblems, solvers, branchCreators,Double.MAX_VALUE,0.6,0.3,0.1);
//		bap.setNodeOrdering(new BFSbapNodeComparator());
		bap.setNodeOrdering(new NodeBoundbapNodeComparator());
		
		//OPTIONAL: Attach a debugger
//		SimpleDebugger debugger=new SimpleDebugger(bap, true);

		//OPTIONAL: Attach a logger to the Branch-and-Price procedure.
//		SimpleBAPLogger logger=new SimpleBAPLogger(bap, new File("./output/SNDRC.log"));
		BapLogger logger=new BapLogger(bap, new File("./output/BAPlogger.log"));

		//Solve the TSP problem through Branch-and-Price
		bap.runBranchAndPrice(System.currentTimeMillis()+36000000L);
		
		
		//Print solution:
		System.out.println("================ Solution ================");
		System.out.println("BAP terminated with objective : "+bap.getObjective());
		System.out.println("Total Number of iterations: "+bap.getTotalNrIterations());
		System.out.println("Total Number of processed nodes: "+bap.getNumberOfProcessedNodes());
		System.out.println("Total Time spent on master problems: "+bap.getMasterSolveTime()+" Total time spent on pricing problems: "+bap.getPricingSolveTime());
		System.out.println("Best bound : "+bap.getBound());
		if(bap.hasSolution()) {
			System.out.println("Solution is optimal: "+bap.isOptimal());
			System.out.println("Columns (only non-zero columns are returned):");
			List<Cycle> solution = bap.getSolution();
			for (Cycle column : solution) {
				System.out.println(column);
				System.out.println(out(column)+":"+bap.GetOptSolutionValueMap().get(column));
				System.out.println();
			}
				
		}
//		System.out.println(bap.getBound());
		
		List<Map<Integer,Double>> optXValues=bap.GetOptXValues();
		//output x variables
		for(int demand=0;demand<dataModel.numDemand;demand++){
		    for(int edgeIndex:optXValues.get(demand).keySet()){
		        if(optXValues.get(demand).get(edgeIndex)>0.01){
		            Edge edge=dataModel.edgeSet.get(edgeIndex);
		            System.out.println("x[" + demand + "]:" + edge.start + "->" + edge.end + "= "+optXValues.get(demand).get(edgeIndex));
		        }
		    }
		    System.out.println();
		}
		
		bap.close();
		cutHandler.close();
		
		
		
	}
	
	public String out(Cycle column) {
		
		Queue<Edge> path=new PriorityQueue<>();
		for(int edgeIndex:column.edgeIndexSet) {
			path.add(dataModel.edgeSet.get(edgeIndex));
		}
		
		StringBuilder pathRecord=new StringBuilder();
		
//		int count=0;
//		for(Edge edge:path) {
//			count++;
//			pathRecord.append(edge.start);
//			
//			if(count!=column.edgeIndexSet.size()) {
//				pathRecord.append("->");
//			}
//
//		}
		
		Edge edge=null;
		int size=path.size();
		for(int i=0;i<size;i++) {

		    edge=path.poll();
			pathRecord.append('(');
			pathRecord.append(edge.u);
			pathRecord.append(',');
			pathRecord.append(edge.t1);
			pathRecord.append(')');
			
			
			pathRecord.append("->");

		}
		
		pathRecord.append('(');
		pathRecord.append(edge.v);
		pathRecord.append(',');
		pathRecord.append(edge.t2);
		pathRecord.append(')');
		
		
		return pathRecord.toString();
		
	}
	
	public static void main(String[] args) throws IOException {

//		SNDRC sndrc=new SNDRC("./data/test2_5_8.txt");
//		SNDRC sndrc=new SNDRC("./data/test4_5_8_5.txt");
//		SNDRC sndrc=new SNDRC("./data/test3_5_15.txt");
//		SNDRC sndrc=new SNDRC("./data/test1_5_15_5.txt");
//		SNDRC sndrc=new SNDRC("./data/change_fixedCost4.txt");
//		SNDRC sndrc=new SNDRC("./data/change_fixedCost.txt");
//		SNDRC sndrc=new SNDRC("./data/test5_25_10.txt");
//		SNDRC sndrc=new SNDRC("./data/testset/test1_5_10_15_20.txt");
		
		SNDRC sndrc;
		for(String arg:args) {
			long time0=System.currentTimeMillis();
			sndrc=new SNDRC(arg);
			sndrc.Output();
			Properties properties=new Properties();
//			properties.setProperty("EXPORT_MODEL", "True");
//			properties.setProperty("MAXTHREADS", "10");
//			properties.setProperty("PRECISION", "0.001");
			Configuration.readFromFile(properties);
			
			new SNDRCSolver(sndrc);
			
			long time1=System.currentTimeMillis();
			System.out.println("Total time= "+(time1-time0));
			
		}
		

		
	}

}
