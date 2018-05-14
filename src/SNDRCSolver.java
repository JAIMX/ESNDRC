import java.io.File;
import java.io.IOException;
import java.util.*;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.bapNodeComparators.BFSbapNodeComparator;
import org.jorlib.frameworks.columnGeneration.io.SimpleDebugger;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.Configuration;

import bap.BranchAndPrice;
import bap.BranchAndPriceA;
import bap.BranchAndPriceA_M;
import bap.BranchAndPriceB;
import bap.BranchAndPriceB_M;
import bap.bapNodeComparators.NodeBoundbapNodeComparator;
import bap.bapNodeComparators.NodeBoundbapNodeComparatorMaxBound;
import bap.branching.BranchOnHoldingEdge;
import bap.branching.BranchOnLocalService;
import bap.branching.BranchOnLocalServiceForAllPricingProblems;
import bap.branching.BranchOnQVarible;
import bap.branching.BranchOnServiceEdge;
import bap.branching.BranchOnServiceEdgeForAllPricingProblems;
import bap.branching.BranchOnTimeService;
import bap.branching.BranchOnTimeServiceForAllPricingProblems;
import cg.Cycle;
import cg.ExactPricingProblemSolver;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import cg.master.SNDRCMasterData;
import cg.master.cuts.StrongInequalityGenerator;
import logger.BapLogger;
import logger.BapLoggerA;
import logger.BapLoggerA_M;
import logger.BapLoggerB;
import logger.BapLoggerB_M;
import model.SNDRC;
import model.SNDRC.Edge;

public class SNDRCSolver {
	SNDRC dataModel;
	boolean ifOptGetFromSubGraph;
	
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
		StrongInequalityGenerator cutGen=new StrongInequalityGenerator(dataModel,pricingProblems,0);
//		cutHandler.addCutGenerator(cutGen);
		
		//Create the Master Problem
		Master master=new Master(dataModel,pricingProblems,cutHandler,cutGen,false);
		
		//Define which solvers to use
		List<Class<?extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers=Collections.singletonList(ExactPricingProblemSolver.class);
		
		//Define one or more Branch creators
//		List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators=Arrays.asList(new BranchOnQVarible(dataModel, pricingProblems,0.5),new BranchOnServiceEdgeForAllPricingProblems(dataModel, pricingProblems, 0.5),new BranchOnServiceEdge(dataModel, pricingProblems, 0.5));
		List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators=Arrays.asList( new BranchOnLocalServiceForAllPricingProblems(dataModel, pricingProblems, 0.5),new BranchOnLocalService(dataModel, pricingProblems, 0.5),new BranchOnServiceEdge(dataModel, pricingProblems, 0.5));
//		List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators=Arrays.asList(new BranchOnTimeServiceForAllPricingProblems(dataModel, pricingProblems, 0.5, 5),new BranchOnServiceEdge(dataModel, pricingProblems, 0.5));
	      
		//Create a Branch-and-Price instance
//		BranchAndPriceA bap=new BranchAndPriceA(dataModel, master, pricingProblems, solvers, branchCreators,Double.MAX_VALUE,0.6,0.2,0.1,10,0.0001,3,false);
//		BranchAndPriceB bap=new BranchAndPriceB(dataModel, master, pricingProblems, solvers, branchCreators,Double.MAX_VALUE,0.65,0.2,0.1,1,0.001,3,0.1,true);
//		BranchAndPriceB_M bap=new BranchAndPriceB_M(dataModel, master, pricingProblems, solvers, branchCreators,Double.MAX_VALUE,0.65,0.3,0.1,1,-0.001,1,0,true,false);
//		BranchAndPriceA_M bap=new BranchAndPriceA_M(dataModel, master, pricingProblems, solvers, branchCreators,Double.MAX_VALUE,0.6,0.3,0.1,10,0.001,10,0.1,false,true);
		BranchAndPriceA_M bap=new BranchAndPriceA_M(dataModel, master, pricingProblems, solvers, branchCreators,5000000,0.6,0.3,0.1,10,0.001,10,0.1,false,true);
//		bap.setNodeOrdering(new BFSbapNodeComparator());
//		bap.setNodeOrdering(new NodeBoundbapNodeComparatorMaxBound());
//		bap.setNodeOrdering(new NodeBoundbapNodeComparator());
		
		//OPTIONAL: Attach a debugger
//		SimpleDebugger debugger=new SimpleDebugger(bap, true);

		//OPTIONAL: Attach a logger to the Branch-and-Price procedure.
//		BapLoggerA logger=new BapLoggerA(bap, new File("./output/BAPlogger.log"));
//		BapLoggerB logger=new BapLoggerB(bap, new File("./output/BAPlogger.log"));
//		BapLoggerB_M logger=new BapLoggerB_M(bap, new File("./output/BAPlogger.log"));
		BapLoggerA_M logger=new BapLoggerA_M(bap, new File("./output/BAPlogger.log"));

		//Solve the TSP problem through Branch-and-Price
//		bap.runBranchAndPrice(System.currentTimeMillis()+18000000L);    //5 hours
		bap.runBranchAndPrice(System.currentTimeMillis()+36000000L);    //10 hours
		
		this.ifOptGetFromSubGraph=bap.GetIfOptGetFromSubGraph();
		
		
		
		
		//Print solution:
		System.out.println("================ Solution ================");
		System.out.println("BAP terminated with objective : "+bap.getObjective());
		System.out.println("Total Number of iterations: "+bap.getTotalNrIterations());
		System.out.println("Total Number of processed nodes: "+bap.getNumberOfProcessedNodes());
		System.out.println("Total Time spent on master problems: "+bap.getMasterSolveTime()+" Total time spent on pricing problems: "+bap.getPricingSolveTime());
		System.out.println("Best bound : "+bap.getBound());
		
		//Count the number of key edges
		Set<Integer> keyServiceEdgeIndexSet=new TreeSet<>();
		if(bap.hasSolution()){
		    List<Cycle> solution = bap.getSolution();
            for (Cycle column : solution) {
                Set<Integer> tempEdgeSet=column.edgeIndexSet;
                for(int edgeIndex:tempEdgeSet){
                    if(((ifOptGetFromSubGraph)&&(dataModel.subEdgeSet.get(edgeIndex).edgeType==0))||((!ifOptGetFromSubGraph)&&(dataModel.edgeSet.get(edgeIndex).edgeType==0))){
                        if(!keyServiceEdgeIndexSet.contains(edgeIndex)){
                            keyServiceEdgeIndexSet.add(edgeIndex);
                        }
                    }
                }
            }
		}
		
		System.out.println(keyServiceEdgeIndexSet.toString());
		System.out.println("The number of service edges used= "+keyServiceEdgeIndexSet.size());
		System.out.println();
		
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
		
//		List<Map<Integer,Double>> optXValues=bap.GetOptXValues();
//		//output x variables
//		for(int demand=0;demand<dataModel.numDemand;demand++){
//		    for(int edgeIndex:optXValues.get(demand).keySet()){
//		        if(optXValues.get(demand).get(edgeIndex)>0.01){
//		            Edge edge=dataModel.edgeSet.get(edgeIndex);
//		            System.out.println("x[" + demand + "]:" + edge.start + "->" + edge.end + "= "+optXValues.get(demand).get(edgeIndex));
//		        }
//		    }
//		    System.out.println();
//		}
		
		bap.close();
		cutHandler.close();
		
		
		
	}
	
	public String out(Cycle column) {
		
		Queue<Edge> path=new PriorityQueue<>();
		
		if(!ifOptGetFromSubGraph){
			for(int edgeIndex:column.edgeIndexSet) {
				path.add(dataModel.edgeSet.get(edgeIndex));
			}
		}else{
			for(int edgeIndex:column.edgeIndexSet) {
				path.add(dataModel.subEdgeSet.get(edgeIndex));
			}
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
		
		SNDRC sndrc;
		for(String arg:args) {
			long time0=System.currentTimeMillis();
			sndrc=new SNDRC(arg);
//			sndrc.Output();
			Properties properties=new Properties();
//			properties.setProperty("EXPORT_MODEL", "True");
//			properties.setProperty("MAXTHREADS", "10");
//			properties.setProperty("PRECISION", "0.001");
//			properties.setProperty("CUTSENABLED", "false");
			Configuration.readFromFile(properties);
			
			new SNDRCSolver(sndrc);
			
			long time1=System.currentTimeMillis();
			System.out.println("Total time= "+(time1-time0));
			
		}
		

		
	}

}
