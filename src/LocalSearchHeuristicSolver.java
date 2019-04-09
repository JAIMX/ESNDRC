import java.awt.Dialog.ModalExclusionType;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Properties;
import java.util.Queue;
import java.util.Set;


import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.colgenMain.ColGen;
import org.jorlib.frameworks.columnGeneration.io.SimpleDebugger;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.Configuration;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;


import bap.BranchAndPriceA;
import bap.bapNodeComparators.NodeBoundbapNodeComparator;
import bap.branching.BranchOnLocalService;
import bap.branching.BranchOnLocalServiceForAllPricingProblems;
import bap.branching.BranchOnServiceEdge;
import bap.primalHeuristic.ColumnGenerationBasedHeuristic;
import cg.Cycle;
import cg.ExactPricingProblemSolver;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import cg.master.MasterForVehicleCover;
import cg.master.SNDRCMasterData;
import cg.master.cuts.StrongInequalityGenerator;
import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import logger.BapLoggerA;
import model.SNDRC;
import model.SNDRC.Demand;
import model.SNDRC.Edge;
import model.SNDRC.Path;
import model.SNDRC.Service;

public class LocalSearchHeuristicSolver {
	private SNDRC modelData;
	private int r;    //r shortest path in the initialization phase
	private double balanceValue1,balanceValue2;
	private int timeZone;
	private List<Set<Integer>> edgesForXRecord;
	private int tabuListLengthLimit;//parameters for tabu list
//	private List<Integer> tabuCommodityList,tabuEdgeIndexList;
//	private boolean[][] tabuEdge;
    private ArrayDeque<Integer> tabuNodeList;
	public LocalSearchHeuristicSolver(String filename, int r,double balanceValue1,double balanceValue2,int timeZone,int tabuListLengthLimit) throws IOException {
		modelData = new SNDRC(filename);
		this.r = r;
		this.balanceValue1=balanceValue1;
		this.balanceValue2=balanceValue2;
		this.timeZone=timeZone;
		this.tabuListLengthLimit=tabuListLengthLimit;
	}

	class FeasibleSolution {
		List<Map<Integer, Double>> optXValues;
		List<Cycle> cycleValues;
		double fixCost,variableCost,flowCost,totalCost;

		public FeasibleSolution(List<Map<Integer, Double>> optXValues, List<Cycle> cycleValues) {
			this.optXValues = optXValues;
			this.cycleValues = cycleValues;
			
			fixCost=0;
			variableCost=0;
			flowCost=0;
			
			for(Map<Integer,Double> map:optXValues){
				for(int edgeIndex:map.keySet()){
					Edge edge=modelData.edgeSet.get(edgeIndex);
					if(edge.edgeType==0){
						flowCost+=modelData.beta*edge.duration*map.get(edgeIndex);
					}
				}
			}
			
			for(Cycle cycle:cycleValues){
				fixCost+=modelData.fixedCost[cycle.associatedPricingProblem.originNodeO][cycle.associatedPricingProblem.capacityTypeS]*cycle.value;
				for(int edgeIndex:cycle.edgeIndexSet){
					Edge edge=modelData.edgeSet.get(edgeIndex);
					if(edge.edgeType==0){
						//dataModel.alpha * totalLength / (dataModel.speed * dataModel.drivingTimePerDay)
						variableCost+=modelData.alpha/(modelData.speed*modelData.drivingTimePerDay)*edge.duration*cycle.value;
					}
				}
			}
			totalCost=fixCost+variableCost+flowCost;
			
		}
	}
	
	class Node implements Comparable<Node>{
		int nodeIndex,index;
		List<Integer> pathRecord;
		
		public Node(int nodeIndex,List<Integer> pathRecord,int index) {
			this.nodeIndex=nodeIndex;
			this.pathRecord=new ArrayList<>(pathRecord);
			this.index=index;
		}
		
		public int compareTo(Node that){
			return this.index-that.index;
		}
		
	}
	
	class CommoditySubPath{
		int commodityIndex,startNodeIndex,endNodeIndex;
		double amount,flowCost;
		List<Integer> pathEdgeIndexList;
//		int tabuEdgeIndex;
		int tabuNodeIndex;
		
		public CommoditySubPath(int commodityIndex,double amount,List<Integer> list,int tabuNodeIndex){
			this.commodityIndex=commodityIndex;
			this.amount=amount;
			this.pathEdgeIndexList=new ArrayList<>(list);
			this.tabuNodeIndex=tabuNodeIndex;
			
			Edge edge=modelData.edgeSet.get(pathEdgeIndexList.get(0));
			startNodeIndex=edge.start;
			edge=modelData.edgeSet.get(pathEdgeIndexList.get(pathEdgeIndexList.size()-1));
			endNodeIndex=edge.end;
			
			//calculate flowCost(para=dataModel.beta * dataModel.edgeSet.get(edgeIndex).duration)
			flowCost=0;
			for(int edgeIndex:pathEdgeIndexList){
				edge=modelData.edgeSet.get(edgeIndex);
				if(edge.edgeType==0){
					flowCost+=modelData.beta*edge.duration;
				}
			}
			flowCost=flowCost*amount;
			
		}
		
		
		public String toString(){
			String str="";
			for(int edgeIndex:this.pathEdgeIndexList){
				str+=modelData.edgeSet.get(edgeIndex).toString()+" ";
			}
			str+=" commodityIndex="+commodityIndex;
			str+=" amount="+amount;
			
			return str;
		}
		
		
		
	}

	public List<FeasibleSolution> InitializationCGHeuristic() throws TimeLimitExceededException, IloException{
		// We first calculate how many shortest paths service edge (i,j) belong to
		// commodity k's r shortest paths
		int[][] serviceEdgeImportance = new int[modelData.numServiceArc][modelData.numDemand];
		int[] serviceEdgeImportanceSum = new int[modelData.numServiceArc];
		ArrayList<ArrayList<Path>> rPathSet = new ArrayList<>();
		for (int k = 0; k < modelData.numDemand; k++) {
			rPathSet.add(modelData.findRShortestPath(r, k));
		}

		for (int k = 0; k < modelData.numDemand; k++) {
			Demand demand = modelData.demandSet.get(k);
			ArrayList<Path> pathList = rPathSet.get(k);
			for (Path path : pathList) {

				for (int i = 0; i < path.serviceIndexList.size(); i++) {

					int serviceIndex = path.serviceIndexList.get(i);
					Set<Integer> timeSet = new HashSet<>();
					int time = path.timeList.get(i);

					for (int t = 0; t <= demand.duration - path.totalDuration; t++) {
						// timeAvalible+time+t
						int cTime = demand.timeAvailable + time + t;
						cTime = cTime % modelData.timePeriod;
						timeSet.add(cTime);
					}

					for (int edgeIndex = 0; edgeIndex < modelData.numServiceArc; edgeIndex++) {
						Edge edge = modelData.edgeSet.get(edgeIndex);
						if (edge.serviceIndex == serviceIndex && timeSet.contains(edge.t1)) {
							serviceEdgeImportance[edgeIndex][k]++;
						}
					}
				}
			}
		}

		for (int edgeIndex = 0; edgeIndex < modelData.numServiceArc; edgeIndex++) {
			for (int k = 0; k < modelData.numDemand; k++) {
				serviceEdgeImportanceSum[edgeIndex] += serviceEdgeImportance[edgeIndex][k];
			}
		}

		// We use dp to calculate subEdgeSet for each commodity
		List<Set<Integer>> edgesForX = new ArrayList<>();
		for (int k = 0; k < modelData.numDemand; k++) {
			Set<Integer> set = new HashSet<Integer>();
			edgesForX.add(set);
		}

		for (int k = 0; k < modelData.numDemand; k++) {
			// System.out.println();
			// System.out.println("k= "+k);
			Set<Integer> edgeSet = edgesForX.get(k);
			Demand demand = modelData.demandSet.get(k);

			for (Path path : rPathSet.get(k)) {
				int numOfService = path.serviceIndexList.size();
				int totalTime = demand.duration;

				int[][] f = new int[numOfService + 1][totalTime + 1];
				int[][] record = new int[numOfService + 1][totalTime + 1];
				for (int i = 0; i <= numOfService; i++) {
					for (int j = 0; j <= totalTime; j++) {
						record[i][j] = -1;
					}
				}

				for (int serviceIndex = 0; serviceIndex < numOfService; serviceIndex++) {
					Service service = modelData.serviceSet.get(path.serviceIndexList.get(serviceIndex));
					int startTime = path.timeList.get(serviceIndex) + service.duration;

					for (int time = startTime; time <= startTime + totalTime - path.totalDuration; time++) {
						int edgeTimeTransfer = demand.timeAvailable + time - service.duration;
						edgeTimeTransfer = edgeTimeTransfer % modelData.timePeriod;

						int serviceEdgeIndex = path.serviceIndexList.get(serviceIndex) * modelData.timePeriod
								+ edgeTimeTransfer;
						Edge e = modelData.edgeSet.get(serviceEdgeIndex);
						if (e.serviceIndex != path.serviceIndexList.get(serviceIndex) || e.t1 != edgeTimeTransfer)
							System.out.println("Wrong serviceEdgeIndex calculation!!!");

						if (time == startTime) {
							f[serviceIndex + 1][time] = f[serviceIndex][time - service.duration]
									+ serviceEdgeImportanceSum[serviceEdgeIndex];
							record[serviceIndex + 1][time] = serviceEdgeIndex;
						} else {
							int value = f[serviceIndex][time - service.duration]
									+ serviceEdgeImportanceSum[serviceEdgeIndex];
							if (value > f[serviceIndex + 1][time - 1]) {
								f[serviceIndex + 1][time] = value;
								record[serviceIndex + 1][time] = serviceEdgeIndex;
							} else {
								f[serviceIndex + 1][time] = f[serviceIndex + 1][time - 1];
							}
						}

					}

				}

				int i = numOfService;
				int j = totalTime;
				List<Integer> resultList = new ArrayList<>();
				while (i > 0) {

					if (record[i][j] >= 0) {
						resultList.add(record[i][j]);
						j = j - modelData.edgeSet.get(record[i][j]).duration;
						i = i - 1;
					} else {
						j = j - 1;
					}
				}

				// System.out.println(resultList);

				// add edges to edgesForX
				for (int edgeIndex : resultList) {
					edgeSet.add(edgeIndex);
					// System.out.println(modelData.edgeSet.get(edgeIndex).toString());
				}
			}
		}

		// Add holding arcs to edgesForX
		List<Set<Integer>> legalEdgesForX = modelData.edgesForX;
		for (int k = 0; k < modelData.numDemand; k++) {
			Set<Integer> edgeSet = edgesForX.get(k);
			Set<Integer> legalEdgeSet = legalEdgesForX.get(k);

			for (int edgeIndex = modelData.numServiceArc; edgeIndex < modelData.numArc; edgeIndex++) {
				if (legalEdgeSet.contains(edgeIndex)) {
					edgeSet.add(edgeIndex);
				}
			}
		}

		edgesForXRecord=modelData.edgesForX;
		modelData.ModifyEdgesForX(edgesForX);
		
		//cg-based heuristic
		ColumnGenerationBasedHeuristic cgHeuristic=new ColumnGenerationBasedHeuristic(modelData, 0.6, false);
		cgHeuristic.Solve();
		List<FeasibleSolution> solutionList=new ArrayList<>();
		FeasibleSolution feasibleSolution=new FeasibleSolution(cgHeuristic.optXValues, cgHeuristic.incumbentSolution);
		solutionList.add(feasibleSolution);
		
		return solutionList;
		
	}
	public List<FeasibleSolution> Initialization() {

		// We first calculate how many shortest paths service edge (i,j) belong to
		// commodity k's r shortest paths
		int[][] serviceEdgeImportance = new int[modelData.numServiceArc][modelData.numDemand];
		int[] serviceEdgeImportanceSum = new int[modelData.numServiceArc];
		ArrayList<ArrayList<Path>> rPathSet = new ArrayList<>();
		for (int k = 0; k < modelData.numDemand; k++) {
			rPathSet.add(modelData.findRShortestPath(r, k));
		}

		for (int k = 0; k < modelData.numDemand; k++) {
			Demand demand = modelData.demandSet.get(k);
			ArrayList<Path> pathList = rPathSet.get(k);
			for (Path path : pathList) {

				for (int i = 0; i < path.serviceIndexList.size(); i++) {

					int serviceIndex = path.serviceIndexList.get(i);
					Set<Integer> timeSet = new HashSet<>();
					int time = path.timeList.get(i);

					for (int t = 0; t <= demand.duration - path.totalDuration; t++) {
						// timeAvalible+time+t
						int cTime = demand.timeAvailable + time + t;
						cTime = cTime % modelData.timePeriod;
						timeSet.add(cTime);
					}

					for (int edgeIndex = 0; edgeIndex < modelData.numServiceArc; edgeIndex++) {
						Edge edge = modelData.edgeSet.get(edgeIndex);
						if (edge.serviceIndex == serviceIndex && timeSet.contains(edge.t1)) {
							serviceEdgeImportance[edgeIndex][k]++;
						}
					}
				}
			}
		}

		for (int edgeIndex = 0; edgeIndex < modelData.numServiceArc; edgeIndex++) {
			for (int k = 0; k < modelData.numDemand; k++) {
				serviceEdgeImportanceSum[edgeIndex] += serviceEdgeImportance[edgeIndex][k];
			}
		}

		// We use dp to calculate subEdgeSet for each commodity
		List<Set<Integer>> edgesForX = new ArrayList<>();
		for (int k = 0; k < modelData.numDemand; k++) {
			Set<Integer> set = new HashSet<Integer>();
			edgesForX.add(set);
		}

		for (int k = 0; k < modelData.numDemand; k++) {
			// System.out.println();
			// System.out.println("k= "+k);
			Set<Integer> edgeSet = edgesForX.get(k);
			Demand demand = modelData.demandSet.get(k);

			for (Path path : rPathSet.get(k)) {
				int numOfService = path.serviceIndexList.size();
				int totalTime = demand.duration;

				int[][] f = new int[numOfService + 1][totalTime + 1];
				int[][] record = new int[numOfService + 1][totalTime + 1];
				for (int i = 0; i <= numOfService; i++) {
					for (int j = 0; j <= totalTime; j++) {
						record[i][j] = -1;
					}
				}

				for (int serviceIndex = 0; serviceIndex < numOfService; serviceIndex++) {
					Service service = modelData.serviceSet.get(path.serviceIndexList.get(serviceIndex));
					int startTime = path.timeList.get(serviceIndex) + service.duration;

					for (int time = startTime; time <= startTime + totalTime - path.totalDuration; time++) {
						int edgeTimeTransfer = demand.timeAvailable + time - service.duration;
						edgeTimeTransfer = edgeTimeTransfer % modelData.timePeriod;

						int serviceEdgeIndex = path.serviceIndexList.get(serviceIndex) * modelData.timePeriod
								+ edgeTimeTransfer;
						Edge e = modelData.edgeSet.get(serviceEdgeIndex);
						if (e.serviceIndex != path.serviceIndexList.get(serviceIndex) || e.t1 != edgeTimeTransfer)
							System.out.println("Wrong serviceEdgeIndex calculation!!!");

						if (time == startTime) {
							f[serviceIndex + 1][time] = f[serviceIndex][time - service.duration]
									+ serviceEdgeImportanceSum[serviceEdgeIndex];
							record[serviceIndex + 1][time] = serviceEdgeIndex;
						} else {
							int value = f[serviceIndex][time - service.duration]
									+ serviceEdgeImportanceSum[serviceEdgeIndex];
							if (value > f[serviceIndex + 1][time - 1]) {
								f[serviceIndex + 1][time] = value;
								record[serviceIndex + 1][time] = serviceEdgeIndex;
							} else {
								f[serviceIndex + 1][time] = f[serviceIndex + 1][time - 1];
							}
						}

					}

				}

				int i = numOfService;
				int j = totalTime;
				List<Integer> resultList = new ArrayList<>();
				while (i > 0) {

					if (record[i][j] >= 0) {
						resultList.add(record[i][j]);
						j = j - modelData.edgeSet.get(record[i][j]).duration;
						i = i - 1;
					} else {
						j = j - 1;
					}
				}

				// System.out.println(resultList);

				// add edges to edgesForX
				for (int edgeIndex : resultList) {
					edgeSet.add(edgeIndex);
					// System.out.println(modelData.edgeSet.get(edgeIndex).toString());
				}
			}
		}

		// Add holding arcs to edgesForX
		List<Set<Integer>> legalEdgesForX = modelData.edgesForX;
		for (int k = 0; k < modelData.numDemand; k++) {
			Set<Integer> edgeSet = edgesForX.get(k);
			Set<Integer> legalEdgeSet = legalEdgesForX.get(k);

			for (int edgeIndex = modelData.numServiceArc; edgeIndex < modelData.numArc; edgeIndex++) {
				if (legalEdgeSet.contains(edgeIndex)) {
					edgeSet.add(edgeIndex);
				}
			}
		}

		edgesForXRecord=modelData.edgesForX;
		modelData.ModifyEdgesForX(edgesForX);

		// try to solve use branchAndPriceA

		// Create the pricing problems
		List<SNDRCPricingProblem> pricingProblems = new LinkedList<SNDRCPricingProblem>();
		for (int capacityType = 0; capacityType < modelData.numOfCapacity; capacityType++) {
			for (int originNode = 0; originNode < modelData.numNode; originNode++) {
				String name = "capacity type: " + capacityType + " origin node: " + originNode;
				SNDRCPricingProblem pricingProblem = new SNDRCPricingProblem(modelData, name, capacityType, originNode);
				pricingProblems.add(pricingProblem);
			}

		}

		// Create a cutHandler
		CutHandler<SNDRC, SNDRCMasterData> cutHandler = new CutHandler<>();
		StrongInequalityGenerator cutGen = new StrongInequalityGenerator(modelData, pricingProblems, 0);
		// cutHandler.addCutGenerator(cutGen);

		// Create the Master Problem
		Master master = new Master(modelData, pricingProblems, cutHandler, cutGen, false);

		// Define which solvers to use
		List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers = Collections
				.singletonList(ExactPricingProblemSolver.class);

		// Define one or more Branch creators
		List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators = Arrays.asList(
				new BranchOnLocalServiceForAllPricingProblems(modelData, pricingProblems, 0.5),
				new BranchOnLocalService(modelData, pricingProblems, 0.5),
				new BranchOnServiceEdge(modelData, pricingProblems, 0.5));

		// Create a Branch-and-Price instance
		BranchAndPriceA bap = new BranchAndPriceA(modelData, master, pricingProblems, solvers, branchCreators,
				Double.MAX_VALUE, 0.6, 1, 0.1, 10, 0.0001, 3, false);
		bap.setNodeOrdering(new NodeBoundbapNodeComparator());

		BapLoggerA logger = new BapLoggerA(bap, new File("./output/BAPlogger.log"));

		bap.runBranchAndPrice(System.currentTimeMillis() + 300000L); // 5 mins

		// Print solution:
		System.out.println("================ Solution ================");
		System.out.println("BAP terminated with objective : " + bap.getObjective());
		System.out.println("Total Number of iterations: " + bap.getTotalNrIterations());
		System.out.println("Total Number of processed nodes: " + bap.getNumberOfProcessedNodes());
		System.out.println("Total Time spent on master problems: " + bap.getMasterSolveTime()
				+ " Total time spent on pricing problems: " + bap.getPricingSolveTime());
		System.out.println("Best bound : " + bap.getBound());

		
		ArrayList<FeasibleSolution> solutionList=new ArrayList<>();
		if (bap.hasSolution()) {
			System.out.println("Solution is optimal: " + bap.isOptimal());
			System.out.println("Columns (only non-zero columns are returned):");
			List<Cycle> solution = bap.getSolution();
			int totalNumVehicle = 0;
			int vehicleFixTotalCost = 0;
			int vehicleVarTotalCost = 0;

			for (Cycle column : solution) {
				System.out.println(column);
				System.out.println(out(column) + ":" + bap.GetOptSolutionValueMap().get(column));
				totalNumVehicle += MathProgrammingUtil.doubleToInt((double) bap.GetOptSolutionValueMap().get(column));

				int fixCost = (int) modelData.fixedCost[column.associatedPricingProblem.originNodeO][column.associatedPricingProblem.capacityTypeS];
				int varCost = (int) (column.cost - fixCost);

				double columnValue = (double) bap.GetOptSolutionValueMap().get(column);
				int value = MathProgrammingUtil.doubleToInt(columnValue);
				column.value=value;

				vehicleFixTotalCost += fixCost * columnValue;
				vehicleVarTotalCost += varCost * columnValue;
				// vehicleVarCost.put(column, varCost);
				System.out.println("Fix cost= " + fixCost + " variable cost= " + varCost);

				System.out.println();

			}

			// calculate cost of commodities
			List<Map<Integer, Double>> optXValues = bap.GetOptXValues();
			List<Double> commodityCost = new ArrayList<>();
			for (int commodity = 0; commodity < modelData.numDemand; commodity++) {
				Map<Integer, Double> xValues = optXValues.get(commodity);
				double cost = 0;
				for (int edgeIndex : xValues.keySet()) {
					Edge edge;
					edge = modelData.edgeSet.get(edgeIndex);

					if (edge.edgeType == 0) {
						cost += modelData.beta * edge.duration * xValues.get(edgeIndex);

						// commodityDowork+=edge.duration*xValues.get(edgeIndex);

						// commodityFlowIntoTerminal[edge.v][edge.t2]+=xValues.get(edgeIndex);
					}
				}

				commodityCost.add(cost);
				// commodityTotalCost += cost;
			}
			// output x variables
			for (int demand = 0; demand < modelData.numDemand; demand++) {
				for (int edgeIndex : optXValues.get(demand).keySet()) {
					if (optXValues.get(demand).get(edgeIndex) > 0.01) {
						Edge edge;

						edge = modelData.edgeSet.get(edgeIndex);

						if (edge.edgeType == 0) {
							System.out.println("x[" + demand + "]:" + edge.u + "," + edge.t1 + "->" + edge.v + ","
									+ edge.t2 + "= " + optXValues.get(demand).get(edgeIndex) + " " + edgeIndex);
						}

					}
				}
				System.out.println("total cost= " + commodityCost.get(demand));
				System.out.println();
			}
			
			
			 FeasibleSolution feasibleSolution=new FeasibleSolution(optXValues, solution);
			 solutionList.add(feasibleSolution);

		}

		bap.close();
		cutHandler.close();
		

		return solutionList;

	}

	
	public void TabuSearch(FeasibleSolution currentSolution0,int maxIteration,int columnManageLimit) throws Exception{
	    long startTime=System.currentTimeMillis();
		FeasibleSolution bestFoundSolution=currentSolution0;
		double bestObjectiveValue=currentSolution0.totalCost;
		FeasibleSolution currentSolution=currentSolution0;
		
		//for intensification
		Map<Cycle,Integer> columnManager=new HashMap<>();
		Comparator<Cycle> columnManageSort;
		columnManageSort=new Comparator<Cycle>() {
			@Override
			public int compare(Cycle cycle1,Cycle cycle2){
				return (columnManager.get(cycle2)-columnManager.get(cycle1));
			}
		};

		PriorityQueue<Cycle> columnQueue=new PriorityQueue<Cycle>(columnManageSort);
		for(Cycle cycle:currentSolution.cycleValues){
			columnManager.put(cycle, 0);
		}

		
		//for tabu list
		this.tabuNodeList=new ArrayDeque();
		
		int nonImprovementCount=0;
		boolean ifIntensifyImprove=false;
		int iterationCount=0;
		while(nonImprovementCount<maxIteration||ifIntensifyImprove){
			iterationCount++;
		    long currentTime=System.currentTimeMillis();
		    if(currentTime-startTime>3600000) break;
			
			//update tabuNodeList
			while(tabuNodeList.size()>tabuListLengthLimit){
				tabuNodeList.remove();
			}		
			
			//check if we should make bestSolution as our current solution
			if(nonImprovementCount>=maxIteration){
				currentSolution=bestFoundSolution;
				nonImprovementCount=0;
				ifIntensifyImprove=false;
				System.out.println("We make bestSolution as our current solution!!!");
			}
			
			System.out.println("===========================new round "+iterationCount+"=========================");

			FeasibleSolution newSolution=Neighbourhood(currentSolution);
			if(newSolution.totalCost<bestObjectiveValue-0.001){
				bestObjectiveValue=newSolution.totalCost;
				bestFoundSolution=newSolution;
				nonImprovementCount=0;
				System.out.println("A new better solution has been found. obj="+bestObjectiveValue);
			}else{
				nonImprovementCount++;
			}
			
			currentSolution=newSolution;
			
			//intensification
			for(Cycle cycle:currentSolution.cycleValues){
				boolean check=false;
				for(Cycle cycle0:columnManager.keySet()){
					if(cycle.equals(cycle0)){
						check=true;
						break;
					}
				}
				if(!check){
					columnManager.put(cycle, 0);
				}
			}

			if(columnManager.keySet().size()>columnManageLimit){
				columnQueue.clear();
				for(Cycle cycle:columnManager.keySet()){
					columnQueue.add(cycle);
				}
				while(columnQueue.size()>columnManageLimit){
					columnManager.remove(columnQueue.poll());
				}
			}

			
			FeasibleSolution intensifySolution=Intensification(columnManager.keySet());
			System.out.println("# columns in columnManager="+columnManager.size());
			for(Cycle cycle:columnManager.keySet()){
				if(!intensifySolution.cycleValues.contains(cycle)){
					columnManager.put(cycle, columnManager.get(cycle)+1);
				}else{
					columnManager.put(cycle, 0);
				}
			}

			if(intensifySolution.totalCost<bestObjectiveValue){
				bestFoundSolution=intensifySolution;
				bestObjectiveValue=bestFoundSolution.totalCost;
				ifIntensifyImprove=true;
				System.out.println("A new better solution has been found. obj="+bestObjectiveValue);
				System.out.println("# columns in columnManager="+columnManager.size());
			}
			
		}
		
		System.out.println("best objective value="+bestObjectiveValue);
		for(Cycle cycle:bestFoundSolution.cycleValues){
			System.out.println(cycle.toString());
			System.out.println(out(cycle)+":"+cycle.value);
		}
		
	}
	/**
	 * We search for the neighbourhoods and return a best estimated flow distribution
	 * @param currentSolution
	 * @throws Exception 
	 */
	public FeasibleSolution Neighbourhood(FeasibleSolution currentSolution) throws Exception {
		
		//1. For start, we need to identify the empty vehicle run edge
		Map<Cycle,Map<Integer,Integer>> emptyVehicleEdgeRecord=new HashMap<>();//inside map record how many edges are empty for edgeIndex(key)
		double[] averageFixCostForCapacityType=new double[modelData.numOfCapacity];
		identifyEmptyVehicleRun(emptyVehicleEdgeRecord, averageFixCostForCapacityType, currentSolution);
		
		
		
		
		//2. We search the neighbourhoods based on each terminal node
//		System.out.println();
//		System.out.println("Check commoditySubpathList:");
		double[][] flowSum=new double[modelData.numNode][modelData.numServiceArc];//[terminalIndex][serviceEdgeIndex]
		double[] costModification=new double[modelData.numNode];//[terminalIndex]
		double[] totalFlowCostArray=new double[modelData.numNode];
		List<List<CommoditySubPath>> commoditySubpathListRecord=new ArrayList<>();
		
		for(int terminalIndex=0;terminalIndex<modelData.numNode;terminalIndex++){
			
//			System.out.println();
//			System.out.println("terminal index="+terminalIndex);
			List<Map<Integer,Double>> copyOptXValues=new ArrayList<>();
			for(Map<Integer,Double> map:currentSolution.optXValues){
				copyOptXValues.add(new HashMap<>(map));
			}
			
			List<CommoditySubPath> commoditySubpathList=new ArrayList<>();
			//2.1. find all the commodity paths related to terminalIndex
			for(int commodityIndex=0;commodityIndex<modelData.numDemand;commodityIndex++){
//				Map<Integer, Double> xValues=currentSolution.optXValues.get(commodityIndex);
				Map<Integer, Double> xValues=copyOptXValues.get(commodityIndex);
				
				for(int edgeIndex:xValues.keySet()){
					Edge edge=modelData.edgeSet.get(edgeIndex);
					if(edge.edgeType==0&&(edge.u==terminalIndex||edge.v==terminalIndex)){
						
						double amount=xValues.get(edgeIndex);
						if(amount>0.001){
							List<Integer> pathEdgeIndexList=SourceSearch(edgeIndex, commodityIndex,copyOptXValues);
							//modify xValues in copyOptXValues
							for(int tempEdgeIndex:pathEdgeIndexList){
								xValues.put(tempEdgeIndex, xValues.get(tempEdgeIndex)-amount);
							}
							
							CommoditySubPath subPath=new CommoditySubPath(commodityIndex,amount,pathEdgeIndexList,edgeIndex);
//							System.out.println(subPath.toString());
							commoditySubpathList.add(subPath);
						}
						
					}
				}
			}
			commoditySubpathListRecord.add(new ArrayList<>(commoditySubpathList));
			
			
			
			
			//2.2 remove the unnecessary vehicle edges
			Map<Cycle,Map<Integer,Integer>> removeVehicleEdgeRecord=new HashMap<>();
			double totalRemoveFixCost=removeEmptyVehicleEdge(removeVehicleEdgeRecord, copyOptXValues, currentSolution.cycleValues,emptyVehicleEdgeRecord);
//			System.out.println("totalRemoveFixCost="+totalRemoveFixCost);
	
			
			
			//2.3 Create the residual network for each subPath in commoditySubpathList and redirect the flow
			
			//Now we can drop the cycle info and only record the vehicle edge info on each service edge.
			int[][] vehicleEdgeCover=new int[modelData.numServiceArc][modelData.numOfCapacity];
			double[] flowEdgeCover=new double[modelData.numServiceArc];
			for(Cycle cycle:currentSolution.cycleValues){
				Map<Integer,Integer> map1=new HashMap<>();
				Map<Integer,Integer> map2=new HashMap<>();
				if(emptyVehicleEdgeRecord.containsKey(cycle)){
					map1=emptyVehicleEdgeRecord.get(cycle);
				}
				if(removeVehicleEdgeRecord.containsKey(cycle)){
					map2=removeVehicleEdgeRecord.get(cycle);
				}
				
				for(int edgeIndex:cycle.edgeIndexSet){
					Edge edge=modelData.edgeSet.get(edgeIndex);
					int capacityType=cycle.associatedPricingProblem.capacityTypeS;
					if(edge.edgeType==0){
						vehicleEdgeCover[edgeIndex][capacityType]+=MathProgrammingUtil.doubleToInt(cycle.value);
						if(map1.containsKey(edgeIndex)){
							vehicleEdgeCover[edgeIndex][capacityType]-=map1.get(edgeIndex);
						}
						if(map2.containsKey(edgeIndex)){
							vehicleEdgeCover[edgeIndex][capacityType]-=map2.get(edgeIndex);
						}
					}
				}
			}
			
			for(Map<Integer,Double> map:copyOptXValues){
				for(int edgeIndex:map.keySet()){
					if(modelData.edgeSet.get(edgeIndex).edgeType==0){
						flowEdgeCover[edgeIndex]+=map.get(edgeIndex);
					}
				}
			}
			
			int[][] vehicleEdgeCoverCopy=new int[modelData.numServiceArc][modelData.numOfCapacity];
			for(int i=0;i<vehicleEdgeCover.length;i++){
				vehicleEdgeCoverCopy[i]=Arrays.copyOf(vehicleEdgeCover[i], vehicleEdgeCover[i].length);
			}
			

			//Currently, we plan the commodities in a sequence in commoditySubpathList
			double residualNetworkCost=rebuildSubPath(flowEdgeCover, vehicleEdgeCover, commoditySubpathList, averageFixCostForCapacityType);
			
			//Calculate total cost modification
			double variableCost1=0;
			for(Cycle cycle:emptyVehicleEdgeRecord.keySet()){
				Map<Integer,Integer> map=emptyVehicleEdgeRecord.get(cycle);
				//cost = dataModel.alpha * totalLength / (dataModel.speed * dataModel.drivingTimePerDay)
                //+ dataModel.fixedCost[pricingProblem.originNodeO][pricingProblem.capacityTypeS];
				int totalLength=0;
				for(int edgeIndex:map.keySet()){
					Edge edge=modelData.edgeSet.get(edgeIndex);
					if(edge.edgeType==1) System.out.println("Error edge type!!!");
					totalLength+=edge.duration*map.get(edgeIndex);
				}
				variableCost1=-modelData.alpha/(modelData.speed*modelData.drivingTimePerDay)*totalLength;
			}
			
			double variableCost2=0;
			for(Cycle cycle:removeVehicleEdgeRecord.keySet()){
				Map<Integer,Integer> map=removeVehicleEdgeRecord.get(cycle);
				//cost = dataModel.alpha * totalLength / (dataModel.speed * dataModel.drivingTimePerDay)
                //+ dataModel.fixedCost[pricingProblem.originNodeO][pricingProblem.capacityTypeS];
				int totalLength=0;
				for(int edgeIndex:map.keySet()){
					Edge edge=modelData.edgeSet.get(edgeIndex);
					if(edge.edgeType==1) System.out.println("Error edge type!!!");
					totalLength+=edge.duration*map.get(edgeIndex);
				}
				variableCost2=-modelData.alpha/(modelData.speed*modelData.drivingTimePerDay)*totalLength;
			}
			double fixCost2=-totalRemoveFixCost;
			double flowCost2=0;
			for(CommoditySubPath subPath:commoditySubpathList){
				flowCost2+=subPath.flowCost;
			}
			flowCost2=-flowCost2;
			
			double variableCost3=0;
			double fixCost3=0;
			int totalLength=0;
			
			for(int edgeIndex=0;edgeIndex<modelData.numServiceArc;edgeIndex++){
				Edge edge=modelData.edgeSet.get(edgeIndex);
				for(int capacityIndex=0;capacityIndex<modelData.numOfCapacity;capacityIndex++){
					totalLength+=(vehicleEdgeCover[edgeIndex][capacityIndex]-vehicleEdgeCoverCopy[edgeIndex][capacityIndex])*edge.duration;
					fixCost3+=(vehicleEdgeCover[edgeIndex][capacityIndex]-vehicleEdgeCoverCopy[edgeIndex][capacityIndex])*edge.duration*averageFixCostForCapacityType[capacityIndex];
				}
			}
			variableCost3=modelData.alpha/(modelData.speed*modelData.drivingTimePerDay)*totalLength;
			
			double totalFlowCost=0;
			for(int edgeIndex=0;edgeIndex<modelData.numServiceArc;edgeIndex++){
				Edge edge=modelData.edgeSet.get(edgeIndex);
//				totalFlowCost+=(flowEdgeCover[edgeIndex]-flowEdgeCoverCopy[edgeIndex])*edge.duration*modelData.beta;
				totalFlowCost+=flowEdgeCover[edgeIndex]*edge.duration*modelData.beta;
			}
			
			
//			double totalCostModification=variableCost1+variableCost2+fixCost2+variableCost3+fixCost3+totalFlowCost;
			double totalCostModification=variableCost1+flowCost2+variableCost2+fixCost2+residualNetworkCost;
//			System.out.println(variableCost1+" "+flowCost2+" "+variableCost2+" "+fixCost2+" "+residualNetworkCost);
			if(commoditySubpathList.size()==0){
				totalCostModification=0;
			}
//			System.out.println("totalCostModification="+totalCostModification);
			
			
			flowSum[terminalIndex]=Arrays.copyOf(flowEdgeCover, flowEdgeCover.length);
			costModification[terminalIndex]=totalCostModification;
			totalFlowCostArray[terminalIndex]=totalFlowCost;
			
		}
		
		
		//3. Using bnp to solve vehicle cover problem
//		System.out.println();
//		System.out.println("Step 3-SolveVehicleCover");
		int minIndex=0;
		double minCost=costModification[0];
		for(int i=1;i<costModification.length;i++){
			if(minCost>costModification[i]){
				minCost=costModification[i];
				minIndex=i;
			}
		}
		
		List<CommoditySubPath> commoditySubPathList=commoditySubpathListRecord.get(minIndex);
		for(CommoditySubPath subpath:commoditySubPathList){
			tabuNodeList.add(subpath.tabuNodeIndex);
		}
		
//		List<Cycle> cycleSolution=SolveVehicleCover(flowSum[minIndex],totalFlowCostArray[minIndex]);
		List<Cycle> cycleSolution=SolveVehicleCoverCGHeuristic(flowSum[minIndex],totalFlowCostArray[minIndex],false);
		
		//4. Adjust flow based on the cycleSolution
		System.out.println("Step 4-AdjustFlow");
		List<Map<Integer, Double>> optXValues=AdjustFlow(cycleSolution);
		
		FeasibleSolution output=new FeasibleSolution(optXValues, cycleSolution);
		
		// output x variables
//		for (int demand = 0; demand < modelData.numDemand; demand++) {
//			for (int edgeIndex : optXValues.get(demand).keySet()) {
//				if (optXValues.get(demand).get(edgeIndex) > 0.01) {
//					Edge edge;
//					edge = modelData.edgeSet.get(edgeIndex);
//					if (edge.edgeType == 0) {
//						System.out.println("x[" + demand + "]:" + edge.u + "," + edge.t1 + "->" + edge.v + ","
//								+ edge.t2 + "= " + optXValues.get(demand).get(edgeIndex) + " " + edgeIndex);
//					}
//
//				}
//			}
//			System.out.println();
//		}
		
		
		
		return output;
	}
	
	public List<Map<Integer, Double>> AdjustFlow(List<Cycle> cycleSolution) throws IloException{
		modelData.ModifyEdgesForX(edgesForXRecord);
		
		double[] totalVehicleCapacity=new double[modelData.numServiceArc];
		//calculate total vehicle capacity for each service edge
		for(Cycle cycle:cycleSolution){
			int capacity=modelData.capacity[cycle.associatedPricingProblem.capacityTypeS];
			for(int edgeIndex:cycle.edgeIndexSet){
				Edge edge=modelData.edgeSet.get(edgeIndex);
				if(edge.edgeType==0){
					totalVehicleCapacity[edgeIndex]+=capacity*cycle.value;
				}
			}
		}
		
		IloCplex cplex = new IloCplex();

        cplex.setOut(null);
//        cplex.setParam(IloCplex.IntParam.Threads, config.MAXTHREADS);
        cplex.setParam(IloCplex.Param.Simplex.Tolerances.Markowitz, 0.1);
        List<Map<Integer, IloNumVar>> x; // map:edgeIndex, x variable
        IloNumVar[][] q;

        // Define variables x
        x = new ArrayList<>();
        for (int p = 0; p < modelData.numDemand; p++) {
            Map<Integer, IloNumVar> initialX = new HashMap<Integer, IloNumVar>();
            x.add(initialX);
        }

        // add x variables
        for (int p = 0; p < modelData.numDemand; p++) {
            for (int edgeIndex : modelData.edgesForX.get(p)) {
                Edge edge = modelData.edgeSet.get(edgeIndex);
                IloNumVar varX = cplex.numVar(0,modelData.demandSet.get(p).volume,
                        "x" + p + "," + edge.start + "," + edge.end);
                x.get(p).put(edgeIndex, varX);
            }
        }

        // Define the objective
        /**
         * Here we assume the cost of edge AT is 0
         */
        IloLinearNumExpr exprObj = cplex.linearNumExpr();


        for (int p = 0; p < modelData.numDemand; p++) {
            Map<Integer, IloNumVar> map = x.get(p);
            for (int edgeIndex : map.keySet()) {
                exprObj.addTerm(modelData.beta *modelData.edgeSet.get(edgeIndex).duration, map.get(edgeIndex));
            }
        }

        IloObjective obj = cplex.addMinimize(exprObj);
        
        // Define flowBalanceConstraints
        IloRange[][] flowBalanceConstraints = new IloRange[modelData.numDemand][modelData.abstractNumNode];

        IloLinearNumExpr expr = cplex.linearNumExpr();
        for (int p = 0; p < modelData.numDemand; p++) {
            Map<Integer, IloNumVar> map = x.get(p);

            for (int i = 0; i < modelData.abstractNumNode; i++) {
                expr.clear();
                // edges which point from i
                for (int edgeIndex : modelData.pointToEdgeSet.get(i)) {
                    if (map.containsKey(edgeIndex)) {
                        expr.addTerm(1, map.get(edgeIndex));
                    }
                }

                // edges which point to i
                for (int edgeIndex : modelData.pointFromEdgeSet.get(i)) {
                    if (map.containsKey(edgeIndex)) {
                        expr.addTerm(-1, map.get(edgeIndex));
                    }
                }
                flowBalanceConstraints[p][i] = cplex.addEq(modelData.b[p][i], expr);

            }
        }

        // Define weakForcingConstraints
        IloRange[] weakForcingConstraints = new IloRange[modelData.numServiceArc];
        for (int arcIndex = 0; arcIndex < modelData.numServiceArc; arcIndex++) {
            expr.clear();
            for (int p = 0; p < modelData.numDemand; p++) {
                if (x.get(p).containsKey(arcIndex)) {
                    expr.addTerm(1, x.get(p).get(arcIndex));
                }
            }

            weakForcingConstraints[arcIndex] = cplex.addGe(totalVehicleCapacity[arcIndex], expr);
        }
        
        cplex.solve();
        System.out.println("After adjust, flowCost="+cplex.getObjValue());
        
        List<Map<Integer, Double>> xValues = new ArrayList<>();
        for (int commodity = 0; commodity < modelData.numDemand; commodity++) {
            Map tempMap = new HashMap<>();
            for (int edgeIndex : x.get(commodity).keySet()) {
                tempMap.put(edgeIndex, cplex.getValue(x.get(commodity).get(edgeIndex)));
            }
            xValues.add(tempMap);
        }
        
        return xValues;
        
	}

	public FeasibleSolution Intensification(Set<Cycle> cycleSolution) throws IloException, FileNotFoundException{
		// firstly, we create edgesForX(including holding arc)		
		Set<Integer> edgeIndexSet=new HashSet<>();
		for(Cycle cycle:cycleSolution){
			for(int edgeIndex:cycle.edgeIndexSet){
				if(modelData.edgeSet.get(edgeIndex).edgeType==0){
					edgeIndexSet.add(edgeIndex);
				}
			}
		}
		
		List<Set<Integer>> edgesForX = new ArrayList<>();
		for (int k = 0; k < modelData.numDemand; k++) {
			Set<Integer> set = new HashSet<Integer>();
			Set<Integer> setRecord = edgesForXRecord.get(k);
			for(int edgeIndex:setRecord){
				if(modelData.edgeSet.get(edgeIndex).edgeType==1){
					set.add(edgeIndex);
				}else{
					if(edgeIndexSet.contains(edgeIndex)){
						set.add(edgeIndex);
					}
				}
			}
			edgesForX.add(set);
		}
		
		
		
		//cplex solving
        IloCplex cplex=new IloCplex();
//		FileOutputStream outputStream=new FileOutputStream(new File("./output/cplexOut.txt"));
        cplex.setParam(IloCplex.DoubleParam.EpGap, 0.005);
		cplex.setOut(null);
		
		cplex.setParam(IloCplex.IntParam.Threads, 4);
		cplex.setParam(IloCplex.Param.Simplex.Tolerances.Markowitz, 0.1);
		cplex.setParam(IloCplex.DoubleParam.TiLim, 180); //3 mins
		List<Map<Integer,IloNumVar>> x; //map:edgeIndex, x variable
		Map<Path,IloNumVar> pathVarMap=new HashMap<>();
		
		// Define variables x
		x=new ArrayList<Map<Integer,IloNumVar>>();
		for(int p=0;p<modelData.numDemand;p++) {
			Map<Integer,IloNumVar> initialX=new HashMap<Integer,IloNumVar>();
			x.add(initialX);
		}
		
		// add x variables
		for(int p=0;p<modelData.numDemand;p++) {
			for(int edgeIndex:edgesForX.get(p)) {
				Edge edge=modelData.edgeSet.get(edgeIndex);
				IloNumVar varX=cplex.numVar(0, modelData.demandSet.get(p).volume,"x" +p+","+ edge.start + "," + edge.end );
				x.get(p).put(edgeIndex, varX);
			}
		}
		
		
		// Define the objective
		/**
		 * Here we assume the cost of edge AT is 0
		 */
		IloLinearNumExpr exprObj = cplex.linearNumExpr();
		
		for(int p=0;p<modelData.numDemand;p++) {
			Map<Integer,IloNumVar> map=x.get(p);
			for(int edgeIndex:map.keySet()) {
				exprObj.addTerm(modelData.beta*modelData.edgeSet.get(edgeIndex).duration, map.get(edgeIndex));
			}
		}
		

		IloObjective obj = cplex.addMinimize(exprObj);
		
		
		
		// Define flowBalanceConstraints
		IloRange[][] flowBalanceConstraints = new IloRange[modelData.numDemand][modelData.abstractNumNode];

		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int p = 0; p < modelData.numDemand; p++) {
			Map<Integer,IloNumVar> map=x.get(p);
			
			for (int i = 0; i < modelData.abstractNumNode; i++) {
				expr.clear();
				// edges which point from i
				for (int edgeIndex : modelData.pointToEdgeSet.get(i)) {
					if(map.containsKey(edgeIndex)) {
						expr.addTerm(1, map.get(edgeIndex));
					}
				}

				// edges which point to i
				for (int edgeIndex : modelData.pointFromEdgeSet.get(i)) {
					if(map.containsKey(edgeIndex)) {
						expr.addTerm(-1, map.get(edgeIndex));
					}
				}
				flowBalanceConstraints[p][i] = cplex.addEq(modelData.b[p][i], expr);

			}
		}
		
		
		// Define weakForcingConstraints
		IloRange[] weakForcingConstraints = new IloRange[modelData.numServiceArc];
		for (int arcIndex = 0; arcIndex < modelData.numServiceArc; arcIndex++) {
			expr.clear();
			for (int p = 0; p < modelData.numDemand; p++) {
				if(x.get(p).containsKey(arcIndex)) {
					expr.addTerm(1, x.get(p).get(arcIndex));
				}
			}

			weakForcingConstraints[arcIndex] = cplex.addGe(0, expr);
		}
		
		
		
		
		// Define resourceBoundConstraints
		IloRange[][] resourceBoundConstraints = new IloRange[modelData.numOfCapacity][modelData.numNode];


		for (int s = 0; s < modelData.numOfCapacity; s++) {
			for (int o = 0; o < modelData.numNode; o++) {
				expr.clear();
				resourceBoundConstraints[s][o]=cplex.addRange(0,modelData.vehicleLimit[s][o]);
			}
		}
		
		
		
		//add columns
		int count=0;
		Map<Cycle,IloNumVar> columnMap=new HashMap<Cycle,IloNumVar>();
		for(Cycle cycle:cycleSolution){			
			// Register column with objective
			IloColumn iloColumn = cplex.column(obj, cycle.cost);
			// weak forcing constraints
			for (int edgeIndex : cycle.edgeIndexSet) {
				if(modelData.edgeSet.get(edgeIndex).edgeType==0) {
					iloColumn = iloColumn.and(cplex.column(weakForcingConstraints[edgeIndex],
							-modelData.capacity[cycle.associatedPricingProblem.capacityTypeS]));
				}
			}
			
			
			// resource bound constraints
			int capacityType=cycle.associatedPricingProblem.capacityTypeS;
			int originNode=cycle.associatedPricingProblem.originNodeO;
			iloColumn = iloColumn.and(cplex.column(resourceBoundConstraints[capacityType][cycle.associatedPricingProblem.originNodeO],1));
			
			
			// Create the variable and store it
			IloNumVar var =cplex.intVar(iloColumn, 0, Integer.MAX_VALUE,"z_" + capacityType + ","+originNode+","+count++);
			columnMap.put(cycle,var);
		}
		
		cplex.solve();
		System.out.println("Intensification optimal objective= "+cplex.getObjValue());
//		System.out.println("# columns= "+cycleSolution.size());
//		//output check
//		
//		for(Cycle cycle:cycleSolution){
//			Double value=cplex.getValue((IloNumVar) columnMap.get(cycle));
//			if(value>0.001){
//				System.out.println(cycle.toString()+":"+value);
//				System.out.println(out(cycle));
//			}
//		}
		List<Map<Integer, Double>> optXValues=new ArrayList<>();
		List<Cycle> cycleValues=new ArrayList<>();
		for(int k=0;k<modelData.numDemand;k++){
			Map<Integer, Double> tempMap=new HashMap<>();
			for(int edgeIndex:edgesForX.get(k)){
				double value=cplex.getValue(x.get(k).get(edgeIndex));
				if(value>0.001){
					tempMap.put(edgeIndex, value);
				}
			}
			optXValues.add(tempMap);
		}
		
		for(Cycle cycle:cycleSolution){
		  Double value=cplex.getValue(columnMap.get(cycle));
		  if(value>0.001){
			  cycle.value=MathProgrammingUtil.doubleToRoundedDouble(value);
			  cycleValues.add(cycle);
		  }
		}
		

		FeasibleSolution feasibleSolution=new FeasibleSolution(optXValues, cycleValues);
		return feasibleSolution;

		
	}
	
 	public List<Cycle> SolveVehicleCoverCGHeuristic(double[] flowCover,double totalFlowCost,boolean ifAccelerationForUB) throws TimeLimitExceededException, IloException{
        List<Cycle> cycleList = new ArrayList<>();
        int bestObjValue=Integer.MAX_VALUE;
//		System.out.println("-------------Start to solve SolveVehicleCoverCGHeuristic()---------------");
        long runTime = System.currentTimeMillis();

        /// -----------------------------step1: solve the root node------------------------------------///
        // Create the pricing problems
        List<SNDRCPricingProblem> pricingProblems = new LinkedList<SNDRCPricingProblem>();
        for (int capacityType = 0; capacityType < modelData.numOfCapacity; capacityType++) {
            for (int originNode = 0; originNode < modelData.numNode; originNode++) {
                String name = "capacity type: " + capacityType + " origin node: " + originNode;
                SNDRCPricingProblem pricingProblem = new SNDRCPricingProblem(modelData, name, capacityType, originNode);
                pricingProblems.add(pricingProblem);
            }
        }

        // Create a cutHandler
        CutHandler<SNDRC, SNDRCMasterData> cutHandler = new CutHandler<>();
        StrongInequalityGenerator cutGen = new StrongInequalityGenerator(modelData, pricingProblems, 0);
        // cutHandler.addCutGenerator(cutGen);

        // Create the Master Problem
        modelData.flowCover=Arrays.copyOf(flowCover, flowCover.length);
        modelData.totalFlowCost=totalFlowCost;
        MasterForVehicleCover master = new MasterForVehicleCover(modelData, pricingProblems, cutHandler, cutGen, false);

        // Define which solvers to use
        List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers = Collections.singletonList(ExactPricingProblemSolver.class);
        
        // ----------------------------------generateInitialFeasibleSolution-------------------------------------------------//
        int[] temp=new int[modelData.numService];
        List<Cycle> artificalVars = new ArrayList<Cycle>();
        // for weak forcing constraints(ifForResourceBoundConstraints=0)
        for (int edgeIndex = 0; edgeIndex < modelData.numServiceArc; edgeIndex++) {
            Set<Integer> set = new HashSet<>();
            set.add(edgeIndex);
            Cycle cycle = new Cycle(pricingProblems.get(0), true, "Artificial", set, 100000000, 0, 0,temp);
            artificalVars.add(cycle);
        }

        // for resource bound constraints(ifForResourceBoundConstraints=1)
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            Set<Integer> set = new HashSet<>();
            Cycle cycle = new Cycle(pricingProblem, true, "Artificial", set, 100000000, 0, 1,temp);
            artificalVars.add(cycle);
        }

        // for holding edge branch constraints(ifForResourceBoundConstraints=2)
        Set<Integer> set = new HashSet<>();
        Cycle cycle0 = new Cycle(pricingProblems.get(0), true, "Artificial", set, 100000000, 0, 2,temp);
        artificalVars.add(cycle0);

        // ----------------------------------generateInitialFeasibleSolution-------------------------------------------------//
        
        
        ColGen<SNDRC, Cycle, SNDRCPricingProblem> cg = new ColGen<SNDRC, Cycle, SNDRCPricingProblem>(modelData, master,
                pricingProblems, solvers, artificalVars, Integer.MAX_VALUE, Double.MIN_VALUE);

        // SimpleDebugger debugger = new SimpleDebugger(cg);

        // SimpleCGLogger logger = new SimpleCGLogger(cg, new
        // File("./output/cgLogger.log"));
        
        
        cg.solve(System.currentTimeMillis() + 36000000L); // 10 hour limit

//        System.out.println("Time of first LP solve= " + (System.currentTimeMillis() - runTime));

        if (ifAccelerationForUB) {
            /// -------------------------------AccelerationForUB------------------------------------///
            List<Cycle> solution = cg.getSolution();

            while (!isIntegerNode(solution)) {

                boolean ifAllBelowThresholdValue = true;
                solution = cg.getSolution();

                for (Cycle cycle : solution) {
                    if (MathProgrammingUtil.isFractional(cycle.value)) {
                        double decimalValue = cycle.value - (int) cycle.value;
                        if (decimalValue > 0.6) {
                            ifAllBelowThresholdValue = false;
                            ((Master) master).addFixVarConstraint(cycle);
                        }
                    }
                }

                // if all cycles' value are below the threshold value, fix the
                // variable with highest fractional decimal value
                double record = 0;
                Cycle cycleRecord = null;
                if (ifAllBelowThresholdValue) {
                    for (Cycle cycle : solution) {
                        if (MathProgrammingUtil.isFractional(cycle.value)) {
                            double decimalValue = cycle.value - (int) cycle.value;
                            if (decimalValue > record) {
                                cycleRecord = cycle;
                                record = decimalValue;
                            }
                        }
                    }

                    if (cycleRecord != null) {
                        ((Master) master).addFixVarConstraint(cycleRecord);
                    } else {
                        break;
                    }
                }

                // here we should check if the master problem is feasible
                if (((Master) master).CheckFeasibility() == false) {
                    break;
                }

                List<Cycle> nullList = new ArrayList<>();
                cg = new ColGen<>(modelData, master, pricingProblems, solvers, nullList, Integer.MAX_VALUE,
                        Double.MIN_VALUE);
                cg.solve(System.currentTimeMillis() + 18000000L); // 5 hour

                if (isInfeasibleNode(cg.getSolution())) {
                    break;
                }

                if (isIntegerNode(cg.getSolution())) {
                    int integerObjective = MathProgrammingUtil.doubleToInt(cg.getObjective());
                    System.out.println("By acceleration,we find a feasible solution: " + integerObjective+"+"+totalFlowCost+"="+(integerObjective+totalFlowCost));

                    bestObjValue=integerObjective;
                    for (Cycle cycle : cg.getSolution()) {
                        cycleList.add(cycle);
                    }
                    
                    
                    break;
                }

            }

            /// -------------------------------AccelerationForUB------------------------------------///
            System.out.println("Time of acceleration LP solve= " + (System.currentTimeMillis() - runTime));
        }
        
        
        // pick up all the cycles in master problem
        int amount = 0;
        Map<SNDRCPricingProblem, Set<Cycle>> cycleSet = new HashMap<>();
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            Set<Cycle> tempSet = master.getColumns(pricingProblem);
            cycleSet.put(pricingProblem, tempSet);
            amount += tempSet.size();
        }
        
        
        
        
        
      /// --------------------------------step2:construct a cplex model------------------------------------///
        IloCplex cplex = new IloCplex();
        cplex.setOut(null);
        cplex.setParam(IloCplex.IntParam.Threads, 4);
        cplex.setParam(IloCplex.Param.Simplex.Tolerances.Markowitz, 0.1);



        // Define the objective
        /**
         * Here we assume the cost of edge AT is 0
         */
        IloLinearNumExpr exprObj = cplex.linearNumExpr();
        IloObjective obj = cplex.addMinimize(exprObj);


        IloLinearNumExpr expr = cplex.linearNumExpr();

        // Define weakForcingConstraints
        IloRange[] weakForcingConstraints = new IloRange[modelData.numServiceArc];
        for (int arcIndex = 0; arcIndex < modelData.numServiceArc; arcIndex++) {
            expr.clear();
//            for (int p = 0; p < modelData.numDemand; p++) {
//                if (x.get(p).containsKey(arcIndex)) {
//                    expr.addTerm(1, x.get(p).get(arcIndex));
//                }
//            }
            weakForcingConstraints[arcIndex] = cplex.addGe(-modelData.flowCover[arcIndex], expr);

//            weakForcingConstraints[arcIndex] = cplex.addGe(0, expr);
        }

        // Define resourceBoundConstraints
        IloRange[][] resourceBoundConstraints = new IloRange[modelData.numOfCapacity][modelData.numNode];

        for (int s = 0; s < modelData.numOfCapacity; s++) {
            for (int o = 0; o < modelData.numNode; o++) {
                expr.clear();
                resourceBoundConstraints[s][o] = cplex.addRange(0, modelData.vehicleLimit[s][o]);
            }
        }

        // add all columns in cycleSet
        Map<Cycle,IloIntVar> cycleVar = new HashMap<>();
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            int count = 0;
            Set<Cycle> tempSet = cycleSet.get(pricingProblem);

            for (Cycle cycle : tempSet) {

                // Register column with objective
                IloColumn iloColumn = cplex.column(obj, cycle.cost);

                // weak forcing constraints
                for (int edgeIndex : cycle.edgeIndexSet) {
                    if (modelData.edgeSet.get(edgeIndex).edgeType == 0) {
                        iloColumn = iloColumn.and(cplex.column(weakForcingConstraints[edgeIndex],
                                -modelData.capacity[cycle.associatedPricingProblem.capacityTypeS]));
                    }
                }

                // resource bound constraints
                iloColumn = iloColumn.and(cplex.column(
                        resourceBoundConstraints[cycle.associatedPricingProblem.capacityTypeS][cycle.associatedPricingProblem.originNodeO],
                        1));

                // Create the variable and store it
                IloIntVar var = cplex.intVar(iloColumn, 0, Integer.MAX_VALUE,
                        "z_" + cycle.associatedPricingProblem.capacityTypeS + ","
                                + cycle.associatedPricingProblem.originNodeO + "," + count);
                cycleVar.put(cycle, var);
                count++;

            }

        }
        
        
        /// -----------------------------step3:cplex solve and output---------------------------------------------///      
        cg.close();
        cutHandler.close();
        cutGen.close();
        master.close();

//        System.out.println();
//        System.out.println("There are " + amount + " columns added to the model.");
//        System.out.println();

        runTime = System.currentTimeMillis() - runTime;
//        long timeLeft = 3600000 - runTime;
        long timeLeft = 180000; //3 mins
//        long timeLeft = 20000 - runTime;
        cplex.setParam(IloCplex.DoubleParam.TiLim, timeLeft / 1000);
        
//        cplex.exportModel("./output/masterLP/" + "check"  + ".lp");
        cplex.solve();
        System.out.println("optimal objective= " + cplex.getObjValue()+"+"+totalFlowCost+"="+(cplex.getObjValue()+totalFlowCost));
        System.out.println();

        
        // record the solution
        if(cplex.getObjValue()<bestObjValue-0.001){
            cycleList=new ArrayList<>();
            for (Cycle cycle : cycleVar.keySet()) {
                IloIntVar var = cycleVar.get(cycle);
                int value = MathProgrammingUtil.doubleToInt(cplex.getValue(var));

                if (value > 0.01) {
                	cycleList.add(cycle);
                    cycle.value=value;
                }
            }

            // output solution
//            for (Cycle cycle : cycleList) {
//                    System.out.println(cycle);
//                    System.out.println(out(cycle) + ":" + cycle.value);
////                    System.out.println();
//            }
        }

        cplex.end();
        return cycleList;
        
        
        
	}
	
	//for SolveVehicleCoverCGHeuristic
    public boolean isIntegerNode(List<Cycle> solution) {
        boolean out = true;
        for (Cycle cycle : solution) {
            if (MathProgrammingUtil.isFractional(cycle.value)) {
                out = false;
                break;
            }
        }

        return out;
    }
    //for SolveVehicleCoverCGHeuristic
    public boolean isInfeasibleNode(List<Cycle> solution) {
        boolean out = false;
        for (Cycle cycle : solution) {
            if (cycle.isArtificialColumn) {
                out = true;
                break;
            }
        }

        return out;

    }
	
	public List<Cycle> SolveVehicleCover(double[] flowCover,double totalFlowCost){
		
        // Create the pricing problems
        List<SNDRCPricingProblem> pricingProblems = new LinkedList<SNDRCPricingProblem>();
        for (int capacityType = 0; capacityType < modelData.numOfCapacity; capacityType++) {
            for (int originNode = 0; originNode < modelData.numNode; originNode++) {
                String name = "capacity type: " + capacityType + " origin node: " + originNode;
                SNDRCPricingProblem pricingProblem = new SNDRCPricingProblem(modelData, name, capacityType, originNode);
                pricingProblems.add(pricingProblem);
            }
        }

        // Create a cutHandler
        CutHandler<SNDRC, SNDRCMasterData> cutHandler = new CutHandler<>();
        StrongInequalityGenerator cutGen = new StrongInequalityGenerator(modelData, pricingProblems, 0);
        // cutHandler.addCutGenerator(cutGen);

        // Create the Master Problem
        modelData.flowCover=Arrays.copyOf(flowCover, flowCover.length);
        modelData.totalFlowCost=totalFlowCost;
        MasterForVehicleCover master = new MasterForVehicleCover(modelData, pricingProblems, cutHandler, cutGen, false);

        // Define which solvers to use
        List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers = Collections
                .singletonList(ExactPricingProblemSolver.class);
        
        List<? extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>> branchCreators = Arrays.asList(
                new BranchOnLocalServiceForAllPricingProblems(modelData, pricingProblems, 0.5),
                new BranchOnLocalService(modelData, pricingProblems, 0.5),
                new BranchOnServiceEdge(modelData, pricingProblems, 0.5));
        
        BranchAndPriceA bap=new BranchAndPriceA(modelData, master,pricingProblems, solvers,branchCreators,Double.MAX_VALUE,0.6,0.2,0.1,10,0.0001,3,false);
        bap.setNodeOrdering(new NodeBoundbapNodeComparator());
        
        BapLoggerA logger = new BapLoggerA(bap, new File("./output/BAPlogger.log"));
        // OPTIONAL: Attach a debugger
//        SimpleDebugger debugger=new SimpleDebugger(bap, true);

        bap.runBranchAndPrice(System.currentTimeMillis() + 180000L); // 3mins
                                                                       
        // Print solution:
        System.out.println("================ Solution ================");
        System.out.println("BAP terminated with objective : " + bap.getObjective()+"+"+totalFlowCost+"="+(bap.getObjective()+totalFlowCost));
        System.out.println("Total Number of iterations: " + bap.getTotalNrIterations());
        System.out.println("Total Number of processed nodes: " + bap.getNumberOfProcessedNodes());
        System.out.println("Total Time spent on master problems: " + bap.getMasterSolveTime()
                + " Total time spent on pricing problems: " + bap.getPricingSolveTime());
        System.out.println("Best bound : " + bap.getBound());
        
        
        List<Cycle> solution=new ArrayList<>();
        if (bap.hasSolution()) {
            System.out.println("Solution is optimal: " + bap.isOptimal());
            System.out.println("Columns (only non-zero columns are returned):");
            solution = bap.getSolution();
            for (Cycle column : solution) {
                System.out.println(column);
                System.out.println(out(column) + ":" + bap.GetOptSolutionValueMap().get(column));

                int fixCost = (int) modelData.fixedCost[column.associatedPricingProblem.originNodeO][column.associatedPricingProblem.capacityTypeS];
                int varCost = (int) (column.cost - fixCost);

                double columnValue= (double) bap.GetOptSolutionValueMap().get(column);
                int value=MathProgrammingUtil.doubleToInt(columnValue);
                column.value=value;
                
                System.out.println("Fix cost= " + fixCost + " variable cost= " + varCost);

                System.out.println();
            }
        }
        
        
        bap.close();
        cutHandler.close();
        
        return solution;
	}
	
	public double rebuildSubPath(double[] flowEdgeCover,int[][] vehicleEdgeCover,List<CommoditySubPath> commoditySubpathList,double[] averageFixCostForCapacityType) throws Exception{
		
		double totalCost=0;
		
		List<Integer> shortestPath=new ArrayList<>();
		
		//flow 安排顺序改为随机
		int[] indexSeq=new int[commoditySubpathList.size()];
		List<Integer> indexList=new LinkedList<>();
		for(int i=0;i<commoditySubpathList.size();i++){
			indexList.add(i);
		}
		for(int i=0;i<commoditySubpathList.size();i++){
			int index=(int)(Math.random()*indexList.size());
			indexSeq[i]=indexList.remove(index);
		}
		
//		for(CommoditySubPath subPath:commoditySubpathList){
//		System.out.println("indexSeq="+Arrays.toString(indexSeq));
		for(int i=0;i<indexSeq.length;i++){
			CommoditySubPath subPath=commoditySubpathList.get(indexSeq[i]);
			
			double amount=subPath.amount;
			int[][] modifyVehicleEdgeCover=new int[modelData.numServiceArc][modelData.numOfCapacity];
			double[] residualNetwork=createResidualNetwork(flowEdgeCover,vehicleEdgeCover,subPath,modifyVehicleEdgeCover,averageFixCostForCapacityType);
			
			
			totalCost+=findShortestPath(subPath,residualNetwork,shortestPath);
			
			//modify flowEdgeCover vehicleEdgeCover based on shortestPath and modifyVehicleEdgeCover
			for(int edgeIndex:shortestPath){
				Edge edge=modelData.edgeSet.get(edgeIndex);
				if(edge.edgeType==0){
					flowEdgeCover[edgeIndex]+=amount;
					for(int capacityIndex=0;capacityIndex<modelData.numOfCapacity;capacityIndex++){
						vehicleEdgeCover[edgeIndex][capacityIndex]+=modifyVehicleEdgeCover[edgeIndex][capacityIndex];
					}
				}
				
			}
			shortestPath.clear();
			
		}
		
//		System.out.println("Check flowEdgeCover:");
//		for(int edgeIndex=0;edgeIndex<modelData.numServiceArc;edgeIndex++){
//			Edge edge=modelData.edgeSet.get(edgeIndex);
//			if(flowEdgeCover[edgeIndex]>0.001){
//				System.out.println(edge.toString()+" "+flowEdgeCover[edgeIndex]);
//			}
//		}
		
		
		
		return totalCost;
		
	}
	
	public double findShortestPath(CommoditySubPath subPath,double[] residualNetwork,List<Integer> shortestPath) throws Exception{
		int startTime=subPath.startNodeIndex%modelData.timePeriod;
		int endTime=subPath.endNodeIndex%modelData.timePeriod;
		int timeLength=endTime-startTime;
		if(timeLength<=0){
			timeLength+=modelData.timePeriod;
		}
		
		
		double[] f=new double[modelData.abstractNumNode];
		for(int i=0;i<f.length;i++){
			f[i]=Integer.MAX_VALUE;
		}
		f[subPath.startNodeIndex]=0;
		int[] record=new int[modelData.abstractNumNode];
		
		for(int time=0;time<timeLength;time++){
			int timeIndex=startTime+time;
			timeIndex=timeIndex%modelData.timePeriod;
			
			for(int terminalIndex=0;terminalIndex<modelData.numNode;terminalIndex++){
				int nodeIndex=terminalIndex*modelData.timePeriod+timeIndex;
				if(f[nodeIndex]<Integer.MAX_VALUE-1000){
					for(int edgeIndex:modelData.pointToEdgeSet.get(nodeIndex)){
						Edge edge=modelData.edgeSet.get(edgeIndex);
						int duration=edge.duration;
						if(edge.edgeType==1){
							duration=1;
						}
						
						if(time+duration<=timeLength){
							if(!tabuNodeList.contains(edge.start)&&!tabuNodeList.contains(edge.end)){
								if(edge.edgeType==0){
										double distance=f[nodeIndex]+residualNetwork[edgeIndex];
										if(f[edge.end]>distance){
											f[edge.end]=distance;
											record[edge.end]=edgeIndex;
										}							
								}else{ //holding edge
									if(f[edge.end]>f[nodeIndex]){
										f[edge.end]=f[nodeIndex];
										record[edge.end]=edgeIndex;
									}
								}
							}

						}
					}
				}
			}
			
		}
		
		
		//闂佹眹鍨硅ぐ澶岃姳椤掑偊鎷峰☉娅亜锕㈤敓锟� tabu list 闂佸憡鐟崹鐢稿礂濮楋拷楠炲秹宕ｉ妷褏鎲归梺鍛婂笧婢ф銆掗崼鏇熷剹妞ゆ搫绱曢悢鍛存煥濞戞鐒风紒缁樼墵閺屽瞼浠﹂幆褏浜繛瀵稿閸撴繄鎹㈤幋锕�鐐婇柣鎰暯閺佸嫰鏌ｉ妸銉ヮ伂缂佹梻鍠栧鍧楀幢濡や礁顏疍ouble.MaxValue,闁哄鏅滈弻銊ッ洪弽顓炶Е閹兼番鍊ら崵濉甽ow婵炲濮寸粔鐢碉拷鍨矒閹ゎ槺缂佺媴缍佸鑽ゅ鐎ｎ剛鍞撮悗鍨緲鐎氼厾绮婇敂鑺ヤ氦闁搞儵顥撶粈澶愭煕濮橆剚婀版俊鐐插�绘禍鎼佸幢濞嗘劗銈查梺绋跨箰閻ゅ洨绮╅悢鍝ヮ浄鐎瑰憡顩竍u list 婵炴垶鎼╅崢鍏兼櫠瑜斿浠嬫晸閿燂拷
		
//		System.out.println(f[subPath.endNodeIndex]);
		if(f[subPath.endNodeIndex]>100000000){
//			shortestPath=subPath.pathEdgeIndexList;
			for(int edgeIndex:subPath.pathEdgeIndexList){
				shortestPath.add(edgeIndex);
			}
			double length=0;
			for(int edgeIndex:shortestPath){
				if(edgeIndex<modelData.numServiceArc){
					length+=residualNetwork[edgeIndex];
				}
				
				//闁诲繐绻愰弫鍝竍u list婵炴垶鎼╅崢鎯р枔閹寸偞缍囬柣鐔告緲缁犵敻鏌熼悮瀛樺
//				if(tabuList.contains(edgeIndex)){
//					tabuList.remove(edgeIndex);
//				}
			}		
			return length;
			
		}else{
			int nodeIndex=subPath.endNodeIndex;
			
			while(nodeIndex!=subPath.startNodeIndex){
				int edgeIndex=record[nodeIndex];
				shortestPath.add(edgeIndex);
//				System.out.print(modelData.edgeSet.get(edgeIndex).toString()+" ");
				nodeIndex=modelData.edgeSet.get(edgeIndex).start;
				
				if(shortestPath.size()>50){
					System.out.println();
					System.out.println("Error!!!");
					System.out.println(f[subPath.endNodeIndex]);
					System.out.println(subPath.toString());
					for(int serviceEdgeIndex=0;serviceEdgeIndex<modelData.numServiceArc;serviceEdgeIndex++){
						Edge edge=modelData.edgeSet.get(serviceEdgeIndex);
						System.out.println(edge.toString()+":"+residualNetwork[serviceEdgeIndex]);
					}
					throw new Exception("闂佸搫鐗為幏鐑芥煟椤撗冨绩闁活厼澧界划璇参旀担鎭掞拷濠囨煕濞嗘劕鐏╅柡浣规崌閺屻劌鈻庨幒婵嗘");
				}
			}
//			System.out.println();

			
			return f[subPath.endNodeIndex];
		}

		
		
		
	}
	
	public double[] createResidualNetwork(double[] flowEdgeCover,int[][] vehicleEdgeCover,CommoditySubPath subPath,int[][] modifyVehicleEdgeCover,double[] averageFixCostForCapacityType){
		int[] totalVehicleCapacityForServiceEdge=new int[modelData.numServiceArc];
		for(int edgeIndex=0;edgeIndex<modelData.numServiceArc;edgeIndex++){
			totalVehicleCapacityForServiceEdge[edgeIndex]=0;
			for(int capacityType=0;capacityType<modelData.numOfCapacity;capacityType++){
				totalVehicleCapacityForServiceEdge[edgeIndex]+=vehicleEdgeCover[edgeIndex][capacityType]*modelData.capacity[capacityType];
			}
		}
		
		double amount=subPath.amount;
		double variableCostPara=modelData.alpha / (modelData.speed * modelData.drivingTimePerDay);
		double[] residualNetworkCost=new double[modelData.numServiceArc];
		for(int serviceEdgeIndex=0;serviceEdgeIndex<modelData.numServiceArc;serviceEdgeIndex++){
			
			Edge edge=modelData.edgeSet.get(serviceEdgeIndex);
			if(totalVehicleCapacityForServiceEdge[serviceEdgeIndex]<0.1){// no vehicle cover
				residualNetworkCost[serviceEdgeIndex]=modelData.beta*edge.duration*amount;
				
				int capacityType=-1;
				double amount_copy=amount;
				while(amount_copy>0.01){
					
					for(int capacityIndex=0;capacityIndex<modelData.numOfCapacity;capacityIndex++){
						if(modelData.capacity[capacityIndex]>amount-0.001){
							capacityType=capacityIndex;
							break;
						}
					}
					if(capacityType>=0){
						modifyVehicleEdgeCover[serviceEdgeIndex][capacityType]+=1;
						residualNetworkCost[serviceEdgeIndex]+=variableCostPara*edge.duration;
						residualNetworkCost[serviceEdgeIndex]+=averageFixCostForCapacityType[capacityType]*edge.duration;
						break;
					}else{
						amount_copy-=modelData.capacity[modelData.numOfCapacity-1];
						modifyVehicleEdgeCover[serviceEdgeIndex][modelData.numOfCapacity-1]+=1;
						residualNetworkCost[serviceEdgeIndex]+=variableCostPara*edge.duration;
						residualNetworkCost[serviceEdgeIndex]+=averageFixCostForCapacityType[modelData.numOfCapacity-1]*edge.duration;
					}
					
				}
			}else{ // there is vehicle cover
				
				double surplusCapacity=totalVehicleCapacityForServiceEdge[serviceEdgeIndex]-flowEdgeCover[serviceEdgeIndex];
				if(surplusCapacity>amount-0.001){ //no new vehicle edges
					residualNetworkCost[serviceEdgeIndex]+=modelData.beta*edge.duration*amount;
				}else{  // we add new vehicle edges or modify the vehicle type
					
					int capacityType=-1;
					double amount_copy=amount-surplusCapacity;
					//if we add new vehicle edges
					double cost1=0;
					int[] modifyCapacityRecord1=new int[modelData.numOfCapacity];
					while(amount_copy>0.01){
						for(int capacityIndex=0;capacityIndex<modelData.numOfCapacity;capacityIndex++){
							if(modelData.capacity[capacityIndex]>amount_copy-0.001){
								capacityType=capacityIndex;
								break;
							}
						}
						if(capacityType>=0){
							modifyCapacityRecord1[capacityType]+=1;
							cost1+=variableCostPara*edge.duration+averageFixCostForCapacityType[capacityType]*edge.duration;
							break;
						}else{
							amount_copy-=modelData.capacity[modelData.numOfCapacity-1];
							modifyCapacityRecord1[modelData.numOfCapacity-1]+=1;
							cost1+=variableCostPara*edge.duration+averageFixCostForCapacityType[modelData.numOfCapacity-1]*edge.duration;
						}
					}
					
					//if we modify the vehicle type
					double cost2=0;
					int[] modifyCapacityRecord2=new int[modelData.numOfCapacity];
					//Firstly, we check if maximal possible vehicle capacity will satisfy the flow
					int count=0;
					for(int value:vehicleEdgeCover[serviceEdgeIndex]){
						count+=value;
					}
					if(count*modelData.capacity[modelData.numOfCapacity-1]>flowEdgeCover[serviceEdgeIndex]+amount-0.001){
						//we modify the vehicle type
						int currentCapacity=totalVehicleCapacityForServiceEdge[serviceEdgeIndex];
						int maxCapacity=modelData.capacity[modelData.numOfCapacity-1];
						boolean check=false;
						
						for(int capacityIndex=0;capacityIndex<modelData.numOfCapacity-1;capacityIndex++){
							int capacity0=modelData.capacity[capacityIndex];
							
							while(vehicleEdgeCover[serviceEdgeIndex][capacityIndex]+modifyCapacityRecord2[capacityIndex]>0){
								if(currentCapacity+maxCapacity-modelData.capacity[capacityIndex]>flowEdgeCover[serviceEdgeIndex]+amount-0.001){
									for(int capacityIndex1=capacityIndex+1;capacityIndex1<modelData.numOfCapacity;capacityIndex1++){
										int capacity1=modelData.capacity[capacityIndex1];
										if(currentCapacity+capacity1-capacity0>flowEdgeCover[serviceEdgeIndex]+amount-0.001){
											modifyCapacityRecord2[capacityIndex]--;
											modifyCapacityRecord2[capacityIndex1]++;
											currentCapacity=currentCapacity+capacity1-capacity0;
											cost2+=(averageFixCostForCapacityType[capacityIndex1]-averageFixCostForCapacityType[capacityIndex])*edge.duration;
											check=true;
											break;
										}
									}
								}
								
								if(!check){ //change the capacityIndex to largest capacity type
									modifyCapacityRecord2[capacityIndex]--;
									modifyCapacityRecord2[modelData.numOfCapacity-1]++;
									currentCapacity=currentCapacity+maxCapacity-modelData.capacity[capacityIndex];
									cost2+=(averageFixCostForCapacityType[modelData.numOfCapacity-1]-averageFixCostForCapacityType[capacityIndex])*edge.duration;
								}else{
									break;
								}
							}
							if(check) break;
						}
						
					}else{
						cost2=Double.MAX_VALUE;
					}
					
					
					
					
					//compare the cost and decide if we add vehicles or modify vehicle types
					if(cost1<cost2){// add vehicles
						residualNetworkCost[serviceEdgeIndex]=modelData.beta*edge.duration*amount+cost1;
						modifyVehicleEdgeCover[serviceEdgeIndex]=modifyCapacityRecord1;
					}else{
						residualNetworkCost[serviceEdgeIndex]=modelData.beta*edge.duration*amount+cost2;
						modifyVehicleEdgeCover[serviceEdgeIndex]=modifyCapacityRecord2;
					}
					
				}
				
				
			}
			
		}
		
		
		//add surplus value on residual network based on flowEdgeCover
		
		//Firstly, we find the nearest node
		List<Set<Integer>> nearestNodeSet=new ArrayList<>();
		int[][] distance=new int[modelData.numNode][modelData.numNode];
		for(int i=0;i<modelData.numNode;i++){
			for(int j=0;j<modelData.numNode;j++){
				distance[i][j]=Integer.MAX_VALUE;
			}
		}
		for(Service service:modelData.serviceSet){
			distance[service.origin][service.destination]=service.duration;
		}
		for(int terminalIndex=0;terminalIndex<modelData.numNode;terminalIndex++){
			int minDistance=Integer.MAX_VALUE;
			for(int i=0;i<modelData.numNode;i++){
				if(minDistance>distance[i][terminalIndex]){
					minDistance=distance[i][terminalIndex];
				}
			}
			
			Set<Integer> set=new HashSet<>();
			if(minDistance<Integer.MAX_VALUE-10){
				for(int i=0;i<modelData.numNode;i++){
					if(distance[i][terminalIndex]==minDistance){
						set.add(i);
					}
				}
			}
			nearestNodeSet.add(set);
		}
		
		//Next, we modify the cost of network
		for(int edgeIndex=0;edgeIndex<modelData.numServiceArc;edgeIndex++){
			if(flowEdgeCover[edgeIndex]>0.01){
				
				Edge edge=modelData.edgeSet.get(edgeIndex);
				int taili=edge.start;
				//part1
				for(int edgeIndex2:modelData.pointFromEdgeSet.get(taili)){
					Edge edge2=modelData.edgeSet.get(edgeIndex2);
					if(edge2.edgeType==1){
						break;
					}
					residualNetworkCost[edgeIndex2]-=balanceValue1;
					double reduce=balanceValue1/timeZone;
					double reduceCost=balanceValue1-reduce;
					Map<Integer,Double> map=new HashMap<>();
					int time=modelData.edgeSet.get(edgeIndex2).t2;
					for(int timeDiff=1;timeDiff<=timeZone-1;timeDiff++){
						int tempTime=time-timeDiff;
						if(tempTime<0){
							tempTime+=modelData.timePeriod;
						}
						map.put(tempTime, reduceCost);
						reduceCost-=reduce;
					}
					
					for(int edgeIndex3=0;edgeIndex3<modelData.numServiceArc;edgeIndex3++){
						Edge edge3=modelData.edgeSet.get(edgeIndex3);
						if(edge3.serviceIndex==edge2.serviceIndex&&map.containsKey(edge3.t2)){
							residualNetworkCost[edgeIndex3]-=map.get(edge3.t2);
						}
					}
				}
				
				//part2
				int terminalIndex=taili/modelData.timePeriod;
				Set<Integer> set=nearestNodeSet.get(terminalIndex);
				for(int edgeIndex2:modelData.pointFromEdgeSet.get(taili)){
					Edge edge2=modelData.edgeSet.get(edgeIndex2);
					if(set.contains(edge2.u)){
						double reduce=balanceValue2/timeZone;
						double reduceCost=balanceValue2;
						
						Map<Integer,Double> map=new HashMap<>();
						int time=edge2.t1;
						int index=edge2.u*modelData.timePeriod+time;
						map.put(index, reduceCost);
						for(int timeDiff=1;timeDiff<=timeZone-1;timeDiff++){
							int tempTime=time-timeDiff;
							if(tempTime<0){
								tempTime+=modelData.timePeriod;
							}
							reduceCost-=reduce;
							index=edge2.u*modelData.timePeriod+tempTime;
							map.put(index, reduceCost);
						}
						
						for(int edgeIndex3=0;edgeIndex3<modelData.numServiceArc;edgeIndex3++){
							Edge edge3=modelData.edgeSet.get(edgeIndex3);
							if(map.containsKey(edge3.end)){
								residualNetworkCost[edgeIndex3]-=map.get(edge3.end);
							}
						}
						
					}
				}
				
			}
		}
		
		
		return residualNetworkCost;
		
	}

	
 	public double removeEmptyVehicleEdge(Map<Cycle,Map<Integer,Integer>> removeVehicleEdgeRecord,List<Map<Integer,Double>> copyOptXValues,List<Cycle> cycleValues,Map<Cycle,Map<Integer,Integer>> emptyVehicleEdgeRecord) throws Exception{
		double totalRemoveFixCost=0;
		
		Map<Cycle,Double> averageFixCostForPartEdges=new HashMap<>();
		for(Cycle cycle:cycleValues){
			double fixCost=modelData.fixedCost[cycle.associatedPricingProblem.originNodeO][cycle.associatedPricingProblem.capacityTypeS]*cycle.value;
			int totalDistance=0;
			for(int i=0;i<cycle.pattern.length;i++){
				totalDistance+=cycle.pattern[i]*modelData.serviceSet.get(i).duration*cycle.value;
			}
			if(emptyVehicleEdgeRecord.containsKey(cycle)){
				Map<Integer,Integer> map=emptyVehicleEdgeRecord.get(cycle);
				for(int edgeIndex:map.keySet()){
					Edge edge=modelData.edgeSet.get(edgeIndex);
					if(edge.edgeType==0){
						totalDistance-=edge.duration*map.get(edgeIndex);// need to check
					}

				}
			}
			
			//闁哄鏅滈悷鈺呭闯缁岀�榯alDistance闂佸搫鐗嗛ˇ顖濄亹閺屻儲鍤勯柟瀛樺笩缁�锟�0闂佹寧绋戦張顒�锕㈡笟锟藉畷锝夘敍濠靛棗骞嬮梺绋款儐閻喚鑺遍埄鍐枖闁跨喕妫勯埢搴ㄣ��閻炵w闁荤姴顑呴崯顖炲汲閿濆瑙﹂幖娣灪濞堣鈽夐幘璺侯嚝ycle闂佸憡鐟崹顖涚閹烘绀嗛柣妯诲絻缁犵敻鏌曢崱顓熷
			double value=fixCost/totalDistance;
			double distance=0;
			if(totalDistance==0){
				for(int i=0;i<cycle.pattern.length;i++){
					distance+=cycle.pattern[i]*modelData.serviceSet.get(i).duration;
				}
				value=fixCost/distance;
			}
			averageFixCostForPartEdges.put(cycle, value);
			
			
//			if(distance==0){
//				System.out.println(cycle.toString()+": "+averageFixCostForPartEdges.get(cycle));
//				double distance1=0;
//				for(int i=0;i<cycle.pattern.length;i++){
//					distance1+=cycle.pattern[i]*modelData.serviceSet.get(i).duration;
//				}
//				System.out.println("distance1="+distance1);
//				
//				double distance2=0;
//				if(emptyVehicleEdgeRecord.containsKey(cycle)){
//					Map<Integer,Integer> map=emptyVehicleEdgeRecord.get(cycle);
//					for(int edgeIndex:map.keySet()){
//						Edge edge=modelData.edgeSet.get(edgeIndex);
//						if(edge.edgeType==0){
//							distance2-=edge.duration*map.get(edgeIndex);// need to check
//						}
//
//					}
//				}
//				System.out.println("distance2="+distance2);
//				throw new Exception("Infinity error!!!");
//			}
		}
		
		
//		System.out.println("Check averageFixCostForPartEdges:");
//		for(Cycle cycle:averageFixCostForPartEdges.keySet()){
//			System.out.println(cycle.toString()+": "+averageFixCostForPartEdges.get(cycle));
//		}
//		System.out.println();
		
		
		//descend sort
		Comparator<Cycle> averageFixCostComparator= new Comparator<Cycle>() {
			@Override
			public int compare(Cycle cycle1,Cycle cycle2){
				return (int) (averageFixCostForPartEdges.get(cycle2)-averageFixCostForPartEdges.get(cycle1));
			}
		};
		
		List<Cycle> cycleSeq=new ArrayList<>();
		for(Cycle cycle:cycleValues){
			cycleSeq.add(cycle);
		}
		Collections.sort(cycleSeq,averageFixCostComparator);

		
//		System.out.println("Check cycleSeq:");
//		for(Cycle cycle:cycleSeq){
//			System.out.println(cycle.toString());
//		}
//		System.out.println();
		
		
		for(int edgeIndex=0;edgeIndex<modelData.numServiceArc;edgeIndex++){
			Edge edge=modelData.edgeSet.get(edgeIndex);
			
			double flowSum=0;
			for(Map<Integer, Double> map:copyOptXValues){
				if(map.containsKey(edgeIndex)){
					flowSum+=map.get(edgeIndex);
				}
			}
			
			double capacitySum=0;
			for(Cycle cycle:cycleValues){
				if(cycle.edgeIndexSet.contains(edgeIndex)){
					capacitySum+=modelData.capacity[cycle.associatedPricingProblem.capacityTypeS]*cycle.value;
					
					if(emptyVehicleEdgeRecord.containsKey(cycle)){
						Map<Integer,Integer> map=emptyVehicleEdgeRecord.get(cycle);
						if(map.containsKey(edgeIndex)){
							capacitySum-=modelData.capacity[cycle.associatedPricingProblem.capacityTypeS]*map.get(edgeIndex);
						}
					}
					
				}
			}
			
			//if no commodity, we remove all the vehicle edges on this service edge
			if(flowSum<0.001){
				for(Cycle cycle:cycleValues){
					if(cycle.edgeIndexSet.contains(edgeIndex)){
						
						int emptyCount=0;
						//check emptyVehicleEdgeRecord
						if(emptyVehicleEdgeRecord.containsKey(cycle)){
							Map<Integer,Integer> map=emptyVehicleEdgeRecord.get(cycle);
							if(map.containsKey(edgeIndex)){
								emptyCount=map.get(edgeIndex);
							}
						}
						
						if(cycle.value>emptyCount+0.01){// we should remove the edge of this vehicle 
							int value=MathProgrammingUtil.doubleToInt(cycle.value-emptyCount);
							if(removeVehicleEdgeRecord.containsKey(cycle)){
								Map<Integer,Integer> map=removeVehicleEdgeRecord.get(cycle);
								map.put(edgeIndex, value);
							}else{
								Map<Integer,Integer> map=new HashMap<>();
								map.put(edgeIndex, value);
								removeVehicleEdgeRecord.put(cycle, map);
							}
//							System.out.println(edge.toString()+" 0:"+edge.duration*averageFixCostForPartEdges.get(cycle)*value);
							totalRemoveFixCost+=edge.duration*averageFixCostForPartEdges.get(cycle)*value;
							
						}
						
					}
				}
			}else{ // there are flow on this service edge

				while(capacitySum>flowSum+0.01){
					//cost = dataModel.alpha * totalLength / (dataModel.speed * dataModel.drivingTimePerDay)
	                //+ dataModel.fixedCost[pricingProblem.originNodeO][pricingProblem.capacityTypeS];
					
					//Here we pick up the removable cycle with highest average fix cost
					boolean ifFindNewOne=false;
					for(Cycle cycle:cycleSeq){
						int cycleCapacity=modelData.capacity[cycle.associatedPricingProblem.capacityTypeS];
						if(cycle.edgeIndexSet.contains(edgeIndex) && capacitySum-flowSum>cycleCapacity-0.01){
							int emptyCount=0;
							//check emptyVehicleEdgeRecord
							if(emptyVehicleEdgeRecord.containsKey(cycle)){
								Map<Integer,Integer> map=emptyVehicleEdgeRecord.get(cycle);
								if(map.containsKey(edgeIndex)){
									emptyCount+=map.get(edgeIndex);
								}
							}
							
							//check removeVehicleEdgeRecord
							if(removeVehicleEdgeRecord.containsKey(cycle)){
								Map<Integer,Integer> map=removeVehicleEdgeRecord.get(cycle);
								if(map.containsKey(edgeIndex)){
									emptyCount+=map.get(edgeIndex);
								}
							}
							int value=MathProgrammingUtil.doubleToInt(cycle.value-emptyCount);
							if(value<0) System.out.println("Error! value should not be negative");
							
							
							if(value==0){
								continue;
							}else{// we can add the edge of cycle to removeVehicleEdgeRecord
								
								if(removeVehicleEdgeRecord.containsKey(cycle)){
									Map<Integer,Integer> map=removeVehicleEdgeRecord.get(cycle);
									if(map.containsKey(edgeIndex)){
										map.put(edgeIndex, map.get(edgeIndex)+1);
									}else{
										map.put(edgeIndex, 1);
									}

								}else{
									Map<Integer,Integer> map=new HashMap<>();
									map.put(edgeIndex, 1);
									removeVehicleEdgeRecord.put(cycle, map);
								}
								ifFindNewOne=true;
								capacitySum=capacitySum-cycleCapacity;
								
//								System.out.println(edge.toString()+" 1:"+edge.duration*averageFixCostForPartEdges.get(cycle)*value);
								totalRemoveFixCost+=averageFixCostForPartEdges.get(cycle)*edge.duration;
							}

						}
						
						if(ifFindNewOne) break;
					}
					
					if(!ifFindNewOne) break;
				}
				
			}
			
		}
		
		return totalRemoveFixCost;
		
	}
	
	public void identifyEmptyVehicleRun(Map<Cycle,Map<Integer,Integer>> emptyVehicleEdgeRecord,double[] averageFixCostForCapacityType,FeasibleSolution currentSolution){
		Map<Cycle,Double> averageFixCostForAllEdges=new HashMap<>();
		for(Cycle cycle:currentSolution.cycleValues){
			double fixCost=modelData.fixedCost[cycle.associatedPricingProblem.originNodeO][cycle.associatedPricingProblem.capacityTypeS];
			int totalDistance=0;
			for(int i=0;i<cycle.pattern.length;i++){
				totalDistance+=cycle.pattern[i]*modelData.serviceSet.get(i).duration;
			}
			double value=fixCost/totalDistance;
			averageFixCostForAllEdges.put(cycle, value);
		}
		
//		System.out.println("Check averageFixCostForAllEdges:");
//		for(Cycle cycle:averageFixCostForAllEdges.keySet()){
//			System.out.println(cycle.toString()+":"+averageFixCostForAllEdges.get(cycle));
//		}
//		System.out.println();
		

		//descend sort
		Comparator<Cycle> averageFixCostComparator= new Comparator<Cycle>() {
			@Override
			public int compare(Cycle cycle1,Cycle cycle2){
				return (int) (averageFixCostForAllEdges.get(cycle2)-averageFixCostForAllEdges.get(cycle1));
			}
		};
		

		List<Cycle> cycleSeq=new ArrayList<>();
		for(Cycle cycle:currentSolution.cycleValues){
			cycleSeq.add(cycle);
		}

		Collections.sort(cycleSeq,averageFixCostComparator);
		
//		System.out.println("Check cycleSeq:");
//		for(Cycle cycle :cycleSeq){
//			System.out.println(cycle.toString()+" "+averageFixCostForAllEdges.get(cycle));
//		}
//		System.out.println();
		
		
		for(int edgeIndex=0;edgeIndex<modelData.numServiceArc;edgeIndex++){
			
			
			double flowSum=0;
			for(Map<Integer, Double> map:currentSolution.optXValues){
				if(map.containsKey(edgeIndex)){
					flowSum+=map.get(edgeIndex);
				}
			}
			
			double capacitySum=0;
			for(Cycle cycle:currentSolution.cycleValues){
				if(cycle.edgeIndexSet.contains(edgeIndex)){
					capacitySum+=modelData.capacity[cycle.associatedPricingProblem.capacityTypeS]*cycle.value;
				}
			}
			
			
			while(capacitySum>flowSum+0.01){
				//cost = dataModel.alpha * totalLength / (dataModel.speed * dataModel.drivingTimePerDay)
                //+ dataModel.fixedCost[pricingProblem.originNodeO][pricingProblem.capacityTypeS];
				
				//Here we pick up the removable cycle with highest average fix cost
				boolean ifFindNewOne=false;
				for(Cycle cycle:cycleSeq){
					int cycleCapacity=modelData.capacity[cycle.associatedPricingProblem.capacityTypeS];
					if(cycle.edgeIndexSet.contains(edgeIndex) && capacitySum-flowSum>cycleCapacity-0.01){
						if(!emptyVehicleEdgeRecord.containsKey(cycle)){
							Map<Integer,Integer> tempMap=new HashMap<>();
							tempMap.put(edgeIndex, 1);
							emptyVehicleEdgeRecord.put(cycle, tempMap);
							ifFindNewOne=true;
							capacitySum=capacitySum-cycleCapacity;
						}else{
							Map<Integer,Integer> tempMap=emptyVehicleEdgeRecord.get(cycle);
							if(!tempMap.containsKey(edgeIndex)){
								tempMap.put(edgeIndex, 1);
								ifFindNewOne=true;
								capacitySum=capacitySum-cycleCapacity;
							}else{
								if(tempMap.get(edgeIndex)<cycle.value-0.01){
									tempMap.put(edgeIndex, tempMap.get(edgeIndex)+1);
									ifFindNewOne=true;
									capacitySum=capacitySum-cycleCapacity;
								}
							}
						}
					}
					
					if(ifFindNewOne) break;
				}
				
				if(!ifFindNewOne) break;
			}
			
		}
		
//		System.out.println("Check emptyVehicleEdgeRecord:");
//		for(Cycle cycle:emptyVehicleEdgeRecord.keySet()){
//			System.out.println(cycle.toString());
//			Map<Integer,Integer> tempMap=emptyVehicleEdgeRecord.get(cycle);
//			for(int edgeIndex:tempMap.keySet()){
//				System.out.println(modelData.edgeSet.get(edgeIndex).toString()+"="+tempMap.get(edgeIndex));
//			}
//			System.out.println(tempMap);
//		}
		
		
		//for each capacity type, we calculate average fix cost
		int[] sumFixCost=new int[modelData.numOfCapacity];
		int[] sumDistance=new int[modelData.numOfCapacity];
		for(int capacityType=0;capacityType<modelData.numOfCapacity;capacityType++){
			for(Cycle cycle:currentSolution.cycleValues){
				if(capacityType==cycle.associatedPricingProblem.capacityTypeS){
					sumFixCost[capacityType]+=cycle.value*modelData.fixedCost[cycle.associatedPricingProblem.originNodeO][cycle.associatedPricingProblem.capacityTypeS];
					
					double distance=0;
					for(int i=0;i<cycle.pattern.length;i++){
						distance+=cycle.pattern[i]*modelData.serviceSet.get(i).duration;
					}
					distance=distance*cycle.value;
					if(emptyVehicleEdgeRecord.containsKey(cycle)){
						Map<Integer,Integer> tempMap=emptyVehicleEdgeRecord.get(cycle);
						for(int edgeIndex:tempMap.keySet()){
							Edge edge=modelData.edgeSet.get(edgeIndex);
							distance-=tempMap.get(edgeIndex)*edge.duration;
						}
					}
					sumDistance[capacityType]+=distance;				
				}
			}
			averageFixCostForCapacityType[capacityType]=(double)sumFixCost[capacityType]/sumDistance[capacityType];
			if(sumFixCost[capacityType]==0&&sumDistance[capacityType]==0){
				if(capacityType>0){
					averageFixCostForCapacityType[capacityType]=averageFixCostForCapacityType[capacityType-1]*(double)modelData.capacity[capacityType]/modelData.capacity[capacityType-1];
				}else{
					averageFixCostForCapacityType[capacityType]=modelData.fixedCost[0][capacityType]/10;
				}

			}
			
		}
		
//		System.out.println("Check averageFixCostForCapacityType:");
//		System.out.println("sumFixCost="+Arrays.toString(sumFixCost));
//		System.out.println("sumDistance="+Arrays.toString(sumDistance));
//		System.out.println("averageFixCostForCapacityType="+Arrays.toString(averageFixCostForCapacityType));
	}
	
	/**
	 * Here we use BFS to search for furthest extension path for commodity k on edge(i,j)
	 * @param edgeIndex
	 * @param commodityIndex
	 * @param xValues
	 * @return
	 */
	public List<Integer> SourceSearch(int edgeIndex,int commodityIndex,List<Map<Integer, Double>> xValues){
		List<Integer> resultPathRecord=new ArrayList<>();
		Edge keyEdge=modelData.edgeSet.get(edgeIndex);
		int i=keyEdge.start;
		int j=keyEdge.end;
		double commodityAmount=xValues.get(commodityIndex).get(edgeIndex);
		
		int count=0;
		//Backward extension
		Queue<Node> searchQueue=new PriorityQueue<Node>();
		List<Integer> tempList=new ArrayList<>();
		Node rootNode=new Node(i, tempList,count);
		searchQueue.add(rootNode);
		
		while(!searchQueue.isEmpty()){
			Node currentNode=searchQueue.poll();
			Set<Integer> pointFromEdgeSet=modelData.pointFromEdgeSet.get(currentNode.nodeIndex);
			
			
			for(int pointFromEdgeIndex:pointFromEdgeSet){
				Edge edge=modelData.edgeSet.get(pointFromEdgeIndex);
				if(xValues.get(commodityIndex).containsKey(pointFromEdgeIndex)){
					if(xValues.get(commodityIndex).get(pointFromEdgeIndex)>=commodityAmount-0.01){
						//add this pointFromEdge to searchQueue
						int sourceNodeIndex=edge.start;
						List<Integer> pathRecord=new ArrayList<>(currentNode.pathRecord);
						pathRecord.add(pointFromEdgeIndex);
						count++;
						Node newNode=new Node(sourceNodeIndex,pathRecord,count);
						searchQueue.add(newNode);
					}
				}
			}
			
			
			
			// we add the first part path to resultPathRecord based on currentNode
			if(searchQueue.isEmpty()){
				int index=currentNode.pathRecord.size()-1;
				while(index>=0){
					resultPathRecord.add(currentNode.pathRecord.get(index));
					index--;
				}
			}
		}
		
		
		resultPathRecord.add(edgeIndex);
		
		count=0;
		//Forward extension
		searchQueue=new PriorityQueue<Node>();
		tempList=new ArrayList<>();
		rootNode=new Node(j, tempList,count);
		searchQueue.add(rootNode);
		
		while(!searchQueue.isEmpty()){
			Node currentNode=searchQueue.poll();
			Set<Integer> pointToEdgeSet=modelData.pointToEdgeSet.get(currentNode.nodeIndex);
			
			for(int pointToEdgeIndex:pointToEdgeSet){
				Edge edge=modelData.edgeSet.get(pointToEdgeIndex);
				if(xValues.get(commodityIndex).containsKey(pointToEdgeIndex)){
					if(xValues.get(commodityIndex).get(pointToEdgeIndex)>=commodityAmount-0.01){
						//add this pointToEdge to searchQueue
						int sourceNodeIndex=edge.end;
						List<Integer> pathRecord=new ArrayList<>(currentNode.pathRecord);
						pathRecord.add(pointToEdgeIndex);
						count++;
						Node newNode=new Node(sourceNodeIndex,pathRecord,count);
						searchQueue.add(newNode);
					}
				}
			}
			
			
			
			// we add the second part path to resultPathRecord based on currentNode
			if(searchQueue.isEmpty()){
				
				for(int index=0;index<currentNode.pathRecord.size();index++){
					resultPathRecord.add(currentNode.pathRecord.get(index));
				}
				
			}
		}
		
	return resultPathRecord;
	
	}
	
	public String out(Cycle column) {

		Queue<Edge> path = new PriorityQueue<>();

		for (int edgeIndex : column.edgeIndexSet) {
			path.add(modelData.edgeSet.get(edgeIndex));
		}

		StringBuilder pathRecord = new StringBuilder();

		Edge edge = null;
		int size = path.size();
		for (int i = 0; i < size; i++) {

			edge = path.poll();
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

	public static void main(String[] args) throws Exception {    
	    
		long startTime = System.currentTimeMillis();
		
        Properties properties = new Properties();
//      // properties.setProperty("EXPORT_MODEL", "True");
//      // properties.setProperty("MAXTHREADS", "10");
        properties.setProperty("PRECISION", "0.001");
//      properties.setProperty("CUTSENABLED", "false");
        Configuration.readFromFile(properties);
		
//		LocalSearchHeuristicSolver solver = new LocalSearchHeuristicSolver("./data/testset/test0_5_10_10_5.txt", 3,5,3,3,20);	
//		LocalSearchHeuristicSolver solver = new LocalSearchHeuristicSolver("./data/testset/test1_5_10_15_20.txt", 3,5,3,3,60);
		LocalSearchHeuristicSolver solver = new LocalSearchHeuristicSolver("./data/testset/test4_5_15_15_100.txt", 3,5,3,3,60);
//		LocalSearchHeuristicSolver solver = new LocalSearchHeuristicSolver("./data/testset/test12_10_50_30_100A.txt", 3,5,3,10,200);

//		List<FeasibleSolution> solutionList=solver.Initialization();
		List<FeasibleSolution> solutionList=solver.InitializationCGHeuristic();
	    long endTime = System.currentTimeMillis();
		solver.TabuSearch(solutionList.get(0),20,300);
		
		endTime = System.currentTimeMillis();
		System.out.println("total run time=" + (endTime - startTime) + "ms");

//		System.setOut(out0);
//		out0.close();
//		mytxt.close();


		
	}

}
