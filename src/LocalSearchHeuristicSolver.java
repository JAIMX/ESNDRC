import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import bap.BranchAndPriceA;
import bap.bapNodeComparators.NodeBoundbapNodeComparator;
import bap.branching.BranchOnLocalService;
import bap.branching.BranchOnLocalServiceForAllPricingProblems;
import bap.branching.BranchOnServiceEdge;
import cg.Cycle;
import cg.ExactPricingProblemSolver;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import cg.master.SNDRCMasterData;
import cg.master.cuts.StrongInequalityGenerator;
import ch.qos.logback.classic.sift.MDCBasedDiscriminator;
import logger.BapLoggerA;
import model.SNDRC;
import model.SNDRC.Demand;
import model.SNDRC.Edge;
import model.SNDRC.Path;
import model.SNDRC.Service;
import sun.awt.AWTAccessor.ToolkitAccessor;

public class LocalSearchHeuristicSolver {
	SNDRC modelData;
	int r;
	// List<Set<Integer>> edgesForX;

	public LocalSearchHeuristicSolver(String filename, int r) throws IOException {
		modelData = new SNDRC(filename);
		this.r = r;

	}

	public void Initialization() {

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
				Double.MAX_VALUE, 0.6, 0.2, 0.1, 10, 0.0001, 3, false);
		bap.setNodeOrdering(new NodeBoundbapNodeComparator());

		BapLoggerA logger = new BapLoggerA(bap, new File("./output/BAPlogger.log"));

		bap.runBranchAndPrice(System.currentTimeMillis() + 7200000L); // 2 hours

		// Print solution:
		System.out.println("================ Solution ================");
		System.out.println("BAP terminated with objective : " + bap.getObjective());
		System.out.println("Total Number of iterations: " + bap.getTotalNrIterations());
		System.out.println("Total Number of processed nodes: " + bap.getNumberOfProcessedNodes());
		System.out.println("Total Time spent on master problems: " + bap.getMasterSolveTime()
				+ " Total time spent on pricing problems: " + bap.getPricingSolveTime());
		System.out.println("Best bound : " + bap.getBound());

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

				vehicleFixTotalCost += fixCost * columnValue;
				vehicleVarTotalCost += varCost * columnValue;
				// vehicleVarCost.put(column, varCost);
				System.out.println("Fix cost= " + fixCost + " variable cost= " + varCost);

				System.out.println();

			}

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
                    
//                    commodityDowork+=edge.duration*xValues.get(edgeIndex);
                    
//                    commodityFlowIntoTerminal[edge.v][edge.t2]+=xValues.get(edgeIndex);
                }
            }

            commodityCost.add(cost);
//            commodityTotalCost += cost;
        }
        // output x variables
        for (int demand = 0; demand < modelData.numDemand; demand++) {
            for (int edgeIndex : optXValues.get(demand).keySet()) {
                if (optXValues.get(demand).get(edgeIndex) > 0.01) {
                    Edge edge;

                    edge = modelData.edgeSet.get(edgeIndex);


                    
                    if(edge.edgeType==0){
                        System.out.println("x[" + demand + "]:" + edge.u + "," + edge.t1 + "->" + edge.v + "," + edge.t2
                                + "= " + optXValues.get(demand).get(edgeIndex) + " " + edge.duration);
                    }

                }
            }
            System.out.println("total cost= " + commodityCost.get(demand));
            System.out.println();
        }

		bap.close();
		cutHandler.close();
	}

	public String out(Cycle column) {

		Queue<Edge> path = new PriorityQueue<>();

		for (int edgeIndex : column.edgeIndexSet) {
			path.add(modelData.edgeSet.get(edgeIndex));
		}

		StringBuilder pathRecord = new StringBuilder();

		// int count=0;
		// for(Edge edge:path) {
		// count++;
		// pathRecord.append(edge.start);
		//
		// if(count!=column.edgeIndexSet.size()) {
		// pathRecord.append("->");
		// }
		//
		// }

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

	public static void main(String[] args) throws IOException {
		LocalSearchHeuristicSolver solver = new LocalSearchHeuristicSolver("./data/testset/test0_5_10_10_5.txt", 3);
		solver.Initialization();
	}

}
