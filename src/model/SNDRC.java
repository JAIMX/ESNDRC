package model;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;

import org.jorlib.frameworks.columnGeneration.model.ModelInterface;

public class SNDRC implements ModelInterface {

	private class Service {
		private int origin;
		private int destination;
		private int LB, UB;
		private int capacity;
		private int duration;
	}

	public class Demand {
		public int origin;
		public int destination;
		public int timeAvailable;
		public int timeDue;
		public int volume;
		public double valueOfTime;
	}

	public class Edge {
		public int start, end;
		public double duration;
		public int u, v, t1, t2;
		public int edgeType;// 0: service arc | 1: holding arc
	}

	public final int fleetSize;
	public final int numNode, abstractNumNode;
	public final int timePeriod;
	public final int numService, numDemand, numOfCapacity;
	public final ArrayList<Service> serviceSet;
	public final ArrayList<Demand> demandSet;
	public final double[][] b;

	public final double alpha, speed, drivingTimePerDay, beta;
	public final double[] fixedCost;
	public final int[] capacity;
	public final int[][] vehicleLimit;
	public final double distanceLimit;
	public final int legLimit;

	// graph parameter
	public final ArrayList<Edge> edgeSet;
	public final ArrayList<HashSet<Integer>> pointToEdgeSet, pointFromEdgeSet; // not
																				// record
																				// holding
																				// arcs
	public final int numServiceArc, numHoldingArc, numArc;

	public SNDRC(String filename) throws IOException {
		// if (readType == 1) {
		// readData1(filename);
		// }
		// if (readType == 2) {
		// readData2(filename);
		// }

		// read data
		Scanner in = new Scanner(Paths.get(filename));

		in.nextLine();
		fleetSize = in.nextInt();
		in.nextLine();
		in.nextLine();
		numNode = in.nextInt();
		in.nextLine();
		in.nextLine();
		timePeriod = in.nextInt();
		in.nextLine();
		in.nextLine();
		numService = in.nextInt();

		in.nextLine();
		in.nextLine();
		serviceSet = new ArrayList<>();
		demandSet = new ArrayList<>();
		for (int i = 0; i < numService; i++) {
			Service service = new Service();
			in.nextInt();
			service.origin = in.nextInt() - 1;
			service.destination = in.nextInt() - 1;
			service.duration = in.nextInt();
			serviceSet.add(service);
		}

		in.nextLine();
		in.nextLine();
		numDemand = in.nextInt();
		in.nextLine();
		in.nextLine();
		for (int i = 0; i < numDemand; i++) {
			Demand demand = new Demand();
			in.nextInt();
			demand.origin = in.nextInt() - 1;
			demand.destination = in.nextInt() - 1;
			demand.timeAvailable = in.nextInt() - 1;
			demand.timeDue = in.nextInt() - 1;
			demand.volume = in.nextInt();
			demandSet.add(demand);
		}

		// read parameter
		in.nextLine();
		in.nextLine();
		alpha = in.nextDouble();
		in.nextLine();
		in.nextLine();
		beta = in.nextDouble();
		in.nextLine();
		in.nextLine();
		speed = in.nextDouble();
		in.nextLine();
		in.nextLine();
		drivingTimePerDay = in.nextDouble();
		in.nextLine();
		in.nextLine();

		numOfCapacity = in.nextInt();
		in.nextLine();
		in.nextLine();
		capacity = new int[numOfCapacity];
		fixedCost = new double[numOfCapacity];
		for (int i = 0; i < numOfCapacity; i++) {
			capacity[i] = in.nextInt();
		}

		in.nextLine();
		in.nextLine();
		for (int i = 0; i < numOfCapacity; i++) {
			fixedCost[i] = in.nextDouble();
		}

		abstractNumNode = numNode * timePeriod;
		// set b[p][i]
		// b = new double[abstractNumNode][numDemand];
		b = new double[numDemand][abstractNumNode];
		for (int i = 0; i < numDemand; i++) {
			Demand demand = demandSet.get(i);
			int start = demand.origin * timePeriod + demand.timeAvailable;
			int end = demand.destination * timePeriod + demand.timeDue;
			// b[start][i] = demand.volume;
			// b[end][i] = -demand.volume;
			b[i][start] = demand.volume;
			b[i][end] = -demand.volume;
		}
		
		//vehicle limit
		in.nextLine();
		in.nextLine();
		vehicleLimit=new int[numOfCapacity][numNode];
		for(int s=0;s<numOfCapacity;s++) {
			for(int o=0;o<numNode;o++) {
				vehicleLimit[s][o]=in.nextInt();
			}
		}
		
		//distance limit
		in.nextLine();
		in.nextLine();
		distanceLimit=in.nextDouble();
		//leg limit
		in.nextLine();
		in.nextLine();
		legLimit=in.nextInt();
		
			
		

		in.close();
		
		
		
		
		///----------------------------------------graph transfer-----------------------------------------///

		// this.graphTransfer();
		edgeSet = new ArrayList<Edge>();
		pointToEdgeSet = new ArrayList<HashSet<Integer>>();
		pointFromEdgeSet = new ArrayList<HashSet<Integer>>();

		for (int i = 0; i < abstractNumNode; i++) {
			HashSet<Integer> templist1 = new HashSet();
			HashSet<Integer> templist2 = new HashSet();
			pointToEdgeSet.add(templist1);
			pointFromEdgeSet.add(templist2);
		}

		// add hat A
		for (int serviceIndex = 0; serviceIndex < numService; serviceIndex++) {
			Service service = serviceSet.get(serviceIndex);
			for (int time = 0; time < timePeriod; time++) {
				int timeEnd = time + service.duration;
				if (timeEnd > timePeriod - 1)
					timeEnd = timeEnd % timePeriod;
				Edge newEdge = new Edge();
				int start = service.origin * timePeriod + time;
				int end = service.destination * timePeriod + timeEnd;
				newEdge.start = start;
				newEdge.end = end;
				newEdge.u = service.origin;
				newEdge.v = service.destination;
				newEdge.t1 = time;
				newEdge.t2 = timeEnd;
				newEdge.duration=service.duration;
				newEdge.edgeType=0;
				edgeSet.add(newEdge);
				pointToEdgeSet.get(start).add(edgeSet.size() - 1);
				pointFromEdgeSet.get(end).add(edgeSet.size() - 1);
				
				//reverse direction arc
				newEdge = new Edge();
				start = service.destination * timePeriod + time;
				end = service.origin * timePeriod + timeEnd;
				newEdge.start = start;
				newEdge.end = end;
				newEdge.u = service.destination;
				newEdge.v = service.origin;
				newEdge.t1 = time;
				newEdge.t2 = timeEnd;
				newEdge.duration=service.duration;
				edgeSet.add(newEdge);
				newEdge.edgeType=0;
				pointToEdgeSet.get(start).add(edgeSet.size() - 1);
				pointFromEdgeSet.get(end).add(edgeSet.size() - 1);
				
			}

		}
		
		numServiceArc = edgeSet.size();
		
		
		// add holding arcs
		for(int localNode=0;localNode<numNode;localNode++) {
			for(int time=0;time<timePeriod;time++) {
				
				Edge newEdge = new Edge();
				int start = localNode*timePeriod+time;
				int end;
				if(time==timePeriod-1) {
					end=localNode*timePeriod;
				}else end=start+1;
				
				newEdge.start = start;
				newEdge.end = end;
				newEdge.u = localNode;
				newEdge.v = localNode;
				newEdge.t1 = time;
				if(time==timePeriod-1) {
					newEdge.t2=0;
				}else newEdge.t2=time+1;
				newEdge.duration=0;
//				newEdge.duration=1;
				newEdge.edgeType=1;
				edgeSet.add(newEdge);
				pointToEdgeSet.get(start).add(edgeSet.size() - 1);
				pointFromEdgeSet.get(end).add(edgeSet.size() - 1);
			}
		}
		

		

		numHoldingArc = abstractNumNode;
		numArc = numServiceArc + numHoldingArc;
		
		

		


	}

	/***
	 * private void readData1(String filename) throws IOException {
	 * 
	 * alpha = 0; speed = 1; drivingTimePerDay = 1; numOfCapacity = 1; beta = 1;
	 * fixedCost = new double[1]; fixedCost[0] = 5000; capacity = new int[1];
	 * capacity[0] = 1;
	 * 
	 * Scanner in = new Scanner(Paths.get(filename));
	 * 
	 * in.nextLine(); fleetSize = in.nextInt(); in.nextLine(); in.nextLine();
	 * numNode = in.nextInt(); in.nextLine(); in.nextLine(); timePeriod =
	 * in.nextInt(); in.nextLine(); in.nextLine(); numService = in.nextInt();
	 * 
	 * in.nextLine(); in.nextLine(); serviceSet = new ArrayList<>(); demandSet = new
	 * ArrayList<>(); for (int i = 0; i < numService; i++) { Service service = new
	 * Service(); in.nextInt(); service.origin = in.nextInt() - 1;
	 * service.destination = in.nextInt() - 1; service.LB = in.nextInt(); service.UB
	 * = in.nextInt(); service.capacity = in.nextInt(); service.duration =
	 * in.nextInt(); serviceSet.add(service); }
	 * 
	 * in.nextLine(); in.nextLine(); numDemand = in.nextInt(); in.nextLine();
	 * in.nextLine(); for (int i = 0; i < numDemand; i++) { Demand demand = new
	 * Demand(); in.nextInt(); demand.origin = in.nextInt() - 1; demand.destination
	 * = in.nextInt() - 1; demand.timeAvailable = in.nextInt() - 1; demand.volume =
	 * in.nextInt(); demand.valueOfTime = in.nextInt(); demandSet.add(demand); }
	 * 
	 * in.close();
	 * 
	 * }
	 * 
	 * private void readData2(String filename) throws IOException {
	 * 
	 * Scanner in = new Scanner(Paths.get(filename));
	 * 
	 * in.nextLine(); fleetSize = in.nextInt(); in.nextLine(); in.nextLine();
	 * numNode = in.nextInt(); in.nextLine(); in.nextLine(); timePeriod =
	 * in.nextInt(); in.nextLine(); in.nextLine(); numService = in.nextInt();
	 * 
	 * in.nextLine(); in.nextLine(); serviceSet = new ArrayList<>(); demandSet = new
	 * ArrayList<>(); for (int i = 0; i < numService; i++) { Service service = new
	 * Service(); in.nextInt(); service.origin = in.nextInt() - 1;
	 * service.destination = in.nextInt() - 1; service.duration = in.nextInt();
	 * serviceSet.add(service); }
	 * 
	 * in.nextLine(); in.nextLine(); numDemand = in.nextInt(); in.nextLine();
	 * in.nextLine(); for (int i = 0; i < numDemand; i++) { Demand demand = new
	 * Demand(); in.nextInt(); demand.origin = in.nextInt() - 1; demand.destination
	 * = in.nextInt() - 1; demand.timeAvailable = in.nextInt() - 1; demand.timeDue =
	 * in.nextInt() - 1; demand.volume = in.nextInt(); demandSet.add(demand); }
	 * 
	 * // read parameter in.nextLine(); in.nextLine(); alpha = in.nextDouble();
	 * in.nextLine(); in.nextLine(); beta = in.nextDouble(); in.nextLine();
	 * in.nextLine(); speed = in.nextDouble(); in.nextLine(); in.nextLine();
	 * drivingTimePerDay = in.nextDouble(); in.nextLine(); in.nextLine();
	 * 
	 * numOfCapacity = in.nextInt(); in.nextLine(); in.nextLine(); capacity = new
	 * int[numOfCapacity]; fixedCost = new double[numOfCapacity]; for (int i = 0; i
	 * < numOfCapacity; i++) { capacity[i] = in.nextInt(); }
	 * 
	 * in.nextLine(); in.nextLine(); for (int i = 0; i < numOfCapacity; i++) {
	 * fixedCost[i] = in.nextDouble(); }
	 * 
	 * in.close(); }
	 * 
	 * private void graphTransfer() {
	 * 
	 * edgeSet = new ArrayList<Edge>(); pointToEdgeSet = new
	 * ArrayList<HashSet<Integer>>(); pointFromEdgeSet = new
	 * ArrayList<HashSet<Integer>>();
	 * 
	 * for (int i = 0; i < abstractNumNode; i++) { HashSet<Integer> templist1 = new
	 * HashSet(); HashSet<Integer> templist2 = new HashSet();
	 * pointToEdgeSet.add(templist1); pointFromEdgeSet.add(templist2); }
	 * 
	 * // add hat A for (int serviceIndex = 0; serviceIndex < numService;
	 * serviceIndex++) { Service service = serviceSet.get(serviceIndex); for (int
	 * time = 0; time < timePeriod; time++) { int timeEnd = time + service.duration;
	 * if (timeEnd > timePeriod - 1) timeEnd = timeEnd % timePeriod; Edge newEdge =
	 * new Edge(); int start = service.origin * timePeriod + time; int end =
	 * service.destination * timePeriod + timeEnd; newEdge.start = start;
	 * newEdge.end = end; newEdge.u = service.origin; newEdge.v =
	 * service.destination; newEdge.t1 = time; newEdge.t2 = timeEnd;
	 * edgeSet.add(newEdge); pointToEdgeSet.get(start).add(edgeSet.size() - 1);
	 * pointFromEdgeSet.get(end).add(edgeSet.size() - 1); }
	 * 
	 * } numServiceArc = edgeSet.size();
	 * 
	 * }
	 ***/
	@Override
	public String getName() {
		return "ServiceNetworkDesignModel";
	}

	public static void main(String[] args) throws IOException {
		SNDRC test = new SNDRC("./data/data0.txt");
	}

}
