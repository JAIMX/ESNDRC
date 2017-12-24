package model;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;

import org.jorlib.frameworks.columnGeneration.model.ModelInterface;

public class SNDRC  implements ModelInterface {

	private class Service {
		private int origin;
		private int destination;
		private int LB, UB;
		private int capacity;
		private int duration;
	}

	private class Demand {
		private int origin;
		private int destination;
		private int timeAvailable;
		private int volume;
		private int valueOfTime;
	}

	private class Edge {
		int start, end;
		double duration;
		int u, v, t1, t2;
	}

	public int fleetSize;
	public int numNode, abstractNumNode;
	public int timePeriod;
	public int numService, numDemand;
	public ArrayList<Service> serviceSet;
	public ArrayList<Demand> demandSet;

	public double alpha, speed, drivingTimePerDay, beta, numOfCapacity;
	public double[] fixedCost;
	public int[] capacity;

	// graph parameter
	public ArrayList<Edge> edgeSet;
	public ArrayList<HashSet<Integer>> pointToEdgeSet, pointFromEdgeSet; // not
																			// record
																			// holding
																			// arcs
	public int numServiceArc, numHoldingArc;

	public SNDRC(int readType, String filename) throws IOException {
		if (readType == 1) {
			readData1(filename);
		}
		if (readType == 2) {
			readData2(filename);
		}

		this.graphTransfer();
	}

	private void readData1(String filename) throws IOException {

		alpha = 0;
		speed = 1;
		drivingTimePerDay = 1;
		numOfCapacity = 1;
		beta = 1;
		fixedCost = new double[1];
		fixedCost[0] = 5000;
		capacity = new int[1];
		capacity[0] = 1;

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
			service.LB = in.nextInt();
			service.UB = in.nextInt();
			service.capacity = in.nextInt();
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
			demand.volume = in.nextInt();
			demand.valueOfTime = in.nextInt();
			demandSet.add(demand);
		}

	}

	private void readData2(String filename) {

	}

	private void graphTransfer() {
		abstractNumNode = numNode * timePeriod;
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
				edgeSet.add(newEdge);
				pointToEdgeSet.get(start).add(edgeSet.size() - 1);
				pointFromEdgeSet.get(end).add(edgeSet.size() - 1);
			}

		}
		numServiceArc = edgeSet.size();

	}
	
	@Override
	public String getName(){
		return "ServiceNetworkDesignModel";
	}

	public static void main(String[] args) throws IOException {
		SNDRC test = new SNDRC(1, "./data/Andersen_1.txt");
	}

}
