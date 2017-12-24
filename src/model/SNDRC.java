package model;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Scanner;

import org.jorlib.frameworks.columnGeneration.model.ModelInterface;

public class SNDRC implements ModelInterface {

	private class service {
		private int origin;
		private int destination;
		private int LB, UB;
		private int capacity;
		private double duration;
	}

	private class demand {
		private int origin;
		private int destination;
		private int timeAvailable;
		private int volume;
		private int valueOfTime;
	}

	public int fleetSize;
	public int numNode;
	public int timePeriod;
	public int numService, numDemand;
	public ArrayList<service> serviceSet;
	public ArrayList<demand> demandSet;

	// public final double alpha,speed,drivingTimePerDay,beta,numOfCapacity;
	// public final double[] fixedCost;
	// public final int[] capacity;

	public SNDRC(int readType,String filename) throws IOException {
		if(readType==1) {
			readData1(filename);
		}
		if(readType==2) {
			readData2(filename);
		}
	}
	
	private void readData1(String filename) throws IOException {
		Scanner in=new Scanner(Paths.get(filename));
		
		in.nextLine();
		fleetSize=in.nextInt();
		
		
	}
	
	private void readData2(String filename) {
		
	}

}
