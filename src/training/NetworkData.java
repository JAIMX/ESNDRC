package training;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import model.SNDRC;

public class NetworkData {
	SNDRC modelData;
	
	int numOfNodes,numOfArcs,totalVolumeOfDemands;
	
	int[] timeSurplus;
	
	int[][] serviceArcFeature;
	
	int[][][] commodityAndArc;
	
	ArrayList<Map<Integer,Double>> lpSolution;
	ArrayList<Set<Integer>> exactSolution;
	
	
	
	public NetworkData(SNDRC modelData){
		this.modelData=modelData;
		timeSurplus=new int[modelData.numDemand];
		serviceArcFeature=new int[modelData.numServiceArc][26];
		commodityAndArc=new int[modelData.numDemand][modelData.numServiceArc][3];
		
		lpSolution=new ArrayList<>();
		for(int k=0;k<modelData.numDemand;k++){
			Map<Integer,Double> temp=new HashMap<>();
			lpSolution.add(temp);
		}
		
		exactSolution=new ArrayList<>();
		for(int k=0;k<modelData.numDemand;k++){
			Set<Integer> temp=new HashSet<>();
			exactSolution.add(temp);
		}
		
	}
	
	public void readData(String filename1,String filename2) throws IOException{
		
		
		//Input part1 data
		//network
		Scanner in = new Scanner(Paths.get(filename1));
		in.nextLine();
		numOfNodes=in.nextInt();
		numOfArcs=in.nextInt();
		totalVolumeOfDemands=in.nextInt();
		
		//commodity
		in.nextLine();
		in.nextLine();
		in.nextLine();
		for(int i=0;i<timeSurplus.length;i++){
			timeSurplus[i]=in.nextInt();
		}
		
		//service arc
		in.nextLine();
		in.nextLine();
		in.nextLine();
		for(int serviceEdgeIndex=0;serviceEdgeIndex<modelData.numServiceArc;serviceEdgeIndex++){
			int edgeIndex=in.nextInt();
			in.nextLine();
			for(int i=0;i<26;i++){
				serviceArcFeature[edgeIndex][i]=in.nextInt();
			}
		}
		
		
		//commodity & service arc
		in.nextLine();
		in.nextLine();
		in.nextLine();
		for(int k=0;k<modelData.numDemand;k++){
			in.nextLine();
			for(int serviceEdgeIndex=0;serviceEdgeIndex<modelData.numServiceArc;serviceEdgeIndex++){
				in.nextLine();
				for(int i=0;i<3;i++){
					commodityAndArc[k][serviceEdgeIndex][i]=in.nextInt();
				}
				in.nextLine();

			}
			
		}
		
		in.close();
		
		
		
		//Input part2 data
		in = new Scanner(Paths.get(filename2));
		
		//lp information
		in.nextLine();
		for(int k=0;k<modelData.numDemand;k++){
			in.nextLine();
			while(in.hasNext("\n")){
				int index=in.nextInt();
				double value=in.nextDouble();
				lpSolution.get(k).put(index, value);
			}
		}
		
		
		
		//Final solution
		in.nextLine();
		in.nextLine();
		in.nextLine();
		for(int k=0;k<modelData.numDemand;k++){
			in.nextLine();
			while(in.hasNext("\n")){
				int index=in.nextInt();
				double value=in.nextDouble();
				exactSolution.get(k).add(index);
			}
		}
		
		in.close();
		
		
		
	}
	
	public static void main(String[] args) throws IOException {
		SNDRC sndrc=new SNDRC("./data/testset/test1_5_10_15_20.txt");
		NetworkData trainingData=new NetworkData(sndrc);
		trainingData.readData("./learningData/result/0_0.txt", "./learningData/result/0_1.txt");
	}
	
}
