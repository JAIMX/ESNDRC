package training;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import com.opencsv.CSVWriter;

import model.SNDRC;
import model.SNDRC.Edge;

public class CreateTrainingData {
	List<NetworkData> networkDataList;

	public CreateTrainingData(List<NetworkData> networkDataList,String filename) throws IOException{
		this.networkDataList=networkDataList;
		
		CSVWriter writer=new CSVWriter(new FileWriter(filename));
		String[] featureName="NumOfNodes,NumOfArcs,VolumeOfDemands,timeSurplus,Arc0,Arc1,Arc2,Arc3,Arc4,Arc5,Arc6,Arc7,Arc8,Arc9,Arc10,Arc11,Arc12,Arc13,Arc14,Arc15,Arc16,Arc17,Arc18,Arc19,Arc20,Arc21,Arc22,Arc23,Arc24,Arc25,iIfOrigin,jIfDestination,rSPCount,iflpCover,lpPercent,timeDiff,y".split(",");
		writer.writeNext(featureName);
		
		for(NetworkData networkData:networkDataList){
			SNDRC modelData=networkData.modelData;
			for(int k=0;k<modelData.numDemand;k++){
				for(int serviceEdgeIndex=0;serviceEdgeIndex<modelData.numServiceArc;serviceEdgeIndex++){
					
					String string=networkData.numOfNodes+","+networkData.numOfArcs+","+networkData.totalVolumeOfDemands+","+networkData.timeSurplus[k];
					for(int i=0;i<networkData.serviceArcFeature[0].length;i++){
						if(i>=2&&i<=13){
							double value=(double)networkData.serviceArcFeature[serviceEdgeIndex][i]/networkData.totalVolumeOfDemands;
							string+=","+String.format("%.2f", value);
						}else string+=","+networkData.serviceArcFeature[serviceEdgeIndex][i];
					}
					
					for(int i=0;i<networkData.commodityAndArc[0][0].length;i++){
						string+=","+networkData.commodityAndArc[k][serviceEdgeIndex][i];
					}
					
					int iflpCover=-1;
					if(networkData.lpSolution.get(k).containsKey(serviceEdgeIndex)){
						iflpCover=1;
					}else iflpCover=0;
					string+=","+iflpCover;
					
					if(iflpCover==0){
						string+=",0.00";
					}else{
						double volume=(double)networkData.lpSolution.get(k).get(serviceEdgeIndex)/modelData.demandSet.get(k).volume;
						string+=","+String.format("%.2f", volume);
					}
					
					
					//calculate timeDiff
					Set<Integer> set=networkData.lpSolution.get(k).keySet();
					Edge edge0=modelData.edgeSet.get(serviceEdgeIndex);
					int minRecord=modelData.timePeriod;
					for(int edgeIndex:set){
						Edge edge=modelData.edgeSet.get(edgeIndex);
						
						if(edge.serviceIndex==edge0.serviceIndex){
							minRecord=Math.min(minRecord, Math.abs(edge0.t1-edge.t1));
						}
					}
					
					string+=","+minRecord;
					
					//lable
					if(networkData.exactSolution.get(k).contains(serviceEdgeIndex)){
						string+=","+1;
					}else string+=","+0;

					
					
					writer.writeNext(string.split(","));
				}
			}
		}
		
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		
//		Scanner in = new Scanner(Paths.get("./learningData/result/nameMap.out"));
//		Map<String,Integer> map=new HashMap<>();
//		int count=0;
//		while(in.hasNextLine()) {
//			String string=in.nextLine();
//			map.put(string, count);
//			count++;
//		}
//		in.close();
//		
//		
//    	String path="./learningData/test";
//    	File file=new File(path);
//    	File[] fs=file.listFiles();
//    	
//    	List<NetworkData> networkDataList=new ArrayList<>();
//    	for(File f:fs){
//    		if(!f.isDirectory()&&!f.isHidden()){
//        		String fileName=f.toString().substring(20);
//        		SNDRC sndrc=new SNDRC("./learningData/test/"+fileName);
//        		NetworkData trainingData=new NetworkData(sndrc);
//        		count=map.get(fileName);
//        		String filename1="./learningData/result/"+count+"_0.txt";
//        		String filename2="./learningData/result/"+count+"_1.txt";
//        		trainingData.readData(filename1, filename2);
//        		networkDataList.add(trainingData);
//    		}
//    	}
//    	CreateTrainingData test=new CreateTrainingData(networkDataList, "./test.csv");

		
	}
}
