import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import model.SNDRC;
import model.SNDRC.Demand;
import model.SNDRC.Edge;
import model.SNDRC.Path;
import model.SNDRC.Service;
import sun.awt.AWTAccessor.ToolkitAccessor;

public class LocalSearchHeuristicSolver {
	SNDRC modelData;
	int r;
//	List<Set<Integer>> edgesForX;

	public LocalSearchHeuristicSolver(String filename,int r) throws IOException {
		modelData=new SNDRC(filename);
		this.r=r;
		
	}
	
	public void Initialization() {
		
		//We first calculate how many shortest paths service edge (i,j) belong to commodity k's r shortest paths
		int[][] serviceEdgeImportance=new int[modelData.numServiceArc][modelData.numDemand];
		int[] serviceEdgeImportanceSum=new int[modelData.numServiceArc];
		ArrayList<ArrayList<Path>> rPathSet=new ArrayList<>();
		for(int k=0;k<modelData.numDemand;k++) {
			rPathSet.add(modelData.findRShortestPath(r, k));
		}
		
		for(int k=0;k<modelData.numDemand;k++) {
			Demand demand=modelData.demandSet.get(k);
			ArrayList<Path> pathList=rPathSet.get(k);
			for(Path path:pathList) {
				
				for(int i=0;i<path.serviceIndexList.size();i++) {
					
					int serviceIndex=path.serviceIndexList.get(i);
					Set<Integer> timeSet=new HashSet<>();
					int time=path.timeList.get(i);
					
					for(int t=0;t<=demand.duration-path.totalDuration;t++) {
						//timeAvalible+time+t
						int cTime=demand.timeAvailable+time+t;
						cTime=cTime%modelData.timePeriod;
						timeSet.add(cTime);
					}
					
					for(int edgeIndex=0;edgeIndex<modelData.numServiceArc;edgeIndex++) {
						Edge edge=modelData.edgeSet.get(edgeIndex);
						if(edge.serviceIndex==serviceIndex&&timeSet.contains(edge.t1)) {
							serviceEdgeImportance[edgeIndex][k]++;
						}
					}
				}
			}
		}
		
		
		
		for(int edgeIndex=0;edgeIndex<modelData.numServiceArc;edgeIndex++) {
			for(int k=0;k<modelData.numDemand;k++) {
				serviceEdgeImportanceSum[edgeIndex]+=serviceEdgeImportance[edgeIndex][k];
			}
		}
		
		
		// We use dp to calculate subEdgeSet for each commodity
		List<Set<Integer>> edgesForX=new ArrayList<>();
		for(int k=0;k<modelData.numDemand;k++) {
			Set<Integer> set=new HashSet<Integer>();
			edgesForX.add(set);
		}
		
		
		for(int k=0;k<modelData.numDemand;k++) {
			System.out.println();
			System.out.println("k= "+k);
			Set<Integer> edgeSet=edgesForX.get(k);
			Demand demand=modelData.demandSet.get(k);
			
			for(Path path:rPathSet.get(k)) {
				int numOfService=path.serviceIndexList.size();
				int totalTime=demand.duration;
				
				int[][] f=new int[numOfService+1][totalTime+1];
				int[][] record=new int[numOfService+1][totalTime+1];
				for(int i=0;i<=numOfService;i++) {
					for(int j=0;j<=totalTime;j++) {
						record[i][j]=-1;
					}
				}
				
				
				for(int serviceIndex=0;serviceIndex<numOfService;serviceIndex++) {
					Service service=modelData.serviceSet.get(path.serviceIndexList.get(serviceIndex));
					int startTime=path.timeList.get(serviceIndex)+service.duration;
					
					for(int time=startTime;time<=startTime+totalTime-path.totalDuration;time++) {
						int edgeTimeTransfer=demand.timeAvailable+time-service.duration;
						edgeTimeTransfer=edgeTimeTransfer%modelData.timePeriod;
						
						int serviceEdgeIndex=path.serviceIndexList.get(serviceIndex)*modelData.timePeriod+edgeTimeTransfer;
						Edge e=modelData.edgeSet.get(serviceEdgeIndex);
						if(e.serviceIndex!=path.serviceIndexList.get(serviceIndex)||e.t1!=edgeTimeTransfer) System.out.println("Wrong serviceEdgeIndex calculation!!!");
						
						if(time==startTime) {
							f[serviceIndex+1][time]=f[serviceIndex][time-service.duration]+serviceEdgeImportanceSum[serviceEdgeIndex];
							record[serviceIndex+1][time]=serviceEdgeIndex;
						}else {
							int value=f[serviceIndex][time-service.duration]+serviceEdgeImportanceSum[serviceEdgeIndex];
							if(value>f[serviceIndex+1][time-1]) {
								f[serviceIndex+1][time]=value;
								record[serviceIndex+1][time]=serviceEdgeIndex;
							}else {
								f[serviceIndex+1][time]=f[serviceIndex+1][time-1];
							}
						}
						
					}
							
				}
				
				int i=numOfService;
				int j=totalTime;
				List<Integer> resultList=new ArrayList<>();
				while(i>0) {
					
					if(record[i][j]>=0) {
						resultList.add(record[i][j]);
						j=j-modelData.edgeSet.get(record[i][j]).duration;
						i=i-1;
					}else {
						j=j-1;
					}
				}
			
				System.out.println(resultList);
				
				//add edges to edgesForX
				for(int edgeIndex:resultList) {
					edgeSet.add(edgeIndex);
					System.out.println(modelData.edgeSet.get(edgeIndex).toString());
				}
			}
		}
				
		
		
		
	}
	
	public static void main(String[] args) throws IOException {
		LocalSearchHeuristicSolver solver=new LocalSearchHeuristicSolver("./data/testset/test0_5_10_10_5.txt", 3);
		solver.Initialization();
	}
	
	
}
