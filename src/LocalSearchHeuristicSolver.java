import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import model.SNDRC;
import model.SNDRC.Demand;
import model.SNDRC.Edge;
import model.SNDRC.Path;

public class LocalSearchHeuristicSolver {
	SNDRC modelData;
	int r;

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
		
//		System.out.println(Arrays.toString(serviceEdgeImportanceSum));
		
		
		
	}
	
	public static void main(String[] args) throws IOException {
		LocalSearchHeuristicSolver solver=new LocalSearchHeuristicSolver("./data/testset/test0_5_10_10_5.txt", 2);
		solver.Initialization();
	}
	
	
}
