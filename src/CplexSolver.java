import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import model.SNDRC;
import model.SNDRC.Edge;

public class CplexSolver {

	SNDRC dataModel;
	String fileName;
	
	private class Node{
		private int nodeIndex;
		private boolean ifUsed;
		private int pathDuration;
		private int leadEdgeIndex;
	}
	
	
	public CplexSolver(SNDRC dataModel,String fileName) {
		this.dataModel=dataModel;
		this.fileName=fileName;
	}
	
	// generate all possible paths for all nodes as start point
	public GeneratePathFile() {
		PrintWriter out=new PrintWriter(fileName);
		
		for(int originNode=0;originNode<dataModel.numNode;originNode++) {
			
			out.println(originNode);
			
			for(int startTime=0;startTime<dataModel.timePeriod;startTime++) {
				int startNodeIndex=originNode*dataModel.timePeriod+startTime;
				
				Stack<Node> stack=new Stack<>();
				List<Integer> pathEdgeRecord=new ArrayList<>();
				
				Node node=new Node();
				node.nodeIndex=startNodeIndex;
				node.ifUsed=false;
				node.pathDuration=0;
				node.leadEdgeIndex=-1;
				stack.add(node);
				
				while(stack.size()>0) {
					Node currentNode=stack.peek();
					if(!currentNode.ifUsed) {   // need to generate new nodes
						
						if(currentNode.leadEdgeIndex>=0) {
							pathEdgeRecord.add(currentNode.leadEdgeIndex);
						}
						
						
						if(currentNode.nodeIndex!=startNodeIndex) {
							for(int edgeIndex:dataModel.pointToEdgeSet.get(node.nodeIndex)) {
								Edge edge=dataModel.edgeSet.get(edgeIndex);
								if((currentNode.pathDuration+edge.duration<dataModel.timePeriod)||(currentNode.pathDuration+edge.duration==dataModel.timePeriod&&edge.end==startNodeIndex)) {
									Node newNode=new Node();
									newNode.nodeIndex=edge.end;
									newNode.ifUsed=false;
									newNode.leadEdgeIndex=edgeIndex;
									newNode.pathDuration=(int) (currentNode.pathDuration+edge.duration);
									
									stack.add(newNode);
								}
						    }
						}else {   // currentNode.nodeIndex==startNodeIndex
							
						}
						

						
						
					}
				}
				
				
				
				
			}
				
		}
		
	}
	
	
}
