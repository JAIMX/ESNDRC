package model;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.util.*;

import javax.xml.crypto.dsig.CanonicalizationMethod;

import org.jgrapht.alg.EulerianCircuit;
import org.jorlib.frameworks.columnGeneration.model.ModelInterface;

import com.google.common.util.concurrent.Service.State;
import com.sun.org.apache.bcel.internal.generic.AASTORE;

import model.SNDRC.Demand;
import model.SNDRC.Edge;

public class SNDRC implements ModelInterface {

    public class Service {
        public int origin;
        public int destination;
//        private int LB, UB;
//        private int capacity;
        public int duration;
        
        public String toString(){
        	String string=origin+"->"+destination+":"+duration;
        	return string;
        }
        
    }

    public class Demand {
        public int origin;
        public int destination;
        public int timeAvailable;
        public int timeDue;
        public int volume;
        public double valueOfTime;
        public int duration;
        
        public Demand(){
        	if(timeAvailable<timeDue){
        		duration=timeDue-timeAvailable;
        	}else{
        		duration=timeDue-timeAvailable+timePeriod;
        	}
        }
        public String toString(){
        	return origin+"->"+destination;
        }
        
    }

    public class Edge implements Comparable<Edge> {
        public int start, end;
        public int duration;
        public int u, v, t1, t2;
        public int edgeType;// 0: service arc | 1: holding arc
        public int serviceIndex;

        public int compareTo(Edge other) {
            return this.t1 - other.t1;
        }
        
        public String toString(){
        	return "("+u+","+t1+")->("+v+","+t2+")";
        }

    }
    
    public class Path implements Comparable<Path>{
    	public List<Integer>  serviceIndexList;
    	public List<Integer> timeList;
    	public int totalDuration;
    	public int origin,destination;
    	
    	public Path(int origin,int destination,List serviceIndexList){
    		this.serviceIndexList=serviceIndexList;
    		this.origin=origin;
    		this.destination=destination;
    		
    		timeList=new ArrayList<Integer>();
    		totalDuration=0;
    		for(int serviceIndex:this.serviceIndexList){
    			timeList.add(totalDuration);
    			Service service=serviceSet.get(serviceIndex);
    			totalDuration+=service.duration;
    		}
    		
    	}
    	
    	
    	public String toString(){
    		String string="";
    		for(int serviceIndex:serviceIndexList){
    			Service service=serviceSet.get(serviceIndex);
    			string=string+" "+service.toString();
    		}
    		
    		return string;
    	}
    	
    	public int compareTo(Path other){
    		return this.totalDuration-other.totalDuration;
    	}
    	
    	public Path merge(List<Integer> serviceIndexSequence){
    		ArrayList<Integer> list=new ArrayList<>();
    		for(int serviceIndex:serviceIndexSequence){
    			list.add(serviceIndex);
    		}
    		for(int serviceIndex:this.serviceIndexList){
    			list.add(serviceIndex);
    		}
    		int origin;
    		if(serviceIndexSequence.size()>0){
        		origin=serviceSet.get(serviceIndexSequence.get(0)).origin;
    		}else{
    			origin=serviceSet.get(this.serviceIndexList.get(0)).origin;
    		}

    		int destination=this.destination;
    		Path newPath=new Path(origin, destination, list);
    		return newPath;
    	}
    	
    }

    public final int fleetSize;
    public final int numNode, abstractNumNode;
    public final int timePeriod;
    public final int numService, numDemand, numOfCapacity;
    public final ArrayList<Service> serviceSet;
    public final ArrayList<Demand> demandSet;
    public final double[][] b;

    public final double alpha, speed, drivingTimePerDay, beta;
    public final double[][] fixedCost;
    public final int[] capacity;
    public final int[][] vehicleLimit;
    public final double distanceLimit;
    public final int legLimit;

    // graph parameter
    public final ArrayList<Edge> edgeSet;
    public final ArrayList<HashSet<Integer>> pointToEdgeSet, pointFromEdgeSet; // record
                                                                               // holding
                                                                               // arcs
    public final int numServiceArc, numHoldingArc, numArc;

    public List<Set<Integer>> edgesForX;
    public boolean isFeasibleForX;
    
    public ArrayList<Edge> subEdgeSet;
    public Map<Integer,Integer> edgeSetIndexMap; //edge index of subEdgeSet to edge index of edgeSet

    private String filename;
    private ArrayList<Set<Integer>> pointToService;
    private boolean[] ifCanBeUsed;
    
    public double[] flowCover;
    public double totalFlowCost;

    public SNDRC(SNDRC sndrcParent, Set<Integer> serviceEdgeSet) {
    	
        this.fleetSize = sndrcParent.fleetSize;
        this.numNode = sndrcParent.numNode;
        this.abstractNumNode = sndrcParent.abstractNumNode;
        this.timePeriod = sndrcParent.timePeriod;
        this.numService = sndrcParent.numService;
        this.numDemand = sndrcParent.numDemand;
        this.numOfCapacity = sndrcParent.numOfCapacity;
        this.serviceSet = sndrcParent.serviceSet;
        this.demandSet = sndrcParent.demandSet;
        this.b = sndrcParent.b;

        this.alpha = sndrcParent.alpha;
        this.speed = sndrcParent.speed;
        this.drivingTimePerDay = sndrcParent.drivingTimePerDay;
        this.beta = sndrcParent.beta;
        this.fixedCost = sndrcParent.fixedCost;
        this.capacity = sndrcParent.capacity;
        this.vehicleLimit = sndrcParent.vehicleLimit;
        this.distanceLimit = sndrcParent.distanceLimit;
        this.legLimit = sndrcParent.legLimit;

        edgeSet = new ArrayList<Edge>();
        pointToEdgeSet = new ArrayList<HashSet<Integer>>();
        pointFromEdgeSet = new ArrayList<HashSet<Integer>>();

        for (int i = 0; i < abstractNumNode; i++) {
            HashSet<Integer> templist1 = new HashSet();
            HashSet<Integer> templist2 = new HashSet();
            pointToEdgeSet.add(templist1);
            pointFromEdgeSet.add(templist2);
        }

        // for service edge
        for (int edgeIndex = 0; edgeIndex < sndrcParent.numServiceArc; edgeIndex++) {
            if (serviceEdgeSet.contains(edgeIndex)) {
                Edge edge = sndrcParent.edgeSet.get(edgeIndex);
                edgeSet.add(edge);
                sndrcParent.edgeSetIndexMap.put(edgeSet.size()-1, edgeIndex);
                
                pointToEdgeSet.get(edge.start).add(edgeSet.size() - 1);
                pointFromEdgeSet.get(edge.end).add(edgeSet.size() - 1);
            }
        }

        numServiceArc = edgeSet.size();
        numHoldingArc = sndrcParent.numHoldingArc;
        numArc = numServiceArc + numHoldingArc;

        // for holding edge
        for (int holdingEdgeIndex = 0; holdingEdgeIndex < sndrcParent.numHoldingArc; holdingEdgeIndex++) {
            int tempIndex = sndrcParent.numServiceArc + holdingEdgeIndex;
            Edge edge = sndrcParent.edgeSet.get(tempIndex);
            edgeSet.add(edge);
            sndrcParent.edgeSetIndexMap.put(edgeSet.size()-1, tempIndex);

            pointToEdgeSet.get(edge.start).add(edgeSet.size() - 1);
            pointFromEdgeSet.get(edge.end).add(edgeSet.size() - 1);
        }

        this.edgesForX = new ArrayList<Set<Integer>>();
        for (int p = 0; p < numDemand; p++) {
            Set<Integer> set = new HashSet<>();
            edgesForX.add(set);
        }
        
        sndrcParent.subEdgeSet=edgeSet;


        
  
///----------------------------------------------reduce size of x variables--------------------------------------------/// 
        
        // add x variables with edges only needed(dp process)
        isFeasibleForX = true;
        for (int p = 0; p < numDemand; p++) {

            boolean ifDestinationAchievable = false;
            boolean[] achieveForward = new boolean[abstractNumNode];
            for (int i = 0; i < achieveForward.length; i++) {
                achieveForward[i] = false;
            }


            Demand demand = demandSet.get(p);
            
            //1st step: forward
            int originNodeIndex = demand.origin * timePeriod + demand.timeAvailable;
            int destinationNodeIndex=demand.destination*timePeriod+demand.timeDue;
            int startTime = demand.timeAvailable;
            int endTime = demand.timeDue;

            int durationLimit;
            achieveForward[originNodeIndex] = true;

            if (endTime > startTime) {
                durationLimit = endTime - startTime;
            } else {
                durationLimit = endTime - startTime + timePeriod;
            }

            int timeDuration = durationLimit;

            for (int t = 0; t < timeDuration; t++) {
                int currentTime = t + startTime;
                currentTime = currentTime % timePeriod;

                for (int localNode = 0; localNode < numNode; localNode++) {
                    int currentNodeIndex = localNode * timePeriod + currentTime;

                    if (achieveForward[currentNodeIndex]) {
                        for (int edgeIndex : pointToEdgeSet.get(currentNodeIndex)) {
                            Edge edge = edgeSet.get(edgeIndex);

                            // here we have some sort of bugs, the duration of holding arcs is 0. But that doesn't affact the correctness of the programs
                            if (edge.duration < durationLimit || (edge.duration == durationLimit
                                    && edge.end == destinationNodeIndex)) {
                                achieveForward[edge.end] = true;

                                if (edge.end == destinationNodeIndex) {
                                    ifDestinationAchievable = true;
                                }

                            }
                        }
                    }

                }

                durationLimit--;

            }

            if (!ifDestinationAchievable) {
                isFeasibleForX = false;
                break;
            }
            
            
            //2st step: backward
            boolean[] achieveBackward = new boolean[abstractNumNode];
            for (int i = 0; i < achieveBackward.length; i++) {
                achieveBackward[i] = false;
            }
            
            achieveBackward[destinationNodeIndex]=true;
            durationLimit=timeDuration;
            
            for(int t=0;t<timeDuration;t++){
                int currentTime=endTime-t;
                if(currentTime<0){
                    currentTime+=timePeriod;
                }
                
                for(int localNode=0;localNode<numNode;localNode++){
                    int currentNodeIndex=localNode*timePeriod+currentTime;
                    
                    if(achieveBackward[currentNodeIndex]){
                        for(int edgeIndex:pointFromEdgeSet.get(currentNodeIndex)){
                            Edge edge=edgeSet.get(edgeIndex);
                            
                            // here we have some sort of bugs, the duration of holding arcs is 0. But that doesn't affact the correctness of the programs
                            if(edge.duration<durationLimit||(edge.duration==durationLimit&&edge.start==originNodeIndex)){
                                achieveBackward[edge.start]=true;
                            }
                        }
                    }
                }
                
                durationLimit--;
            }
            
            
            
            
            //step3:merge the common points
            Set<Integer> commonPointSet=new HashSet<>();
            durationLimit=timeDuration;
            
            for(int t=0;t<timeDuration;t++){
                int currentTime=t+startTime;
                currentTime=currentTime%timePeriod;
                
                if(currentTime!=startTime&&currentTime!=endTime){
                    for(int localNode=0;localNode<numNode;localNode++){
                        int currentNodeIndex=localNode*timePeriod+currentTime;
                        
                        
                        if(achieveForward[currentNodeIndex]&&achieveBackward[currentNodeIndex]){
                            commonPointSet.add(currentNodeIndex);
                        }
                    }
                }
            }
            
            
            commonPointSet.add(originNodeIndex);
            commonPointSet.add(destinationNodeIndex);
            
            
            
            
            for(int edgeIndex=0;edgeIndex<numArc;edgeIndex++){
                Edge edge=edgeSet.get(edgeIndex);
                
                if(commonPointSet.contains(edge.start)&&commonPointSet.contains(edge.end)){
                    edgesForX.get(p).add(edgeIndex);
                }
            }
            
//            System.out.println("Size of commodity "+p+"= "+edgesForX.get(p).size());
            

            
        }
        
        
        
        
        
        
        
        
 ///----------------------------------------------reduce size of x variables--------------------------------------------///      
        
        
        

    }

    public SNDRC(String filename) throws IOException {
        // if (readType == 1) {
        // readData1(filename);
        // }
        // if (readType == 2) {
        // readData2(filename);
        // }

    	this.filename=filename;

    	
    	
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
        fixedCost = new double[numNode][numOfCapacity];
        for (int i = 0; i < numOfCapacity; i++) {
            capacity[i] = in.nextInt();
        }

        in.nextLine();
        in.nextLine();

        for (int city = 0; city < numNode; city++) {
            for (int i = 0; i < numOfCapacity; i++) {
                fixedCost[city][i] = in.nextDouble();
            }
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

        // vehicle limit
        in.nextLine();
        in.nextLine();
        vehicleLimit = new int[numOfCapacity][numNode];
        for (int s = 0; s < numOfCapacity; s++) {
            for (int o = 0; o < numNode; o++) {
//                System.out.println(s+" "+o);
                vehicleLimit[s][o] = in.nextInt();
//                System.out.println(vehicleLimit[s][o]);
            }
        }

        // distance limit
        in.nextLine();
        in.nextLine();
        distanceLimit = in.nextDouble();
        // leg limit
        in.nextLine();
        in.nextLine();
        legLimit = in.nextInt();

        in.close();

        /// ----------------------------------------graph
        /// transfer-----------------------------------------///

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
                newEdge.duration = service.duration;
                newEdge.edgeType = 0;
                newEdge.serviceIndex = serviceIndex;

                edgeSet.add(newEdge);
                pointToEdgeSet.get(start).add(edgeSet.size() - 1);
                pointFromEdgeSet.get(end).add(edgeSet.size() - 1);

                // //reverse direction arc
                // newEdge = new Edge();
                // start = service.destination * timePeriod + time;
                // end = service.origin * timePeriod + timeEnd;
                // newEdge.start = start;
                // newEdge.end = end;
                // newEdge.u = service.destination;
                // newEdge.v = service.origin;
                // newEdge.t1 = time;
                // newEdge.t2 = timeEnd;
                // newEdge.duration=service.duration;
                // edgeSet.add(newEdge);
                // newEdge.edgeType=0;
                // pointToEdgeSet.get(start).add(edgeSet.size() - 1);
                // pointFromEdgeSet.get(end).add(edgeSet.size() - 1);

            }

        }

        numServiceArc = edgeSet.size();

        // add holding arcs
        for (int localNode = 0; localNode < numNode; localNode++) {
            for (int time = 0; time < timePeriod; time++) {

                Edge newEdge = new Edge();
                int start = localNode * timePeriod + time;
                int end;
                if (time == timePeriod - 1) {
                    end = localNode * timePeriod;
                } else
                    end = start + 1;

                newEdge.start = start;
                newEdge.end = end;
                newEdge.u = localNode;
                newEdge.v = localNode;
                newEdge.t1 = time;
                if (time == timePeriod - 1) {
                    newEdge.t2 = 0;
                } else
                    newEdge.t2 = time + 1;
                newEdge.duration = 0;
                // newEdge.duration=1;
                newEdge.edgeType = 1;
                newEdge.serviceIndex = -1;
                edgeSet.add(newEdge);
                pointToEdgeSet.get(start).add(edgeSet.size() - 1);
                pointFromEdgeSet.get(end).add(edgeSet.size() - 1);
            }
        }

        numHoldingArc = abstractNumNode;
        numArc = numServiceArc + numHoldingArc;
        
    	this.pointToService=new ArrayList<>();
    	HashSet<Integer> tempSet;
    	for(int i=0;i<numNode;i++){
    		tempSet=new HashSet<>();
    		this.pointToService.add(tempSet);
    	}
  
    	for(int i=0;i<numService;i++){
    		Service service=serviceSet.get(i);
    		pointToService.get(service.origin).add(i);
    	}
    	
    	
        System.out.println("number of service arcs=" + numServiceArc);
        System.out.println("number of holding arcs=" + numHoldingArc);
        System.out.println();
        

        this.edgesForX = new ArrayList<Set<Integer>>();
        for (int p = 0; p < numDemand; p++) {
            Set<Integer> set = new HashSet<>();
            edgesForX.add(set);
        }

        // add x variables with edges only needed(dp process)
//        for (int p = 0; p < numDemand; p++) {
//            boolean[] achieve = new boolean[abstractNumNode];
//            for (int i = 0; i < achieve.length; i++) {
//                achieve[i] = false;
//            }
//
//            Demand demand = demandSet.get(p);
//            int originNodeIndex = demand.origin * timePeriod + demand.timeAvailable;
//            int startTime = demand.timeAvailable;
//            int endTime = demand.timeDue;
//            int durationLimit;
//            achieve[originNodeIndex] = true;
//
//            if (endTime > startTime) {
//                durationLimit = endTime - startTime;
//            } else {
//                durationLimit = endTime - startTime + timePeriod;
//            }
//
//            int timeDuration = durationLimit;
//
//            for (int t = 0; t < timeDuration; t++) {
//                int currentTime = t + startTime;
//                currentTime = currentTime % timePeriod;
//
//                for (int localNode = 0; localNode < numNode; localNode++) {
//                    int currentNodeIndex = localNode * timePeriod + currentTime;
//
//                    if (achieve[currentNodeIndex]) {
//                        for (int edgeIndex : pointToEdgeSet.get(currentNodeIndex)) {
//                            Edge edge = edgeSet.get(edgeIndex);
//
//                            if (edge.duration < durationLimit || (edge.duration == durationLimit
//                                    && edge.end == demand.destination * timePeriod + endTime)) {
//                                edgesForX.get(p).add(edgeIndex);
//                                achieve[edge.end] = true;
//                            }
//                        }
//                    }
//
//                }
//
//                durationLimit--;
//
//            }
//            
//            
////            System.out.println("Size of commodity "+p+"= "+edgesForX.get(p).size());
//
//        }
        
        
///----------------------------------------------reduce size of x variables--------------------------------------------/// 
        
        // add x variables with edges only needed(dp process)
        isFeasibleForX = true;
        for (int p = 0; p < numDemand; p++) {

            boolean ifDestinationAchievable = false;
            boolean[] achieveForward = new boolean[abstractNumNode];
            for (int i = 0; i < achieveForward.length; i++) {
                achieveForward[i] = false;
            }


            Demand demand = demandSet.get(p);
            
            //1st step: forward
            int originNodeIndex = demand.origin * timePeriod + demand.timeAvailable;
            int destinationNodeIndex=demand.destination*timePeriod+demand.timeDue;
            int startTime = demand.timeAvailable;
            int endTime = demand.timeDue;

            int durationLimit;
            achieveForward[originNodeIndex] = true;

            if (endTime > startTime) {
                durationLimit = endTime - startTime;
            } else {
                durationLimit = endTime - startTime + timePeriod;
            }

            int timeDuration = durationLimit;

            for (int t = 0; t < timeDuration; t++) {
                int currentTime = t + startTime;
                currentTime = currentTime % timePeriod;

                for (int localNode = 0; localNode < numNode; localNode++) {
                    int currentNodeIndex = localNode * timePeriod + currentTime;

                    if (achieveForward[currentNodeIndex]) {
                        for (int edgeIndex : pointToEdgeSet.get(currentNodeIndex)) {
                            Edge edge = edgeSet.get(edgeIndex);

                            // here we have some sort of bugs, the duration of holding arcs is 0. But that doesn't affact the correctness of the programs
                            if (edge.duration < durationLimit || (edge.duration == durationLimit
                                    && edge.end == destinationNodeIndex)) {
                                achieveForward[edge.end] = true;

                                if (edge.end == destinationNodeIndex) {
                                    ifDestinationAchievable = true;
                                }

                            }
                        }
                    }

                }

                durationLimit--;

            }

            if (!ifDestinationAchievable) {
                isFeasibleForX = false;
                break;
            }
            
            
            
            //2st step: backward
            boolean[] achieveBackward = new boolean[abstractNumNode];
            for (int i = 0; i < achieveBackward.length; i++) {
                achieveBackward[i] = false;
            }
            
            achieveBackward[destinationNodeIndex]=true;
            durationLimit=timeDuration;
            
            for(int t=0;t<timeDuration;t++){
                int currentTime=endTime-t;
                if(currentTime<0){
                    currentTime+=timePeriod;
                }
                
                for(int localNode=0;localNode<numNode;localNode++){
                    int currentNodeIndex=localNode*timePeriod+currentTime;
                    
                    if(achieveBackward[currentNodeIndex]){
                        for(int edgeIndex:pointFromEdgeSet.get(currentNodeIndex)){
                            Edge edge=edgeSet.get(edgeIndex);
                            
                            // here we have some sort of bugs, the duration of holding arcs is 0. But that doesn't affact the correctness of the programs
                            if(edge.duration<durationLimit||(edge.duration==durationLimit&&edge.start==originNodeIndex)){
                                achieveBackward[edge.start]=true;
                            }
                        }
                    }
                }
                
                durationLimit--;
            }
            
            
            
            //step3:merge the common points
            Set<Integer> commonPointSet=new HashSet<>();
            durationLimit=timeDuration;
            
            for(int t=0;t<timeDuration;t++){
                int currentTime=t+startTime;
                currentTime=currentTime%timePeriod;
                
                if(currentTime!=startTime&&currentTime!=endTime){
                    for(int localNode=0;localNode<numNode;localNode++){
                        int currentNodeIndex=localNode*timePeriod+currentTime;
                        
                        
                        if(achieveForward[currentNodeIndex]&&achieveBackward[currentNodeIndex]){
                            commonPointSet.add(currentNodeIndex);
                        }
                    }
                }
            }
            
            
            commonPointSet.add(originNodeIndex);
            commonPointSet.add(destinationNodeIndex);
            
            
            
            
            for(int edgeIndex=0;edgeIndex<numArc;edgeIndex++){
                Edge edge=edgeSet.get(edgeIndex);
                
                if(commonPointSet.contains(edge.start)&&commonPointSet.contains(edge.end)){
                    edgesForX.get(p).add(edgeIndex);
                }
            }
            
//            System.out.println("Size of commodity "+p+"= "+edgesForX.get(p).size());
            

            
        }
        
        
        
        
        
        
        
        
 ///----------------------------------------------reduce size of x variables--------------------------------------------/// 
        

        this.edgeSetIndexMap=new HashMap<>();
    }

    public void ModifyEdgesForX(List<Set<Integer>> edgesForX) {
    	this.edgesForX=edgesForX;
    }
    
    @Override
    public String getName() {
        return "ServiceNetworkDesignModel";
    }

    public void Output() {
        for (int edgeIndex = 0; edgeIndex < numArc; edgeIndex++) {
            Edge edge = edgeSet.get(edgeIndex);
            System.out.println("edgeIndex=" + edgeIndex + " " + edge.start + "->" + edge.end + " " + edge.u + ","
                    + edge.t1 + "->" + edge.v + "," + edge.t2);
        }
    }

    public void outputFeature(String filename) throws FileNotFoundException{
    	PrintWriter out=new PrintWriter(filename);
    	
    	
    	//Network related features
    	out.println("Network:");
    	
    	int totalSupply=0;
    	for(Demand demand:demandSet){
    		totalSupply+=demand.volume;
    	}
    	out.println(abstractNumNode+" "+numArc+" "+totalSupply);
    	
    	
    	
    	ArrayList<ArrayList<Path>> rPathSet=new ArrayList<>();
    	for(int k=0;k<numDemand;k++){
    		rPathSet.add(findRShortestPath(3, k));
    		
//        	ArrayList<Path> tempPathSet=rPathSet.get(rPathSet.size()-1);
//        	out.println();
//        	out.println("commodity "+k+":"+demandSet.get(k).toString());
//        	for(Path path:tempPathSet){
//        		out.println(path.toString());
//        	}
    	}
//    	out.println();
    	out.println();

    	
    	out.println("Commodity:");
    	int[] difference=new int[numDemand];
    	for(int k=0;k<numDemand;k++){
        	Demand demand=demandSet.get(k);
        	int duration=demand.timeDue-demand.timeAvailable;
        	if(duration<=0){
        		duration+=timePeriod;
        	}
        	int diff=duration-rPathSet.get(k).get(0).totalDuration;
        	difference[k]=diff;
        	out.print(diff+" ");
    	}
    	out.println();
    	out.println();
    	
    	
    	
    	//Calculate supply and demand information for each node
    	ArrayList<Set<Integer>> supplyNodeInfo,demandNodeInfo;
    	supplyNodeInfo=new ArrayList<>();
    	Set<Integer> tempSet=new HashSet<>();
    	for(int i=0;i<abstractNumNode;i++){
    		tempSet=new HashSet<>();
    		supplyNodeInfo.add(tempSet);
    	}
    	demandNodeInfo=new ArrayList<>();
    	for(int i=0;i<abstractNumNode;i++){
    		tempSet=new HashSet<>();
    		demandNodeInfo.add(tempSet);
    	}
 
    	for(int k=0;k<numDemand;k++){
    		Demand demand=demandSet.get(k);
    		//supply
    		int node=demand.origin;
    		for(int time=demand.timeAvailable;time<=demand.timeAvailable+difference[k];time++){
    			int t=time%timePeriod;
    			int nodeIndex=node*timePeriod+t;
    			supplyNodeInfo.get(nodeIndex).add(k);
    		}
    		//demand
    		node=demand.destination;
    		for(int time=demand.timeDue;time>=demand.timeDue-difference[k];time--){
    			int t=time;
    			if(t<0){
    				t+=timePeriod;
    			}
    			int nodeIndex=node*timePeriod+t;
    			demandNodeInfo.get(nodeIndex).add(k);
    		}
    	}   	
    	
    	out.println("Service Arc:");
    	for(int edgeIndex=0;edgeIndex<numServiceArc;edgeIndex++){
    		Edge edge=edgeSet.get(edgeIndex);
    		out.println(edgeIndex+" "+edge.toString());
    		
    		int count=0;
    		for(int k=0;k<numDemand;k++){
    			ArrayList<Path> pathSet=rPathSet.get(k);
    			Demand demand=demandSet.get(k);
    			boolean check=false;
    			for(Path path:pathSet){
    				for(int i=0;i<path.serviceIndexList.size();i++){
    					int serviceIndex=path.serviceIndexList.get(i);
    					
    					if(edge.serviceIndex==serviceIndex){
    						int time=path.timeList.get(i);
    						for(int t=0;t<=demand.duration-path.totalDuration;t++){
    							//timeAvalible+time+t
    							int cTime=demand.timeAvailable+time+t;
    							cTime=cTime%timePeriod;
    							if(edge.t1==cTime){
    								check=true;
    								break;
    							}
    						}
    					}
    					if(check) break;
    					
    				}
    				if(check) break;
    			}
    			if(check) count++;
    		}
    		out.print(count);
    		
    		
    		
    		out.print(" "+edge.duration);
    		
    		count=0;
    		for(int k:supplyNodeInfo.get(edge.start)){
    			Demand demand=demandSet.get(k);
    			count+=demand.volume;
    		}
    		out.print(" "+count);
    		count=0;
    		for(int k:demandNodeInfo.get(edge.start)){
    			Demand demand=demandSet.get(k);
    			count+=demand.volume;
    		}
    		out.print(" "+count);
  		
    		count=0;
    		for(int k:supplyNodeInfo.get(edge.end)){
    			Demand demand=demandSet.get(k);
    			count+=demand.volume;
    		}
    		out.print(" "+count);
    		count=0;
    		for(int k:demandNodeInfo.get(edge.end)){
    			Demand demand=demandSet.get(k);
    			count+=demand.volume;
    		}
    		out.print(" "+count);
    		
    		
    		
    		
    		
    		
    		int a1,a2,a3,a4,a5,a6,a7,a8;
    		count=0;
    		a1=0;
    		for(int index:pointToEdgeSet.get(edge.start)){
    			Edge tempEdge=edgeSet.get(index);
    			for(int k:supplyNodeInfo.get(tempEdge.end)){
    				count+=demandSet.get(k).volume;
    			}
    			if(supplyNodeInfo.get(tempEdge.end).size()>0){
    				a1++;
    			}
    			
    		}
    		out.print(" "+count);
    		count=0;
    		a2=0;
    		for(int index:pointToEdgeSet.get(edge.start)){
    			Edge tempEdge=edgeSet.get(index);
    			for(int k:demandNodeInfo.get(tempEdge.end)){
    				count+=demandSet.get(k).volume;
    			}
    			if(demandNodeInfo.get(tempEdge.end).size()>0){
    				a2++;
    			}
    		}
    		out.print(" "+count);
    		
    		count=0;
    		a3=0;
    		for(int index:pointFromEdgeSet.get(edge.start)){
    			Edge tempEdge=edgeSet.get(index);
    			for(int k:supplyNodeInfo.get(tempEdge.start)){
    				count+=demandSet.get(k).volume;
    			}
    			if(supplyNodeInfo.get(tempEdge.start).size()>0){
    				a3++;
    			}
    		}
    		out.print(" "+count);
    		count=0;
    		a4=0;
    		for(int index:pointFromEdgeSet.get(edge.start)){
    			Edge tempEdge=edgeSet.get(index);
    			for(int k:demandNodeInfo.get(tempEdge.start)){
    				count+=demandSet.get(k).volume;
    			}
    			if(demandNodeInfo.get(tempEdge.start).size()>0){
    				a4++;
    			}
    		}
    		out.print(" "+count);
    		
    		
    		count=0;
    		a5=0;
    		for(int index:pointToEdgeSet.get(edge.end)){
    			Edge tempEdge=edgeSet.get(index);
    			for(int k:supplyNodeInfo.get(tempEdge.end)){
    				count+=demandSet.get(k).volume;
    			}
    			if(supplyNodeInfo.get(tempEdge.end).size()>0){
    				a5++;
    			}
    		}
    		out.print(" "+count);
    		count=0;
    		a6=0;
    		for(int index:pointToEdgeSet.get(edge.end)){
    			Edge tempEdge=edgeSet.get(index);
    			for(int k:demandNodeInfo.get(tempEdge.end)){
    				count+=demandSet.get(k).volume;
    			}
    			if(demandNodeInfo.get(tempEdge.end).size()>0){
    				a6++;
    			}
    		}
    		out.print(" "+count);
    		
    		count=0;
    		a7=0;
    		for(int index:pointFromEdgeSet.get(edge.end)){
    			Edge tempEdge=edgeSet.get(index);
    			for(int k:supplyNodeInfo.get(tempEdge.start)){
    				count+=demandSet.get(k).volume;
    			}
    			if(supplyNodeInfo.get(tempEdge.start).size()>0){
    				a7++;
    			}
    		}
    		out.print(" "+count);
    		count=0;
    		a8=0;
    		for(int index:pointFromEdgeSet.get(edge.end)){
    			Edge tempEdge=edgeSet.get(index);
    			for(int k:demandNodeInfo.get(tempEdge.start)){
    				count+=demandSet.get(k).volume;
    			}
    			if(demandNodeInfo.get(tempEdge.start).size()>0){
    				a8++;
    			}
    		}
    		out.print(" "+count);
    		
    		
    		
    		
    		
    		
    		
    		
    		//degree information
    		out.print(" "+pointToEdgeSet.get(edge.start).size());
    		out.print(" "+a1);
    		out.print(" "+a2);
    		out.print(" "+pointFromEdgeSet.get(edge.start).size());
    		out.print(" "+a3);
    		out.print(" "+a4);
    		out.print(" "+pointToEdgeSet.get(edge.end).size());
    		out.print(" "+a5);
    		out.print(" "+a6);
    		out.print(" "+pointFromEdgeSet.get(edge.end).size());
    		out.print(" "+a7);
    		out.print(" "+a8);
    		

    		
    		
    		out.println();
    	}
    	
		out.println();
		out.println("Commodity k & Service arc:");
		for(int k=0;k<numDemand;k++){
			Demand demand=demandSet.get(k);
			out.println(k);
			for(int edgeIndex=0;edgeIndex<numServiceArc;edgeIndex++){
				Edge edge=edgeSet.get(edgeIndex);
				out.println(edgeIndex+" "+edge.toString());
				
				int nodeIndex=edge.start;
				if(supplyNodeInfo.get(nodeIndex).contains(k)){
					out.print(1);
				}else{
					out.print(0);
				}
				
				nodeIndex=edge.end;
				if(demandNodeInfo.get(nodeIndex).contains(k)){
					out.print(" "+1);
				}else{
					out.print(" "+0);
				}
				
				
				
				int count=0;
				for(Path path:rPathSet.get(k)){
					boolean check=false;
					
    				for(int i=0;i<path.serviceIndexList.size();i++){
    					int serviceIndex=path.serviceIndexList.get(i);
    					
    					if(edge.serviceIndex==serviceIndex){
    						int time=path.timeList.get(i);
    						for(int t=0;t<=demand.duration-path.totalDuration;t++){
    							//timeAvalible+time+t
    							int cTime=demand.timeAvailable+time+t;
    							cTime=cTime%timePeriod;
    							if(edge.t1==cTime){
    								check=true;
    								break;
    							}
    						}
    					}
    					if(check) break;
    					
    				}
    				if(check) count++;
    				
				}
				out.println(" "+count);
				
				
			}

		}
    		
    	
    	
    	
    	out.close();
    }
    
    /*
     * Based on original network. there are not definitely r paths in the output list when the length of path exceeds
     * time window of certain commodity
     */
    public ArrayList<Path> findRShortestPath(int r,int demandIndex){
    	
    	ArrayList<Path> pathListA=new ArrayList<>();
    	Queue<Path> pathListB=new PriorityQueue<>();
    	Demand demand=demandSet.get(demandIndex);
    	int destination=demand.destination;
    	int duration=demand.timeDue-demand.timeAvailable;
 
    	if(duration<=0){
    		duration+=timePeriod;
    	}
       	
    	ifCanBeUsed=new boolean[numService];
    	for(int i=0;i<numService;i++){
    		ifCanBeUsed[i]=true;
    	}
    	Set<Integer> forbiddenNodeSet=new HashSet<>();
    	
    	Path path0=BellmanFordSP(demand.origin, demand.destination,forbiddenNodeSet);
    	pathListA.add(path0);
    	
    	for(int k=2;k<=r;k++){
    		
    		if(pathListA.size()>=r||pathListA.get(pathListA.size()-1).totalDuration>duration){
    			break;
    		}
    		
    		ifCanBeUsed=new boolean[numService];
    		for(int i=0;i<numService;i++){
    			ifCanBeUsed[i]=true;
    		}
    		
    		Path path=pathListA.get(pathListA.size()-1);
    		List<Integer> serviceIndexSequence=new ArrayList<Integer>();
    		forbiddenNodeSet=new HashSet<>();
    		
    		for(int pathIndex=0;pathIndex<path.serviceIndexList.size();pathIndex++){
    			Service service=serviceSet.get(path.serviceIndexList.get(pathIndex));
    			forbiddenNodeSet.add(service.origin);

    			
    			if(serviceIndexSequence.size()==0){
    				for(Path tempPath:pathListA){
    					int serviceIndex=tempPath.serviceIndexList.get(0);
    					ifCanBeUsed[serviceIndex]=false;
    				}
    			}else{
    				for(Path tempPath:pathListA){
    					boolean check=true;
    					if(serviceIndexSequence.size()<tempPath.serviceIndexList.size()){
        					for(int index=0;index<serviceIndexSequence.size();index++){
        						if(tempPath.serviceIndexList.get(index)!=serviceIndexSequence.get(index)){
        							check=false;
        							break;
        						}
        					}
    					}else{
    						check=false;
    					}
    					
    					if(check){
    						ifCanBeUsed[tempPath.serviceIndexList.get(serviceIndexSequence.size())]=false;
    					}

    					
    				}
    			}
    			
    			Path pathS=BellmanFordSP(service.origin, destination, forbiddenNodeSet);
    			if(pathS!=null){
        			Path newPath=pathS.merge(serviceIndexSequence);
        			pathListB.add(newPath);
    			}
    			
    			serviceIndexSequence.add(path.serviceIndexList.get(pathIndex));
    		}
    		
    		
    		Path tempPath;
    		int length;
    		if(pathListB.size()>0){
        		tempPath=pathListB.poll();
        		length=tempPath.totalDuration;
        		pathListA.add(tempPath);
        		
        		while(pathListB.size()>0 && pathListB.peek().totalDuration==length){
        			tempPath=pathListB.poll();
        			pathListA.add(tempPath);
        		}
        		
    		}else{
    			break;
    		}

    		

    	}
    	
    	
    	if(pathListA.size()>r){
    		for(int index=pathListA.size()-1;index>r-1;index--){
    			pathListA.remove(index);
    		}
    	}
    	
    	Path tempPath=pathListA.get(pathListA.size()-1);
    	while(tempPath.totalDuration>duration){
    		pathListA.remove(pathListA.size()-1);
    		tempPath=pathListA.get(pathListA.size()-1);
    	}
    	
    	
    	
    	return pathListA;
    	
    }
    
    
    public Path BellmanFordSP(int origin,int destination,Set<Integer> forbiddenNodeSet){
    	int[] distTo=new int[numNode];
    	int[] edgeTo=new int[numNode];
    	boolean[] onQ=new boolean[numNode];
    	ArrayDeque<Integer> queue=new ArrayDeque();
    	
    	for(int i=0;i<numNode;i++){
    		distTo[i]=Integer.MAX_VALUE;
    		onQ[i]=false;
    		edgeTo[i]=-1;
    	}
    	
    	distTo[origin]=0;
    	onQ[origin]=true;
    	queue.add(origin);
    	
    	while(!queue.isEmpty()){
    		int v=queue.poll();
    		onQ[v]=false;
    		
    		Set<Integer> serviceIndexSet=pointToService.get(v);
    		for(int serviceIndex:serviceIndexSet){
    			
    			Service service=serviceSet.get(serviceIndex);
    			int w=service.destination;
    			

    			if(ifCanBeUsed[serviceIndex]&&!forbiddenNodeSet.contains(w)){
        			if(distTo[w]>distTo[v]+service.duration){
        				distTo[w]=distTo[v]+service.duration;
        				edgeTo[w]=serviceIndex;
        				
        				if(!onQ[w]){
        					queue.add(w);
        					onQ[w]=true;
        				}
        			}
    			}

    		}
    	}
    	
    	
    	//create the shortest path
    	if(edgeTo[destination]<0){
    		return null;
    	}else{
        	Stack<Integer> pathRecord=new Stack<>();
        	for(int serviceIndex=edgeTo[destination];serviceIndex>=0;serviceIndex=edgeTo[serviceSet.get(serviceIndex).origin]){
        		pathRecord.push(serviceIndex);
        	}
        	ArrayList<Integer> serviceIndexList=new ArrayList<>();
        	while(!pathRecord.isEmpty()){
        		serviceIndexList.add(pathRecord.pop());
        	}
        	
        	Path path=new Path(origin, destination, serviceIndexList);
        	return path;
    	}

    	
    }
    
    public static void main(String[] args) throws IOException {
        SNDRC test = new SNDRC("./data/testset/test0_5_10_10_5.txt");
        test.outputFeature("./featureSet/test0_5_10_10_5.txt");
        
    }

}
