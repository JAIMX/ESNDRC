package model;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Currency;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.model.ModelInterface;

import com.google.common.util.concurrent.Service.State;
import com.sun.org.apache.bcel.internal.generic.AASTORE;

import model.SNDRC.Demand;
import model.SNDRC.Edge;

public class SNDRC implements ModelInterface {

    public class Service {
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

    public class Edge implements Comparable<Edge> {
        public int start, end;
        public double duration;
        public int u, v, t1, t2;
        public int edgeType;// 0: service arc | 1: holding arc
        public int serviceIndex;

        public int compareTo(Edge other) {
            return this.t1 - other.t1;
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
    public final ArrayList<HashSet<Integer>> pointToEdgeSet, pointFromEdgeSet; // not
                                                                               // record
                                                                               // holding
                                                                               // arcs
    public final int numServiceArc, numHoldingArc, numArc;

    public List<Set<Integer>> edgesForX;
    public boolean isFeasibleForX;
    
    public ArrayList<Edge> subEdgeSet;


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

            pointToEdgeSet.get(edge.start).add(edgeSet.size() - 1);
            pointFromEdgeSet.get(edge.end).add(edgeSet.size() - 1);
        }

        this.edgesForX = new ArrayList<Set<Integer>>();
        for (int p = 0; p < numDemand; p++) {
            Set<Integer> set = new HashSet<>();
            edgesForX.add(set);
        }
        
        sndrcParent.subEdgeSet=edgeSet;


        // add x variables with edges only needed(dp process)
//        isFeasibleForX = true;
//        for (int p = 0; p < numDemand; p++) {
//
//            boolean ifDestinationAchievable = false;
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
//                            // here we have some sort of bugs, the duration of holding arcs is 0. But that doesn't affact the correctness of the programs
//                            if (edge.duration < durationLimit || (edge.duration == durationLimit
//                                    && edge.end == demand.destination * timePeriod + endTime)) {
//                                edgesForX.get(p).add(edgeIndex);
//                                achieve[edge.end] = true;
//
//                                if (edge.end == demand.destination * timePeriod + endTime) {
//                                    ifDestinationAchievable = true;
//                                }
//
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
//            if (!ifDestinationAchievable) {
//                isFeasibleForX = false;
//                break;
//            }
//            
//            System.out.println("Size of commodity "+p+"= "+edgesForX.get(p).size());
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
        
        
        

    }

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
     * in.nextLine(); in.nextLine(); serviceSet = new ArrayList<>(); demandSet =
     * new ArrayList<>(); for (int i = 0; i < numService; i++) { Service service
     * = new Service(); in.nextInt(); service.origin = in.nextInt() - 1;
     * service.destination = in.nextInt() - 1; service.LB = in.nextInt();
     * service.UB = in.nextInt(); service.capacity = in.nextInt();
     * service.duration = in.nextInt(); serviceSet.add(service); }
     * 
     * in.nextLine(); in.nextLine(); numDemand = in.nextInt(); in.nextLine();
     * in.nextLine(); for (int i = 0; i < numDemand; i++) { Demand demand = new
     * Demand(); in.nextInt(); demand.origin = in.nextInt() - 1;
     * demand.destination = in.nextInt() - 1; demand.timeAvailable =
     * in.nextInt() - 1; demand.volume = in.nextInt(); demand.valueOfTime =
     * in.nextInt(); demandSet.add(demand); }
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
     * in.nextLine(); in.nextLine(); serviceSet = new ArrayList<>(); demandSet =
     * new ArrayList<>(); for (int i = 0; i < numService; i++) { Service service
     * = new Service(); in.nextInt(); service.origin = in.nextInt() - 1;
     * service.destination = in.nextInt() - 1; service.duration = in.nextInt();
     * serviceSet.add(service); }
     * 
     * in.nextLine(); in.nextLine(); numDemand = in.nextInt(); in.nextLine();
     * in.nextLine(); for (int i = 0; i < numDemand; i++) { Demand demand = new
     * Demand(); in.nextInt(); demand.origin = in.nextInt() - 1;
     * demand.destination = in.nextInt() - 1; demand.timeAvailable =
     * in.nextInt() - 1; demand.timeDue = in.nextInt() - 1; demand.volume =
     * in.nextInt(); demandSet.add(demand); }
     * 
     * // read parameter in.nextLine(); in.nextLine(); alpha = in.nextDouble();
     * in.nextLine(); in.nextLine(); beta = in.nextDouble(); in.nextLine();
     * in.nextLine(); speed = in.nextDouble(); in.nextLine(); in.nextLine();
     * drivingTimePerDay = in.nextDouble(); in.nextLine(); in.nextLine();
     * 
     * numOfCapacity = in.nextInt(); in.nextLine(); in.nextLine(); capacity =
     * new int[numOfCapacity]; fixedCost = new double[numOfCapacity]; for (int i
     * = 0; i < numOfCapacity; i++) { capacity[i] = in.nextInt(); }
     * 
     * in.nextLine(); in.nextLine(); for (int i = 0; i < numOfCapacity; i++) {
     * fixedCost[i] = in.nextDouble(); }
     * 
     * in.close(); }
     * 
     * private void graphTransfer() {
     * 
     * edgeSet = new ArrayList<Edge>(); pointToEdgeSet = new ArrayList<HashSet
     * <Integer>>(); pointFromEdgeSet = new ArrayList<HashSet<Integer>>();
     * 
     * for (int i = 0; i < abstractNumNode; i++) { HashSet<Integer> templist1 =
     * new HashSet(); HashSet<Integer> templist2 = new HashSet();
     * pointToEdgeSet.add(templist1); pointFromEdgeSet.add(templist2); }
     * 
     * // add hat A for (int serviceIndex = 0; serviceIndex < numService;
     * serviceIndex++) { Service service = serviceSet.get(serviceIndex); for
     * (int time = 0; time < timePeriod; time++) { int timeEnd = time +
     * service.duration; if (timeEnd > timePeriod - 1) timeEnd = timeEnd %
     * timePeriod; Edge newEdge = new Edge(); int start = service.origin *
     * timePeriod + time; int end = service.destination * timePeriod + timeEnd;
     * newEdge.start = start; newEdge.end = end; newEdge.u = service.origin;
     * newEdge.v = service.destination; newEdge.t1 = time; newEdge.t2 = timeEnd;
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

    public void Output() {
        for (int edgeIndex = 0; edgeIndex < numArc; edgeIndex++) {
            Edge edge = edgeSet.get(edgeIndex);
            System.out.println("edgeIndex=" + edgeIndex + " " + edge.start + "->" + edge.end + " " + edge.u + ","
                    + edge.t1 + "->" + edge.v + "," + edge.t2);
        }
    }

    public static void main(String[] args) throws IOException {
        SNDRC test = new SNDRC("./data/data0.txt");
    }

}
