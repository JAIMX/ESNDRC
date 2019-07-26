import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Properties;
import java.util.Queue;
import java.util.Scanner;
import java.util.Set;
import java.util.Stack;
import org.jorlib.frameworks.columnGeneration.util.Configuration;

import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import model.SNDRC;
import model.SNDRC.Demand;
import model.SNDRC.Edge;

public class CplexSolver {

    SNDRC dataModel;
    String fileName;

    private class Node {
        private int nodeIndex;
        private boolean ifUsed;
        // private int pathDuration;
        private int powerLeft;
        private int leadEdgeIndex;
    }

    private class Path {
        private int originNode;
        private int capacityType;
        private Set<Integer> edgeIndexSet;
        private Set<Integer> ifChargeSet;
    }

    public CplexSolver(SNDRC dataModel, String fileName) {
        this.dataModel = dataModel;
        this.fileName = fileName;
    }

    public void Solve() throws IloException, IOException {

        IloCplex cplex = new IloCplex();

        FileOutputStream outputStream = new FileOutputStream(new File("./output/cplexOut.txt"));

        // cplex.setOut(outputStream);

        cplex.setParam(IloCplex.IntParam.Threads, 4);
        cplex.setParam(IloCplex.Param.Simplex.Tolerances.Markowitz, 0.1);
        cplex.setParam(IloCplex.DoubleParam.TiLim, 7200); // 2 hours
        List<Map<Integer, IloNumVar>> x; // map:edgeIndex, x variable
        Map<Path, IloNumVar> pathVarMap = new HashMap<>();

        // Define variables x
        x = new ArrayList<Map<Integer, IloNumVar>>();
        for (int p = 0; p < dataModel.numDemand; p++) {
            Map<Integer, IloNumVar> initialX = new HashMap<Integer, IloNumVar>();
            x.add(initialX);
        }

        // add x variables
        for (int p = 0; p < dataModel.numDemand; p++) {
            for (int edgeIndex : dataModel.edgesForX.get(p)) {
                Edge edge = dataModel.edgeSet.get(edgeIndex);
                IloNumVar varX = cplex.numVar(0, dataModel.demandSet.get(p).volume,
                        "x" + p + "," + edge.start + "," + edge.end);
                x.get(p).put(edgeIndex, varX);
            }
        }

        // Define the objective
        /**
         * Here we assume the cost of edge AT is 0
         */
        IloLinearNumExpr exprObj = cplex.linearNumExpr();

        for (int p = 0; p < dataModel.numDemand; p++) {
            Map<Integer, IloNumVar> map = x.get(p);
            for (int edgeIndex : map.keySet()) {
                exprObj.addTerm(dataModel.beta * dataModel.edgeSet.get(edgeIndex).duration, map.get(edgeIndex));
            }
        }

        IloObjective obj = cplex.addMinimize(exprObj);

        // Define flowBalanceConstraints
        IloRange[][] flowBalanceConstraints = new IloRange[dataModel.numDemand][dataModel.abstractNumNode];

        IloLinearNumExpr expr = cplex.linearNumExpr();
        for (int p = 0; p < dataModel.numDemand; p++) {
            Map<Integer, IloNumVar> map = x.get(p);

            for (int i = 0; i < dataModel.abstractNumNode; i++) {
                expr.clear();
                // edges which point from i
                for (int edgeIndex : dataModel.pointToEdgeSet.get(i)) {
                    if (map.containsKey(edgeIndex)) {
                        expr.addTerm(1, map.get(edgeIndex));
                    }
                }

                // edges which point to i
                for (int edgeIndex : dataModel.pointFromEdgeSet.get(i)) {
                    if (map.containsKey(edgeIndex)) {
                        expr.addTerm(-1, map.get(edgeIndex));
                    }
                }
                flowBalanceConstraints[p][i] = cplex.addEq(dataModel.b[p][i], expr);

            }
        }

        // Define weakForcingConstraints
        IloRange[] weakForcingConstraints = new IloRange[dataModel.numServiceArc];
        for (int arcIndex = 0; arcIndex < dataModel.numServiceArc; arcIndex++) {
            expr.clear();
            for (int p = 0; p < dataModel.numDemand; p++) {
                if (x.get(p).containsKey(arcIndex)) {
                    expr.addTerm(1, x.get(p).get(arcIndex));
                }
            }

            weakForcingConstraints[arcIndex] = cplex.addGe(0, expr);
        }

        // Define resourceBoundConstraints
        IloRange[][] resourceBoundConstraints = new IloRange[dataModel.numOfCapacity][dataModel.numNode];
        for (int s = 0; s < dataModel.numOfCapacity; s++) {
            for (int o = 0; o < dataModel.numNode; o++) {
                // expr.clear();
                // expr.addTerm(1, q[s][o]);
                expr.clear();
                // resourceBoundConstraints[s][o] =
                // cplex.addEq(dataModel.vehicleLimit[s][o], 0);
                // resourceBoundConstraints[s][o]=cplex.range(0,dataModel.vehicleLimit[s][o]);
                resourceBoundConstraints[s][o] = cplex.addRange(0, dataModel.vehicleLimit[s][o]);
            }
        }

        // Define storeBoundConstraints
        IloRange[][] storeBoundConstraints = new IloRange[dataModel.numNode][dataModel.timePeriod];
        for (int l = 0; l < dataModel.numNode; l++) {
            for (int t = 0; t < dataModel.timePeriod; t++) {
                expr.clear();
                int nodeIndex = l * dataModel.timePeriod + t;
                // search holding arc index
                int holdingEdgeIndex = -1;
                for (int i = dataModel.numServiceArc; i < dataModel.numArc; i++) {
                    Edge edge = dataModel.edgeSet.get(i);
                    if (edge.start == nodeIndex) {
                        holdingEdgeIndex = i;
                        break;
                    }
                }

                for (int k = 0; k < dataModel.numDemand; k++) {
                    Demand demand = dataModel.demandSet.get(k);
                    if (demand.destination != l) {
                        if (x.get(k).containsKey(holdingEdgeIndex)) {
                            expr.addTerm(1, x.get(k).get(holdingEdgeIndex));
                        }
                    }
                }

                storeBoundConstraints[l][t] = cplex.addGe(dataModel.storeLimit[l], expr);
            }
        }

        // Define chargeBoundConstraints
        IloRange[][] chargeBoundConstraints = new IloRange[dataModel.numNode][dataModel.timePeriod];

        for (int l = 0; l < dataModel.numNode; l++) {
            for (int t = 0; t < dataModel.timePeriod; t++) {
                expr.clear();
                chargeBoundConstraints[l][t] = cplex.addGe(dataModel.chargeLimit[l], expr);
            }
        }

        // add all columns
        Scanner in = new Scanner(Paths.get(fileName));

        int count2 = 0;

        for (int originNode = 0; originNode < dataModel.numNode; originNode++) {
            // System.out.println("originNode="+originNode);
            int count = 0;

            in.nextLine();

            Set<Integer> pathEdgeSet;
            Set<Integer> chargeEdgeSet;
            String path = in.nextLine();

            while (path.charAt(0) != 'e') {
                // System.out.println(path);

                String[] stringSet = path.split(" ");
                pathEdgeSet = new HashSet<Integer>();
                chargeEdgeSet = new HashSet<Integer>();
                for (int i = 0; i < stringSet.length; i++) {
                    int value = Integer.parseInt(stringSet[i]);
                    if (value >= 0) {
                        pathEdgeSet.add(value);
                    } else {
                        pathEdgeSet.add(-value);
                        chargeEdgeSet.add(-value);
                    }
                }

                // System.out.println(pathEdgeSet.toString());

                int pathLength = 0;
                for (int edgeIndex : pathEdgeSet) {
                    pathLength += dataModel.edgeSet.get(edgeIndex).duration;
                }

                for (int capacityType = 0; capacityType < dataModel.numOfCapacity; capacityType++) {

                    double cost = 0;
                    cost += dataModel.alpha * pathLength / (dataModel.speed * dataModel.drivingTimePerDay)
                            + dataModel.fixedCost[originNode][capacityType];
                    cost += dataModel.chargeObjPara * chargeEdgeSet.size();

                    // Register column with objective
                    IloColumn iloColumn = cplex.column(obj, cost);

                    // weak forcing constraints
                    for (int edgeIndex : pathEdgeSet) {
                        if (dataModel.edgeSet.get(edgeIndex).edgeType == 0) {
                            iloColumn = iloColumn.and(
                                    cplex.column(weakForcingConstraints[edgeIndex], -dataModel.capacity[capacityType]));
                        }
                    }

                    // resource bound constraints
                    iloColumn = iloColumn.and(cplex.column(resourceBoundConstraints[capacityType][originNode], 1));

                    // charge bound constraints
                    for (int chargeEdgeIndex : chargeEdgeSet) {
                        int chargeNodeIndex = dataModel.edgeSet.get(chargeEdgeIndex).start;
                        int l = chargeNodeIndex / dataModel.timePeriod;
                        int t = chargeNodeIndex % dataModel.timePeriod;
                        iloColumn = iloColumn.and(cplex.column(chargeBoundConstraints[l][t], 1));
                    }
                    // cplex.exportModel("check.lp");

                    // Create the variable and store it
                    // IloNumVar var =cplex.intVar(iloColumn, 0,
                    // dataModel.vehicleLimit[capacityType][originNode],"z_" +
                    // capacityType + ","+originNode+","+count);
                    IloNumVar var = cplex.intVar(iloColumn, 0, Integer.MAX_VALUE,
                            "z_" + capacityType + "," + originNode + "," + count);
                    count2++;

                    Path newPath = new Path();
                    newPath.capacityType = capacityType;
                    newPath.originNode = originNode;
                    newPath.edgeIndexSet = pathEdgeSet;
                    newPath.ifChargeSet = chargeEdgeSet;

                    pathVarMap.put(newPath, var);
                    count++;
                }

                path = in.nextLine();
            }

            // path.charAt(0)=='-'

        }
        in.close();

        System.out.println("count2= " + count2);

        cplex.solve();
        System.out.println("optimal objective= " + cplex.getObjValue());
        // cplex.exportModel("check.lp");

        // ----------------------------output--------------------------------//
        for (int demand = 0; demand < dataModel.numDemand; demand++) {
            for (int edgeIndex = 0; edgeIndex < dataModel.numArc; edgeIndex++) {

                if (x.get(demand).containsKey(edgeIndex)) {
                    if (cplex.getValue(x.get(demand).get(edgeIndex)) > 0.01) {
                        Edge edge = dataModel.edgeSet.get(edgeIndex);
                        System.out.println("x[" + demand + "]:" + edge.start + "->" + edge.end + "= "
                                + cplex.getValue(x.get(demand).get(edgeIndex)));
                    }
                }

            }

            System.out.println();
        }

        for (Path path : pathVarMap.keySet()) {
            double value = cplex.getValue(pathVarMap.get(path));
            if (Math.abs(value) > 0.01) {

                // ------------------------------------check---------------------------------------//
                Queue<Edge> outPath = new PriorityQueue<>();
                Set<Edge> chargeEdgeSet = new HashSet<>();
                for (int edgeIndex : path.edgeIndexSet) {
                    outPath.add(dataModel.edgeSet.get(edgeIndex));
                    if (path.ifChargeSet.contains(edgeIndex)) {
                        chargeEdgeSet.add(dataModel.edgeSet.get(edgeIndex));
                    }
                }

                System.out.println();

                StringBuilder pathRecord = new StringBuilder();

                Edge edge = null;
                int size = outPath.size();
                for (int i = 0; i < size; i++) {

                    edge = outPath.poll();
                    pathRecord.append('(');
                    pathRecord.append(edge.u);
                    pathRecord.append(',');
                    pathRecord.append(edge.t1);
                    pathRecord.append(')');

                    // pathRecord.append(edge.start);
                    pathRecord.append("->");
                    if (chargeEdgeSet.contains(edge)) {
                        pathRecord.append("charge");
                    }

                }

                pathRecord.append('(');
                pathRecord.append(edge.v);
                pathRecord.append(',');
                pathRecord.append(edge.t2);
                pathRecord.append(')');
                // pathRecord.append(edge.end);
                // out.println(pathRecord.toString());
                System.out.println("originNode= " + path.originNode + " capacityType= " + path.capacityType);
                System.out.println(path.edgeIndexSet.toString());
                System.out.println(pathRecord.toString() + ":" + value);

                // -----------------------------------------------------------------------------//
            }

        }

    }

    // generate all possible paths for all nodes as start point
    public void GeneratePathFile() throws FileNotFoundException {
        PrintWriter out = new PrintWriter(fileName);

        int count = 0;

        for (int originNode = 0; originNode < dataModel.numNode; originNode++) {

            Set<Set<Integer>> checkRepeat = new HashSet<>();// if charge, we set
                                                            // the edge index as
                                                            // negative one

            out.println(originNode);

            for (int startTime = 0; startTime < dataModel.timePeriod; startTime++) {
                int startNodeIndex = originNode * dataModel.timePeriod + startTime;

                Stack<Node> stack = new Stack<>();
                List<Integer> pathEdgeRecord = new ArrayList<>();
                int pathDuration = 0;

                Node node = new Node();
                node.nodeIndex = startNodeIndex;
                node.ifUsed = false;
                // node.pathDuration=0;
                node.leadEdgeIndex = -100000000;
                node.powerLeft = dataModel.powerCapacity;
                stack.add(node);

                while (stack.size() > 0) {
                    Node currentNode = stack.peek();
                    if (!currentNode.ifUsed) { // need to generate new nodes

                        if (currentNode.leadEdgeIndex != -100000000) {
                            pathEdgeRecord.add(currentNode.leadEdgeIndex);
                            Edge edge = dataModel.edgeSet.get(Math.abs(currentNode.leadEdgeIndex));
                            if (edge.edgeType == 0) {
                                pathDuration += edge.duration;
                            } else {
                                pathDuration += 1;
                            }
                        }

                        if (currentNode.leadEdgeIndex == -100000000 || currentNode.nodeIndex != startNodeIndex) {
                            int holdingEdgeIndex = -1;
                            for (int edgeIndex : dataModel.pointToEdgeSet.get(currentNode.nodeIndex)) {
                                Edge edge = dataModel.edgeSet.get(edgeIndex);
                                if (edge.edgeType == 1)
                                    holdingEdgeIndex = edgeIndex;

                                int duration = 0;
                                if (edge.edgeType == 0) {
                                    duration = (int) edge.duration;
                                } else
                                    duration = 1;
                                if ((pathDuration + duration < dataModel.timePeriod)
                                        || (pathDuration + duration == dataModel.timePeriod
                                                && edge.end == startNodeIndex)) {
                                    if (currentNode.powerLeft >= edge.duration) {
                                        Node newNode = new Node();
                                        newNode.nodeIndex = edge.end;
                                        newNode.ifUsed = false;
                                        newNode.leadEdgeIndex = edgeIndex;
                                        newNode.powerLeft = currentNode.powerLeft - edge.duration;
                                        // newNode.pathDuration=(int)
                                        // (currentNode.pathDuration+duration);
                                        stack.add(newNode);
                                    }
                                }
                            }

                            // charge
                            if (currentNode.powerLeft < dataModel.powerCapacity) {
                                if (pathDuration + 1 < dataModel.timePeriod) {
                                    Edge edge = dataModel.edgeSet.get(holdingEdgeIndex);
                                    Node newNode = new Node();
                                    newNode.nodeIndex = edge.end;
                                    newNode.ifUsed = false;
                                    newNode.leadEdgeIndex = -holdingEdgeIndex;
                                    newNode.powerLeft = Math.min(dataModel.powerCapacity,
                                            currentNode.powerLeft + dataModel.chargeOnceDistance);
                                    stack.add(newNode);
                                }
                            }

                        } else { // currentNode.nodeIndex==startNodeIndex

                            Set pathEdgeSet = new HashSet<>();
                            for (int i = 0; i < pathEdgeRecord.size(); i++) {
                                pathEdgeSet.add(pathEdgeRecord.get(i));
                            }
                            if (!checkRepeat.contains(pathEdgeSet)) {
                                checkRepeat.add(pathEdgeSet);

                                count++;

                                for (int i = 0; i < pathEdgeRecord.size(); i++) {
                                    out.print(pathEdgeRecord.get(i));
                                    out.print(" ");
                                }
                                out.println();
                            }

                            // //------------------------------------check---------------------------------------//
                            // Queue<Edge> path=new PriorityQueue<>();
                            // for(int edgeIndex:pathEdgeRecord) {
                            // path.add(dataModel.edgeSet.get(edgeIndex));
                            // }
                            //
                            // out.println();
                            //
                            // StringBuilder pathRecord=new StringBuilder();
                            //
                            // Edge edge=null;
                            // int size=path.size();
                            // for(int i=0;i<size;i++) {
                            //
                            // edge=path.poll();
                            //// pathRecord.append('(');
                            //// pathRecord.append(edge.u);
                            //// pathRecord.append(',');
                            //// pathRecord.append(edge.t1);
                            //// pathRecord.append(')');
                            //
                            // pathRecord.append(edge.start);
                            //
                            // pathRecord.append("->");
                            //
                            // }
                            //
                            //// pathRecord.append('(');
                            //// pathRecord.append(edge.v);
                            //// pathRecord.append(',');
                            //// pathRecord.append(edge.t2);
                            //// pathRecord.append(')');
                            // pathRecord.append(edge.end);
                            // out.println(pathRecord.toString());
                            //
                            // //-----------------------------------------------------------------------------//

                        }

                        currentNode.ifUsed = true;

                    } else { // currentNode.ifUsed==true
                        if (pathEdgeRecord.size() > 0) {
                            int edgeIndex = pathEdgeRecord.remove(pathEdgeRecord.size() - 1);
                            Edge edge = dataModel.edgeSet.get(Math.abs(edgeIndex));
                            if (edge.edgeType == 0) {
                                pathDuration -= edge.duration;
                            } else {
                                pathDuration -= 1;
                            }
                        }

                        stack.remove(currentNode);
                    }

                }

            }

            out.println("end");

        }

        out.close();

        System.out.println("number of cycles= " + count);
    }

    public static void main(String[] args) throws IOException, IloException {

        // for(String arg:args) {
        // long time0=System.currentTimeMillis();
        // SNDRC sndrc=new SNDRC(arg);
        //
        // CplexSolver cplexSolver=new
        // CplexSolver(sndrc,"./output/path/outpath.txt");
        // cplexSolver.GeneratePathFile();
        // cplexSolver.Solve();
        //
        // long time1=System.currentTimeMillis();
        // System.out.println("Total time= "+(time1-time0));
        // }

        String path = "./testdata/";
//        String path = "./data/testset/";
        File file = new File(path);
        File[] fs = file.listFiles();

        for (File f : fs) {
            if (!f.isDirectory() && !f.isHidden()) {

                System.out.println("Solve for "+f.getName());
                long time0 = System.currentTimeMillis();
                SNDRC sndrc = new SNDRC(f.toString());
                CplexSolver cplexSolver = new CplexSolver(sndrc, "./output/path/outpath.txt");
                cplexSolver.GeneratePathFile();
                cplexSolver.Solve();

                long time1 = System.currentTimeMillis();
                System.out.println("Total time= " + (time1 - time0));

            }
        }

    }

}
