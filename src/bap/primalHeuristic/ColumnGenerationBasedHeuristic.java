package bap.primalHeuristic;

import java.io.File;
import java.io.IOException;
import java.util.*;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.colgenMain.ColGen;
import org.jorlib.frameworks.columnGeneration.io.SimpleCGLogger;
import org.jorlib.frameworks.columnGeneration.io.SimpleDebugger;
import org.jorlib.frameworks.columnGeneration.io.TimeLimitExceededException;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.CutHandler;
import org.jorlib.frameworks.columnGeneration.pricing.AbstractPricingProblemSolver;
import org.jorlib.frameworks.columnGeneration.util.Configuration;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import bap.branching.BranchOnLocalService;
import bap.branching.BranchOnLocalServiceForAllPricingProblems;
import bap.branching.BranchOnServiceEdge;
import bap.branching.BranchOnTimeService;
import cg.Cycle;
import cg.ExactPricingProblemSolver;
import cg.SNDRCPricingProblem;
import cg.master.Master;
import cg.master.SNDRCMasterData;
import cg.master.cuts.StrongInequalityGenerator;
import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import model.SNDRC;
import model.SNDRC.Edge;

public class ColumnGenerationBasedHeuristic {

    SNDRC dataModel;
    double thresholdValue;
    boolean ifAccelerationForUB;
    int objectiveIncumbentSolution;
    List<Cycle> incumbentSolution = new ArrayList<>();
    Map<Cycle, Double> optSolutionValueMap;

    Map<Cycle, IloIntVar> cycleVar;
    long runTime = 0;

    public ColumnGenerationBasedHeuristic(SNDRC dataModel, Double thresholdValue, boolean ifAccelerationForUB) {
        this.dataModel = dataModel;
        this.thresholdValue = thresholdValue;
        this.ifAccelerationForUB = ifAccelerationForUB;
    }

    public void Solve() throws TimeLimitExceededException, IloException {

        runTime = System.currentTimeMillis();

        /// -----------------------------------------------------step1: solve
        /// the root node
        /// ------------------------------------------------------------------///
        // Create the pricing problems
        List<SNDRCPricingProblem> pricingProblems = new LinkedList<SNDRCPricingProblem>();
        for (int capacityType = 0; capacityType < dataModel.numOfCapacity; capacityType++) {
            for (int originNode = 0; originNode < dataModel.numNode; originNode++) {
                String name = "capacity type: " + capacityType + " origin node: " + originNode;
                SNDRCPricingProblem pricingProblem = new SNDRCPricingProblem(dataModel, name, capacityType, originNode);
                pricingProblems.add(pricingProblem);
            }
        }

        // Create a cutHandler
        CutHandler<SNDRC, SNDRCMasterData> cutHandler = new CutHandler<>();
        StrongInequalityGenerator cutGen = new StrongInequalityGenerator(dataModel, pricingProblems, 0);
        // cutHandler.addCutGenerator(cutGen);

        // Create the Master Problem
        Master master = new Master(dataModel, pricingProblems, cutHandler, cutGen, false);

        // Define which solvers to use
        List<Class<? extends AbstractPricingProblemSolver<SNDRC, Cycle, SNDRCPricingProblem>>> solvers = Collections
                .singletonList(ExactPricingProblemSolver.class);

        // ----------------------------------generateInitialFeasibleSolution-------------------------------------------------//
        List<Cycle> artificalVars = new ArrayList<Cycle>();
        // for weak forcing constraints(ifForResourceBoundConstraints=0)
        for (int edgeIndex = 0; edgeIndex < dataModel.numServiceArc; edgeIndex++) {
            Set<Integer> set = new HashSet<>();
            set.add(edgeIndex);
            Cycle cycle = new Cycle(pricingProblems.get(0), true, "Artificial", set, 100000000, 0, 0);
            artificalVars.add(cycle);
        }

        // for resource bound constraints(ifForResourceBoundConstraints=1)
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            Set<Integer> set = new HashSet<>();
            Cycle cycle = new Cycle(pricingProblem, true, "Artificial", set, 100000000, 0, 1);
            artificalVars.add(cycle);
        }

        // for holding edge branch constraints(ifForResourceBoundConstraints=2)
        Set<Integer> set = new HashSet<>();
        Cycle cycle0 = new Cycle(pricingProblems.get(0), true, "Artificial", set, 100000000, 0, 2);
        artificalVars.add(cycle0);

        // ----------------------------------generateInitialFeasibleSolution-------------------------------------------------//

        ColGen<SNDRC, Cycle, SNDRCPricingProblem> cg = new ColGen<SNDRC, Cycle, SNDRCPricingProblem>(dataModel, master,
                pricingProblems, solvers, artificalVars, Integer.MAX_VALUE, Double.MIN_VALUE);

        // SimpleDebugger debugger = new SimpleDebugger(cg);

        // SimpleCGLogger logger = new SimpleCGLogger(cg, new
        // File("./output/cgLogger.log"));

        cg.solve(System.currentTimeMillis() + 18000000L); // 5 hour limit

        System.out.println("Time of first LP solve= " + (System.currentTimeMillis() - runTime));

        if (ifAccelerationForUB) {
            /// -------------------------------AccelerationForUB------------------------------------///
            List<Cycle> solution = cg.getSolution();

            while (!isIntegerNode(solution)) {

                boolean ifAllBelowThresholdValue = true;
                solution = cg.getSolution();

                for (Cycle cycle : solution) {
                    if (MathProgrammingUtil.isFractional(cycle.value)) {
                        double decimalValue = cycle.value - (int) cycle.value;
                        if (decimalValue > thresholdValue) {
                            ifAllBelowThresholdValue = false;
                            ((Master) master).addFixVarConstraint(cycle);
                        }
                    }
                }

                // if all cycles' value are below the threshold value, fix the
                // variable with highest fractional decimal value
                double record = 0;
                Cycle cycleRecord = null;
                if (ifAllBelowThresholdValue) {
                    for (Cycle cycle : solution) {
                        if (MathProgrammingUtil.isFractional(cycle.value)) {
                            double decimalValue = cycle.value - (int) cycle.value;
                            if (decimalValue > record) {
                                cycleRecord = cycle;
                                record = decimalValue;
                            }
                        }
                    }

                    if (cycleRecord != null) {
                        ((Master) master).addFixVarConstraint(cycleRecord);
                    } else {
                        break;
                    }
                }

                // here we should check if the master problem is feasible
                if (((Master) master).CheckFeasibility() == false) {
                    break;
                }

                List<Cycle> nullList = new ArrayList<>();
                cg = new ColGen<>(dataModel, master, pricingProblems, solvers, nullList, Integer.MAX_VALUE,
                        Double.MIN_VALUE);
                cg.solve(System.currentTimeMillis() + 18000000L); // 5 hour

                if (isInfeasibleNode(cg.getSolution())) {
                    break;
                }

                if (isIntegerNode(cg.getSolution())) {

                    int integerObjective = MathProgrammingUtil.doubleToInt(cg.getObjective());
                    System.out.println("We have found a feasible solution: " + integerObjective);

                    this.objectiveIncumbentSolution = integerObjective;
                    this.incumbentSolution = new ArrayList<>();
                    for (Cycle cycle : cg.getSolution()) {
                        this.incumbentSolution.add(cycle);
                    }
                    optSolutionValueMap = new HashMap<>();
                    for (Cycle cycle : incumbentSolution) {
                        optSolutionValueMap.put(cycle, cycle.value);
                    }

                    break;
                }

            }

            /// -------------------------------AccelerationForUB------------------------------------///
            System.out.println("Time of acceleration LP solve= " + (System.currentTimeMillis() - runTime));
        }

        // pick up all the cycles in master problem
        int amount = 0;
        Map<SNDRCPricingProblem, Set<Cycle>> cycleSet = new HashMap<>();
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            Set<Cycle> tempSet = master.getColumns(pricingProblem);
            cycleSet.put(pricingProblem, tempSet);
            amount += tempSet.size();
        }

        /// -----------------------------------------------------step2:
        /// construct a cplex
        /// model------------------------------------------------------------------///
        IloCplex cplex = new IloCplex();
        // cplex.setOut(null);
        cplex.setParam(IloCplex.IntParam.Threads, 4);
        cplex.setParam(IloCplex.Param.Simplex.Tolerances.Markowitz, 0.1);

        List<Map<Integer, IloNumVar>> x; // map:edgeIndex, x variable

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
                expr.clear();
                resourceBoundConstraints[s][o] = cplex.addRange(0, dataModel.vehicleLimit[s][o]);
            }
        }

        // add all columns in cycleSet
        this.cycleVar = new HashMap<>();
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            int count = 0;
            Set<Cycle> tempSet = cycleSet.get(pricingProblem);

            for (Cycle cycle : tempSet) {

                // Register column with objective
                IloColumn iloColumn = cplex.column(obj, cycle.cost);

                // weak forcing constraints
                for (int edgeIndex : cycle.edgeIndexSet) {
                    if (dataModel.edgeSet.get(edgeIndex).edgeType == 0) {
                        iloColumn = iloColumn.and(cplex.column(weakForcingConstraints[edgeIndex],
                                -dataModel.capacity[cycle.associatedPricingProblem.capacityTypeS]));
                    }

                }

                // resource bound constraints
                iloColumn = iloColumn.and(cplex.column(
                        resourceBoundConstraints[cycle.associatedPricingProblem.capacityTypeS][cycle.associatedPricingProblem.originNodeO],
                        1));

                // Create the variable and store it
                IloIntVar var = cplex.intVar(iloColumn, 0, Integer.MAX_VALUE,
                        "z_" + cycle.associatedPricingProblem.capacityTypeS + ","
                                + cycle.associatedPricingProblem.originNodeO + "," + count);
                cycleVar.put(cycle, var);
                count++;

            }

        }

        /// -----------------------------------------------------step3:cplex
        /// solve and
        /// output------------------------------------------------------------------///

        cg.close();
        cutHandler.close();

        System.out.println();
        System.out.println("There are " + amount + " columns added to the model.");
        System.out.println();

        runTime = System.currentTimeMillis() - runTime;
        long timeLeft = 36000000 - runTime;
        cplex.setParam(IloCplex.DoubleParam.TiLim, timeLeft / 1000);
        // cplex.setParam(IloCplex.DoubleParam.TiLim, 36000);
        cplex.solve();
        System.out.println("optimal objective= " + cplex.getObjValue());
        System.out.println();

        // output solution
        for (Cycle cycle : cycleVar.keySet()) {
            IloIntVar var = cycleVar.get(cycle);
            int value = (int) cplex.getValue(var);
            if (value > 0.1) {
                System.out.println(cycle);
                System.out.println(out(cycle) + ":" + value);
                System.out.println();
            }
        }
        
        
        

    }

    public String out(Cycle column) {

        Queue<Edge> path = new PriorityQueue<>();

        for (int edgeIndex : column.edgeIndexSet) {
            path.add(dataModel.edgeSet.get(edgeIndex));
        }

        StringBuilder pathRecord = new StringBuilder();

        Edge edge = null;
        int size = path.size();
        for (int i = 0; i < size; i++) {

            edge = path.poll();
            pathRecord.append('(');
            pathRecord.append(edge.u);
            pathRecord.append(',');
            pathRecord.append(edge.t1);
            pathRecord.append(')');

            pathRecord.append("->");

        }

        pathRecord.append('(');
        pathRecord.append(edge.v);
        pathRecord.append(',');
        pathRecord.append(edge.t2);
        pathRecord.append(')');

        return pathRecord.toString();

    }

    public boolean isIntegerNode(List<Cycle> solution) {
        boolean out = true;
        for (Cycle cycle : solution) {
            if (MathProgrammingUtil.isFractional(cycle.value)) {
                out = false;
                break;
            }
        }

        return out;
    }

    public boolean isInfeasibleNode(List<Cycle> solution) {
        boolean out = false;
        for (Cycle cycle : solution) {
            if (cycle.isArtificialColumn) {
                out = true;
                break;
            }
        }

        return out;

    }

    public static void main(String[] args) throws IOException, TimeLimitExceededException, IloException {
        SNDRC sndrc;
        for (String arg : args) {
            long time0 = System.currentTimeMillis();
            sndrc = new SNDRC(arg);
            // sndrc.Output();
            Properties properties = new Properties();
            // properties.setProperty("EXPORT_MODEL", "True");
            // properties.setProperty("MAXTHREADS", "10");
            // properties.setProperty("PRECISION", "0.001");
            Configuration.readFromFile(properties);

            ColumnGenerationBasedHeuristic solver = new ColumnGenerationBasedHeuristic(sndrc, 0.65, true);
            solver.Solve();
            long time1 = System.currentTimeMillis();
            System.out.println("Total time= " + (time1 - time0));
        }

    }
}
