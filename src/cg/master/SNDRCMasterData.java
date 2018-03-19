package cg.master;

import java.util.*;

import org.jorlib.frameworks.columnGeneration.master.MasterData;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;
import org.jorlib.frameworks.columnGeneration.util.OrderedBiMap;

import com.google.common.collect.RangeSet;

import bap.branching.branchingDecisions.RoundHoldingEdge;
import bap.branching.branchingDecisions.RoundLocalService;
import bap.branching.branchingDecisions.RoundLocalServiceForAllPricingProblems;
import bap.branching.branchingDecisions.RoundQ;
import bap.branching.branchingDecisions.RoundServiceEdge;
import bap.branching.branchingDecisions.RoundServiceEdgeForAllPricingProblems;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import cg.master.cuts.StrongInequality;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import model.SNDRC;

/**
 * Container which stores information coming from the master problem. It
 * contains: a reference to the cplex model reference to the pricing problem
 * 
 * @author sxx
 *
 */
public class SNDRCMasterData extends MasterData<SNDRC, Cycle, SNDRCPricingProblem, IloNumVar> {
    public final IloCplex cplex;
    public List<Map<Integer, IloNumVar>> x; // map:edgeIndex, x variable
    public List<Map<Integer, Double>> xValues;// record the value of x variables
    public IloNumVar[][] q;
    public final List<SNDRCPricingProblem> pricingProblems;

    // branch on q variables
    public final Set<RoundQ> qBranchingSet;
    public final Map<RoundQ, IloRange> qBranchingconstraints;
    public final Map<SNDRCPricingProblem, Integer> qVariableLimit;

    // branch on service edge
    public final Set<RoundServiceEdge> serviceEdgeBrachingSet;
    public final Map<RoundServiceEdge, IloRange> ServiceEdgeBranchingConstraints;

    // branch on service edge for all pricing problems
    public Set<RoundServiceEdgeForAllPricingProblems> serviceEdge4AllBrachingSet;
    public Map<RoundServiceEdgeForAllPricingProblems, IloRange> serviceEdge4AllBranchingConstraints;

    // branch on local service for all cycles
    public Set<RoundLocalService> localServiceBranchingSet;
    public Map<RoundLocalService, IloRange> localServiceBranchingConstraints;

    // branch on holding edges
    public Set<RoundHoldingEdge> holdingEdgeBranchingSet;
    public Map<RoundHoldingEdge, IloRange> holdingEdgeBranchingConstraints;

    // branch on local service 4 all
    public Set<RoundLocalServiceForAllPricingProblems> localService4AllBranchingSet;
    public Map<RoundLocalServiceForAllPricingProblems, IloRange> localService4AllBranchingConstraints;

    // for acceleration:
    public Map<Cycle, IloRange> fixVarConstraints;

    // for strong cut
    public Map<Integer, Double> edgeValueMap;
    public final Map<StrongInequality, IloRange> strongInequalities;

    public SNDRCMasterData(IloCplex cplex, List<SNDRCPricingProblem> pricingProblems,
            Map<SNDRCPricingProblem, OrderedBiMap<Cycle, IloNumVar>> varMap,
            Map<RoundQ, IloRange> qBranchingconstraints,
            Map<RoundServiceEdge, IloRange> ServiceEdgeBranchingConstraints, List<Map<Integer, IloNumVar>> x,
            IloNumVar[][] q, Map<RoundServiceEdgeForAllPricingProblems, IloRange> serviceEdge4AllBranchingConstraints,
            Map<RoundLocalService, IloRange> localServiceBranchingConstraints,
            Map<RoundHoldingEdge, IloRange> holdingEdgeBranchingConstraints, Map<RoundLocalServiceForAllPricingProblems, IloRange> localService4AllBranchingConstraints,SNDRC dataModel) {
        super(varMap);
        this.cplex = cplex;
        this.x = x;
        this.q = q;
        this.pricingProblems = pricingProblems;

        this.qBranchingSet = new HashSet<>();
        this.qBranchingconstraints = new HashMap<>();

        if (qBranchingconstraints != null) {
            for (RoundQ roundQ : qBranchingconstraints.keySet()) {
                this.qBranchingSet.add(roundQ);
                this.qBranchingconstraints.put(roundQ, qBranchingconstraints.get(roundQ));
            }
        }

        // calculate qVariableLimit
        qVariableLimit = new HashMap<>();
        HashMap<SNDRCPricingProblem, Integer> temp = new HashMap<>();
        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            temp.put(pricingProblem, 0);
            qVariableLimit.put(pricingProblem,
                    dataModel.vehicleLimit[pricingProblem.capacityTypeS][pricingProblem.originNodeO]);
        }

        for (RoundQ roundQ : qBranchingSet) {
            if (roundQ.roundUpOrDown == 1 && temp.get(roundQ.associatedPricingProblem) < roundQ.qValue) {
                // temp.put(roundQ.associatedPricingProblem,
                // MathProgrammingUtil.doubleToInt(roundQ.qValue));
                temp.put(roundQ.associatedPricingProblem, MathProgrammingUtil.doubleToInt(Math.ceil(roundQ.qValue)));
            }
        }

        for (SNDRCPricingProblem pricingProblem : pricingProblems) {
            qVariableLimit.put(pricingProblem, qVariableLimit.get(pricingProblem) - temp.get(pricingProblem));
        }

        this.serviceEdgeBrachingSet = new HashSet<>();
        this.ServiceEdgeBranchingConstraints = new HashMap<>();
        if (ServiceEdgeBranchingConstraints != null) {
            for (RoundServiceEdge roundServiceEdge : ServiceEdgeBranchingConstraints.keySet()) {
                this.serviceEdgeBrachingSet.add(roundServiceEdge);
                this.ServiceEdgeBranchingConstraints.put(roundServiceEdge,
                        ServiceEdgeBranchingConstraints.get(roundServiceEdge));
            }
        }

        this.serviceEdge4AllBrachingSet = new HashSet<>();
        this.serviceEdge4AllBranchingConstraints = new HashMap<>();
        if (serviceEdge4AllBranchingConstraints != null) {
            for (RoundServiceEdgeForAllPricingProblems roundServiceEdge4All : serviceEdge4AllBranchingConstraints
                    .keySet()) {
                this.serviceEdge4AllBrachingSet.add(roundServiceEdge4All);
                this.serviceEdge4AllBranchingConstraints.put(roundServiceEdge4All,
                        serviceEdge4AllBranchingConstraints.get(roundServiceEdge4All));
            }
        }

        this.localServiceBranchingSet = new HashSet<>();
        this.localServiceBranchingConstraints = new HashMap<>();
        if (localServiceBranchingConstraints != null) {
            for (RoundLocalService roundLocalService : localServiceBranchingConstraints.keySet()) {
                this.localServiceBranchingSet.add(roundLocalService);
                this.localServiceBranchingConstraints.put(roundLocalService,
                        localServiceBranchingConstraints.get(roundLocalService));
            }
        }

        this.holdingEdgeBranchingSet = new HashSet<>();
        this.holdingEdgeBranchingConstraints = new HashMap<>();
        if (holdingEdgeBranchingConstraints != null) {
            for (RoundHoldingEdge branch : holdingEdgeBranchingConstraints.keySet()) {
                this.holdingEdgeBranchingSet.add(branch);
                this.holdingEdgeBranchingConstraints.put(branch, holdingEdgeBranchingConstraints.get(branch));
            }
        }

        this.localService4AllBranchingSet = new HashSet<>();
        this.localService4AllBranchingConstraints = new HashMap<>();
        if (localService4AllBranchingConstraints != null) {
            for (RoundLocalServiceForAllPricingProblems roundLocalService : localService4AllBranchingConstraints.keySet()) {
                this.localService4AllBranchingSet.add(roundLocalService);
                this.localService4AllBranchingConstraints.put(roundLocalService,
                        localService4AllBranchingConstraints.get(roundLocalService));
            }
        }

        fixVarConstraints = new HashMap<>();

        edgeValueMap = new HashMap<Integer, Double>();
        strongInequalities = new HashMap<>();

    }

}
