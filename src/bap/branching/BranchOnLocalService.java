package bap.branching;

import java.util.Arrays;
import java.util.List;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import bap.branching.branchingDecisions.RoundLocalService;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import model.SNDRC;
import model.SNDRC.Edge;

public class BranchOnLocalService extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem> {

    private int branchServiceIndex;
    private double thresholdValue; 
    double bestServiceValue;
    
    public BranchOnLocalService(SNDRC modelData,List<SNDRCPricingProblem> pricingProblems,double thresholdValue){
        super(modelData, pricingProblems);
        this.thresholdValue=thresholdValue;
    }
    
    @Override
    protected boolean canPerformBranching(List<Cycle> solution){
        //reset values
        bestServiceValue=0;
        
        Double[] serviceTotalValue=new Double[dataModel.numService];
        for(int i=0;i<serviceTotalValue.length;i++){
            serviceTotalValue[i]=(double) 0;
        }
        for(Cycle cycle:solution){
            for(int edgeIndex:cycle.edgeIndexSet){
                Edge edge=dataModel.edgeSet.get(edgeIndex);
                if(edge.edgeType==0){
                    serviceTotalValue[edge.serviceIndex]=serviceTotalValue[edge.serviceIndex]+cycle.value;
                }
            }
        }
        
        
        //check if all values are integers
        boolean isAllInteger=true;
        double bestDifference=1;
        
        //Select the local service closest to threshold value
        for(int serviceIndex=0;serviceIndex<dataModel.numService;serviceIndex++){
            double value=serviceTotalValue[serviceIndex];
            
            if(MathProgrammingUtil.isFractional(value)){
                isAllInteger=false;
                double decimalPart=value-Math.floor(value);
                if(Math.abs(thresholdValue-decimalPart)<bestDifference){
                    branchServiceIndex=serviceIndex;
                    bestServiceValue=value;
                    bestDifference=Math.abs(thresholdValue-decimalPart);
                }
            }
        }
        
        return(!isAllInteger);
    }
    
    protected List<BAPNode<SNDRC,Cycle>> getBranches(BAPNode<SNDRC,Cycle> parentNode){
        
        //Branch 1:round q down to the nearest integer
        RoundLocalService branchingDecision1=new RoundLocalService(0,branchServiceIndex,bestServiceValue);
        BAPNode<SNDRC,Cycle> node1=this.createBranch(parentNode, branchingDecision1, parentNode.getSolution(), parentNode.getInequalities());
//      BAPNode<SNDRC,Cycle> node1=this.createBranch(parentNode, branchingDecision1,initialSolution , parentNode.getInequalities());
        
        
//      initialSolution=new ArrayList<Cycle>();
        //Branch 2:round q up to the nearest integer
        RoundLocalService branchingDecision2=new RoundLocalService(1,branchServiceIndex,bestServiceValue);
        BAPNode<SNDRC,Cycle> node2=this.createBranch(parentNode, branchingDecision2, parentNode.getSolution(), parentNode.getInequalities());
//      BAPNode<SNDRC,Cycle> node2=this.createBranch(parentNode, branchingDecision2, initialSolution, parentNode.getInequalities());
        
        return Arrays.asList(node2,node1);
    }
    
    
}
