package bap.branching;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchCreator;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.BAPNode;
import org.jorlib.frameworks.columnGeneration.util.MathProgrammingUtil;

import bap.branching.branchingDecisions.RoundLocalService;
import bap.branching.branchingDecisions.RoundTimeService;
import cg.Cycle;
import cg.SNDRCPricingProblem;
import model.SNDRC;
import model.SNDRC.Edge;

public class BranchOnTimeService extends AbstractBranchCreator<SNDRC, Cycle, SNDRCPricingProblem>{

    private int branchServiceIndex;
    private double thresholdValue; 
    private double bestServiceValue;
    private SNDRCPricingProblem branchPricingProblem;
    private int timePeriodLength;
    private Set branchTimeSet;
    
    public BranchOnTimeService(SNDRC modelData,List<SNDRCPricingProblem> pricingProblems,double thresholdValue,int timePeriodLength){
        super(modelData, pricingProblems);
        this.thresholdValue=thresholdValue;
        this.timePeriodLength=timePeriodLength;
    }
    
    @Override
    protected boolean canPerformBranching(List<Cycle> solution){
        
        Map<SNDRCPricingProblem,Double[]> serviceTotalValue=new HashMap<>();
        for(SNDRCPricingProblem associatedPricingProblem:pricingProblems){
            Double[] temp=new Double[dataModel.numServiceArc];
            
            for (int t = 0; t < dataModel.numServiceArc; t++) {
                temp[t] = (double) 0;
            }
            serviceTotalValue.put(associatedPricingProblem, temp);
        }
        
        
        
        for(Cycle cycle:solution){
            SNDRCPricingProblem associatedPricingProblem=cycle.associatedPricingProblem;
            Double[] temp=serviceTotalValue.get(associatedPricingProblem);
            
            for(int edgeIndex:cycle.edgeIndexSet){
                Edge edge=dataModel.edgeSet.get(edgeIndex);
                if(edge.edgeType==0){
                    temp[edgeIndex]=temp[edgeIndex]+cycle.value;
                }
            }
        }
        
        
        
        //check if all values are integers
        boolean isAllInteger=true;
        double bestDifference=1;
        
        
        
       //Select the service during a time period closest to threshold value
        for(SNDRCPricingProblem pricingProblem:pricingProblems){
            Double[] temp=serviceTotalValue.get(pricingProblem);
            
            for(int serviceIndex=0;serviceIndex<dataModel.numService;serviceIndex++){

                
                for(int t1=0;t1<dataModel.timePeriod;t1++){
                    
                    double value=0;
                    Set timeSet=new HashSet<>();
                    for(int t=0;t<timePeriodLength;t++){
                        int currentTime=t1+t;
                        currentTime=currentTime%dataModel.timePeriod;
                        timeSet.add(currentTime);
                    }
                    
                    // calculate value
                    for(int serviceEdgeIndex=0;serviceEdgeIndex<dataModel.numServiceArc;serviceEdgeIndex++){
                        Edge edge=dataModel.edgeSet.get(serviceEdgeIndex);
                        if(edge.serviceIndex==serviceIndex&&timeSet.contains(edge.t1)){
                            value+=temp[serviceEdgeIndex];
                        }
                    }
                    
                    
                    
                    //check if variable value is integer
                    if(MathProgrammingUtil.isFractional(value)){
                        isAllInteger=false;
                        double decimalPart=value-Math.floor(value);
                        
                        if(Math.abs(thresholdValue-decimalPart)<bestDifference){
                            this.branchPricingProblem=pricingProblem;
                            branchServiceIndex=serviceIndex;
                            branchTimeSet=timeSet;
                            bestServiceValue=value;
                            bestDifference=Math.abs(thresholdValue-decimalPart);
                        }
                    }
                    
                    
                }
            }
        }
        
        return(!isAllInteger);
        
    }
    
    
    protected List<BAPNode<SNDRC,Cycle>> getBranches(BAPNode<SNDRC,Cycle> parentNode){
        
        //Branch 1:round down to the nearest integer
        RoundTimeService branchingDecision1=new RoundTimeService(0,branchServiceIndex,bestServiceValue,branchPricingProblem,branchTimeSet);
        BAPNode<SNDRC,Cycle> node1=this.createBranch(parentNode, branchingDecision1, parentNode.getSolution(), parentNode.getInequalities());
        
        
        //Branch 2:round up to the nearest integer
        RoundTimeService branchingDecision2=new RoundTimeService(1,branchServiceIndex,bestServiceValue,branchPricingProblem,branchTimeSet);
        BAPNode<SNDRC,Cycle> node2=this.createBranch(parentNode, branchingDecision2, parentNode.getSolution(), parentNode.getInequalities());
        
        return Arrays.asList(node2,node1);
    }
    
    
    
}
