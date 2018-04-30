package bap.branching.branchingDecisions;

import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import model.SNDRC;

/**
 * Round on a local service during a period for a specific pricing problem
 * @author DELL
 *
 */
public class RoundTimeService implements BranchingDecision<SNDRC, Cycle>{
    public int roundUpOrDown;//0:round down;1:round up
    public int localServiceIndex; //correspond to serviceSet in SNDRC
    public double branchValue;
    public SNDRCPricingProblem associatedPricingProblem;
    public Set timeSet;
    
    public RoundTimeService(int roundUpOrDown,int localServiceIndex,double branchValue,SNDRCPricingProblem associatedPricingProblem,Set timeSet){
        this.roundUpOrDown=roundUpOrDown;
        this.localServiceIndex=localServiceIndex;
        this.branchValue=branchValue;
        this.associatedPricingProblem=associatedPricingProblem;
        this.timeSet=timeSet;
    }
    
    @Override
    public boolean inEqualityIsCompatibleWithBranchingDecision(AbstractInequality inequality) {
        return true;
    }
    
    @Override
    public boolean columnIsCompatibleWithBranchingDecision(Cycle cycle) {
        return true;
    }
    
    
    @Override
    public String toString() {
        if(roundUpOrDown==0) {
            return "Round service("+localServiceIndex+")"+" capacityType="+associatedPricingProblem.capacityTypeS+" originNode="+associatedPricingProblem.originNodeO+"timeSet "+" down:"+branchValue;
        }else {
            return "Round service("+localServiceIndex+")"+" capacityType="+associatedPricingProblem.capacityTypeS+" originNode="+associatedPricingProblem.originNodeO+"timeSet "+" up:"+branchValue;
        }
        
    }

}
