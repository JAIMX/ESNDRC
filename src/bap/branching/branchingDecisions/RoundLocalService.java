package bap.branching.branchingDecisions;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import model.SNDRC;

/**
 * Round on local movement between cities (in order to remove symmetry), for all cycles
 * @author DELL
 *
 */
public class RoundLocalService implements BranchingDecision<SNDRC, Cycle> {

    public int roundUpOrDown;//0:round down;1:round up
    public int localServiceIndex; //correspond to serviceSet in SNDRC
    public double branchValue;
    public SNDRCPricingProblem associatedPricingProblem;
    
    public RoundLocalService(int roundUpOrDown,int localServiceIndex,double branchValue,SNDRCPricingProblem associatedPricingProblem){
        this.roundUpOrDown=roundUpOrDown;
        this.localServiceIndex=localServiceIndex;
        this.branchValue=branchValue;
        this.associatedPricingProblem=associatedPricingProblem;
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
            return "Round service("+localServiceIndex+")"+" capacityType="+associatedPricingProblem.capacityTypeS+" originNode="+associatedPricingProblem.originNodeO+" down:"+branchValue;
        }else {
        	return "Round service("+localServiceIndex+")"+" capacityType="+associatedPricingProblem.capacityTypeS+" originNode="+associatedPricingProblem.originNodeO+" up:"+branchValue;
        }
        
    }
    
    
    
}
