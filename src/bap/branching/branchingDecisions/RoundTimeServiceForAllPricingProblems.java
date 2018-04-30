package bap.branching.branchingDecisions;

import java.util.Set;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

import cg.Cycle;
import cg.SNDRCPricingProblem;
import model.SNDRC;

/**
 * Round on a local service during a period for all pricing problems
 * @author DELL
 *
 */
public class RoundTimeServiceForAllPricingProblems implements BranchingDecision<SNDRC, Cycle>{
    public int roundUpOrDown;//0:round down;1:round up
    public int localServiceIndex; //correspond to serviceSet in SNDRC
    public double branchValue;
    public Set timeSet;
    public RoundTimeServiceForAllPricingProblems(int roundUpOrDown,int localServiceIndex,double branchValue,Set branchTimeSet){
        this.roundUpOrDown=roundUpOrDown;
        this.localServiceIndex=localServiceIndex;
        this.branchValue=branchValue;
        this.timeSet=branchTimeSet;
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
            return "Round service("+localServiceIndex+")"+"timeSet "+" down:"+branchValue;
        }else {
            return "Round service("+localServiceIndex+")"+"timeSet "+" up:"+branchValue;
        }
        
    }

}
