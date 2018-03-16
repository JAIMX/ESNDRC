package bap.branching.branchingDecisions;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

import cg.Cycle;
import model.SNDRC;

/**
 * for each time t, the number of holding arcs ran by all the cycles from time t should be integer
 * @author DELL
 *
 */
public class RoundHoldingEdge implements BranchingDecision<SNDRC, Cycle>{

    public int roundUpOrDown;//0:round down;1:round up
    public int branchTime;
    public double branchValue;
    
    
    public RoundHoldingEdge(int roundUpOrDown, int branchTime,double branchValue){
        this.roundUpOrDown=roundUpOrDown;
        this.branchTime=branchTime;
        this.branchValue=branchValue;
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
            return "Round time="+branchTime+" down:"+branchValue;
        }else {
            return "Round time="+branchTime+" up:"+branchValue;
        }
        
    }
    
}
