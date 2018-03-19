package bap.branching.branchingDecisions;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.branchingDecisions.BranchingDecision;
import org.jorlib.frameworks.columnGeneration.master.cutGeneration.AbstractInequality;

import cg.Cycle;
import model.SNDRC;

/**
 * Round on local movement between cities (in order to remove symmetry), for all
 * cycles
 * 
 * @author DELL
 *
 */
public class RoundLocalServiceForAllPricingProblems implements BranchingDecision<SNDRC, Cycle> {

    public int roundUpOrDown;// 0:round down;1:round up
    public int localServiceIndex; // correspond to serviceSet in SNDRC
    public double branchValue;

    public RoundLocalServiceForAllPricingProblems(int roundUpOrDown, int localServiceIndex, double branchValue) {
        this.roundUpOrDown = roundUpOrDown;
        this.localServiceIndex = localServiceIndex;
        this.branchValue = branchValue;
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
        if (roundUpOrDown == 0) {
            return "Round service(" + localServiceIndex + ") down:" + branchValue;
        } else {
            return "Round service(" + localServiceIndex + ") up:" + branchValue;
        }

    }
    
    
}
