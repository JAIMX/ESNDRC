package debugger;

import org.jorlib.frameworks.columnGeneration.branchAndPrice.AbstractBranchAndPrice;
import org.jorlib.frameworks.columnGeneration.io.SimpleDebugger;
import org.jorlib.frameworks.columnGeneration.branchAndPrice.EventHandling.FinishProcessingNodeEvent;

import cg.Cycle;

public class FinishProcessingNodeDebugger extends SimpleDebugger {

	public FinishProcessingNodeDebugger(AbstractBranchAndPrice bap, boolean captureColumnGenerationEventsBAP) {
		super(bap,true);
	}
	
	@Override
	public void finishedColumnGenerationForNode(FinishProcessingNodeEvent event) {
		super.finishedColumnGenerationForNode(event);
		for(Cycle c:(List<Cycle>)event.node.getSolution()) {
			System.out.println(c);
		}
	}
	
}
