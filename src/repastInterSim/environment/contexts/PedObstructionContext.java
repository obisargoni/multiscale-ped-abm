package repastInterSim.environment.contexts;

import repast.simphony.context.DefaultContext;
import repastInterSim.environment.PedObstruction;
import repastInterSim.main.GlobalVars;

public class PedObstructionContext extends DefaultContext<PedObstruction> {
	
	public PedObstructionContext() {
		super(GlobalVars.CONTEXT_NAMES.PED_OBSTRUCTION_CONTEXT);
	}

}
