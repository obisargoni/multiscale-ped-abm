package repastInterSim.environment.contexts;

import repast.simphony.context.DefaultContext;
import repastInterSim.environment.PedObstruction;
import repastInterSim.main.GlobalVars;

public class PedObstructionPointsContext extends DefaultContext<PedObstruction> {
	
	public PedObstructionPointsContext() {
		super(GlobalVars.CONTEXT_NAMES.PED_OBSTRUCTION_POINTS_CONTEXT);
	}

}
