package repastInterSim.environment.contexts;

import repast.simphony.context.DefaultContext;
import repastInterSim.environment.Road;
import repastInterSim.main.GlobalVars;

public class RoadContext extends DefaultContext<Road> {
	
	public RoadContext() {
		super(GlobalVars.CONTEXT_NAMES.ROAD_CONTEXT);
	}


}
