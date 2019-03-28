/*
 * Taken from Nick Malleson's RepastCity project
 */

package contexts;

import repast.simphony.context.DefaultContext;
import repastInterSim.Road;
import repastInterSim.UserPanel;

public class RoadContext extends DefaultContext<Road> {
	
	public RoadContext() {
		super(UserPanel.ROAD_CONTEXT);
	}

}
