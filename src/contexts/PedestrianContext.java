package contexts;

import repast.simphony.context.DefaultContext;
import repastInterSim.Ped;
import repastInterSim.UserPanel;

public class PedestrianContext extends DefaultContext<Ped> {
	
	public PedestrianContext() {
		super(UserPanel.PEDESTRIAN_CONTEXT);
	}
}
