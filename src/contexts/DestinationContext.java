package contexts;

import repast.simphony.context.DefaultContext;
import repastInterSim.Destination;
import repastInterSim.UserPanel;

public class DestinationContext extends DefaultContext<Destination> {
	
	public DestinationContext() {
		super(UserPanel.DESTINATION_CONTEXT);
	}

}
