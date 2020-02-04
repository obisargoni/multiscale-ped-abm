package repastInterSim.environment;

import java.util.ArrayList;

import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.Geography;

import repastInterSim.agent.MobileAgent;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class Destination implements FixedGeography{
	
	private Context<Object> context;
	private Geography<Object> geography;
	private Geography<Destination> destinationGeography;
	private Geometry geom;
	private double arrivalDist = GlobalVars.MOBILE_AGENT_PARAMS.destinationArrivalDistance; // Distance from destination at which agents are removed from the simulation
		
	public Destination() {

	}
	
	@ScheduledMethod(start = 1, interval = 1, priority = ScheduleParameters.LAST_PRIORITY)
	public void removeAgent() {
        ArrayList<MobileAgent> AgentsToRemove = new ArrayList<MobileAgent>();
        
        // If there are no agents to remove end the simulation
        if (context.getObjects(MobileAgent.class).size() == 0) {
        	RunEnvironment.getInstance().endRun();
        }
        
        // Iterate over peds and remove them if they have arrived at the destination
        for (Object o :context.getObjects(MobileAgent.class)) {
        	MobileAgent mA  = (MobileAgent) o;
        	
        	// Get the geometries in the CRS used for spatial calculations
        	Geometry dGeom =  mA.getDestination().getGeom();
        	Geometry mAGeom = SpaceBuilder.getAgentGeometry(this.geography, mA);
        	
        	// If the pedestrian agent in within the bounds of the destination then remove it from the context as it has reached its destination
        	if (dGeom.isWithinDistance(mAGeom, this.arrivalDist)) {
        		AgentsToRemove.add(mA);
        		break; // End the iteration, only one pedestrian can be removed at a time
        	}
        }
        // Now iterate over all of the peds to remove and remove them from the context
        // Need to do this separately from iterating over the peds in the context since altering the context whilst iterating over it throws and exception
        for (MobileAgent mA : AgentsToRemove) {
        	mA.tidyForRemoval();
        	context.remove(mA);
        }
    }

	@Override
	public Geometry getGeom() {
		return this.geom;
	}

	@Override
	public void setGeom(Geometry g) {
		this.geom = g;
	}
	
	public void setRootGeography(Geography<Object> G) {
		this.geography = G;
	}
	
	public void setDestinationGeography(Geography<Destination> G) {
		this.destinationGeography = G;
	}
	
	public void setRootContext(Context<Object> C) {
		this.context = C;
	}

}
