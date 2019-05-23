package repastInterSim.environment;

import java.util.ArrayList;

import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.Geography;

import repastInterSim.agent.mobileAgent;
import repastInterSim.main.SpaceBuilder;

public class Destination implements FixedGeography{
	
	private Context<Object> context;
	private Geography<Object> geography;
	private Geography<Destination> destinationGeography;
	private Geometry geom;
	private double arrivalDist = 1; // Distance in metres at which agents are considered to have arrived at the destination and are removed from the space
		
	public Destination() {

	}
	
	@ScheduledMethod(start = 1, interval = 1, priority = ScheduleParameters.LAST_PRIORITY)
	public void removeAgent() {
        ArrayList<mobileAgent> AgentsToRemove = new ArrayList<mobileAgent>();
        
        // If there are no agents to remove end the simulation
        if (context.getObjects(mobileAgent.class).size() == 0) {
        	RunEnvironment.getInstance().endRun();
        }
        
        // Iterate over peds and remove them if they have arrived at the destination
        for (Object o :context.getObjects(mobileAgent.class)) {
        	mobileAgent mA  = (mobileAgent) o;
        	
        	// Get the geometries in the CRS used for spatial calculations
        	Geometry dGeom = SpaceBuilder.getAgentGeometry(this.destinationGeography, mA.getDestination());
        	Geometry coordP = SpaceBuilder.getAgentGeometry(this.geography, mA);
        	
        	// If the pedestrian agent in within the bounds of the destination then remove it from the context as it has reached its destination
        	if (dGeom.isWithinDistance(coordP, this.arrivalDist)) {
        		AgentsToRemove.add(mA);
        		break; // End the iteration, only one pedestrian can be removed at a time
        	}
        }
        // Now iterate over all of the peds to remove and remove them from the context
        // Need to do this separately from iterating over the peds in the context since altering the context whilst iterating over it throws and exception
        for (mobileAgent mA : AgentsToRemove) {
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
	
	public void setObjectGeography(Geography<Object> G) {
		this.geography = G;
	}
	
	public void setDestinationGeography(Geography<Destination> G) {
		this.destinationGeography = G;
	}
	
	public void setObjectContext(Context<Object> C) {
		this.context = C;
	}

}
