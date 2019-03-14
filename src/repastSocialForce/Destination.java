package repastSocialForce;

import java.awt.Color;
import java.util.ArrayList;

import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.gis.Geography;
import repast.simphony.util.ContextUtils;

public class Destination {

	private Geography<Object> geography;
	private Color colour; // Colour of the destination
		
	public Destination(Geography<Object> geography, Color col) {
		this.geography = geography;
		this.colour = col;
	}
	
	@ScheduledMethod(start = 1, interval = 1, priority = ScheduleParameters.LAST_PRIORITY)
	 public void removePeds() {
	        Context<Object> context = ContextUtils.getContext(this);
	        Geography<Object> geography = this.geography;
	        
	        ArrayList<Ped> PedsToRemove = new ArrayList<Ped>();
	        
	        // If there are no pedestrians to remove end the simulation
	        if (context.getObjects(Ped.class).size() == 0) {
	        	RunEnvironment.getInstance().endRun();
	        }
	        
	        // Iterate over peds and remove them if they have arrive at the destination
	        for (Object p :context.getObjects(Ped.class)) {
	        	Ped P = (Ped) p;

	        	Geometry dGeom = geography.getGeometry(P.destination);
	        	Geometry coordP = geography.getGeometry(P);
	        	
	        	// If the pedestrian agent in within the bounds of the destination then remove it from the context as it has reached its destination
	        	if (dGeom.contains(coordP)) {
	        		//Source.removedPeds.add(P);
	        		PedsToRemove.add(P);
	        		break; // End the iteration since iterating over context having modified it can throw an exception
	        	}
	        }
	        // Now iterate over all of the peds to remove and remove them from the context
	        // Need to do this separately from iterating over the peds in the context since altering the context whilst iterating over it throws and exception
	        for (Ped P : PedsToRemove) {
	        	context.remove(P);
	        }
	    }

}
