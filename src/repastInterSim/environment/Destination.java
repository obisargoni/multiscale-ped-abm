package repastInterSim.environment;

import java.util.ArrayList;

import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.Geography;

import repastInterSim.agent.Ped;
import repastInterSim.main.SpaceBuilder;

public class Destination implements FixedGeography{
	
	private Context<Object> context;
	private Geography<Object> geography;
	private Geometry geom;
	private double arrivalDist = 1; // Distance in metres at which agents are considered to have arrived at the destination and are removed from the space
		
	public Destination() {

	}
	
	@ScheduledMethod(start = 1, interval = 1, priority = ScheduleParameters.LAST_PRIORITY)
	 public void removePeds() throws MismatchedDimensionException, NoSuchAuthorityCodeException, FactoryException, TransformException {
		
	        ArrayList<Ped> PedsToRemove = new ArrayList<Ped>();
	        
	        // If there are no pedestrians to remove end the simulation
	        if (context.getObjects(Ped.class).size() == 0) {
	        	RunEnvironment.getInstance().endRun();
	        }
	        
	        // Iterate over peds and remove them if they have arrive at the destination
	        for (Object p :context.getObjects(Ped.class)) {
	        	Ped P = (Ped) p;
	        	
	        	// Get the geometries in the CRS used for spatial calculations
	        	Geometry dGeom = SpaceBuilder.getAgentGeometry(this.geography, P.destination);
	        	Geometry coordP = SpaceBuilder.getAgentGeometry(this.geography, P);
	        	
	        	// If the pedestrian agent in within the bounds of the destination then remove it from the context as it has reached its destination
	        	if (dGeom.isWithinDistance(coordP, this.arrivalDist)) {
	        		PedsToRemove.add(P);
	        		break; // End the iteration, only one pedestrian can be removed at a time
	        	}
	        }
	        // Now iterate over all of the peds to remove and remove them from the context
	        // Need to do this separately from iterating over the peds in the context since altering the context whilst iterating over it throws and exception
	        for (Ped P : PedsToRemove) {
	        	context.remove(P);
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
	
	public void setGeography(Geography<Object> G) {
		this.geography = G;
	}
	
	public void setContext(Context<Object> C) {
		this.context = C;
	}

}
