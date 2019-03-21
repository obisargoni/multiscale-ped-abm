package repastSocialForce;

import java.awt.Color;
import java.util.ArrayList;

import org.geotools.geometry.jts.JTS;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.Geography;
import repast.simphony.util.ContextUtils;

public class Destination {

	private Geography<Object> geography;
	private Color colour; // Colour of the destination
	private double arrivalDist = 5; // Distance in metres at which agents are considered to have arrived at the destination and are removed from the space
	private MathTransform transformtoMetre;
	private MathTransform transformtoDegree;
		
	public Destination(Geography<Object> geography, Color col, MathTransform ttM, MathTransform ttD) {
		this.geography = geography;
		this.colour = col;
		this.transformtoMetre = ttM;
		this.transformtoDegree = ttD;
	}
	
	@ScheduledMethod(start = 1, interval = 1, priority = ScheduleParameters.LAST_PRIORITY)
	 public void removePeds() throws MismatchedDimensionException, NoSuchAuthorityCodeException, FactoryException, TransformException {
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
	        	
	        	// Get the geometries and transform to projection in metres
	        	
	        	Geometry dGeom = JTS.transform(geography.getGeometry(P.destination), this.transformtoMetre);
	        	Geometry coordP = JTS.transform(geography.getGeometry(P), this.transformtoMetre);
	        	
	        	// If the pedestrian agent in within the bounds of the destination then remove it from the context as it has reached its destination
	        	if (dGeom.isWithinDistance(coordP, this.arrivalDist)) {
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
