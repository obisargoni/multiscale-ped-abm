package repastSocialForce;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.gis.util.GeometryUtil;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;
import repast.simphony.space.gis.Geography;
import repast.simphony.util.ContextUtils;

/*
 * An agent which generates pedestrian agents and adds them to the space
 */
public class Source {
	// static means this variable is shared between all instances of this class
	static int countRemovedPeds = 0;
	private int nPeds; // World dimensions and number of pedestrians to add
	private Geometry boundary;
	private Destination d;
	private Color color; // Colour to use for styling the pedestrians
	private GeometryFactory geoFac;
	
	/*
	 * Instance method for Source. Source also contains method calls to move peds.
	 * Source instantiates a number of pedestrians in the space and assigns them
	 * a destination and a colour.
	 */
	public Source(Geometry bndry, GeometryFactory gF, Destination d, int nP, Color col) {
		
		this.boundary = bndry;
		this.d = d;
		this.nPeds = nP;
		this.color = col;
		this.geoFac = gF;
	}

	
	/*
	 * Add multiple pedestrian agents to the space. This is triggered to 
	 * occur once at the start of the simulation.
	 * 
	 *  In future will want to permit pedestrians to be added throughout the 
	 *  simulation.
	 */
    @ScheduledMethod(start = 1, priority = ScheduleParameters.FIRST_PRIORITY) // 
    public void addPeds() {  
    	
    	// Get the number of pedetrian agent to add to the space from the parameters
    	Parameters  params = RunEnvironment.getInstance().getParameters();
    	int nP = (int)params.getInteger("nPeds");

		// Generate random points in the area to create agents.
		List<Coordinate> agentCoords = GeometryUtil.generateRandomPointsInPolygon(boundary, nP);
		
		// Create the agents from the collection of random coords.
		int cnt=0;
		Random randCoord = new Random();
		for (Coordinate coord : agentCoords) {
			
			// Generate a random initial direction for the pedestrian
    		double randBearing = randCoord.nextFloat() * FastMath.PI * 2;
    		double[] dir = {FastMath.sin(randBearing), FastMath.cos(randBearing)};
			
    		Ped newPed = addPed(dir,coord, d); // GisAgent agent = new GisAgent("Site " + cnt);
			
			cnt++;
		}
    }

    public Ped addPed(double[] direction, Coordinate coord, Destination d) {
    	
        Context<Object> context = ContextUtils.getContext(this);
        Geography<Object> geography = (Geography<Object>)context.getProjection("Geography");
        
        // Instantiate a new pedestrian agent and add the agent to the context
        Ped newPed = new Ped(geography, direction, d, color);
        context.add(newPed);
        
        // Create a new point geometry. Move the pedestrian to this point. In doing so this 
        // pedestrian agent becomes associated with this geometry.
		Point geom = geoFac.createPoint(coord);
		geography.move(newPed, geom);
        	
        return newPed;
    }
    
    //@ScheduledMethod(start = 1, interval = 1, priority = ScheduleParameters.LAST_PRIORITY)
    public void removePeds() {
        Context<Object> context = ContextUtils.getContext(this);
        Geography<Object> geography = (Geography<Object>) context.getProjection("geography");
        
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
        	this.countRemovedPeds ++;	
        }
    }
}
