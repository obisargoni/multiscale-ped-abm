package repastInterSim.agent;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.geotools.coverage.grid.GridCoordinates2D;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.operation.distance.DistanceOp;

import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.Network;
import repastInterSim.environment.OD;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.Vector;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.PedPathFinder;

public class Ped extends MobileAgent {
	private static Logger LOGGER = Logger.getLogger(Ped.class.getName());

    private PedPathFinder pathFinder;
        
    private double k; // Constant related to the interaction between agents and the desired velocity of this agent
    
    // Variables related to the pedestrians vision and movements
    private double theta; // Field of vision extends from -theta to + theta from the normal to the agent (ie agent's direction)
    private double m; // Agent's mass
    private double dmax; // Maximum distance within which object impact pedestrian movement, can be though of as horizon of field of vision
    private double v0; // Desired walking speed of pedestrian agent
    private double a0; // Angle to the destination
    private double angres; // Angular resolution used when sampling the field of vision
    private double[] v, newV; // Velocity and direction vectors
    private double rad; // Radius of circle representing pedestrian, metres
    
    // Variables used in accumulator model of crossing choice
	private double lambda; // Used to control effect of salience distance on contribution of option utility to activation
	private double alpha; // Controls sensitivity to traffic exposure
	private double gamma; // Controls the rate at which historic activations decay
	private double epsilon; // Proportion of median activation that ca activation must be to be considered dominant

    private List<Coordinate> pedPrimaryRoute; // The primary route are the coordinates the pedestrian commits to the route when first added to the model
    private List<GridCoordinates2D> nextPathSection;
    
    private boolean enteringCrossing = false; // Indicates whether the pedestrian agent should interact with vehicle agents to determine whether to proceed
    private boolean yieldAtCrossing = false; // Indicates whether the pedestrian agent is in a yield state or not, which determines how they move
    
    private String roadLinkFID = null;
    
    private Double pHorizon; // Tactical planning horizon of ped agent in degrees
        
    private int yieldTime = 0;
    
    private Color col; // Colour of the pedestrian
    
    private String chosenCrossingType = "none";
    
    private int stepsSinceReachedTarget = 0; // Counter used to identify when peds get struck and remove them from the simulation.
    
    /*
     * Instance method for the ped class that sets the ped speed and mass to be the default (average) values
     */
    public Ped(OD o, OD d, Double alpha, Double lambda, Double gamma, Double epsilon, boolean minimiseCrossings, Geography<Junction> paveG, Network<Junction> paveNetwork) {
    	super(o,d);
    	init(GlobalVars.pedVavg, GlobalVars.pedMassAv, alpha, lambda, gamma, epsilon, minimiseCrossings, GlobalVars.deafultTacticalPlanningHorizon, paveG, paveNetwork);
    }
    
    /*
     * Instance method for the Ped class.
     * 
     * @param space the continuous space the Ped exists in
     * @param direction the pedestrian's direction
     */
    public Ped(OD o, OD d, Double s, Double m, Double alpha, Double lambda, Double gamma, Double epsilon, boolean minimiseCrossings, Double pH, Geography<Junction> paveG, Network<Junction> paveNetwork) {
    	super(o, d);
    	init(s, m, alpha, lambda, gamma, epsilon, minimiseCrossings, pH, paveG, paveNetwork);
    }
    
    private void init(Double s, Double m, Double alpha, Double lambda, Double gamma, Double epsilon, boolean minimiseCrossings, Double pH, Geography<Junction> paveG, Network<Junction> paveNetwork) {
        this.v0  = s;
        this.m  = m;
        this.rad = m / 320; // As per Moussaid 2011
        
        // Set the pedestrian velocity - half of max velocity in the direction the pedestrian is facing
        double[] v = {0,0};
        this.v =  v;

        this.tau  = 0.5/GlobalVars.tStep; // as per Moussaid 2011
        this.dmax = 10/GlobalVars.spaceScale; // as per Moussaid 2011
        this.angres = (2*Math.PI) / (36 / 3); // Equivalent to 30 degrees
        this.theta = (2*Math.PI*75) / 360; // 75 degrees, as per Moussaid 2011
        this.k = GlobalVars.interactionForceConstant;
        
        // Set the tactical planning horizon
        this.pHorizon = pH;
        
        // Set parameters for crossing choice model (operational level path finding)
		this.alpha = alpha;
		this.lambda = lambda;
		this.gamma = gamma;
		this.epsilon = epsilon;
		
		this.pathFinder = new PedPathFinder(this, paveG, paveNetwork, minimiseCrossings);
    }

	/*
     * Calculate the pedestrian's acceleration and resulting velocity
     * given its location, north and destination.
     */
    @ScheduledMethod(start = 1, interval = 1, priority = 2)
    public void step() throws Exception {        
    	
    	this.stepsSinceReachedTarget++;
    	
    	// Decide yield process involves checking route for road crossing coordinates. Needs to happen before agents updates
    	// its route coordinate because this involves removing coordinates from the route.
   		//decideYield();
    	
    	// If the ped has reached the end of its tactical path, this update will result in a null current junction
    	// In this case need to update the tactical path
    	Boolean tacticalCoordUpdateRequired = false;
    	if(this.pathFinder.getTacticalPath().getCurrentJunction()==null) {
    		tacticalCoordUpdateRequired = true;
    	}

    	
   		// If agent does not intend to yield, agent walks and, if a route coordinate is reached, updates list of route coordinates
   		if (!this.yieldAtCrossing) {
        	if (tacticalCoordUpdateRequired) {
        		pathFinder.updateTacticalPath();
        	}
        	walk(pathFinder.getTacticalPath().getTargetCoordinate());
        	pathFinder.step();
    	}
   		
   		// If agent does intend to yield, agent walks as usual until the crossing point is reached
   		// Once crossing point is reached agent does not move whilst in yield state
   		// Separation here between intention to yield and performing of yielding action (during which intention could be allowed to change, in principle)
    	else if (this.yieldAtCrossing) {
    		double distanceToCrossing = this.maLoc.distance(pathFinder.getNextCrossingCoord());
        	if (distanceToCrossing > 2) {
            	// Walk towards the next coordinate along the route
            	if (tacticalCoordUpdateRequired) {
            		pathFinder.updateTacticalPath();
            	}
            	walk(pathFinder.getTacticalPath().getTargetCoordinate());
            	pathFinder.step();

        	}
        	else {
        		assert true;
        	}
    	}
   		
    	// Finally update the target coordinate if current target coordinate has been reached
    	if (this.maLoc.distance(this.pathFinder.getTacticalPath().getTargetCoordinate()) < 0.5) {
    		this.stepsSinceReachedTarget=0;
    		this.pathFinder.getTacticalPath().updateTargetCoordiante();
    	}
    	
    	// Remove stuck agents from the simulation
    	if (this.stepsSinceReachedTarget>GlobalVars.stuckPedNSteps) {
    		LOGGER.log(Level.FINE, "Removed stuck ped. Origin ID: " + this.origin.getFID() + " Dest ID: " + this.destination.getFID());
    		this.tidyForRemoval();
    		SpaceBuilder.context.remove(this);
    	}
    }
    
    public void walk(Coordinate dLoc) {
    	
    	// Update pedestrians knowledge of which direction the current destination is
        this.a0 = getBearingToDestinationCoord(dLoc);
    	
        // Assume that velocity is updated instantly, and then position is update with new velocity since (walking takes time but adjusting velocity is quick)
        double[] a = accel();
        double[] dv = {a[0]*this.tau, a[1]*this.tau};
        this.newV  = Vector.sumV(v,dv);
        this.v = newV;

        maLoc.x += this.v[0]*this.tau;
        maLoc.y += this.v[1]*this.tau;
        
        // Now create new geometry at the location of the new centroid
        Point pt = GISFunctions.pointGeometryFromCoordinate(maLoc);
		Geometry pGeomNew = pt.buffer(this.rad);
        
        // Move the agent to the new location. This requires transforming the geometry 
        // back to the geometry used by the geography, which is what this function does.
        GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, pGeomNew, this);
        
        // Update the coordinate after moving the pedestrian
        setLoc();
        
        // Set the direction the pedestrian faces to be the direction of its velocity vector
        setPedestrianBearingFromVelocity(this.v);
    }
   
    /*
     * Calculate the acceleration of the pedestrian.
     * 
     * @param location ndpoint representing the pedestrian's location
     * @param north Vector indicating the direction against which bearings are taken
     * @param endPt ndpoint representing the pedestrian's destination 
     * 
     * @return a double representing the pedestrian's new acceleration
     */
    public double[] accel()  {
        
        double[] totA, fovA, contA;
        
        // Start by finding obstacle objects (Peds, Vehicles, PedObstructions close to the ped
        Polygon fieldOfVisionApprox = getPedestrianFieldOfVisionPolygon(this.a0);
        
        List<Geometry> obstructionGeoms = SpatialIndexManager.findIntersectingGeometries(SpaceBuilder.pedObstructGeography, fieldOfVisionApprox, "intersects");
        HashMap<Ped, Geometry> pedsWithGeoms = getFOVPedsAndGeoms(fieldOfVisionApprox);
        
        // Calculate acceleration due to field of vision consideration
        fovA = motiveAcceleration(obstructionGeoms, pedsWithGeoms.keySet());

        // Calculate acceleration due to avoiding collisions with other agents and walls.
        contA = totalContactAcceleration(obstructionGeoms, pedsWithGeoms);
        
        totA = Vector.sumV(fovA, contA);
        
        return totA;
    }

	// Calculate the acceleration towards the destination accounting for objects in the field of vision but not collisions
    public double[] motiveAcceleration(Iterable<Geometry> obstGeoms, Iterable<Ped> peds)  {
    	
    	double[] desiredVelocity = desiredVelocity(obstGeoms, peds);
    	
    	// Acceleration is set as the acceleration required to reach the desired velocity within tau amount of time
    	double[] a = {0,0};
    	a[0] = (desiredVelocity[0] - this.v[0]) / this.tau;
    	a[1] = (desiredVelocity[1] - this.v[1]) / this.tau;
    	
    	return a;
    }
    
    /* 
     * Calculates the acceleration due to contact with other agents.
     * performs spatial query to identify pedestrian agents that are in 
     * contact with the ego agent. Calculates the force due to each contacting 
     * agents and sums the forces and divides by the ego agent's mass to produce
     * the acceleration. 
     */
    public double[] totalContactAcceleration(List<Geometry> obstrGeoms, HashMap<Ped, Geometry> peds)  {
    	double[] cATotal = {0,0};
    	
    	// Get the geometry  and context of the ego agent
    	Geometry thisGeom = GISFunctions.getAgentGeometry(SpaceBuilder.geography, this);
    	
    	// Iterate over all other pedestrian agents and for those that touch the 
    	// ego agent calculate the interaction force
    	// Check to see if this line intersects with any agents
        for (Entry<Ped, Geometry> entry: peds.entrySet()) {
        	Ped p = entry.getKey();
        	Geometry pGeom = entry.getValue();
        	DistanceOp distOp = new DistanceOp(thisGeom, pGeom);
           	if (Double.compare(distOp.distance(), 0.0)==0) {
           		double[] pCA = pedestrianContactAcceleration(p, pGeom);
           		cATotal = Vector.sumV(cATotal, pCA);
           	}
        }

        for (Geometry obstrGeom :obstrGeoms) {
        	Coordinate[] intersectionCoords = thisGeom.intersection(obstrGeom).getCoordinates();
           	if (intersectionCoords.length>0) {
           		double[] oCA = obstructionContactAcceleration(thisGeom, intersectionCoords);
           		cATotal = Vector.sumV(cATotal, oCA);
           	}
        } 
    	
    	return cATotal;
    }
    
    public double[] pedestrianContactAcceleration(Ped agentPed, Geometry agentGeom) {
    	
    	// Get the radius of the circles representing the pedestrians and the distance between the circles' centroids
    	double r_i = this.rad;
    	double r_j = agentPed.rad;
    	
    	Coordinate agentCoord = agentGeom.getCentroid().getCoordinate();
    	double d_ij = maLoc.distance(agentCoord);
    	
    	// Get the vector that points from centorid of other agent to the ego agent,
    	// this is the direction that the force acts in
    	double[] n = {maLoc.x - agentCoord.x, maLoc.y - agentCoord.y};
    	n = Vector.unitV(n);
    	
    	double magA = this.k * (r_i + r_j - d_ij) / this.m;
    	double[] A  = {magA*n[0], magA*n[1]};
    	
    	return A;
    }
    
    public double[] obstructionContactAcceleration(Geometry egoGeom, Coordinate[] intCoords) {
    	
    	// Get the radius of the circles representing the pedestrians and the distance between the circles' centroids
    	double r_i = this.rad;
    	
    	// Find midpoint between intersecting coords
    	int nCoords = intCoords.length;
    	Coordinate midPoint = GISFunctions.midwayBetweenTwoCoordinates(intCoords[0], intCoords[nCoords-1]);
    	double d_ij = maLoc.distance(midPoint);
    	
    	// Get the vector that points from centroid of other agent to the ego agent,
    	// this is the direction that the force acts in.
    	// This should also be perpendicular to the obstacle 
    	double[] n = {maLoc.x - midPoint.x, maLoc.y - midPoint.y};
    	n = Vector.unitV(n);
    	
    	double magA = this.k * (r_i - d_ij) / this.m;
    	double[] A  = {magA*n[0], magA*n[1]};
    	
    	return A;
    	
    }
    
    public double[] desiredVelocity(Iterable<Geometry> obstGeoms, Iterable<Ped> peds)  {
    	
    	// Get the desired direction of travel and minimum distance to collision in that direction
    	Map<String, Double> desiredDirection = desiredDirection(obstGeoms, peds);
    	
    	// Calculate the desired speed, minimum between desired speed and speed required to avoid colliding
    	double desiredSpeed = Math.min(this.v0, (desiredDirection.get("collision_distance") - this.rad) / this.tau);
    	
    	// Get the desired direction for the pedestrian and use to set velocity
    	double alpha = desiredDirection.get("angle");
    	double[] v = {desiredSpeed*Math.sin(alpha), desiredSpeed*Math.cos(alpha)};
    	
    	return v;
    }
    
    // Wrapper function that identifies the chosen walking direction
    public Map<String, Double> desiredDirection(Iterable<Geometry> obstGeoms, Iterable<Ped> peds)  {
    	
    	// Then sample field of vision and initialise arrays of distances that correspond to each angle 
    	List<Double> sampledAngles = sampleFoV();
    	int nAngleOptions = sampledAngles.size()-1; // Each sector between a sampled angle results in an angle option for direction
    	double[] distances = new double[sampledAngles.size()];
    	double[] displacementDistances = new double[sampledAngles.size()];
    	
    	// Initialise values as dmax since this is limit of peds field of vision
    	for (int i=0; i<distances.length; i++) {
    		distances[i] = this.dmax;
    	}
    	
    	// First find the minimum displacement distance and associated angle for obstruction geometries
    	displacementDistancesToObstacleGeometries(obstGeoms, sampledAngles, distances, displacementDistances);
    	displacementDistancesToPeds(peds, sampledAngles, distances, displacementDistances);
    	
    	Integer lowi = null;
    	double minDD = Double.MAX_VALUE;
    	
    	// Now loop through the remaining angles and calculate displacement distances for any angles that don't have obstacles in
    	for (int i = 0;i<nAngleOptions; i++) {
    		if (displacementDistances[i] < 0) {
    			distances[i] = this.dmax;
    			displacementDistances[i] = displacementDistance(sampledAngles.get(i), distances[i]);
    		}
    		
    		if (displacementDistances[i] < minDD) {
    			minDD = displacementDistances[i];
    			lowi = i;
    		}
    	}
    	
    	Map<String, Double> output = new HashMap<String, Double>();
    	output.put("angle", sampledAngles.get(lowi));
    	output.put("collision_distance", distances[lowi]);
    	
    	return output;    	
    }

	/*
     * Given a set of geometries, calculate displacement distance for each of the geometries. Fill arrays of distances and displacement distances 
     * with the lowest distance at that angle.
     * 
     * This is an alternative approach to finding minimum displacement distance that avoids many intersection computations.
     * 
     * @param Iterable<Geometry> obstGeoms
     * 		The obstacle geometries to calculate the displacement distances for
     * 
     * @return double[]
     * 		Array. First element is distance to the object. Second is displacement distance. Third is angle to geometry that gives least displacement distance. 
     */
    public void dispalcementDistancesToPointGeometries(Iterable<Geometry> obstGeoms, List<Double> fovAngles, double[] ds, double[] dds) {
    	double[] output = new double[3];
    	output[1] = Double.MAX_VALUE;
    	Geometry agentG = GISFunctions.getAgentGeometry(SpaceBuilder.geography, this);
    	for (Geometry g: obstGeoms) {
    		DistanceOp distOp = new DistanceOp(agentG, g);
    		
    		// Calculate angle
    		double alphaToGeom = GISFunctions.bearingBetweenCoordinates(maLoc, g.getCentroid().getCoordinate());
    		
    		// Check whether angle lies within field of vision and find which field of vision angle sample this angle corresponds to
    		
    		// Translate angle to be relative to peds bearing and in range -pi -> pi
    		double alpha = alphaToGeom - this.bearing;
    		if (alpha > Math.PI) {
    			alpha = alpha - 2*Math.PI;
    		}
    		else if(alpha < -Math.PI) {
    			alpha = alpha + 2*Math.PI;
    		}
    		
    		assert (alpha >= -Math.PI) & (alpha <= Math.PI);
    		
    		// Check whether angle lies within field of vision, if not continue to next geom
    		if ( (alpha < -this.theta) | (alpha > this.theta) ) {
    			continue;
    		}
    		
    		// Get sample index of angle
    		int ai = (int) ( (alpha - (- this.theta)) / this.angres);
    		    		
    		// Calculate distance - do I need to limit to dmax?
    		double fAlpha = distOp.distance();
    		
    		// Calculate displacement distance. Use the actual angle to the object.
    		double dAlpha = displacementDistance(alphaToGeom, fAlpha);
    		
    		// Check if displacement distance is lower than previous min
    		if (fAlpha <= this.dmax) { // if object within field of vision
	    		if (dds[ai] < 0) { // if value not yet set for this angle
	    			fovAngles.set(ai, alphaToGeom);
	    			ds[ai] = fAlpha;
	    			dds[ai] = dAlpha;
	    		}
	    		else if (fAlpha < ds[ai]) { // If this geometry is closer to the pedestrian
	    			fovAngles.set(ai, alphaToGeom);
	    			ds[ai] = fAlpha;
	    			dds[ai] = dAlpha;
	    		}
    		}
    	}
    }
    
    /*
     * Calculate the distance and displacement distance to each ped in the input iterable.
     * 
     * Distance are the distance of first collision with the ped if the ego ped moves at it's preferred walking
     * speed at the sample angle and the other ped continues with its current velocity.
     * 
     * @param Iterable<Ped>
     * 		The peds to consider when calculating distances.
     * @param List<Double> sampledAngles
     * 		The bearings to consider the ego ped moving in.
     * @param double[] distances
     * 		The distance of nearest collision for each bearing.
     * @param double[] displacementDistances
     * 		The displacement distances for each bearing.
     */
    private void displacementDistancesToPeds(Iterable<Ped> peds, List<Double> sampledAngles, double[] distances,
			double[] displacementDistances) {
		for (Ped p: peds) {
			// Calculate difference in current position
			double[] dR = {this.maLoc.x - p.getLoc().x, this.maLoc.y - p.getLoc().y};
			
			// Iterate over possible directions of travel
			for (int i=0; i<sampledAngles.size(); i++) {
				double b = sampledAngles.get(i);
				double[] vThis = {this.v0*Math.sin(b), this.v0*Math.cos(b)};
				double[] dV = {vThis[0] - p.getV()[0], vThis[1] - p.getV()[1]};
				
				// Calculate time of closest approach
				double tClosest = -(dR[0] *dV[0] + dR[1] * dV[1]) / ( dV[0]*dV[0] + dV[1]*dV[1]);				
				
				// If tClosest is in the past set time of closest approach to now
				if (tClosest<0) {
					tClosest = 0;
				}
				
				// Calculate distance of closest approach
				double[] futurePLoc = {p.getLoc().x+p.getV()[0]*tClosest, p.getLoc().y+p.getV()[1]*tClosest};
				double[] futureThisLoc = {maLoc.x+this.v0*Math.sin(b)*tClosest, maLoc.y+this.v0*Math.cos(b)*tClosest};
				double dClosest = Math.sqrt( Math.pow(futureThisLoc[0] - futurePLoc[0], 2) +  Math.pow(futureThisLoc[1] - futurePLoc[1], 2) );
				
				// If peds collide on this course find distance this ped can travel in direction b until collision
				// If this is less that the current distance set for this angle, update the distance and displacement distance for this angle
				if (dClosest<0) {
					double fAlpha = this.v0*tClosest;
					
					// Assumes that angles without a detected obstruction have dmax entered as distance
					if (fAlpha < distances[i]) {
						distances[i] = fAlpha;
						displacementDistances[i] = displacementDistance(b, fAlpha);
					}
						
				}
			}
		}
		
	}
    
    /*
     * Function to calculate distance to obstruction geometries in each of the bearing directions given in the input sampled angles.
     * 
     * @param Iterable<Geometry> obstGeoms
     * 		The obstruction geometries to check distances to
     * @param List<Double> sampledAngles
     * 		The bearings to consider the ego ped moving in.
     * @param double[] distances
     * 		The distance of nearest collision for each bearing.
     * @param double[] displacementDistances
     * 		The displacement distances for each bearing.
     */
    public void displacementDistancesToObstacleGeometries(Iterable<Geometry> obstGeoms, List<Double> fovAngles, double[] ds, double[] dds)  {    	
    	
    	for(int i=0; i<fovAngles.size(); i++) {
        	// Get the distance to nearest object for this angle
        	double fAlpha =  distanceToObject(fovAngles.get(i), obstGeoms);
        	double dAlpha = displacementDistance(fovAngles.get(i), fAlpha);
        	
        	ds[i] = fAlpha;
        	dds[i] = dAlpha;
    	}
    }
    
    /*
     * Function to calculate distance to nearest collision with objects passed in as iterables, for a given angle alpha.
     * 
     * Correction needed to account for movements of other peds.
     * 
     * @param double alpha
     * 		The bearing to look for objects along
     * @param Iterable<Object> maObjs
     * 		An iterable containing mobile agents (Peds and Vehicle). Search through these to avoid colliding with Peds and Vehicles
     * @param Iterable<PedObstruction> pedObstObjs
     * 		An iterable containing PedObstruction objects. Search through these to avoid colliding with PedObstructions
     * 
     * @return double
     * 		Distance to the nearest object in the direction of the angle alpha
     */
    public double distanceToObject(double alpha, Iterable<Geometry> obstGeoms)  {
    	
    	// Initialise distance to nearest object as the max distance in the field of vision
    	double d = this.dmax;
    	
    	LineString sampledRay = GISFunctions.linestringRay(maLoc, alpha, dmax);
    	
    	// Check to see if this line intersects with any pedestrian agents
        for (Geometry obstG :obstGeoms) {
               	
           	DistanceOp distOP = new DistanceOp(obstG, sampledRay);
           	// This check is equivalent to agentG.intersects(sampledRay) (tested this using assertions in a simulation run) but is slightly faster.
           	int i=0;
           	while (Double.compare(distOP.distance(), 0.0) == 0) {
           		// DistanceOp can be used to approximate intersection.
           		// When distance is zero geometries overlapp. The nearest point between the two geometries are the intersection.
           		// Distance Op approximately returns one of the intersection points
           		// Since there might be multiple, need to create a new ray that extends just up to the first point distance op found
           		// and re run until the distance is no zero, meaning that the nearest intersecting coord is found.
           		Coordinate intersectingCoord = distOP.nearestPoints()[1];
           		d = maLoc.distance(intersectingCoord);
           		
           		sampledRay = GISFunctions.linestringRay(maLoc, alpha, d*0.9999);
           		distOP = new DistanceOp(obstG, sampledRay);
           		i++;
           		
           		// Add break clause to avoid infinite loop occuring
           		if (i>100) {
           			d = 0;
           			break;
           		}
           	}              	
        }
        
        return d;    	
    }
    
    /*
     * Sample angles in field of vision
     */
    public List<Double> sampleFoV() {
    	
    	// Initialise a list to hole the sampled field of vision vectors
    	List<Double> sampledAngles = new ArrayList<Double>();
    	
    	double sampleAngle = this.bearing-this.theta; // First angle to sample
    	double sampleAnglemax = this.bearing + this.theta;
    	while (sampleAngle <= sampleAnglemax) {
    		sampledAngles.add(sampleAngle);
    		sampleAngle+=this.angres;
    	}
    	
    	return sampledAngles;
    	
    }
    
    /*
     * Calculate the displacement distance in a direction given the direction angle, angle to destination, and distance to nearest object in direction angle 
     */
    public double displacementDistance(double alpha, double fAlpha) {
    	return Math.pow(this.dmax, 2) + Math.pow(fAlpha, 2) - 2*this.dmax*fAlpha*Math.cos(this.a0 - alpha);
    }
    
    
    public double[] limitV(double[] input) {
        double totalV = Vector.mag(input);
        
        if (totalV > v0) {
        	double norm = v0/totalV;
            input[0] = input[0]*norm;
            input[1] = input[1]*norm;}
        return input;
    }
    
    /**
     * Determines whether the agents is close to a crossing point by comparing grid coverage values,
     * used for routing, at the current location and estimates near future location.
     * 
     * If grid cell values are expected to increase the agent is heading towards a lower priority area
     * and is considered to be approaching a crossing point.
     */
    public void decideYield() {
    	
    	checkIfEnteringCrossing();
    	
    	if (this.enteringCrossing) {
    		
    		// Decide whether to yield or not. For now yield for fixed amount of time by default
    		if (this.yieldTime < 10) {
            	this.yieldAtCrossing = true;
            	this.yieldTime++;
        	}
        	
    		// If not yielding allow agent to progress by setting enteringCrossing state to false and removing the crossing coord from list of crossing coords
	    	else {
	    		this.yieldAtCrossing = false;
	    		this.yieldTime = 0;
				this.enteringCrossing = false;
				
				// Set crossing coord to null because agent has stopped yielding so this crossing can be considered to be passed
				pathFinder.setNextCrossingCoord(null);
	    	}
    	}
    }
    
    /*
     * Use ray that extends ahead of pedestrian to check whether pedestrian is about to enter the carriageway.
     * 
     * The method not needed because AccumulatorRoute records when pedestrian is entering the carriadgeway using the route crossing coordinates.
     */
    public void checkIfEnteringCrossing() {    	
    	// Crossing coord could be null if there isn't a crossing coming up on the route
    	if (pathFinder.getNextCrossingCoord() != null) {
        	Point nextCrossing = GISFunctions.pointGeometryFromCoordinate(pathFinder.getNextCrossingCoord());
        	
        	Coordinate lookAhead  = getPedestrianLookAheadCoord(GlobalVars.lookAheadTimeSteps);
        	Coordinate[] lookAheadLineCoords = {this.maLoc, lookAhead};
        	Geometry lookAheadLine = new GeometryFactory().createLineString(lookAheadLineCoords).buffer(GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.SMALLPLUS.dist);
        	
        	// Ped agent considered to be about to enter crossing if crossing point is close to expected future location 
        	// Expected future location is only about 1m ahead, may need to rethink how its calculated
        	// Also the buffer distance also needs considering
        	if (lookAheadLine.intersects(nextCrossing)) {
        		this.enteringCrossing = true;
        	}
        	else {
        		this.enteringCrossing = false;
        	}
    	}
    	else {
    		this.enteringCrossing = false;
    	}
    }	
    
    /**
     * Method to be run when agent is removed from the context.
     */
    public void tidyForRemoval() {
    	//IO.exportPedGridRouteData(this, "final_", false);
    }
    
    public Color getColor() {
    	return this.col;
    }
    
    /*
     * Calculate bearing to destination and convert to a unit vector
     */
    public double getBearingToDestinationCoord(Coordinate dLoc)  {
    	
        double[] dirToEnd = {dLoc.x - maLoc.x, dLoc.y - maLoc.y};        
        dirToEnd = Vector.unitV(dirToEnd);
        
        double a0 = Vector.angleBetweenNorthAndUnitVector(dirToEnd);
        
        return a0;
    }
    
    /*
     * Set the direction of the pedestrian to be the same as the direction of the velocity vector
     */
    public double setPedestrianBearingFromVelocity(double[] v) {
    	
    	// If velocity is 0 then don't update the pedestrian direction
    	if (Vector.mag(v)==0) {
    		return this.bearing;
    	}
    	
    	double[] unitV = Vector.unitV(v);
    	
    	this.bearing = Vector.angleBetweenNorthAndUnitVector(unitV);
    	
    	return this.bearing;
    	
    }
    
    /**
     * Set the estimated expected location of the pedestrian agent in a number of timesteps time
     * 
     * @param nTimeSteps
     * 			Integer number of timesteps to estimate location at
     */
    public Coordinate getPedestrianLookAheadCoord(int nTimeSteps) {
    	
    	// aP is the bearing from north represents pedestrian direction.
    	// To get expected location in n timesteps multiply the distance covered in three timesteps
    	// by the bearing resolved in the x and y directions
    	
    	double dx = this.v0*nTimeSteps*GlobalVars.stepToTimeRatio*Math.sin(this.bearing);
    	double dy = this.v0*nTimeSteps*GlobalVars.stepToTimeRatio*Math.cos(this.bearing);
    	
    	Coordinate newLookAhead = new Coordinate(this.maLoc.x + dx, this.maLoc.y + dy);
    	return newLookAhead;
    }
    
    /*
     * Returns a triangle polygon that approximates pedestrian's field of vision. The triangle 
     * is wider and longer than the actual field of vision to ensure than the triangle fully encompasses the field of vision.
     * 
     * @param double b
     * 		The bearing of the pedestrian (direction pedestrian is facing)
     */
    public Polygon getPedestrianFieldOfVisionPolygon(double b) {
        double envAng = Math.min(Math.PI/2, this.theta * 1.1); // Angle of cone, a bit larger than field of vision, but not larger than pi/2
        double r = Math.abs((this.dmax *1.1) / Math.cos(envAng / 2)); // Length of cone side. Chosen as length needed to ensure polygon extends beyond field of vision.
        double a1 = b - envAng;
        double a2 = b + envAng;
        
        Coordinate c1 = new Coordinate(maLoc.x + r*Math.sin(a1), maLoc.y + r*Math.cos(a1));
        Coordinate c2 = new Coordinate(maLoc.x + r*Math.sin(b), maLoc.y + r*Math.cos(b));
        Coordinate c3 = new Coordinate(maLoc.x + r*Math.sin(a2), maLoc.y + r*Math.cos(a2));
        
        Coordinate[] fovCoords = {maLoc, c1, c2, c3, maLoc};
        Polygon fovArea = SpaceBuilder.fac.createPolygon(fovCoords);
        
        return fovArea;
    }
    
    public List<Geometry> getObstacleGeometries(Polygon fieldOfVisionApprox, Geography<PedObstruction> pedObstGeog) {
        // Get list of all geometries of other pedestrian agents and pedestrian obstructions that intersect field of vision
        List<Geometry> obstacleGeoms = SpatialIndexManager.findIntersectingGeometries(pedObstGeog, fieldOfVisionApprox, "intersects");
        Iterable<Ped> pedsInArea = SpaceBuilder.geography.getObjectsWithin(fieldOfVisionApprox.getEnvelopeInternal(), Ped.class);
        for (Ped p: pedsInArea) {
        	if (p != this) {
        		Geometry pGeom = GISFunctions.getAgentGeometry(SpaceBuilder.geography, p);
        		obstacleGeoms.add(pGeom);
        	}
        }
		return obstacleGeoms;
	}
    
    public HashMap<Ped, Geometry> getFOVPedsAndGeoms(Polygon fieldOfVisionApprox) {
    	HashMap<Ped, Geometry> peds = new HashMap<Ped, Geometry>();
        Iterable<Ped> pedsInArea = SpaceBuilder.geography.getObjectsWithin(fieldOfVisionApprox.getEnvelopeInternal(), Ped.class);
        for (Ped p: pedsInArea) {
        	if (p != this) {
        		Geometry pGeom = GISFunctions.getAgentGeometry(SpaceBuilder.geography, p);
        		peds.put(p, pGeom);
        	}
        }
        return peds;
    }
    
    public void setRoadLinkFID(String rlFID) {
    	this.roadLinkFID = rlFID;
    }
    
    public String getRoadLinkFID() {
    	return this.roadLinkFID;
    }
    
    public double getRad() {
    	return this.rad;
    }
    
    public double getSpeed() {
    	return Vector.mag(this.v);
    }
    
    
    /*
     * Set the location attribute of the agent to be the coordinate of its 
     * centroid, in the coordinate reference frame used by the agent for GIS calculations. 
     */
    @Override
    public void setLoc()  {
    	// Get centroid coordinate of this agent
    	Coordinate pL = GISFunctions.getAgentGeometry(SpaceBuilder.geography, this).getCentroid().getCoordinate();
    	this.maLoc = pL;
    }
    
    /*
     * Get the destination of this pedestrian
     * 
     *  @returns
     *  	The Destination object of this pedestrian
     */
    @Override
    public OD getDestination() {
    	return this.destination;
    }
    
    
    public PedPathFinder getPathFinder() {
    	return this.pathFinder;
    }
    
    @Override
    public Geography<Object> getGeography() {
    	return SpaceBuilder.geography;
    }
    
    public String getPrimaryRouteCoordinatesString() {    	
    	String coordString = IO.getCoordinateListString(this.pedPrimaryRoute);
    	return coordString;
    }
    
    public List<GridCoordinates2D> getNextPathSection(){
    	return this.nextPathSection;
    }


	public Double getpHorizon() {
		return pHorizon;
	}


	public void setpHorizon(Double pHorizon) {
		this.pHorizon = pHorizon;
	}


	public double getLambda() {
		return lambda;
	}


	public double getAlpha() {
		return alpha;
	}


	public double getGamma() {
		return gamma;
	}


	public double getEpsilon() {
		return epsilon;
	}


	public String getChosenCrossingType() {
		return chosenCrossingType;
	}
	
	public String getChosenCrossingTypeString() {
		return chosenCrossingType.toString();
	}


	public void setChosenCrossingType(String chosenCrossingType) {
		this.chosenCrossingType = chosenCrossingType;
	}
	
	/*
	 * Flag that indicates whether this pedestrian has chosen a crossing and is therefore attempting to cross the road
	 */
	public boolean isCrossing() {
		return this.pathFinder.getTacticalPath().getAccumulatorRoute().isCrossing();
	}
	
    public double[] getV() {
		return v;
	}
    
    // Used for testing only
    public void setV(double[] v) {
    	this.v = v;
    }
	
}
