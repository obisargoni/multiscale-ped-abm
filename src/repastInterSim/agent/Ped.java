package repastInterSim.agent;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.geotools.coverage.grid.GridCoordinates2D;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;

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
    		this.pathFinder.getTacticalPath().updateTargetCoordiante();
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
        Iterable<Object> mobileAgentsInArea = SpaceBuilder.geography.getObjectsWithin(fieldOfVisionApprox.getEnvelopeInternal());        
        List<PedObstruction> pedObsInSearchArea = SpatialIndexManager.findIntersectingObjects(SpaceBuilder.pedObstructGeography, fieldOfVisionApprox);
        
        // Calculate acceleration due to field of vision consideration
        fovA = motiveAcceleration(mobileAgentsInArea, pedObsInSearchArea);

        // To Do: Calculate acceleration due to avoiding collisions with other agents and walls.
        contA = totalContactAcceleration();
        
        totA = Vector.sumV(fovA, contA);
        
        return totA;
    }
    
    
    // Calculate the acceleration towards the destination accounting for objects in the field of vision but not collisions
    public double[] motiveAcceleration(Iterable<Object> maObjs, Iterable<PedObstruction> pedObstObjs)  {
    	
    	double[] desiredVelocity = desiredVelocity(maObjs, pedObstObjs);
    	
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
    public double[] totalContactAcceleration()  {
    	double[] cATotal = {0,0};
    	
    	// Get the geometry  and context of the ego agent
    	Geometry thisGeom = GISFunctions.getAgentGeometry(SpaceBuilder.geography, this);
    	
    	// Iterate over all other pedestrian agents and for those that touch the 
    	// ego agent calculate the interaction force
    	// Check to see if this line intersects with any agents
        for (Object agent :SpaceBuilder.context.getObjects(Ped.class)) {
        	Ped P = (Ped)agent;
        	if (P != this) {
               	Geometry agentG = GISFunctions.getAgentGeometry(SpaceBuilder.geography, P);
               	if (agentG.intersects((thisGeom))) {
               		double[] pCA = pedestrianContactAcceleration(this, P, agentG);
               		cATotal = Vector.sumV(cATotal, pCA);
               	}
        	}
        }

        for (Object obstr :SpaceBuilder.pedObstructGeography.getAllObjects()) {
        	PedObstruction Obstr = (PedObstruction)obstr;
           	Geometry obstrGeom = Obstr.getGeom();
           	if (obstrGeom.intersects((thisGeom))) {
           		double[] oCA = obstructionContactAcceleration(this, thisGeom, Obstr, obstrGeom);
           		cATotal = Vector.sumV(cATotal, oCA);
           	}
        } 
    	
    	return cATotal;
    }
    
    public double[] pedestrianContactAcceleration(Ped egoPed, Ped agentPed, Geometry agentGeom) {
    	
    	// Get the radius of the circles representing the pedestrians and the distance between the circles' centroids
    	double r_i = egoPed.rad;
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
    
    public double[] obstructionContactAcceleration(Ped egoPed, Geometry egoGeom, PedObstruction Obstr, Geometry obstrGeom) {
    	
    	// Get the radius of the circles representing the pedestrians and the distance between the circles' centroids
    	double r_i = egoPed.rad;
    	
    	Geometry obstIntersection = egoGeom.intersection(obstrGeom);
    	Coordinate intersectionCoord = obstIntersection.getCentroid().getCoordinate();
    	double d_ij = maLoc.distance(intersectionCoord);
    	
    	// Get the vector that points from centroid of other agent to the ego agent,
    	// this is the direction that the force acts in.
    	// This should also be perpendicular to the obstacle 
    	double[] n = {maLoc.x - intersectionCoord.x, maLoc.y - intersectionCoord.y};
    	n = Vector.unitV(n);
    	
    	double magA = this.k * (r_i - d_ij) / this.m;
    	double[] A  = {magA*n[0], magA*n[1]};
    	
    	return A;
    	
    }
    
    // Function to sample field of vision
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
    
    // Function to calculate distance to nearest collision for a given angle f(a) -  this will need to account for movements of other peds
    public double distanceToObject(double alpha)  {
    	  	
    	double d = distanceToObject(alpha, SpaceBuilder.context.getObjects(Ped.class), SpaceBuilder.pedObstructGeography.getAllObjects());
        
        return d;    	
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
     * 		Distance to the nearest object in the direction of th angle alpha
     */
    public double distanceToObject(double alpha, Iterable<Object> maObjs, Iterable<PedObstruction> pedObstObjs)  {
    	
    	// Initialise distance to nearest object as the max distance in the field of vision
    	double d = this.dmax;
    	
    	LineString sampledRay = GISFunctions.linestringRay(maLoc, alpha, dmax);
    	
    	// Check to see if this line intersects with any pedestrian agents
        for (Object agent :maObjs) {
        	if (agent != this) {
               	Geometry agentG = GISFunctions.getAgentGeometry(SpaceBuilder.geography, agent);
               	if (agentG.intersects(sampledRay)) {
               		// The intersection geometry could be multiple points.
               		// Iterate over them find the distance to the nearest pedestrian
               		Coordinate[] agentIntersectionCoords = agentG.intersection(sampledRay).getCoordinates();
               		for(Coordinate c: agentIntersectionCoords) {
                   		double dAgent = maLoc.distance(c);
                   		if (dAgent < d) {
                   			d = dAgent;
                   		}
               		}
               	}
        	}
        }
        
    	// Check to see if this line intersects with any obstacle agents
        for (PedObstruction obstr : pedObstObjs) {
           	Geometry obstG = obstr.getGeom();
           	if (obstG.intersects(sampledRay)) {
           		// The intersection geometry could be multiple points.
           		// Iterate over them and take the smallest distance - this is the distance to the nearest obstacle
           		Coordinate[] obstIntersectionCoords = obstG.intersection(sampledRay).getCoordinates();
           		for(Coordinate c: obstIntersectionCoords) {
           			double dAgent = maLoc.distance(c);
               		if (dAgent < d) {
               			d = dAgent;
               		}
           		}
           	}
        }
        
        return d;    	
    }
    
    // Function to calculate d(a) using cos rule
    public double displacementDistance(double alpha, Iterable<Object> maObjs, Iterable<PedObstruction> pedObstObjs)  {
    	
    	// Get the distance to nearest object for this angle
    	double fAlpha =  distanceToObject(alpha, maObjs, pedObstObjs);
    	
    	double dAlpha = Math.pow(this.dmax, 2) + Math.pow(fAlpha, 2) - 2*this.dmax*fAlpha*Math.cos(this.a0 - alpha);
    	
    	return dAlpha;
    }
    
    // Wrapper function that identifies the chosen walking direction
    public Map<String, Double> desiredDirection(Iterable<Object> maObjs, Iterable<PedObstruction> pedObstObjs)  {
    	
    	// Sample field of vision
    	List<Double> sampledAngles = sampleFoV();
    	
    	// Initialise the displacement distance (which must be minimised) and the direction of travel
    	// The angle here is relative to the direction of the agent
    	double d = displacementDistance(sampledAngles.get(0), maObjs, pedObstObjs);
    	double alpha = sampledAngles.get(0);    
    	
    	// Loop through the remaining angles and find the angle which minimises the displacement distance
    	for (int i = 1;i<sampledAngles.size(); i++) {
    		
    		double dDist = displacementDistance(sampledAngles.get(i), maObjs, pedObstObjs);
    		
    		if (dDist < d) {
    			d = dDist;
    			alpha = sampledAngles.get(i);
    		}
    	}
    	
    	Map<String, Double> output = new HashMap<String, Double>();
    	output.put("angle", alpha);
    	output.put("collision_distance", distanceToObject(alpha)); // don't need to calculate distanceToObject again here, can use the minimum distance found so far, d
    	
    	return output;    	
    }
    
    public double[] desiredVelocity(Iterable<Object> maObjs, Iterable<PedObstruction> pedObstObjs)  {
    	
    	// Get the desired direction of travel and minimum distance to collision in that direction
    	Map<String, Double> desiredDirection = desiredDirection(maObjs, pedObstObjs);
    	
    	// Calculate the desired speed, minimum between desired speed and speed required to avoid colliding
    	double desiredSpeed = Math.min(this.v0, (desiredDirection.get("collision_distance") - this.rad) / this.tau);
    	
    	// Get the desired direction for the pedestrian and use to set velocity
    	double alpha = desiredDirection.get("angle");
    	double[] v = {desiredSpeed*Math.sin(alpha), desiredSpeed*Math.cos(alpha)};
    	
    	return v;
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
	
}
