package repastInterSim.agent;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.geotools.coverage.grid.GridCoordinates2D;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.Geography;
import repast.simphony.util.ContextUtils;
import repastInterSim.environment.Destination;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.Vector;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class Ped implements MobileAgent {
    private Geography<Object> geography; // Space the agent exists in
    public Destination destination; // The destination agent that this pedestrian agents is heading towards.
    
    private GridRoute route;
    
    private Random rnd = new Random(); // Random seed used to give a distribution of velocities 
    
    private double k; // Constant related to the interaction between agents and the desired velocity of this agent
    
    // Variables related to the pedestrians vision and movements
    private double theta; // Field of vision extends from -theta to + theta from the normal to the agent (ie agent's direction)
    private double m; // Agent's mass
    private double dmax; // Maximum distance within which object impact pedestrian movement, can be though of as horizon of field of vision
    private double v0; // Desired walking speed of pedestrian agent
    private double a0; // Angle to the destination
    private double aP; // Angle of pedestrian direction
    private double tau; // Time period in which pedestrian agent is able to come to a complete stop. Used to set acceleration required to avoid collisions.
    private double angres; // Angular resolution used when sampling the field of vision
    private double[] v, newV; // Velocity and direction vectors
    private double rad; // Radius of circle representing pedestrian, metres
    private Coordinate pLoc; // The coordinate of the centroid of the pedestrian agent.
    private Coordinate routeCoord; // The next coordiante of the agent's route
    
    private boolean enteringCrossing = false; // Indicates whether the pedestrian agent should interact with vehicle agents to determine whether to proceed
    private boolean yieldAtCrossing = false; // Indicates whether the pedestrian agent is in a yield state or not, which determines how they move
    
    private String roadLinkFID = null;

    private HashMap<Integer, Double> gridSummandPriorityMap = new HashMap<Integer, Double>(); // Used to get grid cell summand value when running flood fill algorithm for routing
    private double vehiclePriorityCostRatio; // The ratio of pedestrian priority cell cost to vehicle priority cell cost. Represents pedestrian's perception of cost of moving in vehicle priority space.
    
    private int yieldTime = 0;
    
    private Color col; // Colour of the pedestrian
    
    
    /*
     * Instance method for the Ped class.
     * 
     * @param space the continuous space the Ped exists in
     * @param direction the pedestrian's direction
     */
    public Ped(Geography<Object> geography, Geography<Destination> destinationGeography, Destination d) {
        this.geography = geography;
        this.destination = d;
        this.v0  = rnd.nextGaussian() * GlobalVars.pedVsd + GlobalVars.pedVavg;
        this.m  = rnd.nextGaussian() * GlobalVars.pedMasssd + GlobalVars.pedMassAv;
        this.rad = m / 320; // As per Moussaid
        
        // Set the pedestrian velocity - half of max velocity in the direction the pedestrian is facing
        double[] v = {0,0};
        this.v =  v;

        this.tau  = 0.5/GlobalVars.tStep;
        this.dmax = 10/GlobalVars.spaceScale; 
        this.angres = (2*Math.PI) / 36; // Equivalent to 10 degrees
        this.theta = (2*Math.PI*75) / 360; // 75 degrees
        this.k = GlobalVars.interactionForceConstant;
        
        // Set the cost to the agent of moving in pedestrian and vehicle priority areas. Used when running flood fill for routing
        this.vehiclePriorityCostRatio = GlobalVars.MOBILE_AGENT_PARAMS.cautiousPriorityCostRatio;
        this.gridSummandPriorityMap.put(GlobalVars.GRID_PARAMS.getPriorityValueMap().get("pedestrian"), 1.0);
        this.gridSummandPriorityMap.put(GlobalVars.GRID_PARAMS.getPriorityValueMap().get("vehicle"), this.vehiclePriorityCostRatio);

		// Get the destination coordinate, initialise new route and generate a pedestrian route
		Coordinate dCoord = this.destination.getGeom().getCentroid().getCoordinate(); 
		this.route = new GridRoute(geography, this, this.gridSummandPriorityMap, dCoord);
		
    }
    
    
    /*
     * Calculate the pedestrian's acceleration and resulting velocity
     * given its location, north and destination.
     */
    @ScheduledMethod(start = 1, interval = 1, priority = 2)
    public void step() throws Exception {        
    	// Agent decides whether to yield, in which case it doesn't progress to the next route coord
        // Order here assumes that agent does this (which requires evaluating if they are approaching a crossing) before moving
   		decideYield(); 
    	
   		// If agent does not intend to yield, agent walks and, if a route coordinate is reached, updates list of route coordinates
   		if (!this.yieldAtCrossing) {
        	if (this.pLoc.distance(routeCoord) < 0.5) {
        		this.route.removeNextFromRoute();;
        		this.routeCoord = this.route.getNextRouteCoord();
        	}
        	walk(routeCoord);
    	}
   		
   		// If agent does intend to yield, agent walks as usual until the crossing point is reached
   		// Once crossing point is reached agent does not move whilst in yield state
   		// Separation here between intention to yield and performing of yielding action (during which intention could be allowed to change, in principle)
    	else if (this.yieldAtCrossing) {
    		double distanceToCrossing = this.pLoc.distance(this.route.getNextRouteCrossingCoord());
        	if (distanceToCrossing > 2) {
            	// Walk towards the next coordinate along the route
            	if (this.pLoc.distance(routeCoord) < 0.5) {
            		// Get next coordinate
            		this.route.removeNextFromRoute();;
            		this.routeCoord = this.route.getNextRouteCoord();
            	}
            	walk(routeCoord);

        	}
        	else {
        		assert true;
        	}
    	}
    }
    
    public void walk(Coordinate dLoc) {
    	
    	// Update pedestrians knowledge of which direction the current destination is
        setBearingToDestinationCoord(dLoc);
    	
        double[] a = accel();
        double[] dv = {a[0]*this.tau, a[1]*this.tau};
        this.newV  = Vector.sumV(v,dv);
        this.v = newV;

        pLoc.x += this.v[0]*this.tau;
        pLoc.y += this.v[1]*this.tau;
        
        // Now create new geometry at the location of the new centroid
        Point pt = GISFunctions.pointGeometryFromCoordinate(pLoc);
		Geometry pGeomNew = pt.buffer(this.rad);
        
        // Move the agent to the new location. This requires transforming the geometry 
        // back to the geometry used by the geography, which is what this function does.
        SpaceBuilder.moveAgentToGeometry(this.geography, pGeomNew, this);
        
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
        
        // Calculate acceleration due to field of vision consideration
        fovA = motiveAcceleration();


        // To Do: Calculate acceleration due to avoiding collisions with other agents and walls.
        contA = totalContactAcceleration();
        
        totA = Vector.sumV(fovA, contA);
        
        return totA;
    }
    
    
    // Calculate the acceleration towards the destination accounting for objects in the field of vision but not collisions
    public double[] motiveAcceleration()  {
    	
    	double[] desiredVelocity = desiredVelocity();
    	
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
    	Geometry thisGeom = SpaceBuilder.getAgentGeometry(geography, this);
        Context<Object> context = ContextUtils.getContext(this);
    	
    	
    	// Iterate over all other pedestrian agents and for those that touch the 
    	// ego agent calculate the interaction force
    	// Check to see if this line intersects with any agents
        for (Object agent :context.getObjects(Ped.class)) {
        	Ped P = (Ped)agent;
        	if (P != this) {
               	Geometry agentG = SpaceBuilder.getAgentGeometry(geography, P);
               	if (agentG.intersects((thisGeom))) {
               		double[] pCA = pedestrianContactAcceleration(this, P, agentG);
               		cATotal = Vector.sumV(cATotal, pCA);
               	}
        	}
        }
        
        for (Object obstr :context.getObjects(PedObstruction.class)) {
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
    	double d_ij = pLoc.distance(agentCoord);
    	
    	// Get the vector that points from centorid of other agent to the ego agent,
    	// this is the direction that the force acts in
    	double[] n = {pLoc.x - agentCoord.x, pLoc.y - agentCoord.y};
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
    	double d_ij = pLoc.distance(intersectionCoord);
    	
    	// Get the vector that points from centroid of other agent to the ego agent,
    	// this is the direction that the force acts in.
    	// This should also be perpendicular to the obstacle 
    	double[] n = {pLoc.x - intersectionCoord.x, pLoc.y - intersectionCoord.y};
    	n = Vector.unitV(n);
    	
    	double magA = this.k * (r_i - d_ij) / this.m;
    	double[] A  = {magA*n[0], magA*n[1]};
    	
    	return A;
    	
    }
    
    // Function to sample field of vision
    public List<Double> sampleFoV() {
    	
    	// Initialise a list to hole the sampled field of vision vectors
    	List<Double> sampledAngles = new ArrayList<Double>();
    	
    	double sampleAngle = this.aP-this.theta; // First angle to sample
    	double sampleAnglemax = this.aP + this.theta;
    	while (sampleAngle <= sampleAnglemax) {
    		sampledAngles.add(sampleAngle);
    		sampleAngle+=this.angres;
    	}
    	
    	return sampledAngles;
    	
    }
    
    // Function to calculate distance to nearest collision for a given angle f(a) -  this will need to account for movements of other peds
    public double distanceToObject(double alpha)  {
    	
    	// Initialise distance to nearest object as the max distance in the field of vision
    	double d = this.dmax;
    	
    	// Get unit vector in the direction of the sampled angle
    	double[] rayVector = {Math.sin(alpha), Math.cos(alpha)};
    	
    	// Get the coordinate of the end of the field of vision in this direction
    	Coordinate rayEnd = new Coordinate(pLoc.x + rayVector[0]*this.dmax, pLoc.y + rayVector[1]*this.dmax);
    	
    	Coordinate[] lineCoords = {pLoc, rayEnd};
    	// Create a line from the pedestrian to the end of the field of vision in this direction
    	LineString sampledRay = new GeometryFactory().createLineString(lineCoords);
    	
    	// Check to see if this line intersects with any pedestrian agents
        Context<Object> context = ContextUtils.getContext(this);
        for (Object agent :context.getObjects(Ped.class)) {
        	Ped P = (Ped)agent;
        	if (P != this) {
               	Geometry agentG = SpaceBuilder.getAgentGeometry(geography, P);
               	if (agentG.intersects(sampledRay)) {
               		// The intersection geometry could be multiple points.
               		// Iterate over them find the distance to the nearest pedestrian
               		Coordinate[] agentIntersectionCoords = agentG.intersection(sampledRay).getCoordinates();
               		for(Coordinate c: agentIntersectionCoords) {
                   		double dAgent = pLoc.distance(c);
                   		if (dAgent < d) {
                   			d = dAgent;
                   		}
               		}
               	}
        	}
        }
        
    	// Check to see if this line intersects with any obstacle agents
        for (Object obstr :context.getObjects(PedObstruction.class)) {
        	PedObstruction Obstr = (PedObstruction)obstr;
           	Geometry obstG = Obstr.getGeom();
           	if (obstG.intersects(sampledRay)) {
           		// The intersection geometry could be multiple points.
           		// Iterate over them and take the smallest distance - this is the distance to the nearest obstacle
           		Coordinate[] obstIntersectionCoords = obstG.intersection(sampledRay).getCoordinates();
           		for(Coordinate c: obstIntersectionCoords) {
           			double dAgent = pLoc.distance(c);
               		if (dAgent < d) {
               			d = dAgent;
               		}
           		}
           	}
        }
        
        return d;    	
    }
    
    // Function to calculate d(a) using cos rule
    public double displacementDistance(double alpha)  {
    	
    	// Get the distance to nearest object for this angle
    	double fAlpha =  distanceToObject(alpha);
    	
    	double dAlpha = Math.pow(this.dmax, 2) + Math.pow(fAlpha, 2) - 2*this.dmax*fAlpha*Math.cos(this.a0 - alpha);
    	
    	return dAlpha;
    }
    
    // Wrapper function that identifies the chosen walking direction
    public Map<String, Double> desiredDirection()  {
    	
    	// Sample field of vision
    	List<Double> sampledAngles = sampleFoV();
    	
    	// Initialise the displacement distance (which must be minimised) and the direction of travel
    	// The angle here is relative to the direction of the agent
    	double d = displacementDistance(sampledAngles.get(0));
    	double alpha = sampledAngles.get(0);    
    	
    	// Loop through the remaining angles and find the angle which minimises the displacement distance
    	for (int i = 1;i<sampledAngles.size(); i++) {
    		
    		double dDist = displacementDistance(sampledAngles.get(i));
    		
    		if (dDist < d) {
    			d = dDist;
    			alpha = sampledAngles.get(i);
    		}
    	}
    	
    	Map<String, Double> output = new HashMap<String, Double>();
    	output.put("angle", alpha);
    	output.put("collision_distance", distanceToObject(alpha));
    	
    	return output;    	
    }
    
    public double[] desiredVelocity()  {
    	
    	// Get the desired direction of travel and minimum distance to collision in that direction
    	Map<String, Double> desiredDirection = desiredDirection();
    	
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
				
				// Change its description so that not recognised as a crossing coord
		    	Coordinate nextCrossingCoord = this.route.getNextRouteCrossingCoord();
				this.route.updateRouteCoordDescription(nextCrossingCoord, GlobalVars.TRANSPORT_PARAMS.routeDefaultDescription);
	    	}
    	}
    }
    
    public void checkIfEnteringCrossing() {
    	Coordinate nextCrossingCoord = this.route.getNextRouteCrossingCoord();
    	if (nextCrossingCoord != null) {
        	Point nextCrossing = GISFunctions.pointGeometryFromCoordinate(nextCrossingCoord);
        	
        	Coordinate lookAhead  = getPedestrianLookAheadCoord(GlobalVars.lookAheadTimeSteps);
        	Coordinate[] lookAheadLineCoords = {this.pLoc, lookAhead};
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
    
    private void updateRoadLinkRouteSection() {
    	
    	Road nextRoad = null;
		try {
			nextRoad = GISFunctions.getCoordinateRoad(this.routeCoord);
		} catch (RoutingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    	// perceive the space taken up by vehicles on the road links that pass by/though this road
    	double vehicleRoadSpace = estimateVehicleRoadSpace(nextRoad);
    	double updatedVehicleGridCellCost = 1 + GlobalVars.MOBILE_AGENT_PARAMS.gridCellCostParam * vehicleRoadSpace;
    	
    	// Using vehicle dominance figure, update pedestrian perception of costs of moving in vehicle priority areas
    	// Use this updated perception of costs when calculating updated Route
    	HashMap<Integer, Double> updatedGridSummandPriorityMap = this.gridSummandPriorityMap;
    	updatedGridSummandPriorityMap.put(GlobalVars.GRID_PARAMS.getPriorityValueMap().get("vehicle"), this.vehiclePriorityCostRatio * updatedVehicleGridCellCost);
    	
    	// Create new Route object, that evaluates flood fill values over a partial section of the grid, with the coordinate where the road changes as the destination
    	Coordinate nextRoadLinkCoord = this.route.getRouteRoadLinkX().get(0);
    	GridRoute partialRoute = new GridRoute(this.geography, this, updatedGridSummandPriorityMap, nextRoadLinkCoord, true);
    	
    	// Get updated set of route coords to follow to next road link coordinate
    	partialRoute.setGroupedGridPath();
    	
    	// Update the next section of the pedestrian's route with the path produced by this partial route (which only goes up to the end of the current road link)
    	updateGroupedGridPath(partialRoute, thisRoadLinkCoord);
    }
    
	// Consider moving to ped, confusing to be in route - a bit meta
	private void updateGroupedGridPath(GridRoute updatedRoute, Coordinate routeSectionCoord) {
		GridCoordinates2D routeSectionCell = this.route.getRouteCoordMap().get(routeSectionCoord);
		
		// Iterate over the route section cells in the partial route to find the one that matches the one we are replacing
		List<GridCoordinates2D> updatedPathSection = null;
		for (GridCoordinates2D cell: updatedRoute.getGroupedGridPath().keySet()) {
			if((cell.x==routeSectionCell.x)&(cell.y==routeSectionCell.y)) {
				updatedPathSection = updatedRoute.getGroupedGridPath().get(cell);
			}
		}
		this.route.getGroupedGridPath().replace(routeSectionCell, updatedPathSection);
	}
	
    
    private double estimateVehicleRoadSpace(Road r) {
    	int vehicleNumber = r.getRoadLinksVehicleCount();
    	int nRoadLinks = r.getRoadLinks().size();
    	double roadLinkLength = r.getRoadLinks().get(0).getGeom().getLength();
    	
    	double roadArea = roadLinkLength * GlobalVars.MOBILE_AGENT_PARAMS.laneWidth * nRoadLinks;
    	double vehicleDensity = vehicleNumber/roadArea;
    	
    	// Proportion of road taken up by physical presence of vehicles
    	double stationaryVehicleSpace = vehicleDensity * GlobalVars.MOBILE_AGENT_PARAMS.vehicleWidth * GlobalVars.MOBILE_AGENT_PARAMS.vehicleLength;
    	
    	// Proportion of road taken up by spacing between vehicles (assuming time separated by driver reaction time)
    	double vehicleSeparationSpace = vehicleDensity * GlobalVars.MOBILE_AGENT_PARAMS.vehicleWidth * GlobalVars.MOBILE_AGENT_PARAMS.vehicleSpeed * GlobalVars.MOBILE_AGENT_PARAMS.vehicleReactionTime;
    	
    	// Proportion of road space taken up by space in front of lead vehicle required for pedestrian to cross
    	double pedestrianGapSpace;
        if (vehicleDensity == 0) {
        	pedestrianGapSpace = 0;
        } 
        else {
        	// Total width is given by number of road links * lane width (assuming each link is only one lane)
        	double tGap = nRoadLinks*GlobalVars.MOBILE_AGENT_PARAMS.laneWidth / this.v0;
        	
        	// A pedestrian gap space is added for each leader vehicle (1 if only one road link has vehicles on it. 2 if two road likes have vehicles on them)
        	pedestrianGapSpace = r.getNumberLeadVehicles() * GlobalVars.MOBILE_AGENT_PARAMS.vehicleWidth * GlobalVars.MOBILE_AGENT_PARAMS.vehicleSpeed * tGap /roadArea;
        }
        
        return stationaryVehicleSpace + vehicleSeparationSpace + pedestrianGapSpace;
    }
    
    public Color getColor() {
    	return this.col;
    }
    
    /*
     * Calculate bearing to destination and convert to a unit vector
     */
    public double setBearingToDestinationCoord(Coordinate dLoc)  {
    	
        double[] dirToEnd = {dLoc.x - pLoc.x, dLoc.y - pLoc.y};        
        dirToEnd = Vector.unitV(dirToEnd);
        
        this.a0 = Vector.angleBetweenNorthAndUnitVector(dirToEnd);
        
        return this.a0;
    }
    
    /*
     * Set the direction of the pedestrian to be the same as the direction of the velocity vector
     */
    public double setPedestrianBearingFromVelocity(double[] v) {
    	
    	// If velocity is 0 then don't update the pedestrian direction
    	if (Vector.mag(v)==0) {
    		return this.aP;
    	}
    	
    	double[] unitV = Vector.unitV(v);
    	
    	this.aP = Vector.angleBetweenNorthAndUnitVector(unitV);
    	
    	return this.aP;
    	
    }
    
    public void setPedestrianBearing(double aP) {
    	this.aP = aP;
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
    	
    	double dx = this.v0*nTimeSteps*GlobalVars.stepToTimeRatio*Math.sin(this.aP);
    	double dy = this.v0*nTimeSteps*GlobalVars.stepToTimeRatio*Math.cos(this.aP);
    	
    	Coordinate newLookAhead = new Coordinate(this.pLoc.x + dx, this.pLoc.y + dy);
    	return newLookAhead;
    }
    
    public void setRoadLinkFID(String rlFID) {
    	this.roadLinkFID = rlFID;
    }
    
    public String getRoadLinkFID() {
    	return this.roadLinkFID;
    }
    
    /**
     * Method to be run when agent is removed from the context.
     * 
     * Currently nothing needs to be tidied for pedestrian agents.
     */
    public void tidyForRemoval() {
    	;
    }
    
    public double getRad() {
    	return this.rad;
    }
    
    public double getSpeed() {
    	return Vector.mag(this.v);
    }
    
    /*
     * Get the coordinate of the agents centroid in the references frame used 
     * by the agent class for GIS calculations. Note that this attribute is updated
     * as part of the algorithm for producing pedestrian movement.
     * 
     * @returns Coordinate. The coordinate of the centroid of the pedestrian agent.
     */
    @Override
    public Coordinate getLoc() {
    	return this.pLoc;
    }
    
    /*
     * Set the location attribute of the agent to be the coordinate of its 
     * centroid, in the coordinate reference frame used by the agent for GIS calculations. 
     */
    @Override
    public void setLoc()  {
    	// Get centroid coordinate of this agent
    	Coordinate pL = SpaceBuilder.getAgentGeometry(geography, this).getCentroid().getCoordinate();
    	this.pLoc = pL;
    }
    
    /*
     * Get the destination of this pedestrian
     * 
     *  @returns
     *  	The Destination object of this pedestrian
     */
    @Override
    public Destination getDestination() {
    	return this.destination;
    }
    
    
    /*
     * Getter for the route
     * 
     * @returns Route of the pedestrian
     * 
     */
    public GridRoute getRoute() {
    	return this.route;
    }
    
    @Override
    public Geography<Object> getGeography() {
    	return this.geography;
    }
    
    
    public void setRouteCoord(Coordinate rC) {
    	this.routeCoord = rC;
    }
}
