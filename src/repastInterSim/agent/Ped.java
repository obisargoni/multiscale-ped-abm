package repastInterSim.agent;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateSequence;
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
import repastInterSim.environment.Route;
import repastInterSim.environment.Vector;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class Ped implements mobileAgent {
    private Geography<Object> geography; // Space the agent exists in
    public Destination destination; // The destination agent that this pedestrian agents is heading towards.
    
    private Route route;
    
    private Random rnd = new Random(); // Random seed used to give a distribution of velocities 
    
    private double A, B, r, k; // Constants related to the interaction between agents and the desired velocity of this agent
    
    // Variables related to the pedestrians vision and movements
    private double theta; // Field of vision extends from -theta to + theta from the normal to the agent (ie agent's direction)
    private double m; // Agent's mass
    private double dmax; // Maximum distance within which object impact pedestrian movement, can be though of as horizon of field of vision
    private double v0; // Desired walking spped of pedestrian agent
    private double a0; // Angle to the destination
    private double aP; // Angle of pedestrian direction
    private double tau; // Time period in which pedestrian agent is able to come to a complete stop. Used to set acceleration required to avoid collisions.
    private double angres; // Angular resolution used when sampling the field of vision
    private double[] v, newV; // Velocity and direction vectors
    private double rad; // Radius of circle representing pedestrian, metres
    private GeometryFactory gF;
    private Coordinate pLoc; // The coordinate of the centroid of the pedestrian agent.
    
    private Color col; // Colour of the pedestrian
    
    
    /*
     * Instance method for the Ped class.
     * 
     * @param space the continuousspace the Ped exists in
     * @param direction the pedestrian's direction
     */
    public Ped(Geography<Object> geography, Geography<Destination> destinationGeography, GeometryFactory gF, Destination d) {
        this.geography = geography;
        this.destination = d;
        this.v0  = rnd.nextGaussian() * GlobalVars.pedVsd + GlobalVars.pedVavg;
        this.m  = rnd.nextGaussian() * GlobalVars.pedMasssd + GlobalVars.pedMassAv;
        this.rad = m / 320; // As per Moussaid
        this.gF = gF;
        
        // Set the pedestrian velocity - half of max velocity in the direction the pedestrian is facing
        double[] v = {0,0};
        this.v =  v;

        this.tau  = 0.5/GlobalVars.tStep;
        this.m     = 80;
        this.dmax = 10/GlobalVars.spaceScale; 
        this.angres = (2*Math.PI) / 36; // Equivalent to 10 degrees
        this.theta = (2*Math.PI*75) / 360; // 75 degrees
        this.k = GlobalVars.interactionForceConstant;
        
        this.A     = 2000*GlobalVars.tStep*GlobalVars.tStep/GlobalVars.spaceScale;
        this.B     = 0.08/GlobalVars.spaceScale;
        this.r     = 0.275/GlobalVars.spaceScale;
        
		// Get the destination coordinate, initialise new route and generate a pedestrian route
		Coordinate dCoord = this.destination.getGeom().getCentroid().getCoordinate(); 
		this.route = new Route(geography, this, dCoord);
		
    }
    
    
    /*
     * Calculate the pedestrian's acceleration and resulting velocity
     * given its location, north and destination.
     */
    @ScheduledMethod(start = 1, interval = 1, priority = 2)
    public void step() throws Exception {
    	
    	// Walk towards the next coordinate along the route
        Coordinate routeCoord = this.route.getRouteXCoordinate(0);
    	walk(routeCoord);
    	
    	// Check if destination reached within 1m of route coordinate, if true remove that coordinate from the route
    	if (this.pLoc.distance(routeCoord) < 1) {
    		this.route.removeRouteXCoordinate(routeCoord);
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
    @Override
    public Route getRoute() {
    	return this.route;
    }
}
