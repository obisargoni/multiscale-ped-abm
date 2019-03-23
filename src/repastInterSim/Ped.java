package repastInterSim;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.operation.TransformException;

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

public class Ped {
    private Geography<Object> geography; // Space the agent exists in
    public Destination destination; // The destination agent that this pedestrian agents is heading towards.
    
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
    private double angres; // Angular resolution used when ampling the field of vision
    private double[] v, newV; // Velocity and direction vectors
    private double rad; // Radius of circle representing pedestrian, metres
    private GeometryFactory gF;
    
    private Color col; // Colour of the pedestrian
    
    
    /*
     * Instance method for the Ped class.
     * 
     * @param space the continuousspace the Ped exists in
     * @param direction the pedestrian's direction
     */
    public Ped(Geography<Object> geography, GeometryFactory gF, double aP, Destination d, Color col) {
        this.geography = geography;
        this.destination = d;
        this.v0  = rnd.nextGaussian() * UserPanel.pedVsd + UserPanel.pedVavg;
        this.m  = rnd.nextGaussian() * UserPanel.pedMasssd + UserPanel.pedMassAv;
        this.rad = m / 320; // As per Moussaid
        this.aP = aP;
        this.gF = gF;
        
        this.col = col;
        
        // Set the pedestrian velocity - half of max velocity in the direction the pedestrian is facing
        double[] v = {0,0};
        this.v =  v;

        this.tau  = 0.5/UserPanel.tStep;
        this.m     = 80;
        this.dmax = 10/SpaceBuilder.spaceScale; 
        this.angres = (2*Math.PI) / 36; // Equivalent to 10 degrees
        this.theta = (2*Math.PI*75) / 360; // 75 degrees
        
        this.A     = 2000*UserPanel.tStep*UserPanel.tStep/SpaceBuilder.spaceScale;
        this.B     = 0.08/SpaceBuilder.spaceScale;
        this.r     = 0.275/SpaceBuilder.spaceScale;

    }
    
    
    /*
     * Calculate the pedestrian's acceleration and resulting velocity
     * given its location, north and destination.
     */
    @ScheduledMethod(start = 1, interval = 1, priority = 2)
    public void walk() throws MismatchedDimensionException, NoSuchAuthorityCodeException, FactoryException, TransformException {
    	
    	// Get the geometry of the pedestrian agent and its destination in the CRS used for spatial calculations
    	Geometry pGeom = SpaceBuilder.getGeometryForCalculation(this.geography, this);
     	Coordinate pLoc = pGeom.getCentroid().getCoordinate();
    	
        double[] a = accel();
        double[] dv = {a[0]*this.tau, a[1]*this.tau};
        this.newV  = Vector.sumV(v,dv);
        this.v = newV;

        pLoc.x += this.v[0]*this.tau;
        pLoc.y += this.v[1]*this.tau;
        
        // Now create new geometry at the location of the new centroid
        Coordinate[] pLocArray = {pLoc};
        CoordinateSequence cs = gF.getCoordinateSequenceFactory().create(pLocArray);
        Point pt = new Point(cs, this.gF);
		Geometry pGeomNew = pt.buffer(this.rad);
        
        // Move the agent to the new location. This requires transforming the geometry 
        // back to the geometry used by the geography, which is what this function does.
        SpaceBuilder.moveAgentToCalculationGeometry(this.geography, pGeomNew, this);
        
     	Coordinate pLoc_ = pGeomNew.getCentroid().getCoordinate();
        
        seta0FromDestinationCoord();
        setPedAngleFromVelocity(this.v);
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
    public double[] accel() throws MismatchedDimensionException, TransformException {
        
        double[] fovA;
        
        // Calculate acceleration due to field of vision consideration
        fovA = motiveAcceleration();


        // To Do: Calculate acceleration due to avoiding collisions with other agents and walls.
        return fovA;
    }
    
    
    // Calculate the acceleration towards the destination accounting for objects in the field of vision but not collisions
    public double[] motiveAcceleration() throws MismatchedDimensionException, TransformException {
    	
    	double[] desiredVelocity = desiredVelocity();
    	
    	// Acceleration is set as the acceleration required to reach the desired velocity within tau amount of time
    	double[] a = {0,0};
    	a[0] = (desiredVelocity[0] - this.v[0]) / this.tau;
    	a[1] = (desiredVelocity[1] - this.v[1]) / this.tau;
    	
    	return a;
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
    public double distanceToObject(double alpha) throws MismatchedDimensionException, TransformException {
    	
    	// Initialise distance to nearest object as the max distance in the field of vision
    	double d = this.dmax;
    	
    	// Get coordinate of this agent
    	Coordinate pLoc = SpaceBuilder.getGeometryForCalculation(geography, this).getCentroid().getCoordinate();
    	
    	// Get unit vector in the direction of the sampled angle
    	double[] rayVector = {Math.sin(alpha), Math.cos(alpha)};
    	
    	// Get the coordinate of the end of the field of vision in this direction
    	Coordinate rayEnd = new Coordinate(pLoc.x + rayVector[0]*this.dmax, pLoc.y + rayVector[1]*this.dmax);
    	
    	Coordinate[] lineCoords = {pLoc, rayEnd};
    	// Create a line from the pedestrian to the end of the field of vision in this direction
    	LineString sampledRay = new GeometryFactory().createLineString(lineCoords);
    	
    	// Check to see if this line intersects with any agents
        Context<Object> context = ContextUtils.getContext(this);
        for (Object agent :context.getObjects(Ped.class)) {
        	Ped P = (Ped)agent;
        	if (P != this) {
               	Geometry agentG = SpaceBuilder.getGeometryForCalculation(geography, P);
               	if (agentG.intersects(sampledRay)) {
               		double dAgent = agentG.distance(sampledRay);
               		
               		if (dAgent < d) {
               			d = dAgent;
               		}
               	}
        	}
        }
        
        return d;    	
    }
    
    // Function to calculate d(a) using cos rule
    public double displacementDistance(double alpha) throws MismatchedDimensionException, TransformException {
    	
    	// Get the distance to nearest object for this angle
    	double fAlpha =  distanceToObject(alpha);
    	
    	double dAlpha = Math.pow(this.dmax, 2) + Math.pow(fAlpha, 2) - 2*this.dmax*fAlpha*Math.cos(this.a0 - alpha);
    	
    	return dAlpha;
    }
    
    // Wrapper function that identifies the chosen walking direction
    public Map<String, Double> desiredDirection() throws MismatchedDimensionException, TransformException {
    	
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
    
    public double[] desiredVelocity() throws MismatchedDimensionException, TransformException {
    	
    	// Get the desired direction of travel and minimum distance to collision in that direction
    	Map<String, Double> desiredDirection = desiredDirection();
    	
    	// Calculate the desired speed, minimum between desired speed and speed required to avoid colliding
    	double desiredSpeed = Math.min(this.v0, desiredDirection.get("collision_distance") / this.tau);
    	
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
    
    public double seta0FromDestinationCoord() throws MismatchedDimensionException, TransformException {
        // Calculate bearing to destination and convert to a unit vector
        Coordinate dLoc = SpaceBuilder.getGeometryForCalculation(geography, destination).getCoordinate();
        Coordinate pLoc = SpaceBuilder.getGeometryForCalculation(geography, this).getCentroid().getCoordinate();

        double[] dirToEnd = {dLoc.x - pLoc.x, dLoc.y - pLoc.y};        
        dirToEnd = Vector.unitV(dirToEnd);
        
        this.a0 = Vector.angleBetweenNorthAndUnitVector(dirToEnd);
        
        return this.a0;
    }
    
    /*
     * Set the direction of the pedestrian to be the same as the direction of the velocity vector
     */
    public double setPedAngleFromVelocity(double[] v) {
    	
    	// If velocity is 0 then don't update the pedestrian direction
    	if (Vector.mag(v)==0) {
    		return this.aP;
    	}
    	
    	double[] unitV = Vector.unitV(v);
    	
    	this.aP = Vector.angleBetweenNorthAndUnitVector(unitV);
    	
    	return this.aP;
    	
    }
    
    public void setaP(double aP) {
    	this.aP = aP;
    }
    
    public double getRad() {
    	return this.rad;
    }
    
    public double getSpeed() {
    	return Vector.mag(this.v);
    }
    
}
