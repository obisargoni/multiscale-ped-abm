package repastInterSim;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.context.Context;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.Geography;
import repast.simphony.util.ContextUtils;

public class Ped {
    private Geography<Object> geography;
    private List<Double> forcesX, forcesY; 
    private Random rnd = new Random();
    private int age;
    public Destination destination; // Distance from extent to exit simulation
    private double endPtDist, endPtTheta, critGap;

    private double wS, etaS, wV, etaV, sigR;    //errors
    private double m, horiz, A, B, k, r;        //interactive force constants (accT is also)

    private double[] v, dv, newV, dir; 
    private double xTime, accT, maxV;
    private Color col; // Colour of the pedestrian
    
    private MathTransform transformtoMetre;
	private MathTransform transformtoDegree;
    
    /*
     * Instance method for the Ped class.
     * 
     * @param space the continuousspace the Ped exists in
     * @param direction the pedestrian's direction
     */
    public Ped(Geography<Object> geography, double[] direction, Destination d, Color col) {
        this.geography = geography;
        this.destination = d;
        this.maxV  = rnd.nextGaussian() * UserPanel.pedVsd + UserPanel.pedVavg;
        this.dir   = direction; // a 2D unit vector indicating the position
        this.col = col;
        
        // Set the pedestrian velocity - half of max velocity in the direction the pedestrian is facing
        this.v = new double[this.dir.length];
        for (int i = 0; i < dir.length;i++) {
        	this.v[i] = dir[i]*0.5*this.maxV;
        }
        
        this.age   = 0; // Placeholder 

        //3-circle variables - from Helbing, et al (2000) [r from Rouphail et al 1998]
        this.accT  = 0.5/UserPanel.tStep;                        //acceleration time, also termed 'relaxation time'. The time over which the pedestrian regains its desired velocity
        this.m     = 80;                                         //avg ped mass in kg
        this.horiz = 5/SpaceBuilder.spaceScale;                   //distance at which peds affect each other
        this.A     = 2000*UserPanel.tStep*UserPanel.tStep/SpaceBuilder.spaceScale;    //ped interaction constant (kg*space units/time units^2)
        this.B     = 0.08/SpaceBuilder.spaceScale;                    //ped distance interaction constant (space units)
        this.k     = 120000*UserPanel.tStep*UserPanel.tStep;         //wall force constant, no currently used
        this.r     = 0.275/SpaceBuilder.spaceScale; //ped radius (space units)

    }
    
    
    /*
     * Calculate the pedestrian's acceleration and resulting velocity
     * given its location, north and destination.
     */
    @ScheduledMethod(start = 1, interval = 1, priority = 2)
    public void walk() throws MismatchedDimensionException, NoSuchAuthorityCodeException, FactoryException, TransformException {
    	
    	// Get the geometry of the pedestrian agent and its destination in the CRS used for spatial calculations
    	Geometry pGeom = SpaceBuilder.getGeometryForCalculation(this.geography, this);
     	Coordinate pLoc = pGeom.getCoordinate();
    	Geometry dGeom = SpaceBuilder.getGeometryForCalculation(this.geography, this.destination);
    	Coordinate dLoc = dGeom.getCoordinate();
    	
        this.dv    = accel(pLoc,dLoc);
        this.newV  = Vector.sumV(v,dv);
        this.newV  = limitV(newV);
        
    	// Update the direction and velocity of pedestrian
    	this.dir = Vector.unitV(newV);
        this.v = newV;
                
        pLoc.x += this.v[0];
        pLoc.y += this.v[1];
        
        // Move the agent to the new location. This requires transforming the geometry 
        // back to the geometry used by the geography, which is what this function does.
        SpaceBuilder.moveAgentToCalculationGeometry(this.geography, pGeom, this);
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
    public double[] accel(Coordinate pedLocation, Coordinate destLocation) throws MismatchedDimensionException, TransformException {
        
        double[] drivingForce, interactionForce, acc;
        
        // Calculate driving force
        drivingForce = sfDrivingForce(pedLocation, destLocation);


        //calculate interactive forces        
        interactionForce = sfTotalInteractiveForce(pedLocation);
        
        acc = new double[] {drivingForce[0] + interactionForce[0], drivingForce[1] + interactionForce[1]};
        return acc;
    }
    
    public double[] sfDrivingForce(Coordinate pLoc, Coordinate dLoc) {
    	
        //calculate heading to endpoint and convert to a unit vector
        double[] dirToEnd = {dLoc.x - pLoc.x, dLoc.y - pLoc.y};        
        dirToEnd = Vector.unitV(dirToEnd);
        
        // Calculate the bearing to the end point
        double dpEnd = Vector.dotProd(SpaceBuilder.north, dirToEnd);
        double endPtTheta = FastMath.acos(dpEnd);
        
        //calculate motive force 
        double motFx = Math.signum(dirToEnd[0])*Math.abs(Math.abs(maxV*Math.sin(endPtTheta)) - v[0])/accT;
        double motFy = Math.signum(dirToEnd[1])*Math.abs(Math.abs(maxV*Math.cos(endPtTheta)) - v[1])/accT;
        
        double[] motF = {motFx, motFy};

        return motF;
    }
    
    public double[] sfSingleInteractiveForce(Coordinate pLoc, Coordinate otherPLoc, double rij) {
    	
    	double[] F = {0,0};
		// Is other pedestrian in front or behind
        // Get displacement to other pedestrian and use to calculate dot product
        // If dot product is between 0-1 then other pedestrian is visible
        double[] dirToPed = {otherPLoc.x - pLoc.x, otherPLoc.y - pLoc.y};
        dirToPed = Vector.unitV(dirToPed);
        
        double dpPed  = Vector.dotProd(this.dir, dirToPed);
        double dpNPed = Vector.dotProd(SpaceBuilder.north, dirToPed);
        
        if (0 <= dpPed & dpPed <= 1) {     //peds only affected by those in front of them
            double absDist = pLoc.distance(otherPLoc);
            if (absDist < horiz) { // peds only affected by others within their horizon
                double signFx  = Math.signum(dirToPed[0]); 
                double signFy  = Math.signum(dirToPed[1]);
                double theta   = FastMath.acos(dpNPed);
                double interFx = signFx*A*Math.exp((rij-absDist)/B)*Math.abs(Math.sin(theta))/m;
                double interFy = signFy*A*Math.exp((rij-absDist)/B)*Math.abs(Math.cos(theta))/m;
                
                F[0] = interFx;
                F[1] = interFy;
                }
            }
        return F;
    }
    
    public double[] sfTotalInteractiveForce(Coordinate pLoc) throws MismatchedDimensionException, TransformException {
        Context<Object> context = ContextUtils.getContext(this);

    	double Fx, Fy;
    	Fx = Fy = 0;
    	
    	for (Object p :context.getObjects(Ped.class)) {
        	Ped P = (Ped) p;
        	if (P != this) {
        		// Get the geometry of the other pedestrian agent in the CRS used for spatial calculations
        		Geometry otherGeom = SpaceBuilder.getGeometryForCalculation(this.geography, P);
                Coordinate otherLoc = otherGeom.getCoordinate();
                double rij     = this.r + P.r;
                
                double[] interF = sfSingleInteractiveForce(pLoc, otherLoc, rij);
                Fx += interF[0];
                Fy += interF[1];
                }
        	}
    	double[] F = {Fx, Fy};
    	
    	return F;
    	
    }
    
    public double[] limitV(double[] input) {
        double totalV = Vector.mag(input);
        
        if (totalV > maxV) {
        	double norm = maxV/totalV;
            input[0] = input[0]*norm;
            input[1] = input[1]*norm;}
        return input;
    }
    
    public Color getColor() {
    	return this.col;
    }
    
}
