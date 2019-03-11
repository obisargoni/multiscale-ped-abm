package repastSocialForce;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;

import repast.simphony.context.Context;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;
import repast.simphony.util.ContextUtils;

public class Ped {
    private ContinuousSpace<Object> space;
    private List<Double> forcesX, forcesY; 
    private Random rnd = new Random();
    private int age;
    public int destExtent; // Distance from extent to exit simulation
    private double endPtDist, endPtTheta, critGap;

    private double wS, etaS, wV, etaV, sigR;    //errors
    private double m, horiz, A, B, k, r;        //interactive force constants (accT is also)
    private NdPoint myLoc;
    public NdPoint endPt;
    private double[] v, dv, newV, dir; 
    private double xTime, accT, maxV;
    private Color col; // Colour of the pedestrian
    
    /*
     * Instance method for the Ped class.
     * 
     * @param space the continuousspace the Ped exists in
     * @param direction the pedestrian's direction
     */
    public Ped(ContinuousSpace<Object> space, double[] direction, Destination d, Color col) {
        this.space = space;
        this.endPt = space.getLocation(d);
        this.destExtent = d.getExtent();
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
        this.r     = 0.275/SpaceBuilder.spaceScale;                   //ped radius (space units)
    }
    
    /* 
     * The pedestrian method to be scheduled. This method updates the pedestrian's
     * acceleration and velocity and walks the pedestrians following this update
     */
    //@ScheduledMethod(start = 1, interval = 1)
    public void step() {
    	calc(); // Calculate the pedestrian's new acceleration and velocity
    	walk(); // Walk the pedestrian
    }
    
    /*
     * Calculate the pedestrian's acceleration and resulting velocity
     * given its location, north and destination.
     */
    @ScheduledMethod(start = 1, interval = 1, priority = 3)
    public void calc() {
        this.myLoc = space.getLocation(this);
        this.dv    = accel(myLoc,endPt);
        this.newV  = Vector.sumV(v,dv);
        this.newV  = limitV(newV);
    }

    /*
     * Move the pedestrian to a new location.
     */
    @ScheduledMethod(start = 1, interval = 1, priority = 2)
    public void walk() { 
    	// Update the direction and velocity of pedestrian
    	this.dir = Vector.unitV(newV);
        this.v = newV;
        move(v);
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
    public double[] accel(NdPoint location, NdPoint endPt) {
        forcesX = new ArrayList<Double>();
        forcesY = new ArrayList<Double>();
        double xF, yF;
        double[] acc;
        xF = yF = 0;

        //calculate heading to endpoint and convert to a unit vector
        double[] dirToEnd = this.space.getDisplacement(location, endPt);        
        dirToEnd = Vector.unitV(dirToEnd);
        
        // Calculate the bearing to the end point
        double dpEnd = Vector.dotProd(SpaceBuilder.north, dirToEnd);
        double endPtTheta = FastMath.acos(dpEnd);
        
        //calculate motive force 
        double motFx = (maxV*Math.sin(endPtTheta) - v[0])/accT;
        double motFy = (maxV*Math.cos(endPtTheta) - v[1])/accT;
        forcesX.add(motFx);
        forcesY.add(motFy);

        //calculate interactive forces
        //TODO: write code to make a threshold for interaction instead of the arbitrary horizon
        
        // Iterate over peds and remove them if they have arrive at the destination
        Context<Object> context = ContextUtils.getContext(this);
        for (Object p :context.getObjects(Ped.class)) {
        	Ped P = (Ped) p;
        	if (P != this) {
                NdPoint otherLoc = space.getLocation(P);

        		// Is other pedestrian in front or behind
                // Get displacement to other pedestrian and use to calculate dot product
                // If dot product is between 0-1 then other pedestrian is visible
                double[] dirToPed = this.space.getDisplacement(location, otherLoc);
                dirToPed = Vector.unitV(dirToPed);
                
                double dpPed  = Vector.dotProd(this.dir, dirToPed);
                double dpNPed = Vector.dotProd(SpaceBuilder.north, dirToPed);
                
                if (0 <= dpPed & dpPed <= 1) {     //peds only affected by those in front of them
                    double absDist = space.getDistance(location, otherLoc);
                    if (absDist < horiz) { // peds only affected by others within their horizon
                        double signFx  = Math.signum(dirToPed[0]); 
                        double signFy  = Math.signum(dirToPed[1]);
                        double theta   = FastMath.acos(dpNPed);
                        double rij     = this.r + P.r;
                        double interFx = A*Math.exp((rij-absDist)/B)*Math.sin(theta)/m;
                        double interFy = A*Math.exp((rij-absDist)/B)*Math.cos(theta)/m;
                        forcesX.add(interFx);
                        forcesY.add(interFy);
                        }
                    }
                }
        	}
        	

        //sum all forces
        for (Double b : forcesX) {
            xF += b;}
        for (Double c : forcesY) {
            yF += c;}
        acc = new double[] {xF, yF};
        return acc;
    }

    public void move(double[] displacement) {
        double[] zero = new double[] {0,0};

        if (displacement != zero) { 
            space.moveByDisplacement(this,displacement);
            setLoc(space.getLocation(this)); // Update the local location variable
        } 
    }
    
    public void setLoc(NdPoint newLoc) {
    	this.myLoc = newLoc;
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
