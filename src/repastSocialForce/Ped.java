package repastSocialForce;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;

import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;

public class Ped {
    private ContinuousSpace<Object> space;
    private List<Double> forcesX, forcesY; 
    private Random rnd = new Random();
    private int age;
    public int destExtent; // Distance from extent to exit simulation
    private double endPtDist, endPtTheta, critGap;
    //private double side = roadBuilder.sidewalk;
    private double wS, etaS, wV, etaV, sigR;    //errors
    private double m, horiz, A, B, k, r;        //interactive force constants (accT is also)
    private NdPoint myLoc;
    public NdPoint endPt;
    private double[] v, dv, newV;
    private double xTime, accT, maxV;
    private int dir;         // dir = 1 walks up, -1 walks down
    
    /*
     * Instance method for the Ped class.
     * 
     * @param space the continuousspace the Ped exists in
     * @param direction the pedestrian's direction
     */
    public Ped(ContinuousSpace<Object> space, int direction, Destination d) {
        this.space = space;
        this.endPt = space.getLocation(d);
        this.maxV  = rnd.nextGaussian() * UserPanel.pedVsd + UserPanel.pedVavg;
        this.dir   = direction; // 1 moves up, -1 moves down - need to make direction an angle ting
        this.v     = new double[] {0,(double)dir*.5*maxV}; // Pedestrians initialised with some velocity
        this.age   = 0; // Placeholder 

        //3-circle variables - from Helbing, et al (2000) [r from Rouphail et al 1998]
        this.accT  = 0.5/UserPanel.tStep;                        //acceleration time
        this.m     = 80;                                         //avg ped mass in kg
        this.horiz = 5/roadBuilder.spaceScale;                   //distance at which peds affect each other
        this.A     = 2000*UserPanel.tStep*UserPanel.tStep/roadBuilder.spaceScale;    //ped interaction constant (kg*space units/time units^2)
        this.B     = 0.08/roadBuilder.spaceScale;                    //ped distance interaction constant (space units)
        this.k     = 120000*UserPanel.tStep*UserPanel.tStep;         //wall force constant, no currently used
        this.r     = 0.275/roadBuilder.spaceScale;                   //ped radius (space units)
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
     * given its location, direction and destination.
     */
    @ScheduledMethod(start = 1, interval = 1, priority = 2)
    public void calc() {
        this.myLoc = space.getLocation(this);
        this.dv    = accel(myLoc,dir,endPt);
        this.newV  = sumV(v,dv);
        this.newV  = limitV(newV);
    }

    /*
     * Move the pedestrian to a new location.
     */
    @ScheduledMethod(start = 1, interval = 1, priority = 1)
    public void walk() {
        this.v = newV;
        move(v);
    }
    
    public void exit() {
    	if (this.space.getDistanceSq(this.space.getLocation(this), this.endPt) < this.endPtDist) {
    		
    	}
    }
    
    /*
     * Calculate the acceleration of the pedestrian.
     * 
     * @param location ndpoint representing the pedestrian's location
     * @param direct integer giving the pedestrian's direction
     * @param endPt ndpoint representing the pedestrian's destination 
     * 
     * @return a double representing the pedestrian's new acceleration
     */
    public double[] accel(NdPoint location, int direct, NdPoint endPt) {
        forcesX = new ArrayList<Double>();
        forcesY = new ArrayList<Double>();
        double xF, yF;
        double[] acc;
        xF = yF = 0;

        //calculate heading to endpoint
        double endPtDist  = space.getDistance(location, endPt); 
        double endPtDelX  = endPt.getX()-location.getX();
        double endPtTheta = FastMath.asin((double)direct*endPtDelX/endPtDist);
        if (direct == -1) {
            endPtTheta += Math.PI;}

        //calculate motive force - drive towards end point
        double motFx = (maxV*Math.sin(endPtTheta) - v[0])/accT;
        double motFy = (maxV*Math.cos(endPtTheta) - v[1])/accT;
        forcesX.add(motFx);
        forcesY.add(motFy);

        //calculate interactive forces
        //TODO: write code to make a threshold for interaction instead of the arbitrary horizon
        for (Ped a : Source.allPeds) {
            if (a != this) {
                NdPoint otherLoc = space.getLocation(a);

                // Is other pedestrian in front or behind
                // Not sure where dir get's set
                // Why isn't x taken into account?
                double  visible  = Math.signum((double)dir*(otherLoc.getY()-location.getY()));
                
                if (visible == 1) {     //peds only affected by those in front of them
                    double absDist = space.getDistance(location, otherLoc);
                    if (absDist < horiz) {
                        double delX    = location.getX()-otherLoc.getX();
                        double delY    = location.getY()-otherLoc.getY();
                        double delXabs = Math.abs(delX);
                        double signFx  = Math.signum(delX);
                        double signFy  = Math.signum(delY);
                        double theta   = FastMath.asin(delXabs/absDist);
                        double rij     = this.r + a.r;
                        double interFx = signFx*A*Math.exp((rij-absDist)/B)*Math.sin(theta)/m;
                        double interFy = signFy*A*Math.exp((rij-absDist)/B)*Math.cos(theta)/m;
                        forcesX.add(interFx);
                        forcesY.add(interFy);}}}}

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
        double totalV, norm;
        if (this.dir == 1) {
            if (input[1] < 0) {
                input[1] = 0;}}
        else {
            if (input[1] > 0) {
                input[1] = 0;}}
        totalV = Math.sqrt(input[0]*input[0] + input[1]*input[1]);
        if (totalV > maxV) {
            norm = maxV/totalV;
            input[0] = input[0]*norm;
            input[1] = input[1]*norm;}
        return input;
    }
    
    /*
     * Generic function to sum two vectors that are the same size
     * 
     * @param a vector of doubles
     * @param b vector of doubles
     * 
     * @return vector of doubles that is the sum of the two input vectors
     */
    public double[] sumV(double[] a, double[] b) {
    	
    	if (a.length != b.length) {
    		return null;
    	}
    	else {
    		
    	}
        double[] c = new double[a.length];
            for (int i = 0; i < a.length; i++) {
            c[i] = a[i] + b[i];}
        return c;
    }    
    
}
