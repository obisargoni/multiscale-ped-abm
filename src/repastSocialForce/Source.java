package repastSocialForce;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;
import repast.simphony.util.ContextUtils;

/*
 * An agent which generates pedestrian agents and adds them to the space
 */
public class Source {
	// static means this variable is shared between all instances of this class
	static int countRemovedPeds = 0;
	private int worldL, worldW, nPeds; // World dimensions and number of pedestrians to add
	private Destination d;
	private Color col; // Colour to use for styling the pedestrians
	
	/*
	 * Instance method for Source. Source also contains method calls to move peds.
	 */
	public Source(int worldL, int worldW, Destination d, int nPeds, Color col) {
		
		this.worldL = worldL;
		this.worldW = worldW;
		this.d = d;
		this.nPeds = nPeds;
		this.col = col;

	}

	
	// Only want to add a new ped infrequently
    @ScheduledMethod(start = 1, priority = ScheduleParameters.FIRST_PRIORITY)
    public void addPeds() {    	
    	for (int i =0;i<nPeds;i++) {
        	// These should be random
    		Random randCoord = new Random();
    		int xCoord = randCoord.nextInt(this.worldW);
    		int yCoord = randCoord.nextInt(this.worldL);
    		
    		// Generate a random initial direction for the pedestrian
    		double randBearing = randCoord.nextFloat() * FastMath.PI * 2;
    		double[] dir = {FastMath.sin(randBearing), FastMath.cos(randBearing)};
    		
            Ped addedPed = addPed(dir, xCoord, yCoord, d);
    	}
    }

    public Ped addPed(double[] direction, int x, int y, Destination d) {
        Context<Object> context = ContextUtils.getContext(this);
        ContinuousSpace<Object> space = (ContinuousSpace<Object>) context.getProjection("space");
        Ped newPed = new Ped(space,direction, d, col);
        context.add(newPed);
        space.moveTo(newPed,x,y);
        newPed.setLoc(space.getLocation(newPed));
        return newPed;
    }
    
    @ScheduledMethod(start = 1, interval = 1, priority = ScheduleParameters.LAST_PRIORITY)
    public void removePeds() {
        Context<Object> context = ContextUtils.getContext(this);
        ContinuousSpace<Object> space = (ContinuousSpace<Object>) context.getProjection("space");
        
        ArrayList<Ped> PedsToRemove = new ArrayList<Ped>();
        
        // If there are no pedestrians to remove end the simulation
        if (context.getObjects(Ped.class).size() == 0) {
        	RunEnvironment.getInstance().endRun();
        }
        
        // Iterate over peds and remove them if they have arrive at the destination
        for (Object p :context.getObjects(Ped.class)) {
        	Ped P = (Ped) p;
        	NdPoint endPt = P.endPt;
        	int destExtent = P.destExtent;
        	NdPoint locP = space.getLocation(P);
        	if (space.getDistance(locP, endPt) < destExtent) {
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
