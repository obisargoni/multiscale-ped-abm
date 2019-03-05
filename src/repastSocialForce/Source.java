package repastSocialForce;

import java.util.ArrayList;
import java.util.Random;

import repast.simphony.context.Context;
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
	private int worldL, worldW;
	private Destination d;
	
	/*
	 * Instance method for Source. Source also contains method calls to move peds.
	 */
	public Source(int worldL, int worldW, Destination d) {
		
		this.worldL = worldL;
		this.worldW = worldW;
		this.d = d;

	}
	
	
	// Only want to add a new ped infrequently
    @ScheduledMethod(start = 1, interval = 100, priority = ScheduleParameters.FIRST_PRIORITY)
    public void addPeds() {
    	int pedDirection = 1; // currently treated as being either 1 or -1 (up or down), don't think direction changes either
    	
    	// These should be random
		Random randCoord = new Random();
		int xCoord = randCoord.nextInt(this.worldW);
		int yCoord = randCoord.nextInt(this.worldL);
        Ped addedPed = addPed(pedDirection, xCoord, yCoord, d);
    }

    public Ped addPed(int direction, int x, int y, Destination d) {
        Context<Object> context = ContextUtils.getContext(this);
        ContinuousSpace<Object> space = (ContinuousSpace<Object>) context.getProjection("space");
        Ped newPed = new Ped(space,direction, d);
        context.add(newPed);
        space.moveTo(newPed,x,y);
        newPed.setLoc(space.getLocation(newPed));
        return newPed;
    }
    
    @ScheduledMethod(start = 1, interval = 100, priority = ScheduleParameters.LAST_PRIORITY)
    public void removePeds() {
        Context<Object> context = ContextUtils.getContext(this);
        ContinuousSpace<Object> space = (ContinuousSpace<Object>) context.getProjection("space");
        
        ArrayList<Ped> PedsToRemove = new ArrayList<Ped>();
        
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
