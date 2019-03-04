package repastSocialForce;

import java.util.ArrayList;
import java.util.Random;

import repast.simphony.context.Context;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;
import repast.simphony.util.ContextUtils;

/*
 * An agent which generates pedestrian agents and adds them to the space
 */
public class Source {
	// static means this variable is shared between all instances of this class
	static ArrayList<Ped> allPeds = new ArrayList<Ped>();
	static NdPoint endPoint;
	private int worldL, worldW;
	
	/*
	 * Instance method for Source. Source also contains method calls to move peds.
	 */
	public Source(int worldL, int worldW) {
		
		this.worldL = worldL;
		this.worldW = worldW;
		
		// Initialise the endPoint for all pedestrians - this needs to be a separate agent so can be visuallised
		Random randCoord = new Random();
		int xCoord = randCoord.nextInt(worldW);
		int yCoord = randCoord.nextInt(worldL);
		NdPoint pt = new NdPoint(xCoord, yCoord);
		Source.endPoint = pt;
	}
	
	
	// Only want to add a new ped infrequently
    @ScheduledMethod(start = 1, interval = 100, priority = 1)
    public void doStuff() {
    	int pedDirection = 1; // currently treated as being either 1 or -1 (up or down), don't think direction changes either
    	
    	// These should be random
		Random randCoord = new Random();
		int xCoord = randCoord.nextInt(this.worldW);
		int yCoord = randCoord.nextInt(this.worldL);
        Ped addedPed = addPed(pedDirection, xCoord, yCoord, Source.endPoint);
        Source.allPeds.add(addedPed);

        for (Ped a : allPeds) {
            a.calc();
        }
        for (Ped b : allPeds) {
            b.walk();
        }
        
    }

    public Ped addPed(int direction, int x, int y, NdPoint endPoint) {
        Context<Object> context = ContextUtils.getContext(this);
        ContinuousSpace<Object> space = (ContinuousSpace<Object>) context.getProjection("space");
        Ped newPed = new Ped(space,direction, endPoint);
        context.add(newPed);
        space.moveTo(newPed,x,y);
        NdPoint myLoc = space.getLocation(newPed);
        newPed.setLoc(myLoc);
        return newPed;
    }
    
}
