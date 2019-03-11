package repastSocialForce;

import java.util.Random;

import java.awt.Color;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.continuous.ContinuousSpaceFactory;
import repast.simphony.context.space.continuous.ContinuousSpaceFactoryFinder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.SimpleCartesianAdder;
import repast.simphony.space.continuous.WrapAroundBorders;

public class SpaceBuilder extends DefaultContext<Object> implements ContextBuilder<Object> {
	
	static double spaceScale = 1;
	static double[] north = {0,1}; // Defines north, against which bearings are taken
	
	    /* (non-Javadoc)
	 * @see repast.simphony.dataLoader.ContextBuilder#build(repast.simphony.context.Context)
	 * 

	 */
	@Override
	public Context<Object> build(Context<Object> context) {
	    
		context.setId("repastSocialForce");
		int worldL = 500;
		int worldW = 500;
		
		ContinuousSpaceFactory spaceFactory = 
	            ContinuousSpaceFactoryFinder.createContinuousSpaceFactory(null);
	    ContinuousSpace<Object> space = 
	            spaceFactory.createContinuousSpace("space",context, new SimpleCartesianAdder<Object>(),
	                                               new WrapAroundBorders(), worldL, worldW);
	    ISchedule clock = RunEnvironment.getInstance().getCurrentSchedule();
	    context.add(space);
	    context.add(clock);
	    
	    // A separate class is used to handle the creation of pedestrians
	    //Destination d = addRandomDestination(context, space, worldL, worldW, 5);
	    Destination d1 = addUserDestination(context, space, "destX1", "destY1", 5);
	    Destination d2 = addUserDestination(context, space, "destX2", "destY2", 5);
	    Source flowSource1 = new Source(worldL, worldW, d1, 100, Color.RED);
	    Source flowSource2 = new Source(worldL, worldW, d2, 100, Color.BLUE);
	    context.add(flowSource1);
	    context.add(flowSource2);
		return context;
	}
	
	public Destination addRandomDestination(Context<Object> context, ContinuousSpace<Object> space, int worldL, int worldW, int destExtent) {
		// Initialise random coordinates for the destination
		Random randCoord = new Random();
		double xCoord = (double)randCoord.nextInt(worldW);
		double yCoord = (double)randCoord.nextInt(worldL);
		
		Destination d = new Destination(space, destExtent);
		context.add(d);
		space.moveTo(d,  xCoord, yCoord);
		
		return d;
	}
	
	public Destination addUserDestination(Context<Object> context, ContinuousSpace<Object> space, String paramX, String paramY, int destExtent) {
		// Get the x&y coords for the destination set by the user
		Parameters  params = RunEnvironment.getInstance().getParameters();
		double xCoord = (double)params.getInteger(paramX);
		double yCoord = (double)params.getInteger(paramY);
		
		Destination d = new Destination(space, destExtent);
		context.add(d);
		space.moveTo(d,  xCoord, yCoord);
		
		return d;
		
	}

}
