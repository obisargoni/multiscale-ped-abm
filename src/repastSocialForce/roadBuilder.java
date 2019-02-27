package repastSocialForce;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.continuous.ContinuousSpaceFactory;
import repast.simphony.context.space.continuous.ContinuousSpaceFactoryFinder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.SimpleCartesianAdder;
import repast.simphony.space.continuous.StrictBorders;

public class roadBuilder extends DefaultContext<Object> implements ContextBuilder<Object> {
	
	static double spaceScale = 1;
	
	    /* (non-Javadoc)
	 * @see repast.simphony.dataLoader.ContextBuilder#build(repast.simphony.context.Context)
	 * 

	 */
	@Override
	public Context<Object> build(Context<Object> context) {
	    
		context.setId("repastSocialForce");
		int roadL = 500;
		int worldW = 500;
		
		ContinuousSpaceFactory spaceFactory = 
	            ContinuousSpaceFactoryFinder.createContinuousSpaceFactory(null);
	    ContinuousSpace<Object> space = 
	            spaceFactory.createContinuousSpace("space",context, new SimpleCartesianAdder<Object>(),
	                                               new StrictBorders(), roadL, worldW);
	    ISchedule clock = RunEnvironment.getInstance().getCurrentSchedule();
	    
	    // A separate class is used to handle the creation of pedestrians
	    Source flowSource = new Source();
	    
	    context.add(space);
	    context.add(clock);
	    context.add(flowSource);
		return context;
	}

}
