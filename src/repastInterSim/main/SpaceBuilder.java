package repastInterSim.main;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import org.opengis.geometry.MismatchedDimensionException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedulableAction;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.gis.util.GeometryUtil;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.space.graph.Network;
import repast.simphony.util.collections.IndexedIterable;
import repastInterSim.agent.MobileAgent;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.OD;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.VehicleDestinationContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.CAContext;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.PedObstructionContext;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadLinkContext;

public class SpaceBuilder extends DefaultContext<Object> implements ContextBuilder<Object> {
	
	private static Boolean isDirected = true; // Indicates whether the vehicle road network is directed ot not. 
	
	private static Integer pDI = 0; // Pedestrian destination index. Used to select which destination to assign to pedestrians
	private static Integer vDI = 0; // Vehicle destination index. Used to select which destination to assign to pedestrians
		
	public static Context<Object> context;
	public static Geography<Object> geography; 
	
	public static Context<Road> roadContext;
	public static Geography<Road> roadGeography;
	
	public static Context<PedObstruction> pedObstructContext;
	public static Geography<PedObstruction> pedObstructGeography;
	
	public static Context<OD> vehicleDestinationContext;
	public static Geography<OD> vehicleDestinationGeography;
	
	public static Context<OD> pedestrianDestinationContext;
	public static Geography<OD> pedestrianDestinationGeography;
	
	public static Context<RoadLink> roadLinkContext;
	public static Geography<RoadLink> roadLinkGeography;
	
	public static Context<Junction> junctionContext;
	public static Geography<Junction> junctionGeography;
	public static Network<Junction> roadNetwork;
	
	public static Context<RoadLink> orRoadLinkContext;
	public static Geography<RoadLink> orRoadLinkGeography;
	
	public static Context<CrossingAlternative> caContext;
	public static Geography<CrossingAlternative> caGeography;
	
	public static Context<Junction> orJunctionContext;
	public static Geography<Junction> orJunctionGeography;
	public static Network<Junction> orRoadNetwork;
	
	public static Context<Junction> pavementJunctionContext;
	public static Geography<Junction> pavementJunctionGeography;
	public static Context<RoadLink> pavementLinkContext;
	public static Geography<RoadLink> pavementLinkGeography;
	public static Network<Junction> pavementNetwork;
	
	
	private static ArrayList<Geography> fixedGeographies = new ArrayList<Geography>();
	
	public static GeometryFactory fac;
	
	private static ISchedulableAction addVehicleAction;
	private static ISchedulableAction addPedAction;
	private static ISchedulableAction removeMAgentAction;
	
	/*
	 * A logger for this class. Note that there is a static block that is used to configure all logging for the model
	 * (at the bottom of this file).
	 */
	//private static Logger LOGGER = Logger.getLogger(SpaceBuilder.class.getName());
	
	    /* (non-Javadoc)
	 * @see repast.simphony.dataLoader.ContextBuilder#build(repast.simphony.context.Context)
	 * 
	 */
	@Override
	public Context<Object> build(Context<Object> c) {
		
		//RepastInterSimLogging.init();
	    
		context = c;
		context.setId(GlobalVars.CONTEXT_NAMES.MAIN_CONTEXT);
		
		fac = new GeometryFactory();
		
		// Read in the model properties
		try {
			IO.readProperties();
		} catch (IOException ex) {
			throw new RuntimeException("Could not read model properties,  reason: " + ex.toString(), ex);
		}
	   
		// Initiate geographic spaces
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, context, geoParams);
		geography.setCRS(GlobalVars.geographyCRSString);
		context.add(geography);
		
		// Road Geography stores polygons representing road and pavement surfaces
		roadContext = new RoadContext();
		roadGeography = createTypedGeography(Road.class, roadContext, GlobalVars.CONTEXT_NAMES.ROAD_GEOGRAPHY);
		context.addSubContext(roadContext);
		fixedGeographies.add(roadGeography);
		
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		pedObstructContext = new PedObstructionContext();
		pedObstructGeography = createTypedGeography(PedObstruction.class, pedObstructContext, GlobalVars.CONTEXT_NAMES.PED_OBSTRUCTION_GEOGRAPHY);
		context.addSubContext(pedObstructContext);
		fixedGeographies.add(pedObstructGeography);

		// Road link geography is used to create the road network projection
		roadLinkContext = new RoadLinkContext();
		roadLinkGeography = createTypedGeography(RoadLink.class, roadLinkContext, GlobalVars.CONTEXT_NAMES.ROAD_LINK_GEOGRAPHY);
		context.addSubContext(roadLinkContext);
		fixedGeographies.add(roadLinkGeography);
		
		orRoadLinkContext = new RoadLinkContext(GlobalVars.CONTEXT_NAMES.OR_ROAD_LINK_CONTEXT);
		orRoadLinkGeography = createTypedGeography(RoadLink.class, orRoadLinkContext, GlobalVars.CONTEXT_NAMES.OR_ROAD_LINK_GEOGRAPHY);
		context.addSubContext(orRoadLinkContext);
		fixedGeographies.add(orRoadLinkGeography);
		
		pavementLinkContext = new RoadLinkContext(GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_CONTEXT);
		pavementLinkGeography = createTypedGeography(RoadLink.class, pavementLinkContext, GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_GEOGRAPHY);
		context.addSubContext(pavementLinkContext);
		fixedGeographies.add(pavementLinkGeography);

		// Junction geography also used to create the road network
		junctionContext = new JunctionContext();
		junctionGeography = createTypedGeography(Junction.class, junctionContext, GlobalVars.CONTEXT_NAMES.JUNCTION_GEOGRAPHY);
		context.addSubContext(junctionContext);
		fixedGeographies.add(junctionGeography);
		
		orJunctionContext = new JunctionContext(GlobalVars.CONTEXT_NAMES.OR_JUNCTION_CONTEXT);
		orJunctionGeography = createTypedGeography(Junction.class, orJunctionContext, GlobalVars.CONTEXT_NAMES.OR_JUNCTION_GEOGRAPHY);
		context.addSubContext(orJunctionContext);
		fixedGeographies.add(orJunctionGeography);
		
		pavementJunctionContext = new JunctionContext(GlobalVars.CONTEXT_NAMES.PAVEMENT_JUNCTION_CONTEXT);
		pavementJunctionGeography = createTypedGeography(Junction.class, pavementJunctionContext, GlobalVars.CONTEXT_NAMES.PAVEMENT_JUNCTION_GEOGRAPHY);
		context.addSubContext(pavementJunctionContext);
		fixedGeographies.add(pavementJunctionGeography);
		
		// Destinations geography used for creating cache of destinations and their nearest road coordinates		
		vehicleDestinationContext = new VehicleDestinationContext();
		vehicleDestinationGeography = createTypedGeography(OD.class, vehicleDestinationContext, GlobalVars.CONTEXT_NAMES.VEHICLE_DESTINATION_GEOGRAPHY);
		context.addSubContext(vehicleDestinationContext);
		fixedGeographies.add(vehicleDestinationGeography);

		pedestrianDestinationContext = new PedestrianDestinationContext();
		pedestrianDestinationGeography = createTypedGeography(OD.class, pedestrianDestinationContext, GlobalVars.CONTEXT_NAMES.PEDESTRIAN_DESTINATION_GEOGRAPHY);
		context.addSubContext(pedestrianDestinationContext);
		fixedGeographies.add(pedestrianDestinationGeography);
		
		caContext = new CAContext();
		caGeography = createTypedGeography(CrossingAlternative.class, caContext, GlobalVars.CONTEXT_NAMES.CA_CONTEXT);
		context.addSubContext(caContext);
		fixedGeographies.add(caGeography);
		
	    // Load agents from shapefiles
		String GISDataDir = IO.getProperty("GISDataDir");
		try {
			
			// Build the road network and pavement networks
			
			// 1a. Load the vehicle road links
			String roadLinkFile = GISDataDir + IO.getProperty("VehicleRoadLinkShapefile");
			GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, roadLinkContext);
			SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
			
			// 1b. Load the open road road links
			String orRoadLinkFile = GISDataDir + IO.getProperty("ORRoadLinkShapefile");
			GISFunctions.readShapefile(RoadLink.class, orRoadLinkFile, orRoadLinkGeography, orRoadLinkContext);
			SpatialIndexManager.createIndex(orRoadLinkGeography, RoadLink.class);
			
			// 1c. Load pavement nodes and links used for the pavement network
			String pedJPath = GISDataDir + IO.getProperty("PavementJunctionsShapefile");					
			GISFunctions.readShapefile(Junction.class, pedJPath, pavementJunctionGeography, pavementJunctionContext);
			SpatialIndexManager.createIndex(pavementJunctionGeography, Junction.class);
			
			String pavementLinkFile = GISDataDir + IO.getProperty("PavementLinkShapefile");
			GISFunctions.readShapefile(RoadLink.class, pavementLinkFile, pavementLinkGeography, pavementLinkContext);
			SpatialIndexManager.createIndex(pavementLinkGeography, RoadLink.class);
			
			// 2a. vehicle roadNetwork
			NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, isDirected);
			builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
			roadNetwork = builder.buildNetwork();
			GISFunctions.buildGISRoadNetwork(roadLinkGeography, junctionContext,junctionGeography, roadNetwork);
			
			// 2b. open road road network (use by pedestrian agents)
			NetworkBuilder<Junction> orBuilder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.OR_ROAD_NETWORK,orJunctionContext, false);
			orBuilder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
			orRoadNetwork = orBuilder.buildNetwork();
			GISFunctions.buildGISRoadNetwork(orRoadLinkGeography, orJunctionContext,orJunctionGeography, orRoadNetwork);
			
			// 2c pavement network
			NetworkBuilder<Junction> pavementBuilder = new NetworkBuilder<Junction>("PAVEMENT_NETWORK", pavementJunctionContext, false);
			pavementBuilder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
			pavementNetwork = pavementBuilder.buildNetwork();
			GISFunctions.buildGISRoadNetwork(pavementLinkGeography, pavementJunctionContext, pavementJunctionGeography, pavementNetwork);
			
			
			// Build the fixed environment
			
			// 1. Load vehicle and pedestrian destinations
			String vehicleDestinationsFile = GISDataDir + IO.getProperty("VehicleDestinationsFile");
			GISFunctions.readShapefile(OD.class, vehicleDestinationsFile, vehicleDestinationGeography, vehicleDestinationContext);
			
			String pedestrianDestinationsFile = GISDataDir + IO.getProperty("PedestrianDestinationsFile");
			GISFunctions.readShapefile(OD.class, pedestrianDestinationsFile, pedestrianDestinationGeography, pedestrianDestinationContext);
			
			// 2. Load roads
			String vehicleRoadFile = GISDataDir + IO.getProperty("VehicleRoadShapefile");
			String pedestrianRoadFile = GISDataDir + IO.getProperty("PedestrianRoadShapefile");
			GISFunctions.readShapefile(Road.class, vehicleRoadFile, roadGeography, roadContext);
			GISFunctions.readShapefile(Road.class, pedestrianRoadFile, roadGeography, roadContext);
			SpatialIndexManager.createIndex(roadGeography, Road.class);
			
			// Link road with itn and OR road links
			// Also assigns the Road objects to the road links. This enables lookups between OR and ITN road links, through the road objects.
			for (Road r: SpaceBuilder.roadGeography.getAllObjects()) {
				List<RoadLink> roadLinks = new ArrayList<RoadLink>();
				for(RoadLink rl: SpaceBuilder.roadLinkGeography.getAllObjects()) {
					// Iterating over the vehicle road links (ITN) but using their corresponding ped road link (open road) id to check whether they belong to this vehicle polygon
					if (rl.getPedRLID().contentEquals(r.getRoadLinkID())) {
						roadLinks.add(rl);
						rl.getRoads().add(r);
					}
				}
				
				RoadLink orLink = null;
				for(RoadLink rl: SpaceBuilder.orRoadLinkGeography.getAllObjects()) {
					// Iterating over the vehicle road links (ITN) but using their corresponding ped road link (open road) id to check whether they belong to this vehicle polygon
					if (rl.getFID().contentEquals(r.getRoadLinkID())) {
						orLink = rl;
						orLink.getRoads().add(r);
						break;
					}
				}
				
				r.setRoadLinks(roadLinks);
				r.setORRoadLink(orLink);
			}
						
			// 3. Load pedestrian obstruction boundaries
			String pedObstructionFile = GISDataDir + IO.getProperty("PedestrianObstructionShapefile");
			GISFunctions.readShapefile(PedObstruction.class, pedObstructionFile, pedObstructGeography, pedObstructContext);
			SpatialIndexManager.createIndex(pedObstructGeography, PedObstruction.class);
			
			// 4. Load crossing alternatives
			String caFile = GISDataDir + IO.getProperty("CAShapefile");
			GISFunctions.readShapefile(CrossingAlternative.class, caFile, caGeography, caContext);
			SpatialIndexManager.createIndex(caGeography, CrossingAlternative.class);
			
		} catch (MalformedURLException | FileNotFoundException | MismatchedDimensionException e1 ) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		// Set the internal context and geography attributes of the destination agents
		IndexedIterable<OD> vehicleDestinations = vehicleDestinationContext.getObjects(OD.class);
		for (OD d : vehicleDestinations) {
			d.setRootContext(context);
			d.setRootGeography(geography);
			d.setDestinationGeography(vehicleDestinationGeography);
		}
		
		IndexedIterable<OD> pedestrianDestinations = pedestrianDestinationContext.getObjects(OD.class);
		for (OD d : pedestrianDestinations) {
			d.setRootContext(context);
			d.setRootGeography(geography);
			d.setDestinationGeography(pedestrianDestinationGeography);
		}		
		
		// Read in OD matrix data for vehicles from CSV
		List<String[]> vehicleFlows = IO.readCSV(GISDataDir + IO.getProperty("vehicleODFlowsFile"));
		
		// Read in OD matrix data for pedestrians from CSV
		List<String[]> pedestrianFlows = IO.readCSV(GISDataDir + IO.getProperty("pedestrianODFlowsFile"));

		// Schedule the creation of vehicle agents - tried doing this with annotations but it didnt work
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
		Parameters params = RunEnvironment.getInstance ().getParameters();
		
		int  addVehicleTicks = params.getInteger("addVehicleTicks");
	    ScheduleParameters vehicleScheduleParams = ScheduleParameters.createRepeating(1, addVehicleTicks, ScheduleParameters.FIRST_PRIORITY);
	    addVehicleAction = schedule.schedule(vehicleScheduleParams, this, "addVehicleAgents", vehicleFlows);
	    
		// Schedule the creation of pedestrian agents
		int  addPedTicks = params.getInteger("addPedTicks");
	    ScheduleParameters pedestrianScheduleParams = ScheduleParameters.createRepeating(1,addPedTicks,ScheduleParameters.FIRST_PRIORITY);
	    addPedAction = schedule.schedule(pedestrianScheduleParams, this, "addPedestrianAgents", pedestrianFlows);
	    
	    // Schedule method that removes agents
		ScheduleParameters removeMAgentScheduleParameters = ScheduleParameters.createRepeating(1, 1, ScheduleParameters.LAST_PRIORITY);
		removeMAgentAction = schedule.schedule(removeMAgentScheduleParameters, this, "removeAgent");
	    
	    // Stop adding agents to the simulation at 1500 ticks
	    int endTick = 1000;
	    ScheduleParameters stopAddingAgentsScheduleParams = ScheduleParameters.createOneTime(endTick, ScheduleParameters.LAST_PRIORITY);
	    schedule.schedule(stopAddingAgentsScheduleParams, this, "stopAddingAgents");
	    
	    ScheduleParameters endRunScheduleParams = ScheduleParameters.createRepeating(endTick,10,ScheduleParameters.LAST_PRIORITY);
	    schedule.schedule(endRunScheduleParams, this, "endSimulation");
	    
	    //IO.exportGridCoverageData(baseGrid);
		return context;
		
	}
	
	/**
	 * Stop adding mobile agents to the context after a certain number of ticks
	 * @param endTick
	 * 		Number of ticks at which to stop adding mobile agents
	 */
	public void stopAddingAgents() {
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();

		// Remove add agents methods from the scedule
		schedule.removeAction(addVehicleAction);
		schedule.removeAction(addPedAction);
	}
	
	/**
	 * End the simulation if there are no mobile agents in the context.
	 */
	public void endSimulation() {
		if (context.getObjects(MobileAgent.class).size() == 0) {
			RunEnvironment.getInstance().endRun();
		}
	}

	/*
	 * Scheduled method that adds vehicle agents to the simulation. Each method call vehicle agents
	 * are initialised with origins and destinations taken from an OD matrix. The OD matrix values 
	 * are used to control the frequency of vehicles initialised with each OD pair, therefore controlling
	 * the flow of vehicles 
	 * 
	 */
	public void addVehicleAgents(List<String[]> odData) {
		
		// Get number of possible origins/destinations
		int nOD = vehicleDestinationGeography.size();
		
		// Use to generate random number and use to determine which OD pairs to use when creating a vehicle agent - check this is valid
		Random rn = new Random();
		
		// Iterate through all OD pairs and initialise vehicle moving between these two if prob is above threshold
		for (int iO = 0; iO<nOD; iO++) {
			int iD = vDI % nOD;
			
			// First row of data is the IDs of the ODs
			String[] ids = odData.get(0);
			String idO = ids[iO];
			String idD = ids[iD];
			
			// Get the OD matrix entry. Add one to row index to skip header of ids
			Float flow = Float.parseFloat(odData.get(iO + 1)[iD]);
			float threshold = rn.nextFloat();

			
			// Create vehicle instance probabilistically according to flow rates
			if (flow > threshold) {
				OD o = null;
				OD d = null;
				for (OD j: vehicleDestinationContext.getObjects(OD.class)) {
					if (j.getFID().contentEquals(idO)) {
						o = j;
					}
					else if (j.getFID().contentEquals(idD)) {
						d = j;
					}
				}
				addVehicle(o, d);
			}
		}
		
		// Increment the vehicle destination index so that a different destination is considered next time.
		vDI++;
	}
	
    /*
     * Scheduled method that adds pedestrian agents to the simulation. Origin and destination of the pedestrian
     * is set probabilistically based on OD matrix flow values.
     * 
     * @param odData The OD flow data used to create pedestrian agents with origins and destinations that match the flow data
     */
	public void addPedestrianAgents(List<String[]> odData) {
		
		// Get number of possible origins/destinations
		int nOD = pedestrianDestinationGeography.size();
		
		// Use to generate random number and use to determine which OD pairs to use when creating a vehicle agent - check this is valid
		Random rn = new Random();
		
		// Iterate through all OD pairs and initialise vehicle moving between these two if prob is above threshold
		for (int iO = 0; iO<nOD; iO++) {
						
			int iD = pDI % nOD;
			
			// Get the OD matrix entry
			Float flow = Float.parseFloat(odData.get(iO)[iD]);
			float threshold = rn.nextFloat();
			
			// Create vehicle instance probabilistically according to flow rates
			if (flow > threshold) {
				OD o = pedestrianDestinationContext.getObjects(OD.class).get(iO);
				OD d = pedestrianDestinationContext.getObjects(OD.class).get(iD);
				addPed(o, d);
			}
		}
		
		// Increment the vehicle destination index so that a different destination is considered next time.
		pDI++; 
	}
	
	/*
	 * Create a destination agent and add the agent to the context. Generate a random coordinate that lies within a 
	 * boundary in the geography and moves the agent to that coordinate.
	 * 
	 * @param context
	 * 			The context to add the destination agent to
	 * @param geography
	 * 			The geography to move the destination agent to
	 * @param gF
	 * 			The geometry factory used to generate the geometry to move the destination agent to
	 * @param bndry
	 * 			The geometry to select a random coordinate from within
	 * @returns d
	 * 			The destination agent that was added to the context and geography
	 */
	public OD addRandomDestination(Geometry bndry) {
		
		OD d = new OD();
		context.add(d);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = GeometryUtil.generateRandomPointsInPolygon(bndry, 1).get(0);
		
		Geometry destGeom = fac.createPoint(destCoord);
		geography.move(d, destGeom);
				
		return d;
	}
	
	/*
	 * Create a destination agent and add the agent to the context. Creates a coordinate from input x and y values
	 * and moves the agent to that coordinate.
	 * 
	 * @param context
	 * 			The context to add the destination agent to
	 * @param geography
	 * 			The geography to move the destination agent to
	 * @param gF
	 * 			The geometry factory used to generate the geometry to move the destination agent to
	 * @param paramX
	 * 			The x coordinate for the destination
	 * @param paramY
	 * 			The y coordinate for the destination
	 * @returns d
	 * 			The destination agent that was added to the context and geography
	 */
	public OD addUserDestination( String paramX, String paramY) {
		
		OD d = new OD();
		context.add(d);
		
		// Get the x&y coords for the destination set by the user
		Parameters  params = RunEnvironment.getInstance().getParameters();
		double xCoord = (double)params.getInteger(paramX);
		double yCoord = (double)params.getInteger(paramY);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = new Coordinate(xCoord, yCoord);
		
		Geometry destGeom = fac.createPoint(destCoord);
		
		geography.move(d, destGeom);
		
		return d;
		
	}
	
	/*
	 * Adds pedestrian agents to the context and geography. Creates a new pedestrian agent and adds the pedestrian to the context.
	 * Creates a circle geometry in the geography and moves the pedestrian agent to this geometry. Assigns the pedestrian agent
	 * its destination and initial direction.
	 * 
	 * @param context
	 * 			The context to add the pedestrian agent to
	 * @param geography
	 * 			The geography to move the pedestrian agent to
	 * @param 
	 *			The geometry factory used to generate the geometry to move the pedestrian agent to
	 * @param coord
	 * 			The coordinate to move the centroid of the pedestrian to in the geography
	 */
    private Ped addPed(OD o, OD d)  {
        
        // Instantiate a new pedestrian agent and add the agent to the context
		Parameters  params = RunEnvironment.getInstance().getParameters();
		Random rn = new Random();
		boolean minimiseCrossing = false;
		if (params.getDouble("minCrossingProp") > rn.nextFloat()) {
			minimiseCrossing = true;
		}
		
    	Ped newPed = new Ped(o, d, params.getDouble("alpha"), params.getDouble("lambda"), params.getDouble("gamma"), params.getDouble("epsilon"), minimiseCrossing, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementNetwork);
        context.add(newPed);
        
        // Create a new point geometry.
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		
		// Transform the coordinate so that the circle can be created using a radius in metres
		Geometry circle = pt.buffer(newPed.getRad());
		
		// Move the pedestrian to this geometry
		GISFunctions.moveAgentToGeometry(geography, circle, newPed);
		
		// Set the location attribute of the pedestrian agent to be its current location. Simplifies future calculations
		newPed.setLoc();
		
        return newPed;
    }

	/*
     * Initialise a vehicle agent and add to to the context and projection
     */
    private Vehicle addVehicle(OD o, OD d) {
		Vehicle V = new Vehicle(GlobalVars.maxVehicleSpeed, GlobalVars.defaultVehicleAcceleration, GlobalVars.initialVehicleSpeed, o, d);
		context.add(V);
		Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = fac.createPoint(oCoord);
		Geometry vehicleCircle = pt.buffer(2);
		GISFunctions.moveAgentToGeometry(geography, vehicleCircle, V);
		V.setLoc();
		
		return V;
    }
    
	public void removeAgent() {
        ArrayList<MobileAgent> AgentsToRemove = new ArrayList<MobileAgent>();
        
        // Iterate over peds and remove them if they have arrived at the destination
        for (Object o :context.getObjects(MobileAgent.class)) {
        	MobileAgent mA  = (MobileAgent) o;
        	
        	// Get the geometries in the CRS used for spatial calculations
        	Geometry dGeom =  mA.getDestination().getGeom();
        	Geometry mAGeom = GISFunctions.getAgentGeometry(SpaceBuilder.geography, mA);
        	
        	// If the pedestrian agent in within the bounds of the destination then remove it from the context as it has reached its destination
        	if (dGeom.isWithinDistance(mAGeom, GlobalVars.MOBILE_AGENT_PARAMS.destinationArrivalDistance)) {
        		AgentsToRemove.add(mA);
        		continue;
        	}
        	
        	// Also check whether the mobile agent has reached its default destination
        	if (mA.getDefaultDestination().distance(mAGeom.getCoordinate()) < GlobalVars.MOBILE_AGENT_PARAMS.destinationArrivalDistance) {
        		AgentsToRemove.add(mA);
        		continue;
        	}
        }
        // Now iterate over all of the peds to remove and remove them from the context
        // Need to do this separately from iterating over the peds in the context since altering the context whilst iterating over it throws and exception
        for (MobileAgent mA : AgentsToRemove) {
        	mA.tidyForRemoval();
        	context.remove(mA);
        }
    }
	
	
	/*
	 * Calculate the size of an iterable
	 * 
	 * @param i
	 * 			The iterable to calculate the size of
	 * 
	 * @returns 
	 * 			The size of the iterable
	 */
	public static int sizeOfIterable(Iterable<?> i) {
		int size = 0;
		Iterator<?> it = i.iterator();
		while (it.hasNext()) {
			size++;
			it.next();
		}
		return size;
	}
	
	
	/**
	 * Initialise a geography projection. Set the geography name and coordinate reference system.
	 * 
	 * @param <T> 
	 * 			The type of object to be contained in the context and geography
	 * @param cl
	 * 			The class of object to be contained in the context and geography
	 * @param context
	 * 			Context to use when creating geography projection. Must be type Context<T>
	 * @param geographyName
	 * 			The name to use for the geography
	 * @return
	 * 			Geography project of type T
	 */
	private <T> Geography<T> createTypedGeography(Class<T> cl, Context<T> context, String geographyName) {
		GeographyParameters<T> GeoParams = new GeographyParameters<T>();
		Geography<T> g = GeographyFactoryFinder.createGeographyFactory(null).createGeography(geographyName, context, GeoParams);
		g.setCRS(GlobalVars.geographyCRSString);
		
		return g;
	}
}
