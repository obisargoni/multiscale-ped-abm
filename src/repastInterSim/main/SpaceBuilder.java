package repastInterSim.main;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.opengis.geometry.MismatchedDimensionException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

import cern.jet.random.Normal;
import cern.jet.random.Uniform;
import cern.jet.random.engine.RandomEngine;
import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.environment.RunState;
import repast.simphony.engine.schedule.ISchedulableAction;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.gis.util.GeometryUtil;
import repast.simphony.parameter.Parameters;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.space.graph.Network;
import repast.simphony.space.projection.Projection;
import repast.simphony.util.collections.IndexedIterable;
import repastInterSim.agent.MobileAgent;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Route;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.OD;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.FixedGeography;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.VehicleDestinationContext;
import repastInterSim.pathfinding.RoadNetworkRoute;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.CAContext;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.PedObstructionContext;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadLinkContext;

public class SpaceBuilder extends DefaultContext<Object> implements ContextBuilder<Object> {
	
	private static Logger LOGGER = Logger.getLogger(SpaceBuilder.class.getName());
	
	private Boolean isDirected = true; // Indicates whether the vehicle road network is directed ot not. s
	private ArrayList<Geography> fixedGeographies;
	
	// Lookups between or and itn road links
	public static HashMap<RoadLink, List<RoadLink>> orToITN = new HashMap<RoadLink, List<RoadLink>>();
	public static HashMap<RoadLink, RoadLink> itnToOR = new HashMap<RoadLink, RoadLink>();
	public static HashMap<Junction, List<Junction>> orJuncToPaveJunc = new HashMap<Junction, List<Junction>>();
	
	private GeometryFactory fac = new GeometryFactory();
	
	private ISchedulableAction addVehicleAction;
	private ISchedulableAction addPedAction;
	private ISchedulableAction stopAddingPedsAction;
	private ISchedulableAction removeMAgentAction;
	private int nPedsCreated;
	
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
	public Context<Object> build(Context<Object> context) {
		
		//RepastInterSimLogging.init();
		Parameters params = RunEnvironment.getInstance ().getParameters();
		
		// Clear caches before starting
		RoadNetworkRoute.clearCaches();
		Route.clearCaches();
		SpatialIndexManager.clearCaches();
		
		fixedGeographies = new ArrayList<Geography>();
		SpaceBuilder.orToITN = new HashMap<RoadLink, List<RoadLink>>();
		SpaceBuilder.itnToOR = new HashMap<RoadLink, RoadLink>();
		SpaceBuilder.orJuncToPaveJunc = new HashMap<Junction, List<Junction>>();
	    
		context.setId(GlobalVars.CONTEXT_NAMES.MAIN_CONTEXT);	
		this.nPedsCreated=0;
		Ped.resetID();
		Vehicle.resetID();
		
		// Correct way to register multiple random number streams
		
		RandomEngine engPedOD = RandomHelper.registerGenerator("pedODThresholds", params.getInteger("pedODSeed"));
		Uniform pedODUniform = new Uniform(0, 1, engPedOD);
		RandomHelper.registerDistribution("pedODThresholds", pedODUniform);
		
		RandomEngine engVehOD = RandomHelper.registerGenerator("vehODThresholds", params.getInteger("vehODSeed"));
		Uniform vehODUniform = new Uniform(0, 1, engVehOD);
		RandomHelper.registerDistribution("vehODThresholds", vehODUniform);
   
		RandomEngine engCASample = RandomHelper.registerGenerator("caSampleDistribution", params.getInteger("caSampleSeed"));
		Uniform caSampleUniform = new Uniform(0, 1, engCASample);
		RandomHelper.registerDistribution("caSampleDistribution", caSampleUniform);
		
		RandomEngine engPedSpeeds = RandomHelper.registerGenerator("pedSpeeds", params.getInteger("pedSpeedSeed"));
		Normal pedSpeedsNorm= new Normal(GlobalVars.pedVavg, GlobalVars.pedVsd, engPedSpeeds);
		RandomHelper.registerDistribution("pedSpeeds", pedSpeedsNorm);
		
		RandomEngine engPedMasses = RandomHelper.registerGenerator("pedMasses", params.getInteger("pedMassSeed"));
		Normal pedMassesNorm= new Normal(GlobalVars.pedMassAv, GlobalVars.pedMasssd, engPedMasses);
		RandomHelper.registerDistribution("pedMasses", pedMassesNorm);
		
		// Read in the model properties
		try {
			IO.readProperties();
		} catch (IOException ex) {
			throw new RuntimeException("Could not read model properties,  reason: " + ex.toString(), ex);
		}
	   
		// Initiate geographic spaces
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		Geography geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, context, geoParams);
		geography.setCRS(GlobalVars.geographyCRSString);
		context.add(geography);
		
		// Road Geography stores polygons representing road and pavement surfaces
		RoadContext roadContext = new RoadContext();
		Geography<Road> roadGeography = createTypedGeography(Road.class, roadContext, GlobalVars.CONTEXT_NAMES.ROAD_GEOGRAPHY);
		context.addSubContext(roadContext);
		fixedGeographies.add(roadGeography);
		
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		PedObstructionContext pedObstructContext = new PedObstructionContext();
		Geography<PedObstruction> pedObstructGeography = createTypedGeography(PedObstruction.class, pedObstructContext, GlobalVars.CONTEXT_NAMES.PED_OBSTRUCTION_GEOGRAPHY);
		context.addSubContext(pedObstructContext);
		fixedGeographies.add(pedObstructGeography);

		// Road link geography is used to create the road network projection
		RoadLinkContext roadLinkContext = new RoadLinkContext();
		Geography<RoadLink> roadLinkGeography = createTypedGeography(RoadLink.class, roadLinkContext, GlobalVars.CONTEXT_NAMES.ROAD_LINK_GEOGRAPHY);
		context.addSubContext(roadLinkContext);
		fixedGeographies.add(roadLinkGeography);
		
		RoadLinkContext orRoadLinkContext = new RoadLinkContext(GlobalVars.CONTEXT_NAMES.OR_ROAD_LINK_CONTEXT);
		Geography<RoadLink> orRoadLinkGeography = createTypedGeography(RoadLink.class, orRoadLinkContext, GlobalVars.CONTEXT_NAMES.OR_ROAD_LINK_GEOGRAPHY);
		context.addSubContext(orRoadLinkContext);
		fixedGeographies.add(orRoadLinkGeography);
		
		RoadLinkContext pavementLinkContext = new RoadLinkContext(GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_CONTEXT);
		Geography<RoadLink> pavementLinkGeography = createTypedGeography(RoadLink.class, pavementLinkContext, GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_GEOGRAPHY);
		context.addSubContext(pavementLinkContext);
		fixedGeographies.add(pavementLinkGeography);

		// Junction geography also used to create the road network
		JunctionContext junctionContext = new JunctionContext();
		Geography<Junction> junctionGeography = createTypedGeography(Junction.class, junctionContext, GlobalVars.CONTEXT_NAMES.JUNCTION_GEOGRAPHY);
		context.addSubContext(junctionContext);
		fixedGeographies.add(junctionGeography);
		
		JunctionContext orJunctionContext = new JunctionContext(GlobalVars.CONTEXT_NAMES.OR_JUNCTION_CONTEXT);
		Geography<Junction> orJunctionGeography = createTypedGeography(Junction.class, orJunctionContext, GlobalVars.CONTEXT_NAMES.OR_JUNCTION_GEOGRAPHY);
		context.addSubContext(orJunctionContext);
		fixedGeographies.add(orJunctionGeography);
		
		JunctionContext pavementJunctionContext = new JunctionContext(GlobalVars.CONTEXT_NAMES.PAVEMENT_JUNCTION_CONTEXT);
		Geography<Junction> pavementJunctionGeography = createTypedGeography(Junction.class, pavementJunctionContext, GlobalVars.CONTEXT_NAMES.PAVEMENT_JUNCTION_GEOGRAPHY);
		context.addSubContext(pavementJunctionContext);
		fixedGeographies.add(pavementJunctionGeography);
		
		// Destinations geography used for creating cache of destinations and their nearest road coordinates		
		VehicleDestinationContext vehicleDestinationContext = new VehicleDestinationContext();
		Geography<OD> vehicleDestinationGeography = createTypedGeography(OD.class, vehicleDestinationContext, GlobalVars.CONTEXT_NAMES.VEHICLE_DESTINATION_GEOGRAPHY);
		context.addSubContext(vehicleDestinationContext);
		fixedGeographies.add(vehicleDestinationGeography);

		PedestrianDestinationContext pedestrianDestinationContext = new PedestrianDestinationContext();
		Geography<OD> pedestrianDestinationGeography = createTypedGeography(OD.class, pedestrianDestinationContext, GlobalVars.CONTEXT_NAMES.PEDESTRIAN_DESTINATION_GEOGRAPHY);
		context.addSubContext(pedestrianDestinationContext);
		fixedGeographies.add(pedestrianDestinationGeography);
		
		CAContext caContext = new CAContext();
		Geography<CrossingAlternative> caGeography = createTypedGeography(CrossingAlternative.class, caContext, GlobalVars.CONTEXT_NAMES.CA_GEOGRAPHY);
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
			String pedJPath = GISDataDir + IO.getProperty(GlobalVars.PavementJunctionShapeFile);					
			GISFunctions.readShapefile(Junction.class, pedJPath, pavementJunctionGeography, pavementJunctionContext);
			SpatialIndexManager.createIndex(pavementJunctionGeography, Junction.class);
			
			String pavementLinkFile = GISDataDir + IO.getProperty("PavementLinkShapefile");
			GISFunctions.readShapefile(RoadLink.class, pavementLinkFile, pavementLinkGeography, pavementLinkContext);
			SpatialIndexManager.createIndex(pavementLinkGeography, RoadLink.class);
			
			// 2a. vehicle roadNetwork
			NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, isDirected);
			builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
			Network<Junction> roadNetwork = builder.buildNetwork();
			GISFunctions.buildGISRoadNetwork(roadLinkGeography, junctionContext,junctionGeography, roadNetwork);
			
			// 2b. open road road network (use by pedestrian agents)
			NetworkBuilder<Junction> orBuilder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.OR_ROAD_NETWORK,orJunctionContext, false);
			orBuilder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
			Network<Junction> orRoadNetwork = orBuilder.buildNetwork();
			GISFunctions.buildGISRoadNetwork(orRoadLinkGeography, orJunctionContext,orJunctionGeography, orRoadNetwork);
			
			// 2c pavement network
			NetworkBuilder<Junction> pavementBuilder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.PAVEMENT_NETWORK, pavementJunctionContext, false);
			pavementBuilder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
			Network<Junction> pavementNetwork = pavementBuilder.buildNetwork();
			GISFunctions.buildGISRoadNetwork(pavementLinkGeography, pavementJunctionContext, pavementJunctionGeography, pavementNetwork);
			
			// Create lookup from OR road junctions to pavement junctions
			for (Junction orJ : orJunctionGeography.getAllObjects()) {
				if (!SpaceBuilder.orJuncToPaveJunc.containsKey(orJ)) {
					SpaceBuilder.orJuncToPaveJunc.put(orJ, new ArrayList<Junction>());
				}
				for (Junction paveJ : pavementJunctionGeography.getAllObjects()) {
					if(paveJ.getjuncNodeID().contentEquals(orJ.getFID())) {
						SpaceBuilder.orJuncToPaveJunc.get(orJ).add(paveJ);
					}
				}
			}			
			
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
			
			// Create lookups from OR road links to ITN links and visa versa
			for (RoadLink orRL: orRoadLinkGeography.getAllObjects()) {
				List<RoadLink> itnLinks = new ArrayList<RoadLink>();
				
				for(RoadLink itnLink: roadLinkGeography.getAllObjects()) {
					if (itnLink.getPedRLID().contentEquals(orRL.getFID())) {
						SpaceBuilder.itnToOR.put(itnLink, orRL);
						itnLinks.add(itnLink);
					}
				}
				SpaceBuilder.orToITN.put(orRL, itnLinks);
				
				// Also assign Road objects to OR road links, so pedestrians can identify pavement polygons nearby.
				for (Road r: roadGeography.getAllObjects()) {
					if(r.getRoadLinkID().contentEquals(orRL.getFID())) {
						orRL.getRoads().add(r);
					}
				}
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
		
		// Read in OD matrix data for vehicles from CSV
		List<String[]> vehicleFlows = IO.readCSV(GISDataDir + IO.getProperty("vehicleODFlowsFile"));
		
		// Read in OD matrix data for pedestrians from CSV
		List<String[]> pedestrianFlows = IO.readCSV(GISDataDir + IO.getProperty("pedestrianODFlowsFile"));

		// Schedule the creation of vehicle agents - tried doing this with annotations but it didnt work
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
		
		int  addVehicleTicks = params.getInteger("addVehicleTicks");
	    ScheduleParameters vehicleScheduleParams = ScheduleParameters.createRepeating(1, addVehicleTicks, ScheduleParameters.FIRST_PRIORITY);
	    addVehicleAction = schedule.schedule(vehicleScheduleParams, this, "addVehicleAgents", vehicleFlows);
	    
		// Schedule the creation of pedestrian agents
	    int startPedsTick = 3*addVehicleTicks; // Add peds to model after three round of adding vehicles.
		int  addPedTicks = params.getInteger("addPedTicks");
	    ScheduleParameters pedestrianScheduleParams = ScheduleParameters.createRepeating(startPedsTick,addPedTicks,ScheduleParameters.FIRST_PRIORITY);
	    addPedAction = schedule.schedule(pedestrianScheduleParams, this, "addPedestrianAgents", pedestrianFlows);
	    
	    // Schedule method that removes agents
		ScheduleParameters removeMAgentScheduleParameters = ScheduleParameters.createRepeating(1, 1, ScheduleParameters.LAST_PRIORITY);
		removeMAgentAction = schedule.schedule(removeMAgentScheduleParameters, this, "removeAgentsAtDestinations");
	    
	    // Stop adding agents to the simulation at endTick ticks
	    ScheduleParameters stopAddingAgentsScheduleParams = ScheduleParameters.createRepeating(startPedsTick+1, addPedTicks, ScheduleParameters.LAST_PRIORITY);
	    stopAddingPedsAction = schedule.schedule(stopAddingAgentsScheduleParams, this, "stopAddingPedAgents");
	    
	    //IO.exportGridCoverageData(baseGrid);
		return context;
		
	}
	
	/**
	 * Stop adding Ped agents to the simulation
	 */
	public void stopAddingPedAgents() {
		
		// Check number of peds in model
		int nPedsToCreate = RunEnvironment.getInstance().getParameters().getInteger("nPeds");
		
		if (this.nPedsCreated>=nPedsToCreate) {
			ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();

			// Remove add agents methods from the scedule
			boolean success = schedule.removeAction(addPedAction);
			
			// if action not removed, reschedule this method for the following tick
			if (success==false) {
			    ScheduleParameters stopAddingAgentsScheduleParams = ScheduleParameters.createOneTime(schedule.getTickCount()+1, ScheduleParameters.LAST_PRIORITY);
			    schedule.schedule(stopAddingAgentsScheduleParams, this, "stopAddingPedAgents");
			}
			else {
				// If successfully stop adding peds to model, schedule methods that monitors when ti end the simulation
			    ScheduleParameters endRunScheduleParams = ScheduleParameters.createRepeating(schedule.getTickCount()+1,10,ScheduleParameters.LAST_PRIORITY);
			    schedule.schedule(endRunScheduleParams, this, "endSimulation");
			    
			    // Also schedule action that removes the action that triggers this function
			    ScheduleParameters sp = ScheduleParameters.createOneTime(schedule.getTickCount()+1, ScheduleParameters.FIRST_PRIORITY);
			    schedule.schedule(sp, this, "stopStopAddingPedAgents");			    
			}
		}
	}
	
	public void stopStopAddingPedAgents() {
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
		boolean sucess = schedule.removeAction(stopAddingPedsAction);
	}
	
	/**
	 * Stop adding Vehicle agents to the simulation
	 */
	public void stopAddingVehicleAgents() {
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();

		// Remove add agents methods from the scedule
		boolean success = schedule.removeAction(addVehicleAction);
		
		// if action not removed, reschedule this method for the following tick
		if (success==false) {
		    ScheduleParameters stopAddingAgentsScheduleParams = ScheduleParameters.createOneTime(schedule.getTickCount()+1, ScheduleParameters.LAST_PRIORITY);
		    schedule.schedule(stopAddingAgentsScheduleParams, this, "stopAddingVehicleAgents");
		}
		else {
			// If successfully stop adding vehicle agents, can schedule method remove all vehicle agents from the simualtion, to trigger the end of the run.
		    ScheduleParameters removeAllVehicleAgentsParams = ScheduleParameters.createOneTime(schedule.getTickCount()+1, ScheduleParameters.LAST_PRIORITY);
		    schedule.schedule(removeAllVehicleAgentsParams, this, "removeAllVehicleAgents");
		}
	}
	
	public void removeAllVehicleAgents() {
		Context context = RunState.getInstance().getMasterContext();
	        
        // Iterate over vehicles and remove them
		List<Vehicle> vehiclesToRemove = new ArrayList<Vehicle>();
        for (Object o :context.getObjects(Vehicle.class)) {
        	Vehicle v  = (Vehicle) o;
        	vehiclesToRemove.add(v);
        }
		
		for (Vehicle v: vehiclesToRemove) {
			removeMobileAgent(v, null);
		}

	}
	
	
	/**
	 * End the simulation if there are no mobile agents in the context.
	 */
	public void endSimulation() {
		
		// First stop adding vehicle agents if no more pedestrian agents
		Context context = RunState.getInstance().getMasterContext();
		if (context.getObjects(Ped.class).size()==0) {
			stopAddingVehicleAgents();
		}
		
		// Then if all agents have completed trips end simulation
		if (context.getObjects(MobileAgent.class).size() == 0) {
			runCleanUP();
			RunEnvironment.getInstance().endRun();
		}
	}
	
	/*
	 * Method to remove agents from the simulation. Requires some care because of links created between agents.
	 */
	private void runCleanUP() {
		for (Geography g: this.fixedGeographies) {
			for (Object o : g.getAllObjects()) {
				FixedGeography fg = (FixedGeography)o;
				fg.clear();
			}
		}
		
		RoadNetworkRoute.clearCaches();
		Route.clearCaches();
		SpatialIndexManager.clearCaches();
		
		fixedGeographies = new ArrayList<Geography>();
		SpaceBuilder.orToITN = new HashMap<RoadLink, List<RoadLink>>();
		SpaceBuilder.itnToOR = new HashMap<RoadLink, RoadLink>();
		SpaceBuilder.orJuncToPaveJunc = new HashMap<Junction, List<Junction>>();
		
		Context mc = RunState.getInstance().getMasterContext();
		for (Object o: mc.getSubContexts()) {
			Context sc  = (Context) o;
			sc.clear();
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
		Geography<OD> vehicleDestinationGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.VEHICLE_DESTINATION_GEOGRAPHY);
		List<OD[]> ods = mobileAgentODs(vehicleDestinationGeography, odData, "vehODThresholds");
		
		for (int i=0; i< ods.size(); i++) {
			OD[] od = ods.get(i);
			addVehicle(od[0], od[1]);
		}
	}
	
    /*
     * Scheduled method that adds pedestrian agents to the simulation. Origin and destination of the pedestrian
     * is set probabilistically based on OD matrix flow values.
     * 
     * @param odData The OD flow data used to create pedestrian agents with origins and destinations that match the flow data
     */
	public void addPedestrianAgents(List<String[]> odData) {
		
		Geography<OD> pedestrianDestinationGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.PEDESTRIAN_DESTINATION_GEOGRAPHY);
		List<OD[]> ods = mobileAgentODs(pedestrianDestinationGeography, odData, "pedODThresholds");
		
		int nPedsToCreate = RunEnvironment.getInstance().getParameters().getInteger("nPeds");
		
		for (int i=0; i< ods.size(); i++) {
			OD[] od = ods.get(i);
			
			if (this.nPedsCreated<nPedsToCreate) {
				addPed(od[0], od[1]);
				this.nPedsCreated++;
			}
		}
	}
	
	/*
	 * Produce list of origin and destination pairs that are used to initialise a group of mobile agents.
	 */
	private List<OD[]> mobileAgentODs(Geography<OD> odGeography, List<String[]> odFlows, String randomDistribution) {
		
		List<OD[]> ods = new ArrayList<OD[]>();
		
		List<Integer> originsThisStep = new ArrayList<Integer>();
		
		// Iterate through all OD pairs and initialise ped moving between these two if prob is above threshold
		// And ped hasn't been created at this origin already
		int nOD = odGeography.size();
		for (int iD = 0; iD<nOD; iD++) {			
			for (int iO = 0; iO<nOD; iO++) {
				// If ped already going to be created at this origin skip it, cant have two peds created at same location
				boolean inList = originsThisStep.contains(iO);
				if (inList) {
					continue;
				}
			
				// First row of data is the IDs of the ODs
				String[] ids = odFlows.get(0);
				String idO = ids[iO];
				String idD = ids[iD];
				
				// Get the OD matrix entry
				Float flow = Float.parseFloat(odFlows.get(iO+1)[iD]);
				double threshold = RandomHelper.getDistribution(randomDistribution).nextDouble();
				
				// According to flow rate, record this od pair as a pair to create mobile agent moving between.
				if (flow > threshold) {
					originsThisStep.add(iO);
					OD o = null;
					OD d = null;
					for (OD j: odGeography.getAllObjects()) {
						if (j.getFID().contentEquals(idO)) {
							o = j;
						}
						else if (j.getFID().contentEquals(idD)) {
							d = j;
						}
					}
					OD[] od = {o,d};
					ods.add(od);
				}
			}
		}
		
		return ods;
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
		Context context = RunState.getInstance().getMasterContext();
		OD d = new OD();
		context.add(d);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = GeometryUtil.generateRandomPointsInPolygon(bndry, 1).get(0);
		
		Geometry destGeom = fac.createPoint(destCoord);
		
		Geography geography = (Geography) context.getProjection(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY);
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
		Context context = RunState.getInstance().getMasterContext();

		OD d = new OD();
		context.add(d);
		
		// Get the x&y coords for the destination set by the user
		Parameters  params = RunEnvironment.getInstance().getParameters();
		double xCoord = (double)params.getInteger(paramX);
		double yCoord = (double)params.getInteger(paramY);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = new Coordinate(xCoord, yCoord);
		
		Geometry destGeom = fac.createPoint(destCoord);
		
		Geography geography = (Geography) context.getProjection(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY);
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
		Context context = RunState.getInstance().getMasterContext();
		Geography geography = (Geography) context.getProjection(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY);
        
        // Instantiate a new pedestrian agent and add the agent to the context
		Parameters  params = RunEnvironment.getInstance().getParameters();
		
		// Draw velocity and mass from random distribution
		Double v = GlobalVars.pedVavg + 3*GlobalVars.pedVsd; // Initialises as a value far from mean
		while ( (v < GlobalVars.pedVavg - 2*GlobalVars.pedVsd) | (v > GlobalVars.pedVavg + 2*GlobalVars.pedVsd) ){ // Exclude extreme values
			v = RandomHelper.getDistribution("pedSpeeds").nextDouble();
		}
		
		
		Double m = GlobalVars.pedMassAv + 3*GlobalVars.pedMasssd; // Initialises as a value far from mean
		while ( (m < GlobalVars.pedMassAv - 2*GlobalVars.pedMasssd) | (m > GlobalVars.pedMassAv + 2*GlobalVars.pedMasssd) ){ // Exclude extreme values
			m = RandomHelper.getDistribution("pedMasses").nextDouble();
		}
		
		Network<Junction> pavementNetwork = SpaceBuilder.getNetwork(GlobalVars.CONTEXT_NAMES.PAVEMENT_NETWORK);
		Geography<Junction> pavementJunctionGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.PAVEMENT_JUNCTION_GEOGRAPHY);
    	Ped newPed = new Ped(o, d, v, m, params.getDouble("alpha"), params.getDouble("lambda"), params.getDouble("gamma"), params.getDouble("epsilon"), params.getInteger("timeThreshold"), params.getBoolean("minCrossing"), params.getDouble("tacticalPlanHorizon"), pavementJunctionGeography, pavementNetwork);
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
		Context context = RunState.getInstance().getMasterContext();
		Geography geography = (Geography) context.getProjection(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY);
		
		Vehicle V = new Vehicle(GlobalVars.maxVehicleSpeed, GlobalVars.defaultVehicleAcceleration, GlobalVars.initialVehicleSpeed, o, d);
		context.add(V);
		Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = fac.createPoint(oCoord);
		Geometry vehicleCircle = pt.buffer(2);
		GISFunctions.moveAgentToGeometry(geography, vehicleCircle, V);
		V.setLoc();
		
		return V;
    }
    
	public void removeAgentsAtDestinations() {
		Context context = RunState.getInstance().getMasterContext();
		Geography geography = (Geography) context.getProjection(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY);
		
        ArrayList<MobileAgent> AgentsToRemove = new ArrayList<MobileAgent>();
        
        // Iterate over peds and remove them if they have arrived at the destination
        for (Object o :context.getObjects(MobileAgent.class)) {
        	MobileAgent mA  = (MobileAgent) o;
        	
        	// Get the geometries in the CRS used for spatial calculations
        	Geometry dGeom =  mA.getDestination().getGeom();
        	Geometry mAGeom = GISFunctions.getAgentGeometry(geography, mA);
        	
        	// If the pedestrian agent in within the bounds of the destination then remove it from the context as it has reached its destination
        	if (dGeom.isWithinDistance(mAGeom, GlobalVars.MOBILE_AGENT_PARAMS.destinationArrivalDistance)) {
        		AgentsToRemove.add(mA);
        		continue;
        	}
        	
        	// Also check whether the mobile agent has reached its default destination
        	if (mA.getDestination().getGeom().getCoordinate().distance(mAGeom.getCoordinate()) < GlobalVars.MOBILE_AGENT_PARAMS.destinationArrivalDistance) {
        		AgentsToRemove.add(mA);
        		continue;
        	}
        }
        // Now iterate over all of the peds to remove and remove them from the context
        // Need to do this separately from iterating over the peds in the context since altering the context whilst iterating over it throws and exception
        for (MobileAgent mA : AgentsToRemove) {
        	removeMobileAgent(mA, null);
        }
    }
	
	public static void removeMobileAgent(MobileAgent mA, String msg) {
		Context context = RunState.getInstance().getMasterContext();
		
		if (msg!=null) {
    		LOGGER.log(Level.FINE, msg);
		}
    	mA.tidyForRemoval();
    	context.remove(mA);		
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
	
	private static <T> Projection<T> getSpaceBuilderProjection(String name){
		Context c = RunState.getInstance().getMasterContext();
		for (Object o : c.getSubContexts()) {
			Context sc = (Context) o;
			Collection<Projection> ps = sc.getProjections();
			for (Projection p: ps) {
				if (p.getName().contentEquals(name)) {
					return p;
				}
			}
		}
		return null;
	}
	
	public static <T> Geography<T> getGeography(String geogName) {
		Projection<T> p = getSpaceBuilderProjection(geogName);
		Geography<T> g = (Geography<T>) p;
		return g;
	}
	
	public static <T> Network<T> getNetwork(String netName) {
		Projection<T> p = getSpaceBuilderProjection(netName);
		Network<T> n = (Network<T>) p;
		return n;
	}
	
	public static Context getContext(String contextName) {
		Context mc = RunState.getInstance().getMasterContext();
		for (Object o : mc.getSubContexts()) {
			Context sc = (Context) o;
			String name = (String) sc.getId();
			if (name.contentEquals(contextName)) {
				return sc;
			}
		}
		return null;
	}
}
