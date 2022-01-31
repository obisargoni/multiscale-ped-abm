package repastInterSim.tests;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

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
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.environment.RunState;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.OD;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.CAContext;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.PedObstructionContext;
import repastInterSim.environment.contexts.PedObstructionPointsContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.RoadNetworkRoute;

public class EnvironmentSetup {
	
	public static Context<Object> context;
	public static Geography<Object> geography; 
	
	public static Context<Road> roadContext;
	public static Geography<Road> roadGeography;
	
	public static Context<PedObstruction> pedObstructContext;
	public static Geography<PedObstruction> pedObstructGeography;
	
	public static Context<PedObstruction> pedObstructionPointsContext;
	public static Geography<PedObstruction> pedObstructionPointsGeography;
	
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
	
	// Lookups between or and itn road links
	public static HashMap<RoadLink, List<RoadLink>> orToITN = new HashMap<RoadLink, List<RoadLink>>();
	public static HashMap<RoadLink, RoadLink> itnToOR = new HashMap<RoadLink, RoadLink>();
	public static HashMap<Junction, List<Junction>> orJuncToPaveJunc = new HashMap<Junction, List<Junction>>();
	
	static String testGISDir = ".//data//test_gis_data//";
	static String serialisedLookupPath = testGISDir + "road_link_roads_cache.serialised";
	
	public EnvironmentSetup() {
		// TODO Auto-generated constructor stub
	}
	
	static void setUpProperties() throws IOException {
		context = new DefaultContext<Object>();
		IO.readProperties();
		EnvironmentSetup.clearCaches();
		GeometryFactory fac = new GeometryFactory();
		
		setUpObjectGeography();
	}
	
	static void setUpRandomDistributions() {
		// Correct way to register multiple random number streams
		RandomEngine eng = RandomHelper.registerGenerator("maODThresholds", RandomHelper.getSeed()+1);
		Uniform maODUniform = new Uniform(0, 1, eng);
		RandomHelper.registerDistribution("maODThresholds", maODUniform);
		
		RandomEngine engPedMinCross = RandomHelper.registerGenerator("pedMinCrossThresholds", RandomHelper.getSeed()+2);
		Uniform pedMinCrossUniform = new Uniform(0, 1, engPedMinCross);
		RandomHelper.registerDistribution("pedMinCrossThresholds", pedMinCrossUniform);
   
		RandomEngine engCASample = RandomHelper.registerGenerator("caSampleDistribution", RandomHelper.getSeed()+3);
		Uniform caSampleUniform = new Uniform(0, 1, engCASample);
		RandomHelper.registerDistribution("caSampleDistribution", caSampleUniform);
		
		RandomEngine engPedSpeeds = RandomHelper.registerGenerator("pedSpeeds", RandomHelper.getSeed()+4);
		Normal pedSpeedsNorm= new Normal(GlobalVars.pedVavg, GlobalVars.pedVsd, engPedSpeeds);
		RandomHelper.registerDistribution("pedSpeeds", pedSpeedsNorm);
		
		RandomEngine engPedMasses = RandomHelper.registerGenerator("pedMasses", RandomHelper.getSeed()+5);
		Normal pedMassesNorm= new Normal(GlobalVars.pedMassAv, GlobalVars.pedMasssd, engPedMasses);
		RandomHelper.registerDistribution("pedMasses", pedMassesNorm);
	}
	
	static void setUpObjectGeography() {
		RunState.init();
		RunState.getInstance().setMasterContext(context);
		
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		EnvironmentSetup.geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, EnvironmentSetup.context, geoParams);
		EnvironmentSetup.geography.setCRS(GlobalVars.geographyCRSString);
		EnvironmentSetup.context.add(EnvironmentSetup.geography);
	}
	
	static void setUpRoads() throws Exception {
		// Get road geography
		Context<Road> testRoadContext = new RoadContext();
		GeographyParameters<Road> GeoParamsRoad = new GeographyParameters<Road>();
		EnvironmentSetup.roadGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.ROAD_GEOGRAPHY, testRoadContext, GeoParamsRoad);
		EnvironmentSetup.roadGeography.setCRS(GlobalVars.geographyCRSString);
		
		String pedestrianRoadsPath = testGISDir + "topographicAreaPedestrian.shp";
		String vehicleRoadsPath = testGISDir + "topographicAreaVehicle.shp";
		
		// Load vehicle origins and destinations
		try {
			GISFunctions.readShapefile(Road.class, vehicleRoadsPath, EnvironmentSetup.roadGeography, testRoadContext);
			GISFunctions.readShapefile(Road.class, pedestrianRoadsPath, EnvironmentSetup.roadGeography, testRoadContext);
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		SpatialIndexManager.createIndex(EnvironmentSetup.roadGeography, Road.class);
		
		context.addSubContext(testRoadContext);
	}
	
	static Geography<RoadLink> setUpRoadLinks(String roadLinkFile, String contextName, String geogName) throws MalformedURLException, FileNotFoundException {
		String roadLinkPath = testGISDir + roadLinkFile;
		
		// Initialise test road link geography and context
		Context<RoadLink> roadLinkContext = new RoadLinkContext(contextName);
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> rlG = GeographyFactoryFinder.createGeographyFactory(null).createGeography(geogName, roadLinkContext, GeoParams);
		rlG.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(RoadLink.class, roadLinkPath, rlG, roadLinkContext);
		SpatialIndexManager.createIndex(rlG, RoadLink.class);
		
		context.addSubContext(roadLinkContext);
		
		return rlG;
	}
	
	static void setUpORRoadLinks() throws Exception {
		EnvironmentSetup.orRoadLinkGeography = setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp", GlobalVars.CONTEXT_NAMES.OR_ROAD_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.OR_ROAD_LINK_GEOGRAPHY);
	}
	
	static void setUpITNRoadLinks() throws Exception {
		EnvironmentSetup.roadLinkGeography = setUpRoadLinks("mastermap-itn RoadLink Intersect Within with orientation.shp", GlobalVars.CONTEXT_NAMES.ROAD_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.ROAD_LINK_GEOGRAPHY);
	}
	
	static void setUpPavementLinks(String linkFile, String contextName, String geogName) throws MalformedURLException, FileNotFoundException {
		EnvironmentSetup.pavementLinkGeography = setUpRoadLinks(linkFile, contextName, geogName);
	}
		
	static Geography<OD> setUpODs(String odFile, String geogName, String contextName) throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> ODContext = new DefaultContext<OD>(contextName);
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> odG = GeographyFactoryFinder.createGeographyFactory(null).createGeography(geogName, ODContext, GeoParamsOD);
		odG.setCRS(GlobalVars.geographyCRSString);
		
		context.addSubContext(ODContext);
		
		// Load vehicle origins and destinations
		String testODFile = testGISDir + odFile;
		GISFunctions.readShapefile(OD.class, testODFile, odG, ODContext);
		SpatialIndexManager.createIndex(odG, OD.class);
		return odG;		
	}
	
	static void setUpPedODs(String odFileName) throws MalformedURLException, FileNotFoundException {
		EnvironmentSetup.pedestrianDestinationGeography = setUpODs(odFileName, GlobalVars.CONTEXT_NAMES.PEDESTRIAN_DESTINATION_GEOGRAPHY, GlobalVars.CONTEXT_NAMES.PEDESTRIAN_DESTINATION_CONTEXT);
	}
	
	static void setUpPedODs() throws MalformedURLException, FileNotFoundException {
		setUpPedODs("OD_pedestrian_nodes_test.shp");
	}
	
	static void setUpVehicleODs() throws MalformedURLException, FileNotFoundException {
		EnvironmentSetup.vehicleDestinationGeography = setUpODs("OD_vehicle_nodes_intersect_within.shp", "vehicleODGeography", GlobalVars.CONTEXT_NAMES.VEHICLE_DESTINATION_CONTEXT);		
	}
	
	static void setUpVehicleODs(String file) throws MalformedURLException, FileNotFoundException {
		EnvironmentSetup.vehicleDestinationGeography = setUpODs(file, "vehicleODGeography", GlobalVars.CONTEXT_NAMES.VEHICLE_DESTINATION_CONTEXT);		
	}
	
	static void setUpPedObstructions() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<PedObstruction> pedObstructContext = new PedObstructionContext();
		GeographyParameters<PedObstruction> GeoParams = new GeographyParameters<PedObstruction>();
		EnvironmentSetup.pedObstructGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedObstructGeography", pedObstructContext, GeoParams);
		EnvironmentSetup.pedObstructGeography.setCRS(GlobalVars.geographyCRSString);
		context.addSubContext(pedObstructContext);
		
		// Load ped obstructions data
		String testPedObstructFile = testGISDir + IO.getProperty("PedestrianObstructionShapefile");
		GISFunctions.readShapefile(PedObstruction.class, testPedObstructFile, EnvironmentSetup.pedObstructGeography, pedObstructContext);
		SpatialIndexManager.createIndex(EnvironmentSetup.pedObstructGeography, PedObstruction.class);
	}
	
	static void setUpPedObstructionPoints() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<PedObstruction> pedObstructionPointsContext = new PedObstructionPointsContext();
		GeographyParameters<PedObstruction> GeoParams = new GeographyParameters<PedObstruction>();
		EnvironmentSetup.pedObstructionPointsGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedObstructionPointsGeography", pedObstructionPointsContext, GeoParams);
		EnvironmentSetup.pedObstructionPointsGeography.setCRS(GlobalVars.geographyCRSString);
		context.addSubContext(pedObstructionPointsContext);
		
		
		// Load ped obstructions data
		String testPedObstructFile = testGISDir + IO.getProperty("PedestrianObstructionPointsShapefile");
		GISFunctions.readShapefile(PedObstruction.class, testPedObstructFile, EnvironmentSetup.pedObstructionPointsGeography, pedObstructionPointsContext);
		SpatialIndexManager.createIndex(EnvironmentSetup.pedObstructionPointsGeography, PedObstruction.class);
	}
	
	static void setUpCrossingAlternatives(String caFile) throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<CrossingAlternative> caContext = new CAContext();
		GeographyParameters<CrossingAlternative> GeoParams = new GeographyParameters<CrossingAlternative>();
		EnvironmentSetup.caGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.CA_GEOGRAPHY, caContext, GeoParams);
		EnvironmentSetup.caGeography.setCRS(GlobalVars.geographyCRSString);
		context.addSubContext(caContext);
		
		
		// Load ped obstructions data
		String testCAFile = testGISDir + caFile;
		GISFunctions.readShapefile(CrossingAlternative.class, testCAFile, EnvironmentSetup.caGeography, caContext);
		SpatialIndexManager.createIndex(EnvironmentSetup.caGeography, CrossingAlternative.class);
	}
	
	static Network<Junction> setUpRoadNetwork(boolean isDirected, Geography<RoadLink> rlG, String name) {
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		Geography<Junction> junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.JUNCTION_GEOGRAPHY, junctionContext, GeoParamsJunc);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);
		context.addSubContext(junctionContext);
		
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(name,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		Network<Junction> rN = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(rlG, junctionContext, junctionGeography, rN);
		
		return rN;
	}
	
	static void setUpORRoadNetwork(boolean isDirected) {
		// Shouldnt use new context - use existing context
		Context<Junction> junctionContext = new JunctionContext(GlobalVars.CONTEXT_NAMES.OR_JUNCTION_CONTEXT);
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		EnvironmentSetup.orJunctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.OR_JUNCTION_GEOGRAPHY, junctionContext, GeoParamsJunc);
		EnvironmentSetup.orJunctionGeography.setCRS(GlobalVars.geographyCRSString);
		context.addSubContext(junctionContext);
		
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.OR_ROAD_NETWORK,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		EnvironmentSetup.orRoadNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(EnvironmentSetup.orRoadLinkGeography, junctionContext, EnvironmentSetup.orJunctionGeography, EnvironmentSetup.orRoadNetwork);
	}
	
	static void setUpITNRoadNetwork(boolean isDirected) {
		EnvironmentSetup.roadNetwork = setUpRoadNetwork(isDirected, EnvironmentSetup.roadLinkGeography, GlobalVars.CONTEXT_NAMES.ROAD_NETWORK);
	}
	
	static void setUpPavementNetwork() {
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.PAVEMENT_NETWORK, EnvironmentSetup.pavementJunctionContext, false);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		EnvironmentSetup.pavementNetwork = builder.buildNetwork();
		context.addSubContext(pavementJunctionContext);
		
		GISFunctions.buildGISRoadNetwork(EnvironmentSetup.pavementLinkGeography, EnvironmentSetup.pavementJunctionContext, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
	}
	
	static void setUpPedJunctions() throws Exception {
		String pedJPath = testGISDir + IO.getProperty("PavementJunctionsShapefile");
		
		// Initialise test road link geography and context
		EnvironmentSetup.pavementJunctionContext = new JunctionContext(GlobalVars.CONTEXT_NAMES.PAVEMENT_JUNCTION_CONTEXT);
		GeographyParameters<Junction> GeoParams = new GeographyParameters<Junction>();
		EnvironmentSetup.pavementJunctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.PAVEMENT_JUNCTION_GEOGRAPHY, EnvironmentSetup.pavementJunctionContext, GeoParams);
		EnvironmentSetup.pavementJunctionGeography.setCRS(GlobalVars.geographyCRSString);
		context.addSubContext(pavementJunctionContext);
				
		GISFunctions.readShapefile(Junction.class, pedJPath, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementJunctionContext);
		SpatialIndexManager.createIndex(EnvironmentSetup.pavementJunctionGeography, Junction.class);
	}
	
	List<RoadLink> planStrategicPath(Coordinate o, Coordinate d, String j1ID, String j2ID){
		
		// Plan the strategic path
		List<RoadLink> sP = new ArrayList<RoadLink>();		
		RoadNetworkRoute rnr = new RoadNetworkRoute(o , d);
		
		// Define my starting and ending junctions to test
		// Need to to this because although the origin and destination were selected above, this was just to initialise RoadNetworkRoute with
		// The route is actually calculated using junctions. 
		List<Junction> currentJunctions = new ArrayList<Junction>();
		List<Junction> destJunctions = new ArrayList<Junction>();
		for(Junction j: EnvironmentSetup.orRoadNetwork.getNodes()) {
			
			// Set the test current junctions 
			if (j.getFID().contentEquals(j1ID)) {
				currentJunctions.add(j);
			}
			
			// Set the test destination junctions
			if (j.getFID().contentEquals(j2ID)) {
				destJunctions.add(j);
			}
		}
		
		Junction[] routeEndpoints = new Junction[2];
		
		// Get shortest Route according to the Route class
		List<RepastEdge<Junction>> shortestRoute = null;
		try {
			shortestRoute = rnr.getShortestRoute(EnvironmentSetup.orRoadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		for(int i=0; i< shortestRoute.size(); i++) {
			NetworkEdge<Junction> e = (NetworkEdge<Junction>)shortestRoute.get(i);
			sP.add(e.getRoadLink());
		}
		
		return sP;
	}
	
	/*
	 * Code taken from EnvironmentSetup that connects Roads with Road Links and visa versa.
	 */
	static void assocaiteRoadsWithRoadLinks() {
		// Create lookups from OR road links to ITN links and visa versa
		for (RoadLink orRL: EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			List<RoadLink> itnLinks = new ArrayList<RoadLink>();
			
			for(RoadLink itnLink: EnvironmentSetup.roadLinkGeography.getAllObjects()) {
				if (itnLink.getPedRLID().contentEquals(orRL.getFID())) {
					SpaceBuilder.itnToOR.put(itnLink, orRL);
					itnLinks.add(itnLink);
				}
			}
			SpaceBuilder.orToITN.put(orRL, itnLinks);
			
			// Also assign Road objects to OR road links, so pedestrians can identify pavement polygons nearby.
			for (Road r: EnvironmentSetup.roadGeography.getAllObjects()) {
				if(r.getRoadLinkID().contentEquals(orRL.getFID())) {
					orRL.getRoads().add(r);
				}
			}
		}
		
		// Create lookup from OR road junctions to pavement junctions
		for (Junction orJ : EnvironmentSetup.orJunctionGeography.getAllObjects()) {
			if (!SpaceBuilder.orJuncToPaveJunc.containsKey(orJ)) {
				SpaceBuilder.orJuncToPaveJunc.put(orJ, new ArrayList<Junction>());
			}
			for (Junction paveJ : EnvironmentSetup.pavementJunctionGeography.getAllObjects()) {
				if(paveJ.getjuncNodeID().contentEquals(orJ.getFID())) {
					SpaceBuilder.orJuncToPaveJunc.get(orJ).add(paveJ);
				}
			}
		}
	}
	
	static Vehicle createVehicle(OD o, OD d) {
		Vehicle V = new Vehicle(GlobalVars.maxVehicleSpeed, GlobalVars.defaultVehicleAcceleration, GlobalVars.initialVehicleSpeed, o, d);
		EnvironmentSetup.context.add(V);
		Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		GeometryFactory fac = new GeometryFactory();
		Point pt = fac.createPoint(oCoord);
		Geometry vehicleCircle = pt.buffer(2);
		GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, vehicleCircle, V);
		V.setLoc();
		return V;
	}
	
	static Vehicle createVehicle(String oID, String dID) {
		OD o = null;
		OD d = null;
		for (OD od: EnvironmentSetup.vehicleDestinationGeography.getAllObjects()) {
			if (od.getFID().contentEquals(oID)) {
				o = od;
			}
			else if (od.getFID().contentEquals(dID)) {
				d = od;
			}
		}
		return createVehicle(o,d);
	}
	
	static Ped createPedestrian(Integer oID, Integer dID, String oFID, String dFID, boolean minimisesCrossing) {
		double s = GlobalVars.pedVavg;
		double m = GlobalVars.pedMassAv;
		double alpha = 0.5;
		double lambda = 1.0;
		double gamma = 0.9;
		double epsilon = 3.0;
		int tt = 30;
		double pH = 20.0;
		
		if ( (oFID==null) | (dFID==null) ) {
			return createPedestrian(oID, dID, s, m, alpha, lambda, gamma, epsilon, tt, minimisesCrossing, pH);
		}
		else {
			return createPedestrian(oFID, dFID, s, m, alpha, lambda, gamma, epsilon, tt, minimisesCrossing, pH);
		}
		
	}
	
	static Ped createPedestrian(int oID, int dID, Double s, Double m, Double alpha, Double lambda, Double gamma, Double epsilon, Integer tt, boolean minimiseCrossings, Double pH) {
		
		OD o = null;
		OD d = null;
		
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == oID) {
				o = i;
			}
			else if (i.getId() == dID) {
				d = i;
			}
		}
		
		Ped p = createPedestrian(o, d, s, m, alpha, lambda, gamma, epsilon, tt, minimiseCrossings, pH);
		return p;
	}
	
	static Ped createPedestrian(String oFID, String dFID, Double s, Double m, Double alpha, Double lambda, Double gamma, Double epsilon, Integer tt, boolean minimiseCrossings, Double pH) {
		
		OD o = null;
		OD d = null;
		
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getFID().contentEquals(oFID)) {
				o = i;
			}
			else if (i.getFID().contentEquals(dFID)) {
				d = i;
			}
		}
		
		Ped p = createPedestrian(o, d, s, m, alpha, lambda, gamma, epsilon, tt, minimiseCrossings, pH);
		return p;
	}
	
	static Ped createPedestrian(OD o, OD d, Double s, Double m, Double alpha, Double lambda, Double gamma, Double epsilon, Integer tt, boolean minimiseCrossings, Double pH) {

		Ped p = new Ped(o, d, s, m, alpha, lambda, gamma, epsilon, tt, minimiseCrossings, pH, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);

        EnvironmentSetup.context.add(p);        
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(p.getRad());		
		GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, circle, p);
		p.setLoc();
		return p;
	}
	
	public static void clearCaches() {
		RoadNetworkRoute.clearCaches();
		SpatialIndexManager.clearCaches();
	}
	
	public static Ped createPedAtLocation(Integer oID, Integer dID, String oFID, String dFID, boolean minimisesDistance, Coordinate c, double b) {
		// Create pedestrian and point it towards it's destination (which in this case is just across the road)
		Ped ped = EnvironmentSetup.createPedestrian(oID, dID, oFID, dFID, minimisesDistance);
		
		// Move ped to position and bearing that has caused an error in the simulation
        Point pt = GISFunctions.pointGeometryFromCoordinate(c);
		Geometry pGeomNew = pt.buffer(ped.getRad());
        GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, pGeomNew, ped);
		ped.setLoc();
		ped.setBearing(b);
		return ped;
	}

}
