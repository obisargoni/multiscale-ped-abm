package repastInterSim.tests;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
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
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.RoadNetworkRoute;

public class EnvironmentSetup {
	
	public static Context<Road> roadContext;
	public static Geography<Road> roadGeography;
	static String testGISDir = ".//data//test_gis_data//";
	static String serialisedLookupPath = testGISDir + "road_link_roads_cache.serialised";
	
	public EnvironmentSetup() {
		// TODO Auto-generated constructor stub
	}
	
	static void setUpProperties() throws IOException {
		IO.readProperties();
		EnvironmentSetup.clearCaches();
		SpaceBuilder.fac = new GeometryFactory();
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
		SpaceBuilder.context = new DefaultContext<Object>();
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		SpaceBuilder.geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, SpaceBuilder.context, geoParams);
		SpaceBuilder.geography.setCRS(GlobalVars.geographyCRSString);
		SpaceBuilder.context.add(SpaceBuilder.geography);
	}
	
	static void setUpRoads() throws Exception {
		// Get road geography
		Context<Road> testRoadContext = new RoadContext();
		GeographyParameters<Road> GeoParamsRoad = new GeographyParameters<Road>();
		EnvironmentSetup.roadGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testRoadGeography", testRoadContext, GeoParamsRoad);
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
	}
	
	static Geography<RoadLink> setUpRoadLinks(String roadLinkFile) throws MalformedURLException, FileNotFoundException {
		String roadLinkPath = testGISDir + roadLinkFile;
		
		// Initialise test road link geography and context
		Context<RoadLink> roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> rlG = GeographyFactoryFinder.createGeographyFactory(null).createGeography("orRoadLinkGeography", roadLinkContext, GeoParams);
		rlG.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(RoadLink.class, roadLinkPath, rlG, roadLinkContext);
		SpatialIndexManager.createIndex(rlG, RoadLink.class);
		
		return rlG;
	}
	
	static void setUpORRoadLinks() throws Exception {
		SpaceBuilder.orRoadLinkGeography = setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
	}
	
	static void setUpITNRoadLinks() throws Exception {
		SpaceBuilder.roadLinkGeography = setUpRoadLinks("mastermap-itn RoadLink Intersect Within with orientation.shp");
	}
	
	static void setUpPavementLinks(String linkFile) throws MalformedURLException, FileNotFoundException {
		SpaceBuilder.pavementLinkGeography = setUpRoadLinks(linkFile);
	}
		
	static Geography<OD> setUpODs(String odFile, String geogName) throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> ODContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> odG = GeographyFactoryFinder.createGeographyFactory(null).createGeography(geogName, ODContext, GeoParamsOD);
		odG.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String testODFile = testGISDir + odFile;
		GISFunctions.readShapefile(OD.class, testODFile, odG, ODContext);
		SpatialIndexManager.createIndex(odG, OD.class);
		return odG;		
	}
	
	static void setUpPedODs(String odFileName) throws MalformedURLException, FileNotFoundException {
		SpaceBuilder.pedestrianDestinationGeography = setUpODs(odFileName, "pedODGeography");
	}
	
	static void setUpPedODs() throws MalformedURLException, FileNotFoundException {
		setUpPedODs("OD_pedestrian_nodes_test.shp");
	}
	
	static void setUpVehicleODs() throws MalformedURLException, FileNotFoundException {
		SpaceBuilder.vehicleDestinationGeography = setUpODs("OD_vehicle_nodes_intersect_within.shp", "vehicleODGeography");		
	}
	
	static void setUpVehicleODs(String file) throws MalformedURLException, FileNotFoundException {
		SpaceBuilder.vehicleDestinationGeography = setUpODs(file, "vehicleODGeography");		
	}
	
	static void setUpPedObstructions() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<PedObstruction> pedObstructContext = new PedObstructionContext();
		GeographyParameters<PedObstruction> GeoParams = new GeographyParameters<PedObstruction>();
		SpaceBuilder.pedObstructGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedObstructGeography", pedObstructContext, GeoParams);
		SpaceBuilder.pedObstructGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testPedObstructFile = testGISDir + IO.getProperty("PedestrianObstructionShapefile");
		GISFunctions.readShapefile(PedObstruction.class, testPedObstructFile, SpaceBuilder.pedObstructGeography, pedObstructContext);
		SpatialIndexManager.createIndex(SpaceBuilder.pedObstructGeography, PedObstruction.class);
	}
	
	static void setUpPedObstructionPoints() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<PedObstruction> pedObstructionPointsContext = new PedObstructionPointsContext();
		GeographyParameters<PedObstruction> GeoParams = new GeographyParameters<PedObstruction>();
		SpaceBuilder.pedObstructionPointsGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedObstructionPointsGeography", pedObstructionPointsContext, GeoParams);
		SpaceBuilder.pedObstructionPointsGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testPedObstructFile = testGISDir + IO.getProperty("PedestrianObstructionPointsShapefile");
		GISFunctions.readShapefile(PedObstruction.class, testPedObstructFile, SpaceBuilder.pedObstructionPointsGeography, pedObstructionPointsContext);
		SpatialIndexManager.createIndex(SpaceBuilder.pedObstructionPointsGeography, PedObstruction.class);
	}
	
	static void setUpCrossingAlternatives(String caFile) throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<CrossingAlternative> caContext = new CAContext();
		GeographyParameters<CrossingAlternative> GeoParams = new GeographyParameters<CrossingAlternative>();
		SpaceBuilder.caGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("caGeography", caContext, GeoParams);
		SpaceBuilder.caGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testCAFile = testGISDir + caFile;
		GISFunctions.readShapefile(CrossingAlternative.class, testCAFile, SpaceBuilder.caGeography, caContext);
		SpatialIndexManager.createIndex(SpaceBuilder.caGeography, CrossingAlternative.class);
	}
	
	static Network<Junction> setUpRoadNetwork(boolean isDirected, Geography<RoadLink> rlG, String name) {
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		Geography<Junction> junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", junctionContext, GeoParamsJunc);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);
		
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(name,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		Network<Junction> rN = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(rlG, junctionContext, junctionGeography, rN);
		
		return rN;
	}
	
	static void setUpORRoadNetwork(boolean isDirected) {
		
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		SpaceBuilder.orJunctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.OR_JUNCTION_GEOGRAPHY, junctionContext, GeoParamsJunc);
		SpaceBuilder.orJunctionGeography.setCRS(GlobalVars.geographyCRSString);
		
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.OR_ROAD_NETWORK,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		SpaceBuilder.orRoadNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(SpaceBuilder.orRoadLinkGeography, junctionContext, SpaceBuilder.orJunctionGeography, SpaceBuilder.orRoadNetwork);
	}
	
	static void setUpITNRoadNetwork(boolean isDirected) {
		SpaceBuilder.roadNetwork = setUpRoadNetwork(isDirected, SpaceBuilder.roadLinkGeography, GlobalVars.CONTEXT_NAMES.ROAD_NETWORK);
	}
	
	static void setUpPavementNetwork() {
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>("PAVEMENT_NETWORK", SpaceBuilder.pavementJunctionContext, false);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		SpaceBuilder.pavementNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(SpaceBuilder.pavementLinkGeography, SpaceBuilder.pavementJunctionContext, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementNetwork);
	}
	
	static void setUpPedJunctions() throws Exception {
		String pedJPath = testGISDir + IO.getProperty("PavementJunctionsShapefile");
		
		// Initialise test road link geography and context
		SpaceBuilder.pavementJunctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParams = new GeographyParameters<Junction>();
		SpaceBuilder.pavementJunctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pavementJunctionGeography", SpaceBuilder.pavementJunctionContext, GeoParams);
		SpaceBuilder.pavementJunctionGeography.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(Junction.class, pedJPath, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementJunctionContext);
		SpatialIndexManager.createIndex(SpaceBuilder.pavementJunctionGeography, Junction.class);
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
		for(Junction j: SpaceBuilder.orRoadNetwork.getNodes()) {
			
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
			shortestRoute = rnr.getShortestRoute(SpaceBuilder.orRoadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
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
	 * Code taken from SpaceBuilder that connects Roads with Road Links and visa versa.
	 */
	static void assocaiteRoadsWithRoadLinks() {
		// Create lookups from OR road links to ITN links and visa versa
		for (RoadLink orRL: SpaceBuilder.orRoadLinkGeography.getAllObjects()) {
			List<RoadLink> itnLinks = new ArrayList<RoadLink>();
			
			for(RoadLink itnLink: SpaceBuilder.roadLinkGeography.getAllObjects()) {
				if (itnLink.getPedRLID().contentEquals(orRL.getFID())) {
					SpaceBuilder.itnToOR.put(itnLink, orRL);
					itnLinks.add(itnLink);
				}
			}
			SpaceBuilder.orToITN.put(orRL, itnLinks);
			
			// Also assign Road objects to OR road links, so pedestrians can identify pavement polygons nearby.
			for (Road r: SpaceBuilder.roadGeography.getAllObjects()) {
				if(r.getRoadLinkID().contentEquals(orRL.getFID())) {
					orRL.getRoads().add(r);
				}
			}
		}
		
		// Create lookup from OR road junctions to pavement junctions
		for (Junction orJ : SpaceBuilder.orJunctionGeography.getAllObjects()) {
			if (!SpaceBuilder.orJuncToPaveJunc.containsKey(orJ)) {
				SpaceBuilder.orJuncToPaveJunc.put(orJ, new ArrayList<Junction>());
			}
			for (Junction paveJ : SpaceBuilder.pavementJunctionGeography.getAllObjects()) {
				if(paveJ.getjuncNodeID().contentEquals(orJ.getFID())) {
					SpaceBuilder.orJuncToPaveJunc.get(orJ).add(paveJ);
				}
			}
		}
	}
	
	static Vehicle createVehicle(OD o, OD d) {
		Vehicle V = new Vehicle(GlobalVars.maxVehicleSpeed, GlobalVars.defaultVehicleAcceleration, GlobalVars.initialVehicleSpeed, o, d);
		SpaceBuilder.context.add(V);
		Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		GeometryFactory fac = new GeometryFactory();
		Point pt = fac.createPoint(oCoord);
		Geometry vehicleCircle = pt.buffer(2);
		GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, vehicleCircle, V);
		V.setLoc();
		return V;
	}
	
	static Vehicle createVehicle(String oID, String dID) {
		OD o = null;
		OD d = null;
		for (OD od: SpaceBuilder.vehicleDestinationGeography.getAllObjects()) {
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
		double pH = 20.0;
		
		if ( (oFID==null) | (dFID==null) ) {
			return createPedestrian(oID, dID, s, m, alpha, lambda, gamma, epsilon, minimisesCrossing, pH);
		}
		else {
			return createPedestrian(oFID, dFID, s, m, alpha, lambda, gamma, epsilon, minimisesCrossing, pH);
		}
		
	}
	
	static Ped createPedestrian(int oID, int dID, Double s, Double m, Double alpha, Double lambda, Double gamma, Double epsilon, boolean minimiseCrossings, Double pH) {
		
		OD o = null;
		OD d = null;
		
		for (OD i : SpaceBuilder.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == oID) {
				o = i;
			}
			else if (i.getId() == dID) {
				d = i;
			}
		}
		
		Ped p = createPedestrian(o, d, s, m, alpha, lambda, gamma, epsilon, minimiseCrossings, pH);
		return p;
	}
	
	static Ped createPedestrian(String oFID, String dFID, Double s, Double m, Double alpha, Double lambda, Double gamma, Double epsilon, boolean minimiseCrossings, Double pH) {
		
		OD o = null;
		OD d = null;
		
		for (OD i : SpaceBuilder.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getFID().contentEquals(oFID)) {
				o = i;
			}
			else if (i.getFID().contentEquals(dFID)) {
				d = i;
			}
		}
		
		Ped p = createPedestrian(o, d, s, m, alpha, lambda, gamma, epsilon, minimiseCrossings, pH);
		return p;
	}
	
	static Ped createPedestrian(OD o, OD d, Double s, Double m, Double alpha, Double lambda, Double gamma, Double epsilon, boolean minimiseCrossings, Double pH) {

		Ped p = new Ped(o, d, s, m, alpha, lambda, gamma, epsilon, minimiseCrossings, pH, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementNetwork);

        SpaceBuilder.context.add(p);        
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(p.getRad());		
		GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, circle, p);
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
        GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, pGeomNew, ped);
		ped.setLoc();
		ped.setBearing(b);
		return ped;
	}

}
