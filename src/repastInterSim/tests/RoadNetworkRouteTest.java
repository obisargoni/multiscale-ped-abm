package repastInterSim.tests;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.context.Context;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.OD;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.environment.contexts.VehicleDestinationContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.pathfinding.RoadNetworkRoute;

public class RoadNetworkRouteTest {
	
	private String TestDataDir = ".//data//test_gis_data//";
		
	private Boolean isDirected = true;
	
	private Context<RoadLink> roadLinkContext;
	
	private Context<Junction> junctionContext;
	private Network<Junction> roadNetwork;
	
	private Geography<Road> testRoadGeography = null;
	
	@BeforeEach
    public void setUp() throws Exception {
		
	    // Initialise contexts and geographies used by all tests	
		roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
		
		junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		Geography<Junction> junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", junctionContext, GeoParamsJunc);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);		
		
		// 1. Load road network data
		String roadLinkFile = TestDataDir + "mastermap-itn RoadLink Intersect Within with orientation.shp";
		GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
		// 2. Build road network
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		roadNetwork = builder.buildNetwork();
		GISFunctions.buildGISRoadNetwork(roadLinkGeography, junctionContext,junctionGeography, roadNetwork);
		
    }

	
	@Test
	public void testGetShortestRoute() throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> testODContext = new VehicleDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> testODGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", testODContext, GeoParamsOD);
		testODGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String vehicleDestinationsFile = TestDataDir + "OD_vehicle_nodes_intersect_within.shp";
		GISFunctions.readShapefile(OD.class, vehicleDestinationsFile, testODGeography, testODContext);
		
		// Select vehicle origins and destinations to test
		Coordinate o = testODContext.getObjects(OD.class).get(10).getGeom().getCoordinate();
		Coordinate d = testODContext.getObjects(OD.class).get(6).getGeom().getCoordinate();

		RoadNetworkRoute rnr = new RoadNetworkRoute(o , d);
		
		// Define my starting and ending junctions to test
		// Need to to this because although the origin and destination were selected above, this was just to initialise RoadNetworkRoute with
		// The route is actually calculated using junctions. 
		List<Junction> currentJunctions = new ArrayList<Junction>();
		List<Junction> destJunctions = new ArrayList<Junction>();
		for(Junction j: junctionContext.getObjects(Junction.class)) {
			
			// Set the test current junctions 
			if (j.getFID().contentEquals("osgb4000000029971605")) {
				currentJunctions.add(j);
			}
			
			// Set the test destination junctions
			if (j.getFID().contentEquals("osgb4000000029970678")) {
				destJunctions.add(j);
			}
		}
		
		Junction[] routeEndpoints = new Junction[2];
		
		// Get shortest Route according to the Route class
		List<RepastEdge<Junction>> shortestRoute = null;
		try {
			shortestRoute = rnr.getShortestRoute(roadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
		// Assert this matches what I expect - what I expect
		String[] expectedRoadLinkRouteFIDs = {"osgb4000000030419069", "osgb4000000030344373", "osgb5000005123276401",
												"osgb5000005123276405", "osgb4000000030343876", "osgb5000005156082205", 
												"osgb4000000030343920"};
		
		for(int i=0; i< shortestRoute.size(); i++) {
			NetworkEdge<Junction> e = (NetworkEdge<Junction>)shortestRoute.get(i);
			String expectedFID = expectedRoadLinkRouteFIDs[i];
			assert e.getRoadLink().getFID().contentEquals(expectedFID);
		}
		
	}
	
	@Test
	public void testGetShortestRoutePedOD1() throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> testODContext = new VehicleDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> testODGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", testODContext, GeoParamsOD);
		testODGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String vehicleDestinationsFile = TestDataDir + "OD_pedestrian_nodes.shp";
		GISFunctions.readShapefile(OD.class, vehicleDestinationsFile, testODGeography, testODContext);
		
		// Select pedestrian origins and destinations to test
		Coordinate o = testODContext.getObjects(OD.class).get(0).getGeom().getCoordinate();
		Coordinate d = testODContext.getObjects(OD.class).get(1).getGeom().getCoordinate();

		RoadNetworkRoute rnr = new RoadNetworkRoute(o , d);
		
		// Define my starting and ending junctions to test
		List<Junction> currentJunctions = new ArrayList<Junction>();
		List<Junction> destJunctions = new ArrayList<Junction>();
		String currentJunctionFID = "osgb4000000029970684";
		String destJunctionFID = "osgb4000000029970684";
		for(Junction j: junctionContext.getObjects(Junction.class)) {
			
			// Set the test current junctions 
			if (j.getFID().contentEquals(currentJunctionFID)) {
				currentJunctions.add(j);
			}
			
			// Set the test destination junctions
			if (j.getFID().contentEquals(destJunctionFID)) {
				destJunctions.add(j);
			}
		}
		
		Junction[] routeEndpoints = new Junction[2];
		
		// Get shortest Route according to the Route class
		List<RepastEdge<Junction>> shortestRoute = null;
		try {
			shortestRoute = rnr.getShortestRoute(roadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
		// Route is empty when start and end points are the same
		assert shortestRoute.isEmpty() == true;

	}
	
	@Test
	public void testGetShortestRoutePedOD2() throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> testODContext = new VehicleDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> testODGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", testODContext, GeoParamsOD);
		testODGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String vehicleDestinationsFile = TestDataDir + "OD_pedestrian_nodes.shp";
		GISFunctions.readShapefile(OD.class, vehicleDestinationsFile, testODGeography, testODContext);
		
		// Select pedestrian origins and destinations to test
		Coordinate o = testODContext.getObjects(OD.class).get(0).getGeom().getCoordinate();
		Coordinate d = testODContext.getObjects(OD.class).get(1).getGeom().getCoordinate();

		RoadNetworkRoute rnr = new RoadNetworkRoute(o , d);
		
		// Define my starting and ending junctions to test
		List<Junction> currentJunctions = new ArrayList<Junction>();
		List<Junction> destJunctions = new ArrayList<Junction>();
		String currentJunctionFID = "osgb4000000029970684";
		String destJunctionFID = "osgb4000000029970447";
		for(Junction j: junctionContext.getObjects(Junction.class)) {
			
			// Set the test current junctions 
			if (j.getFID().contentEquals(currentJunctionFID)) {
				currentJunctions.add(j);
			}
			
			// Set the test destination junctions
			if (j.getFID().contentEquals(destJunctionFID)) {
				destJunctions.add(j);
			}
		}
		
		Junction[] routeEndpoints = new Junction[2];
		
		// Get shortest Route according to the Route class
		List<RepastEdge<Junction>> shortestRoute = null;
		try {
			shortestRoute = rnr.getShortestRoute(roadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
		// Check route is made up on single expected road link when these start and end points given
		assert shortestRoute.size() == 1;
		NetworkEdge<Junction> e = (NetworkEdge<Junction>)shortestRoute.get(0);
		assert e.getRoadLink().getFID().contentEquals("osgb4000000030238946");

	}
	
	@Test
	public void testCalculateRouteParity() throws MalformedURLException, FileNotFoundException {
		// Initialise OD context and geography
		Context<OD> testODContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> testODGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", testODContext, GeoParamsOD);
		testODGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String testODFile = TestDataDir + "parity_test_OD.shp";
		GISFunctions.readShapefile(OD.class, testODFile, testODGeography, testODContext);
		
		// Select pedestrian origins and destinations to test
		Coordinate o = testODContext.getObjects(OD.class).get(0).getGeom().getCoordinate();
		Coordinate d = testODContext.getObjects(OD.class).get(1).getGeom().getCoordinate();
				
		// Initialise test road link geography and context
		roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
		
		String roadLinkFile = TestDataDir + "parity_test_lines3.shp";
		GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
		List<RoadLink> roads = new ArrayList<RoadLink>();
		roadLinkGeography.getAllObjects().forEach(roads::add);
		
		// Now perform parity calculation
		int p = RoadNetworkRoute.calculateRouteParity(o, d, roads);
		assert p == 1;
	}
	
	void setUpRoads() {
		
		TestDataDir = ".//data//test_gis_data//";
		
		// Get road geography
		Context<Road> testRoadContext = new RoadContext();
		GeographyParameters<Road> GeoParamsRoad = new GeographyParameters<Road>();
		testRoadGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testRoadGeography", testRoadContext, GeoParamsRoad);
		testRoadGeography.setCRS(GlobalVars.geographyCRSString);		
		
		// Load vehicle origins and destinations
		try {
			GISFunctions.readShapefile(Road.class, TestDataDir + "topographicAreaVehicle.shp", testRoadGeography, testRoadContext);
			GISFunctions.readShapefile(Road.class, TestDataDir + "topographicAreaPedestrian.shp", testRoadGeography, testRoadContext);
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	@Test
	public void testGetRoadLinkRoadCacheInstance() throws Exception {
		
		setUpRoads();
		
		File vehcileRoadsFile = new File(TestDataDir + "topographicAreaVehicle.shp");
		File pedestrianRoadsFile = new File(TestDataDir + "topographicAreaPedestrian.shp");
		File serialisedLoc = new File(TestDataDir + "road_link_rodas_cache.serialised");
		
		String roadLinkID = "osgb4000000030343774";
		
		List<Road> roads = RoadNetworkRoute.getRoadLinkRoads(testRoadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Check the expected number of roads have been returned		
		assert roads.size() == 3;
	}
	
	@Test
	public void testGetRoadLinkPedestrianRoads() throws Exception {
		
		setUpRoads();
		
		File vehcileRoadsFile = new File(TestDataDir + "topographicAreaVehicle.shp");
		File pedestrianRoadsFile = new File(TestDataDir + "topographicAreaPedestrian.shp");
		File serialisedLoc = new File(TestDataDir + "road_link_rodas_cache.serialised");
		
		String roadLinkID = "osgb4000000030343774";
		
		
		List<Road> pedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(testRoadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Check the expected number of roads have been returned		
		assert pedRoads.size() == 2;
	}

}
