package repastInterSim.tests;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
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
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.RoadNetworkRoute;

public class RoadNetworkRouteTest {
	
	private String TestDataDir = ".//data//test_gis_data//";
		
	private Boolean isDirected = true;
	
	void setUpProperties() throws IOException {
		EnvironmentSetup.setUpProperties();
	}
	
	@BeforeEach
    public void setUp() throws Exception {
		
		setUpProperties();
		
	    // Initialise contexts and geographies used by all tests	
		EnvironmentSetup.roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", EnvironmentSetup.roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
		
		EnvironmentSetup.junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		Geography<Junction> junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", EnvironmentSetup.junctionContext, GeoParamsJunc);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);		
		
		// 1. Load road network data
		String roadLinkFile = TestDataDir + "mastermap-itn RoadLink Intersect Within with orientation.shp";
		GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, EnvironmentSetup.roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
		// 2. Build road network
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,EnvironmentSetup.junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		EnvironmentSetup.roadNetwork = builder.buildNetwork();
		GISFunctions.buildGISRoadNetwork(roadLinkGeography, EnvironmentSetup.junctionContext,junctionGeography, EnvironmentSetup.roadNetwork, false);
		
    }

	
	@Test
	public void testGetShortestRoute() throws MalformedURLException, FileNotFoundException {
		
		double start = System.currentTimeMillis();
		
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
		for(Junction j: EnvironmentSetup.junctionContext.getObjects(Junction.class)) {
			
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
			shortestRoute = rnr.getShortestRoute(EnvironmentSetup.roadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
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
		
		double duration = System.currentTimeMillis() - start;
		System.out.print("testGetShortestRoute() duration: \n" + duration + "\n");
	}
	
	@Test
	public void testGetShortestRoutePedOD1() throws MalformedURLException, FileNotFoundException {
		double start = System.currentTimeMillis();
		
		// Initialise OD context and geography
		Context<OD> testODContext = new VehicleDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> testODGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", testODContext, GeoParamsOD);
		testODGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String vehicleDestinationsFile = TestDataDir + "OD_pedestrian_nodes_test.shp";
		GISFunctions.readShapefile(OD.class, vehicleDestinationsFile, testODGeography, testODContext);
		
		// Select pedestrian origins and destinations to test
		Coordinate o = null;
		Coordinate d = null;
		for ( OD od: testODContext.getObjects(OD.class)) {
			if (od.getId() == 1) {
				o = od.getGeom().getCoordinate();
			}
			else if (od.getId() == 2) {
				d = od.getGeom().getCoordinate();
			}
		}

		RoadNetworkRoute rnr = new RoadNetworkRoute(o , d);
		
		// Define my starting and ending junctions to test
		List<Junction> currentJunctions = new ArrayList<Junction>();
		List<Junction> destJunctions = new ArrayList<Junction>();
		String currentJunctionFID = "osgb4000000029970684";
		String destJunctionFID = "osgb4000000029970684";
		for(Junction j: EnvironmentSetup.junctionContext.getObjects(Junction.class)) {
			
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
			shortestRoute = rnr.getShortestRoute(EnvironmentSetup.roadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
		// Route is empty when start and end points are the same
		assert shortestRoute.isEmpty() == true;
		
		double duration = System.currentTimeMillis() - start;
		System.out.print("testGetShortestRoutePedOD1() duration: \n" + duration + "\n");
	}
	
	@Test
	public void testGetShortestRoutePedOD2() throws MalformedURLException, FileNotFoundException {
		double start = System.currentTimeMillis();
		
		// Initialise OD context and geography
		Context<OD> testODContext = new VehicleDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> testODGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", testODContext, GeoParamsOD);
		testODGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String vehicleDestinationsFile = TestDataDir + "OD_pedestrian_nodes_test.shp";
		GISFunctions.readShapefile(OD.class, vehicleDestinationsFile, testODGeography, testODContext);
		
		// Select pedestrian origins and destinations to test
		Coordinate o = null;
		Coordinate d = null;
		for ( OD od: testODContext.getObjects(OD.class)) {
			if (od.getId() == 1) {
				o = od.getGeom().getCoordinate();
			}
			else if (od.getId() == 10) {
				d = od.getGeom().getCoordinate();
			}
		}

		RoadNetworkRoute rnr = new RoadNetworkRoute(o , d);
		
		// Define my starting and ending junctions to test
		List<Junction> currentJunctions = new ArrayList<Junction>();
		List<Junction> destJunctions = new ArrayList<Junction>();
		String currentJunctionFID = "osgb4000000029970684";
		String destJunctionFID = "osgb4000000029970447";
		for(Junction j: EnvironmentSetup.junctionContext.getObjects(Junction.class)) {
			
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
			shortestRoute = rnr.getShortestRoute(EnvironmentSetup.roadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
		// Check route is made up on single expected road link when these start and end points given
		assert shortestRoute.size() == 1;
		NetworkEdge<Junction> e = (NetworkEdge<Junction>)shortestRoute.get(0);
		assert e.getRoadLink().getFID().contentEquals("osgb4000000030238946");
		
		double duration = System.currentTimeMillis() - start;
		System.out.print("testGetShortestRoutePedOD2() duration: \n" + duration + "\n");
	}
	
	/*
	 * Calculates shortest route using the OR road Network. Also gets the start junctions using the OD to pavement junction cache.
	 */
	@Test
	public void testGetShortestRoutePedOD3() throws MalformedURLException, FileNotFoundException {
		
		try {
			IO.readProperties();
			EnvironmentSetup.setUpPedODs();
			EnvironmentSetup.setUpPedJunctions();
			
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
		} catch (Exception e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		
		double start = System.currentTimeMillis();

		
		// Select pedestrian origins and destinations to test
		Coordinate o = null;
		Coordinate d = null;
		for ( OD od: EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (od.getId() == 1) {
				o = od.getGeom().getCoordinate();
			}
			else if (od.getId() == 10) {
				d = od.getGeom().getCoordinate();
			}
		}

		RoadNetworkRoute rnr = new RoadNetworkRoute(o , d);
		
		// Get current and destination junction using the cache
		File odFile = new File(TestDataDir + "OD_pedestrian_nodes_test.shp");
		File pavementJunctionFile = new File(TestDataDir + IO.getProperty(GlobalVars.PavementJunctionShapeFile));
		File paveJuncSeriealizedLoc = new File(TestDataDir + IO.getProperty(GlobalVars.ODPavementJunctionCache));
		Junction oPavementJ = null;
		Junction dPavementJ = null;
		try {
			oPavementJ = RoadNetworkRoute.getNearestpavementJunctionToOD(o, odFile, pavementJunctionFile, paveJuncSeriealizedLoc);
			dPavementJ = RoadNetworkRoute.getNearestpavementJunctionToOD(d, odFile, pavementJunctionFile, paveJuncSeriealizedLoc);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		// Define my starting and ending junctions to test
		String currentJunctionFID = "node_id_73";
		String destJunctionFID = "node_id_60";
		
		// Test that the junctions returned are as expected
		assert oPavementJ.getFID().contentEquals("pave_node_111");
		assert dPavementJ.getFID().contentEquals("pave_node_91");
		assert oPavementJ.getjuncNodeID().contentEquals(currentJunctionFID);
		assert dPavementJ.getjuncNodeID().contentEquals(destJunctionFID);
		
		List<Junction> currentJunctions = new ArrayList<Junction>();
		List<Junction> destJunctions = new ArrayList<Junction>();
		for(Junction j: EnvironmentSetup.orRoadNetwork.getNodes()) {
			
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
			shortestRoute = rnr.getShortestRoute(EnvironmentSetup.orRoadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
		// Check route is made up on single expected road link when these start and end points given
		assert shortestRoute.size() == 1;
		NetworkEdge<Junction> e = (NetworkEdge<Junction>)shortestRoute.get(0);
		assert e.getRoadLink().getFID().contentEquals("A8675945-DE94-4E22-9905-B0623A326221_0");
		
		double duration = System.currentTimeMillis() - start;
		System.out.print("testGetShortestRoutePedOD3() duration: \n" + duration + "\n");
	}

	
	@Test
	public void testCalculateRouteParity() throws MalformedURLException, FileNotFoundException {
		double start = System.currentTimeMillis();

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
		EnvironmentSetup.roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", EnvironmentSetup.roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
		
		String roadLinkFile = TestDataDir + "parity_test_lines3.shp";
		GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, EnvironmentSetup.roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
		List<RoadLink> roads = new ArrayList<RoadLink>();
		roadLinkGeography.getAllObjects().forEach(roads::add);
		
		// Now perform parity calculation
		int p = RoadNetworkRoute.calculateRouteParity(o, d, roads);
		assert p == 1;
		
		double duration = System.currentTimeMillis() - start;
		System.out.print("testCalculateRouteParity() duration: \n" + duration + "\n");
	}
	
	void setUpRoads() {
		
		TestDataDir = ".//data//test_gis_data//";
		
		// Get road geography
		Context<Road> testRoadContext = new RoadContext();
		GeographyParameters<Road> GeoParamsRoad = new GeographyParameters<Road>();
		EnvironmentSetup.roadGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadGeography", testRoadContext, GeoParamsRoad);
		EnvironmentSetup.roadGeography.setCRS(GlobalVars.geographyCRSString);		
		
		// Load vehicle origins and destinations
		try {
			GISFunctions.readShapefile(Road.class, TestDataDir + "topographicAreaVehicle.shp", EnvironmentSetup.roadGeography, testRoadContext);
			GISFunctions.readShapefile(Road.class, TestDataDir + "topographicAreaPedestrian.shp", EnvironmentSetup.roadGeography, testRoadContext);
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	@Test
	public void testGetRoadLinkRoadCacheInstance() throws Exception {

		setUpRoads();
		
		double start = System.currentTimeMillis();
		
		File vehcileRoadsFile = new File(TestDataDir + "topographicAreaVehicle.shp");
		File pedestrianRoadsFile = new File(TestDataDir + "topographicAreaPedestrian.shp");
		File serialisedLoc = new File(TestDataDir + "road_link_roads_cache.serialised");
		
		
		// Road link pedestrian roads cache requires OR link id
		String roadLinkID = "F4C0B1FB-762C-4492-BB0D-673CC4950CBE_0";
		List<Road> roads = RoadNetworkRoute.getRoadLinkRoads(EnvironmentSetup.roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Check the expected number of roads have been returned		
		assert roads.size() == 3;
		
		double duration = System.currentTimeMillis() - start;
		System.out.print("testGetRoadLinkRoadCacheInstance() duration: \n" + duration + "\n");
	}
	
	@Test
	public void testGetRoadLinkPedestrianRoads() throws Exception {
		setUpRoads();
		
		double start = System.currentTimeMillis();
		
		File vehcileRoadsFile = new File(TestDataDir + "topographicAreaVehicle.shp");
		File pedestrianRoadsFile = new File(TestDataDir + "topographicAreaPedestrian.shp");
		File serialisedLoc = new File(TestDataDir + "road_link_roads_cache.serialised");
		
		// Road link pedestrian roads cache requires OR link id
		String roadLinkID = "F4C0B1FB-762C-4492-BB0D-673CC4950CBE_0";
		List<Road> pedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(EnvironmentSetup.roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Check the expected number of roads have been returned		
		assert pedRoads.size() == 2;
		
		double duration = System.currentTimeMillis() - start;
		System.out.print("testGetRoadLinkPedestrianRoads() duration: \n" + duration + "\n");
	}
	
	/*
	 * Test is suppose to replicate problem of different shortest path being produced between runs but doesn't do this because problem
	 * is caused by road links appearing in a different order in the context and geography each simulation run.
	 */
	@Test
	public void testGetShortestRouteQuadGrid1() throws MalformedURLException, FileNotFoundException {
		
		try {
			String origTestDir = EnvironmentSetup.testGISDir;
			EnvironmentSetup.testGISDir += "/QuadGrid145NodesBuffer/";
			
			EnvironmentSetup.setUpProperties();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpPedODs("OD_pedestrian_nodes.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
			EnvironmentSetup.testGISDir = origTestDir;
		} catch (Exception e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		
		// Select pedestrian origins and destinations to test
		Coordinate o = null;
		Coordinate d = null;
		for ( OD od: EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (od.getFID().contentEquals("od_132")) {
				o = od.getGeom().getCoordinate();
			}
			else if (od.getFID().contentEquals("od_0")) {
				d = od.getGeom().getCoordinate();
			}
		}
		
		RoadNetworkRoute rnr = new RoadNetworkRoute(o, d);
		
		// initialise object to record the start and end pavement junctions of the route
		Junction[] routeEnds = null;
		
		// Find shortest path using road network route
		try {
			routeEnds = rnr.setRoadLinkRoute(EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get path of road links and set this as the strategic path
		List<RoadLink> strategicPath = rnr.getRoadsX();
		String fullStrategicPathString="";
		for (RoadLink rl: strategicPath) {
			fullStrategicPathString = fullStrategicPathString + ":" + rl.getFID();
		}
		
		for (int i=0; i<10; i++) {
			// Repeatedly calculate shortest path and check that it is the same to the first one calculated
			RoadNetworkRoute rnri = new RoadNetworkRoute(o, d);
			Junction[] routeEndsi = null;
			try {
				routeEndsi = rnri.setRoadLinkRoute(EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			List<RoadLink> strategicPathi = rnri.getRoadsX();
			String fullStrategicPathStringi="";
			for (RoadLink rl: strategicPathi) {
				fullStrategicPathStringi = fullStrategicPathStringi + ":" + rl.getFID();
			}
			
			assert fullStrategicPathString.contentEquals(fullStrategicPathStringi);
		}	
		
	}

}
