package repastInterSim.tests;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;

import repast.simphony.context.Context;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.space.graph.Network;
import repast.simphony.util.collections.IndexedIterable;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.OD;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.PedObstructionContext;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.PedPathFinder;
import repastInterSim.pathfinding.RoadNetworkRoute;

class GISFunctionsTest {
	
	String testGISDir = ".//data//test_gis_data//";
	String pedestrianRoadsPath = null;
	String vehicleRoadsPath = null;
	String roadLinkPath = null;
	String serialisedLookupPath = null;
	
	private String TestDataDir = ".//data//test_gis_data//";

	void setUp(String lineDataFile) throws Exception {
		
	    // Initialise contexts and geographies used by all tests	
		EnvironmentSetup.roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", EnvironmentSetup.roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
		
		EnvironmentSetup.pedestrianDestinationContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> pedestrianDestinationGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedestrianDestinationGeography", EnvironmentSetup.pedestrianDestinationContext, GeoParamsOD);
		pedestrianDestinationGeography.setCRS(GlobalVars.geographyCRSString);
		
		// 1. Load road network data
		String roadLinkFile = TestDataDir + lineDataFile;
		GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, EnvironmentSetup.roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
		String testODFile = TestDataDir + "parity_test_OD.shp";
		GISFunctions.readShapefile(OD.class, testODFile, pedestrianDestinationGeography, EnvironmentSetup.pedestrianDestinationContext);
		SpatialIndexManager.createIndex(pedestrianDestinationGeography, OD.class);
	}
	
	void setUpProperties() throws IOException {
		IO.readProperties();
	}
	
	void setUpRoads() throws Exception {
		pedestrianRoadsPath = testGISDir + "topographicAreaPedestrian.shp";
		vehicleRoadsPath = testGISDir + "topographicAreaVehicle.shp";
		serialisedLookupPath = testGISDir + "road_link_roads_cache.serialised";
		
		// Get road geography
		Context<Road> testRoadContext = new RoadContext();
		GeographyParameters<Road> GeoParamsRoad = new GeographyParameters<Road>();
		EnvironmentSetup.roadGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testRoadGeography", testRoadContext, GeoParamsRoad);
		EnvironmentSetup.roadGeography.setCRS(GlobalVars.geographyCRSString);

		
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
	
	void setUpRoadLinks() throws Exception {
		setUpRoadLinks("mastermap-itn RoadLink Intersect Within with orientation.shp");
	}
	
	void setUpRoadLinks(String roadLinkFile) throws Exception {
		
		roadLinkPath = testGISDir + roadLinkFile;
		
		// Initialise test road link geography and context
		Context<RoadLink> roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		EnvironmentSetup.roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", roadLinkContext, GeoParams);
		EnvironmentSetup.roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(RoadLink.class, roadLinkPath, EnvironmentSetup.roadLinkGeography, roadLinkContext);
		SpatialIndexManager.createIndex(EnvironmentSetup.roadLinkGeography, RoadLink.class);
		
	}
		
	void setUpODs(String odFile) throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> ODContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		EnvironmentSetup.pedestrianDestinationGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", ODContext, GeoParamsOD);
		EnvironmentSetup.pedestrianDestinationGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String testODFile = testGISDir + odFile;
		GISFunctions.readShapefile(OD.class, testODFile, EnvironmentSetup.pedestrianDestinationGeography, ODContext);
		SpatialIndexManager.createIndex(EnvironmentSetup.pedestrianDestinationGeography, OD.class);
	}
	
	void setUpPedObstructions() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<PedObstruction> pedObstructContext = new PedObstructionContext();
		GeographyParameters<PedObstruction> GeoParams = new GeographyParameters<PedObstruction>();
		EnvironmentSetup.pedObstructGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedObstructGeography", pedObstructContext, GeoParams);
		EnvironmentSetup.pedObstructGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testPedObstructFile = testGISDir + "boundaryPedestrianVehicleArea.shp";
		GISFunctions.readShapefile(PedObstruction.class, testPedObstructFile, EnvironmentSetup.pedObstructGeography, pedObstructContext);
		SpatialIndexManager.createIndex(EnvironmentSetup.pedObstructGeography, PedObstruction.class);
	}
	
	void setUpRoadNetwork() {
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		EnvironmentSetup.junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", junctionContext, GeoParamsJunc);
		EnvironmentSetup.junctionGeography.setCRS(GlobalVars.geographyCRSString);		
		
		// 2. Build road network
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, true);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		EnvironmentSetup.roadNetwork = builder.buildNetwork();
		GISFunctions.buildGISRoadNetwork(EnvironmentSetup.roadLinkGeography, junctionContext,EnvironmentSetup.junctionGeography, EnvironmentSetup.roadNetwork, false);
	}
	
	int getNumberIntersectingCoords(String lineDataFile) {
		try {
			setUp(lineDataFile);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get number of intersecting points and compare to expected value
		Coordinate[] odCoords = new Coordinate[EnvironmentSetup.pedestrianDestinationContext.getObjects(OD.class).size()];
		int i = 0;
		for (OD od : EnvironmentSetup.pedestrianDestinationContext.getObjects(OD.class)) {
			odCoords[i] = od.getGeom().getCoordinate();
			i++;
		}
		LineString ODLine = GISFunctions.lineStringGeometryFromCoordinates(odCoords);
		Geometry[] ODLineGeom = {ODLine};

		
		// Loop through road links in the rout and count number of times the ODLine intersects
		Geometry[] rlGeoms = new Geometry[EnvironmentSetup.roadLinkContext.getObjects(RoadLink.class).size()];
		i = 0;
		for (RoadLink rl: EnvironmentSetup.roadLinkContext.getObjects(RoadLink.class)) {
			rlGeoms[i] = rl.getGeom();
			i++;
		}
		
		// Now calculate number of intersecting coordinates
		int nIntersections = GISFunctions.calculateNIntersectionCoords(ODLineGeom, rlGeoms);
		
		return nIntersections;
	}


	@Test
	void testCalculateNIntersectionCoords() {
		
		// Now calculate number of intersecting coordinates
		int nIntersections1 = getNumberIntersectingCoords("parity_test_lines1.shp");
		int nIntersections2 = getNumberIntersectingCoords("parity_test_lines2.shp");
		int nIntersections3 = getNumberIntersectingCoords("parity_test_lines3.shp");
		
		assert nIntersections1 == 1;
		assert nIntersections2 == 2;
		assert nIntersections3 == 3;
	}
	
	@Test
	void testAngleBetweenConnectedLineStrings() throws Exception {
		LineString l1 = null;
		LineString l2 = null;
		
		setUp("parity_test_lines1.shp");
		
		IndexedIterable<RoadLink> lines = EnvironmentSetup.roadLinkContext.getObjects(RoadLink.class);
		
		l1 = (LineString) lines.get(0).getGeom();
		l2 = (LineString) lines.get(1).getGeom();
		
		Double ang = GISFunctions.angleBetweenConnectedLineStrings(l1, l2);
		
		assert ang == 0.0;
		
		Coordinate c1 = new Coordinate(0,0);
		Coordinate c2 = new Coordinate(0,1);
		Coordinate c3 = new Coordinate(1,2);
		Coordinate c4 = new Coordinate(-1,0);
		
		Coordinate[] l1Coords = {c1,c2};
		Coordinate[] l2Coords = {c2,c3};
		l1 = GISFunctions.lineStringGeometryFromCoordinates(l1Coords);
		l2 = GISFunctions.lineStringGeometryFromCoordinates(l2Coords);
		
		ang = GISFunctions.angleBetweenConnectedLineStrings(l1, l2);
		assert Math.round(ang) == 45.0;
		
		l2Coords[1] = c4;
		l2 = GISFunctions.lineStringGeometryFromCoordinates(l2Coords);
		
		ang = GISFunctions.angleBetweenConnectedLineStrings(l1, l2);
		assert Math.round(ang) == 135.0;
		
		l1Coords[0] = c2;
		l1Coords[1] = c1;
		l1 = GISFunctions.lineStringGeometryFromCoordinates(l1Coords);
		
		ang = GISFunctions.angleBetweenConnectedLineStrings(l1, l2);
		assert Math.round(ang) == 135.0;		
	}
	
	@Test
	void testFarthestUnobstructedRoadCoordinate() throws Exception {
		
		setUpRoads();
		setUpODs("test_ped_OD1.shp");
		setUpPedObstructions();
		
		File vehcileRoadsFile = new File(vehicleRoadsPath);
		File pedestrianRoadsFile = new File(pedestrianRoadsPath);
		File serialisedLoc = new File(serialisedLookupPath);
		
		// Select origin coordinate and pedestrian road geometry
		List<OD> ods = new ArrayList<OD>();
		EnvironmentSetup.pedestrianDestinationGeography.getAllObjects().iterator().forEachRemaining(ods::add);
		Coordinate o = ods.get(0).getGeom().getCoordinate();
		
		String roadLinkID = "A8675945-DE94-4E22-9905-B0623A326221_0";
		
		List<Road> currentPedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(EnvironmentSetup.roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Expected coords
		Coordinate c1 = new Coordinate(530507.95, 180893.35);
		Coordinate c3 = new Coordinate(530468.0087569832, 180871.8784368495);
		Coordinate c4 = new Coordinate(530482.8182132206, 180870.19519803385);
		Coordinate[] carray = {c1,null,c3,c4};
		
		List<Coordinate> expectedCoords = new ArrayList<Coordinate>(Arrays.asList(carray));
		
		
		Coordinate destCoord = null;
		for (int i=0;i<currentPedRoads.size(); i++) {
			Road r = currentPedRoads.get(i);
			destCoord = GISFunctions.farthestUnobstructedGeomCoordinate(o, r.getGeom(), EnvironmentSetup.pedObstructGeography);			
			assert expectedCoords.contains(destCoord);
		}
	}
	
	@Test
	void testNearestUnobstructedRoadCoordinate() throws Exception {
		
		setUpRoads();
		setUpODs("test_ped_OD1.shp");
		setUpPedObstructions();
		
		File vehcileRoadsFile = new File(vehicleRoadsPath);
		File pedestrianRoadsFile = new File(pedestrianRoadsPath);
		File serialisedLoc = new File(serialisedLookupPath);
		
		// Select origin coordinate and pedestrian road geometry
		List<OD> ods = new ArrayList<OD>();
		EnvironmentSetup.pedestrianDestinationGeography.getAllObjects().iterator().forEachRemaining(ods::add);
		Coordinate o = ods.get(0).getGeom().getCoordinate();
		
		String roadLinkID = "A8675945-DE94-4E22-9905-B0623A326221_0";
		
		List<Road> currentPedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(EnvironmentSetup.roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Expected coords
		Coordinate c1 = new Coordinate(530506.8,180891.45);
		Coordinate c2 = new Coordinate(530512.5,180907.6);
		Coordinate c4 = new Coordinate(530521.6192518127,180903.04127475937);
		Coordinate[] carray = {c1,c2,null,c4};
		
		List<Coordinate> expectedCoords = new ArrayList<Coordinate>(Arrays.asList(carray));
		
		
		Coordinate destCoord = null;
		for (int i=0;i<currentPedRoads.size(); i++) {
			Road r = currentPedRoads.get(i);
			destCoord = GISFunctions.nearestUnobstructedGeomCoordinate(o, r.getGeom(), EnvironmentSetup.pedObstructGeography);
			assert expectedCoords.contains(destCoord);
		}
	}
	
	/*
	 * Test method for calculating bearing of vector between two coordinates
	 */
	@Test
	void testBearingBetweenCoordinates() {
		// Create some coordinates and calculate bearing between them
		Coordinate c1 = new Coordinate(530683.9458871038, 181012.96436357914);
		Coordinate c2 = new Coordinate(530706, 181029);
		Coordinate c3 = new Coordinate(530612, 181056);
		
		long start = System.currentTimeMillis();
		
		double bearing = GISFunctions.bearingBetweenCoordinates(c1, c2);
		assert Double.compare(bearing, 0.9421103253862521) == 0;
		
		bearing = GISFunctions.bearingBetweenCoordinates(c2, c1);
		assert Double.compare(bearing, 4.083702978976045) == 0;
		
		bearing = GISFunctions.bearingBetweenCoordinates(c1, c3);
		assert Double.compare(bearing, 5.251459401737451) == 0;
		
		bearing = GISFunctions.bearingBetweenCoordinates(c3, c2);
		assert Double.compare(bearing, 1.8505004792585276)==0;
		
		long duration = System.currentTimeMillis() - start;
		System.out.print(duration);
	}

}
