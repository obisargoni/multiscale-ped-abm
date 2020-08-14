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
import repastInterSim.pathfinding.PedPathFinder;
import repastInterSim.pathfinding.RoadNetworkRoute;

class GISFunctionsTest {
	
	Geography<Road> roadGeography = null;
	Geography<RoadLink> roadLinkGeography = null;
	Geography<OD> odGeography = null;
	Geography<PedObstruction> pedObstructGeography = null;
	Network<Junction> roadNetwork = null;
	Geography<Junction> junctionGeography = null;
	private Context<RoadLink> roadLinkContext;
	private Context<OD> pedestrianDestinationContext;
	
	String testGISDir = ".//data//test_gis_data//";
	String pedestrianRoadsPath = null;
	String vehicleRoadsPath = null;
	String roadLinkPath = null;
	String serialisedLookupPath = null;
	
	private String TestDataDir = ".//data//test_gis_data//";

	void setUp(String lineDataFile) throws Exception {
		
	    // Initialise contexts and geographies used by all tests	
		roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
		
		pedestrianDestinationContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> pedestrianDestinationGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedestrianDestinationGeography", pedestrianDestinationContext, GeoParamsOD);
		pedestrianDestinationGeography.setCRS(GlobalVars.geographyCRSString);
		
		// 1. Load road network data
		String roadLinkFile = TestDataDir + lineDataFile;
		GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
		String testODFile = TestDataDir + "parity_test_OD.shp";
		GISFunctions.readShapefile(OD.class, testODFile, pedestrianDestinationGeography, pedestrianDestinationContext);
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
		roadGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testRoadGeography", testRoadContext, GeoParamsRoad);
		roadGeography.setCRS(GlobalVars.geographyCRSString);

		
		// Load vehicle origins and destinations
		try {
			GISFunctions.readShapefile(Road.class, vehicleRoadsPath, roadGeography, testRoadContext);
			GISFunctions.readShapefile(Road.class, pedestrianRoadsPath, roadGeography, testRoadContext);
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		SpatialIndexManager.createIndex(roadGeography, Road.class);

	}
	
	void setUpRoadLinks() throws Exception {
		setUpRoadLinks("mastermap-itn RoadLink Intersect Within with orientation.shp");
	}
	
	void setUpRoadLinks(String roadLinkFile) throws Exception {
		
		roadLinkPath = testGISDir + roadLinkFile;
		
		// Initialise test road link geography and context
		Context<RoadLink> roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(RoadLink.class, roadLinkPath, roadLinkGeography, roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
	}
		
	void setUpODs(String odFile) throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> ODContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		odGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", ODContext, GeoParamsOD);
		odGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String testODFile = testGISDir + odFile;
		GISFunctions.readShapefile(OD.class, testODFile, odGeography, ODContext);
		SpatialIndexManager.createIndex(odGeography, OD.class);
	}
	
	void setUpPedObstructions() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<PedObstruction> pedObstructContext = new PedObstructionContext();
		GeographyParameters<PedObstruction> GeoParams = new GeographyParameters<PedObstruction>();
		pedObstructGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedObstructGeography", pedObstructContext, GeoParams);
		pedObstructGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testPedObstructFile = testGISDir + "boundaryPedestrianVehicleArea.shp";
		GISFunctions.readShapefile(PedObstruction.class, testPedObstructFile, pedObstructGeography, pedObstructContext);
		SpatialIndexManager.createIndex(pedObstructGeography, PedObstruction.class);
	}
	
	void setUpRoadNetwork() {
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", junctionContext, GeoParamsJunc);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);		
		
		// 2. Build road network
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, true);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		roadNetwork = builder.buildNetwork();
		GISFunctions.buildGISRoadNetwork(roadLinkGeography, junctionContext,junctionGeography, roadNetwork);
	}
	
	int getNumberIntersectingCoords(String lineDataFile) {
		try {
			setUp(lineDataFile);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get number of intersecting points and compare to expected value
		Coordinate[] odCoords = new Coordinate[pedestrianDestinationContext.getObjects(OD.class).size()];
		int i = 0;
		for (OD od : pedestrianDestinationContext.getObjects(OD.class)) {
			odCoords[i] = od.getGeom().getCoordinate();
			i++;
		}
		LineString ODLine = GISFunctions.lineStringGeometryFromCoordinates(odCoords);
		Geometry[] ODLineGeom = {ODLine};

		
		// Loop through road links in the rout and count number of times the ODLine intersects
		Geometry[] rlGeoms = new Geometry[roadLinkContext.getObjects(RoadLink.class).size()];
		i = 0;
		for (RoadLink rl: roadLinkContext.getObjects(RoadLink.class)) {
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
		
		IndexedIterable<RoadLink> lines = roadLinkContext.getObjects(RoadLink.class);
		
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
		odGeography.getAllObjects().iterator().forEachRemaining(ods::add);
		Coordinate o = ods.get(0).getGeom().getCoordinate();
		
		String roadLinkID = "osgb4000000030238946";
		
		List<Road> currentPedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Expected coords
		Coordinate c1 = new Coordinate(530507.95, 180893.35);
		Coordinate c3 = new Coordinate(530468.0087569832, 180871.8784368495);
		Coordinate c4 = new Coordinate(530482.8182132206, 180870.19519803385);
		Coordinate[] carray = {c1,null,c3,c4};
		
		List<Coordinate> expectedCoords = new ArrayList<Coordinate>(Arrays.asList(carray));
		
		
		Coordinate destCoord = null;
		for (int i=0;i<currentPedRoads.size(); i++) {
			Road r = currentPedRoads.get(i);
			destCoord = GISFunctions.farthestUnobstructedRoadCoordinate(o, r.getGeom(), pedObstructGeography);			
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
		odGeography.getAllObjects().iterator().forEachRemaining(ods::add);
		Coordinate o = ods.get(0).getGeom().getCoordinate();
		
		String roadLinkID = "osgb4000000030238946";
		
		List<Road> currentPedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Expected coords
		Coordinate c1 = new Coordinate(530506.8,180891.45);
		Coordinate c2 = new Coordinate(530512.5,180907.6);
		Coordinate c4 = new Coordinate(530521.6192518127,180903.04127475937);
		Coordinate[] carray = {c1,c2,null,c4};
		
		List<Coordinate> expectedCoords = new ArrayList<Coordinate>(Arrays.asList(carray));
		
		
		Coordinate destCoord = null;
		for (int i=0;i<currentPedRoads.size(); i++) {
			Road r = currentPedRoads.get(i);
			destCoord = GISFunctions.nearestUnobstructedRoadCoordinate(o, r.getGeom(), pedObstructGeography);
			assert expectedCoords.contains(destCoord);
		}
	}

}
