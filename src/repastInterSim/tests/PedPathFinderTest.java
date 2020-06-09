package repastInterSim.tests;

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.context.Context;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.OD;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.TacticalAlternative;
import repastInterSim.environment.contexts.PedObstructionContext;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.pathfinding.PedPathFinder;
import repastInterSim.pathfinding.RoadNetworkRoute;

class PedPathFinderTest {
	
	Geography<Road> roadGeography = null;
	Geography<RoadLink> roadLinkGeography = null;
	Geography<OD> odGeography = null;
	Geography<PedObstruction> pedObstructGeography = null;
	
	String testGISDir = ".//data//test_gis_data//";
	String pedestrianRoadsPath = null;
	String vehicleRoadsPath = null;
	String roadLinkPath = null;
	String serialisedLookupPath = null;

	void setUpRoads() throws Exception {
		pedestrianRoadsPath = testGISDir + "topographicAreaVehicle.shp";
		vehicleRoadsPath = testGISDir + "topographicAreaPedestrian.shp";
		serialisedLookupPath = testGISDir + "road_link_rodas_cache.serialised";
		
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
		
		roadLinkPath = testGISDir + "mastermap-itn RoadLink Intersect Within with orientation.shp";
		
		// Initialise test road link geography and context
		Context<RoadLink> roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(RoadLink.class, roadLinkPath, roadLinkGeography, roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
	}
	
	void setUpODs() throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> ODContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		odGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", ODContext, GeoParamsOD);
		odGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String testODFile = testGISDir + "test_ped_OD1.shp";
		GISFunctions.readShapefile(OD.class, testODFile, odGeography, ODContext);
		SpatialIndexManager.createIndex(odGeography, OD.class);
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
			destCoord = PedPathFinder.farthestUnobstructedRoadCoordinate(o, r.getGeom(), pedObstructGeography);			
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
			destCoord = PedPathFinder.nearestUnobstructedRoadCoordinate(o, r.getGeom(), pedObstructGeography);
			assert expectedCoords.contains(destCoord);
		}
	}

	@Test
	void testGetTacticalDestinationCoodinateOptionsFar() throws Exception {
		setUpRoads();
		setUpRoadLinks();
		setUpODs("test_ped_OD1.shp");
		setUpPedObstructions();
		
		File vehcileRoadsFile = new File(vehicleRoadsPath);
		File pedestrianRoadsFile = new File(pedestrianRoadsPath);
		File serialisedLoc = new File(serialisedLookupPath);
		
		// Select pedestrian origins and destinations to test
		List<OD> ods = new ArrayList<OD>();
		odGeography.getAllObjects().iterator().forEachRemaining(ods::add);
		Coordinate o = ods.get(0).getGeom().getCoordinate();
		
		// Select section of road link network we are testing
		String roadLinkID = "osgb4000000030238946";
		List<RoadLink> rls = new ArrayList<RoadLink>();
		for (RoadLink rl: roadLinkGeography.getAllObjects()) {
			if (rl.getFID().contentEquals(roadLinkID)) {
				rls.add(rl);
			}
		}
		
		// Select pedestrian roads used in testing
		List<Road> currentPedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Get the tactical desinations for the origin coord and this road link
		HashMap<String, List<Coordinate>> destOptions = PedPathFinder.getTacticalDestinationCoodinateOptions(o, currentPedRoads, rls, pedObstructGeography, true);
		
		// Check that two keys and two coordinates returned
		Set<String> expectedKeys = new HashSet<String>();
		expectedKeys.add("cross");
		expectedKeys.add("nocross");
		assert destOptions.keySet().containsAll(expectedKeys);
		
		// Check the coorrdinates are as expected
		assert destOptions.get("cross").size() == 1;
		assert destOptions.get("nocross").size() == 2;
		
		assert destOptions.get("cross").get(0).equals(new Coordinate(530468.0087569832,180871.8784368495));
		
		assert destOptions.get("nocross").get(0).equals(new Coordinate(530507.95,180893.35)); // nearest
		assert destOptions.get("nocross").get(1).equals(new Coordinate(530482.8182132206, 180870.19519803385)); // farthest
	}
	
	@Test
	void testGetTacticalDestinationAlternativesFar() throws Exception {
		setUpRoads();
		setUpRoadLinks();
		setUpODs("test_ped_OD1.shp");
		setUpPedObstructions();
		
		File vehcileRoadsFile = new File(vehicleRoadsPath);
		File pedestrianRoadsFile = new File(pedestrianRoadsPath);
		File serialisedLoc = new File(serialisedLookupPath);
		
		// Select pedestrian origins and destinations to test
		List<OD> ods = new ArrayList<OD>();
		odGeography.getAllObjects().iterator().forEachRemaining(ods::add);
		Coordinate o = ods.get(0).getGeom().getCoordinate();
		
		Coordinate d = new Coordinate();
		
		// Select section of road link network we are testing
		String roadLinkID = "osgb4000000030238946";
		List<RoadLink> rls = new ArrayList<RoadLink>();
		for (RoadLink rl: roadLinkGeography.getAllObjects()) {
			if (rl.getFID().contentEquals(roadLinkID)) {
				rls.add(rl);
			}
		}
		
		// Select pedestrian roads used in testing
		List<Road> currentPedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Get the tactical desinations for the origin coord and this road link
		int pH = 3;
		ArrayList<TacticalAlternative> alternatives = PedPathFinder.getTacticalDestinationAlternatives(o, currentPedRoads, rls, d, pH, pedObstructGeography, true);
		
		assert alternatives.stream().filter(ta -> ta.parityT == 1).count()==1;
		assert alternatives.stream().filter(ta -> ta.parityT == 0).count()==2;
		
		assert alternatives.stream().filter(ta -> ta.parityT == 1).collect(Collectors.toList()).get(0).c.equals(new Coordinate(530468.0087569832,180871.8784368495));
		assert alternatives.stream().filter(ta -> ta.parityT == 0).sorted( (ta1,ta2) -> ta1.costT.compareTo(ta2.costT)).collect(Collectors.toList()).get(0).c.equals(new Coordinate(530507.95,180893.35)); // nearest
		assert alternatives.stream().filter(ta -> ta.parityT == 0).sorted( (ta1,ta2) -> ta1.costT.compareTo(ta2.costT)).collect(Collectors.toList()).get(1).c.equals(new Coordinate(530482.8182132206, 180870.19519803385)); // farthest
	}
	
	@Test
	void testGetTacticalDestinationCoodinateOptionsNear() throws Exception {
		setUpRoads();
		setUpRoadLinks();
		setUpODs("test_ped_OD1.shp");
		setUpPedObstructions();
		
		File vehcileRoadsFile = new File(vehicleRoadsPath);
		File pedestrianRoadsFile = new File(pedestrianRoadsPath);
		File serialisedLoc = new File(serialisedLookupPath);
		
		// Select pedestrian origins and destinations to test
		List<OD> ods = new ArrayList<OD>();
		odGeography.getAllObjects().iterator().forEachRemaining(ods::add);
		Coordinate o = ods.get(0).getGeom().getCoordinate();
		
		// Select section of road link network we are testing
		String roadLinkID = "osgb4000000030238946";
		List<RoadLink> rls = new ArrayList<RoadLink>();
		for (RoadLink rl: roadLinkGeography.getAllObjects()) {
			if (rl.getFID().contentEquals(roadLinkID)) {
				rls.add(rl);
			}
		}
		
		// Select pedestrian roads used in testing
		List<Road> currentPedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Get the tactical desinations for the origin coord and this road link
		HashMap<String, List<Coordinate>> destOptions = PedPathFinder.getTacticalDestinationCoodinateOptions(o, currentPedRoads, rls, pedObstructGeography, false);
		
		// Check that two keys and two coordinates returned
		Set<String> expectedKeys = new HashSet<String>();
		expectedKeys.add("cross");
		expectedKeys.add("nocross");
		assert destOptions.keySet().containsAll(expectedKeys);
		
		// Check the coorrdinates are as expected
		assert destOptions.get("cross").size() == 1;
		assert destOptions.get("nocross").size() == 2;
		
		assert destOptions.get("cross").get(0).equals(new Coordinate(530512.5,180907.6));
		
		assert destOptions.get("nocross").get(0).equals(new Coordinate(530521.6192518127,180903.04127475937)); // nearest
		assert destOptions.get("nocross").get(1).equals(new Coordinate(530506.8,180891.45)); // farthest
	}
	
	@Test
	void testGetTacticalDestinationAlternativesNear() throws Exception {
		setUpRoads();
		setUpRoadLinks();
		setUpODs("test_ped_OD1.shp");
		setUpPedObstructions();
		
		File vehcileRoadsFile = new File(vehicleRoadsPath);
		File pedestrianRoadsFile = new File(pedestrianRoadsPath);
		File serialisedLoc = new File(serialisedLookupPath);
		
		// Select pedestrian origins and destinations to test
		List<OD> ods = new ArrayList<OD>();
		odGeography.getAllObjects().iterator().forEachRemaining(ods::add);
		Coordinate o = ods.get(0).getGeom().getCoordinate();
		
		Coordinate d = new Coordinate();
		
		// Select section of road link network we are testing
		String roadLinkID = "osgb4000000030238946";
		List<RoadLink> rls = new ArrayList<RoadLink>();
		for (RoadLink rl: roadLinkGeography.getAllObjects()) {
			if (rl.getFID().contentEquals(roadLinkID)) {
				rls.add(rl);
			}
		}
		
		// Select pedestrian roads used in testing
		List<Road> currentPedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc, roadLinkID);
		
		// Get the tactical desinations for the origin coord and this road link
		int pH = 3;
		ArrayList<TacticalAlternative> alternatives = PedPathFinder.getTacticalDestinationAlternatives(o, currentPedRoads, rls, d, pH, pedObstructGeography, false);
		
		// Check the coorrdinates are as expected
		assert alternatives.stream().filter(ta -> ta.parityT == 1).count()==1;
		assert alternatives.stream().filter(ta -> ta.parityT == 0).count()==2;
		
		assert alternatives.stream().filter(ta -> ta.parityT == 1).collect(Collectors.toList()).get(0).c.equals(new Coordinate(530512.5,180907.6));
		assert alternatives.stream().filter(ta -> ta.parityT == 0).sorted( (ta1,ta2) -> ta1.costT.compareTo(ta2.costT)).collect(Collectors.toList()).get(0).c.equals(new Coordinate(530521.6192518127,180903.04127475937)); // nearest
		assert alternatives.stream().filter(ta -> ta.parityT == 0).sorted( (ta1,ta2) -> ta1.costT.compareTo(ta2.costT)).collect(Collectors.toList()).get(1).c.equals(new Coordinate(530506.8,180891.45)); // farthest
		
	}
	
	@Test
	void testCheckCoordinateIntersectingRoads() throws Exception {
		setUpRoads();
		
		Coordinate c = new Coordinate(530506.8,180891.45);
		String currentRoadLinkID = "osgb4000000030238946";
		String nextRoadLinkID = "osgb4000000030238838";
		
		Coordinate midpoly = new Coordinate(530496.982,180882.880);
		List<Road> inters = SpatialIndexManager.findIntersectingObjects(this.roadGeography, midpoly);
		
		String coordType = PedPathFinder.checkCoordinateIntersectingRoads(c, this.roadGeography, currentRoadLinkID, nextRoadLinkID);
		assert coordType.contentEquals("not_at_end");
		
		c = new Coordinate(530468.0087569832,180871.8784368495);
		
		coordType = PedPathFinder.checkCoordinateIntersectingRoads(c, this.roadGeography, currentRoadLinkID, nextRoadLinkID);
		assert coordType.contentEquals("intersects_next");
		
		c = new Coordinate(530482.8182132206, 180870.19519803385);
		
		coordType = PedPathFinder.checkCoordinateIntersectingRoads(c, this.roadGeography, currentRoadLinkID, nextRoadLinkID);
		assert coordType.contentEquals("not_intersects_next");
		
		// Check for null road link entries
		coordType = PedPathFinder.checkCoordinateIntersectingRoads(c, this.roadGeography, "", "");
		assert coordType.contentEquals("not_intersects_next");
	}

}
