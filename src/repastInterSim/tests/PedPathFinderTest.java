package repastInterSim.tests;

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

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
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.TacticalAlternative;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.PedObstructionContext;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.pathfinding.PedPathFinder;
import repastInterSim.pathfinding.RoadNetworkRoute;

class PedPathFinderTest {
	
	Geography<Road> roadGeography = null;
	Geography<RoadLink> roadLinkGeography = null;
	Geography<OD> odGeography = null;
	Geography<PedObstruction> pedObstructGeography = null;
	Network<Junction> roadNetwork = null;
	Geography<Junction> junctionGeography = null;
	
	String testGISDir = ".//data//test_gis_data//";
	String pedestrianRoadsPath = null;
	String vehicleRoadsPath = null;
	String roadLinkPath = null;
	String serialisedLookupPath = null;
	
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
	
	List<RoadLink> planStrategicPath(Coordinate o, Coordinate d){
		
		// Plan the strategic path
		List<RoadLink> sP = new ArrayList<RoadLink>();		
		RoadNetworkRoute rnr = new RoadNetworkRoute(o , d);
		
		// Define my starting and ending junctions to test
		// Need to to this because although the origin and destination were selected above, this was just to initialise RoadNetworkRoute with
		// The route is actually calculated using junctions. 
		List<Junction> currentJunctions = new ArrayList<Junction>();
		List<Junction> destJunctions = new ArrayList<Junction>();
		for(Junction j: junctionGeography.getAllObjects()) {
			
			// Set the test current junctions 
			if (j.getFID().contentEquals("osgb4000000029970684")) {
				currentJunctions.add(j);
			}
			
			// Set the test destination junctions
			if (j.getFID().contentEquals("osgb4000000029970446")) {
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
		
		for(int i=0; i< shortestRoute.size(); i++) {
			NetworkEdge<Junction> e = (NetworkEdge<Junction>)shortestRoute.get(i);
			sP.add(e.getRoadLink());
		}
		
		return sP;
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
		
		assert alternatives.stream().filter(ta -> ta.getParityT() == 1).count()==1;
		assert alternatives.stream().filter(ta -> ta.getParityT() == 0).count()==2;
		
		assert alternatives.stream().filter(ta -> ta.getParityT() == 1).collect(Collectors.toList()).get(0).getC().equals(new Coordinate(530468.0087569832,180871.8784368495));
		assert alternatives.stream().filter(ta -> ta.getParityT() == 0).sorted( (ta1,ta2) -> ta1.getCostT().compareTo(ta2.getCostT())).collect(Collectors.toList()).get(0).getC().equals(new Coordinate(530507.95,180893.35)); // nearest
		assert alternatives.stream().filter(ta -> ta.getParityT() == 0).sorted( (ta1,ta2) -> ta1.getCostT().compareTo(ta2.getCostT())).collect(Collectors.toList()).get(1).getC().equals(new Coordinate(530482.8182132206, 180870.19519803385)); // farthest
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
		assert alternatives.stream().filter(ta -> ta.getParityT() == 1).count()==1;
		assert alternatives.stream().filter(ta -> ta.getParityT() == 0).count()==2;
		
		assert alternatives.stream().filter(ta -> ta.getParityT() == 1).collect(Collectors.toList()).get(0).getC().equals(new Coordinate(530512.5,180907.6));
		assert alternatives.stream().filter(ta -> ta.getParityT() == 0).sorted( (ta1,ta2) -> ta1.getCostT().compareTo(ta2.getCostT())).collect(Collectors.toList()).get(0).getC().equals(new Coordinate(530521.6192518127,180903.04127475937)); // nearest
		assert alternatives.stream().filter(ta -> ta.getParityT() == 0).sorted( (ta1,ta2) -> ta1.getCostT().compareTo(ta2.getCostT())).collect(Collectors.toList()).get(1).getC().equals(new Coordinate(530506.8,180891.45)); // farthest
		
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
	
	@Test
	void testChooseTacticalDestinationCoordinate1() throws Exception {
		
		// Set up GIS data
		setUpProperties();
		setUpRoads();
		setUpRoadLinks();
		setUpPedObstructions();
		setUpODs("OD_pedestrian_nodes.shp");
		setUpRoadNetwork();
		
		// Setup origin and destination coordinates
		List<OD> ods = new ArrayList<OD>();
		odGeography.getAllObjects().iterator().forEachRemaining(ods::add);
		Coordinate o = ods.stream().filter(od -> od.getId() == 1).collect(Collectors.toList()).get(0).getGeom().getCoordinate();
		Coordinate d = ods.stream().filter(od -> od.getId() == 3).collect(Collectors.toList()).get(0).getGeom().getCoordinate();
		
		List<RoadLink> sP = planStrategicPath(o,d);
		
		// Now set the planning horizon and whether secondary crossing is being performed or not and test the tactical destination coordinate returned
		int pH = 1;
		Boolean secondaryCrossing = false;
		Coordinate tacticalDestCoord = PedPathFinder.chooseTacticalDestinationCoordinate(o, d, this.roadGeography, this.pedObstructGeography, pH, sP, secondaryCrossing);
		
		// od parity = 1 and d beyond planning horizon. Therefore chooses to not cross
		assert tacticalDestCoord.equals(new Coordinate(530482.8182132206, 180870.19519803385));
		
		pH = 3;
		tacticalDestCoord = PedPathFinder.chooseTacticalDestinationCoordinate(o, d, this.roadGeography, this.pedObstructGeography, pH, sP, secondaryCrossing);
		
		// Now destination within planning horizon. Choose to cross because realise crossing is required
		assert tacticalDestCoord.equals(new Coordinate(530468.0087569832,180871.8784368495));
		
		
		// Now test a different od (origin the same but destination different)
		// This time no primary crossing required to reach destination
		o = ods.stream().filter(od -> od.getId() == 1).collect(Collectors.toList()).get(0).getGeom().getCoordinate();
		d = ods.stream().filter(od -> od.getId() == 4).collect(Collectors.toList()).get(0).getGeom().getCoordinate();
		sP = planStrategicPath(o,d);
		
		// tactical dest coord should be the same as pH=1 case above since agent doesn't have knowledge of where destination is
		pH = 1;
		tacticalDestCoord = PedPathFinder.chooseTacticalDestinationCoordinate(o, d, this.roadGeography, this.pedObstructGeography, pH, sP, secondaryCrossing);
		assert tacticalDestCoord.equals(new Coordinate(530482.8182132206, 180870.19519803385));
		
		// This time still no crossing since no primary crossing required to reach end destination
		pH=3;
		tacticalDestCoord = PedPathFinder.chooseTacticalDestinationCoordinate(o, d, this.roadGeography, this.pedObstructGeography, pH, sP, secondaryCrossing);
		assert tacticalDestCoord.equals(new Coordinate(530482.8182132206, 180870.19519803385));
	}
	
	@Test
	void testGetLinksWithinAngularDistance() throws Exception {
		
		// Load links
		setUpRoadLinks("test_strategic_path1.shp");
		
		List<RoadLink> sP = new ArrayList<RoadLink>();
		roadLinkGeography.getAllObjects().forEach(sP::add);
		//Collections.reverse(sP);
		
		Collections.sort(sP, (rl1, rl2) -> rl1.getFID().compareTo(rl2.getFID()));
		
		// Initialise list to contain the links that are within an angular dstance threshold
		List<RoadLink> linksInRange = new ArrayList<RoadLink>();
		
		linksInRange = PedPathFinder.getLinksWithinAngularDistance(sP, 0.0);
		
		assert linksInRange.size() == 1;
		
		linksInRange = PedPathFinder.getLinksWithinAngularDistance(sP, 10.0);
		
		assert linksInRange.size() == 2;
		
		linksInRange = PedPathFinder.getLinksWithinAngularDistance(sP, 50.0);
		
		assert linksInRange.size() == 5;
		
		linksInRange = PedPathFinder.getLinksWithinAngularDistance(sP, 90.0);
		
		assert linksInRange.size() == 6;
		
		linksInRange = PedPathFinder.getLinksWithinAngularDistance(sP, 250.0);
		
		assert linksInRange.size() == 7;
		
	}
	
	@Test
	void testGetLinksWithinAngularDistance2() throws Exception {
		
		// Load links
		setUpRoadLinks("test_strategic_path2.shp");
		
		List<RoadLink> sP = new ArrayList<RoadLink>();
		roadLinkGeography.getAllObjects().forEach(sP::add);
		//Collections.reverse(sP);
		
		Collections.sort(sP, (rl1, rl2) -> rl1.getFID().compareTo(rl2.getFID()));
		
		// Initialise list to contain the links that are within an angular dstance threshold
		List<RoadLink> linksInRange = new ArrayList<RoadLink>();
		
		linksInRange = PedPathFinder.getLinksWithinAngularDistance(sP, 0.0);
		
		assert linksInRange.size() == 1;
		
		linksInRange = PedPathFinder.getLinksWithinAngularDistance(sP, 90.0);
		
		assert linksInRange.size() == 3;
	}
	
	@Test
	void testGetLinksWithinAngularDistance3() throws Exception {
		// Difference between "test_strategic_path3.shp" and "test_strategic_path2.shp"
		// is that "test_strategic_path3.shp" reverses the order of coords for one of the line strings 
		// compared to the others. This tests that angle is still correctly calculated.
		
		// Load links
		setUpRoadLinks("test_strategic_path3.shp");
		
		List<RoadLink> sP = new ArrayList<RoadLink>();
		roadLinkGeography.getAllObjects().forEach(sP::add);
		//Collections.reverse(sP);
		
		Collections.sort(sP, (rl1, rl2) -> rl1.getFID().compareTo(rl2.getFID()));
		
		// Initialise list to contain the links that are within an angular dstance threshold
		List<RoadLink> linksInRange = new ArrayList<RoadLink>();
		
		linksInRange = PedPathFinder.getLinksWithinAngularDistance(sP, 0.0);
		
		assert linksInRange.size() == 1;
		
		linksInRange = PedPathFinder.getLinksWithinAngularDistance(sP, 90.0);
		
		assert linksInRange.size() == 3;
	}

}
