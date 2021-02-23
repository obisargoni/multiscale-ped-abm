package repastInterSim.tests;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
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
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.PedPathFinder;
import repastInterSim.pathfinding.RoadNetworkRoute;
import repastInterSim.pathfinding.NetworkPathFinder;
import repastInterSim.pathfinding.TacticalRoute;

class PedPathFinderTest {
	
	Context<Object> context = new DefaultContext<Object>();;
	Geography<Object> geography; 
	Geography<Road> roadGeography = null;
	Geography<OD> odGeography = null;
	Geography<CrossingAlternative> caGeography = null;
	Geography<PedObstruction> pedObstructGeography = null;
	Geography<RoadLink> roadLinkGeography = null;
	Geography<RoadLink> pavementLinkGeography = null;
	Network<Junction> roadNetwork = null;
	Context<Junction> pavementJunctionContext = null;
	Geography<Junction> pavementJunctionGeography = null;
	Network<Junction> pavementNetwork = null;
	
	
	String testGISDir = ".//data//test_gis_data//";
	String pedestrianRoadsPath = null;
	String vehicleRoadsPath = null;
	String roadLinkPath = null;
	String pavementLinkPath = null;
	String pedJPath = null;
	String serialisedLookupPath = null;
	
	void setUpProperties() throws IOException {
		IO.readProperties();
	}
	
	void setUpObjectGeography() {
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, context, geoParams);
		geography.setCRS(GlobalVars.geographyCRSString);
		context.add(geography);
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
		
		SpaceBuilder.roadGeography = roadGeography;

	}
	
	Geography<RoadLink> setUpLinks(String roadLinkFile) throws MalformedURLException, FileNotFoundException {
		roadLinkPath = testGISDir + roadLinkFile;
		
		// Initialise test road link geography and context
		Context<RoadLink> roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> rlG = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", roadLinkContext, GeoParams);
		rlG.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(RoadLink.class, roadLinkPath, rlG, roadLinkContext);
		SpatialIndexManager.createIndex(rlG, RoadLink.class);
		
		return rlG;
	}
	
	void setUpRoadLinks(String roadLinkFile) throws Exception {
		roadLinkGeography = setUpLinks(roadLinkFile);
	}
	
	void setUpRoadLinks() throws Exception {
		setUpRoadLinks("mastermap-itn RoadLink Intersect Within with orientation.shp");
	}
	
	void setUpPavementLinks(String linkFile) throws MalformedURLException, FileNotFoundException {
		pavementLinkGeography = setUpLinks(linkFile);
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
	
	void setUpCrossingAlternatives() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<CrossingAlternative> caContext = new CAContext();
		GeographyParameters<CrossingAlternative> GeoParams = new GeographyParameters<CrossingAlternative>();
		caGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("caGeography", caContext, GeoParams);
		caGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testCAFile = testGISDir + "crossing_lines.shp";
		GISFunctions.readShapefile(CrossingAlternative.class, testCAFile, caGeography, caContext);
		SpatialIndexManager.createIndex(caGeography, CrossingAlternative.class);
		
		SpaceBuilder.caGeography = caGeography;
	}
	
	void setUpRoadNetwork(boolean isDirected) {
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		Geography<Junction> junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", junctionContext, GeoParamsJunc);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);
		
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		roadNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(roadLinkGeography, junctionContext, junctionGeography, roadNetwork);
	}
	
	void setUpPavementNetwork() {
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>("PAVEMENT_NETWORK", pavementJunctionContext, false);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		pavementNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(pavementLinkGeography, pavementJunctionContext, pavementJunctionGeography, pavementNetwork);
		
		SpaceBuilder.pavementNetwork = pavementNetwork;
	}
	
	void setUpPedJunctions() throws Exception {
		setUpProperties();
		pedJPath = testGISDir + IO.getProperty("PavementJunctionsShapefile");
		
		// Initialise test road link geography and context
		pavementJunctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParams = new GeographyParameters<Junction>();
		pavementJunctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pavementJunctionGeography", pavementJunctionContext, GeoParams);
		pavementJunctionGeography.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(Junction.class, pedJPath, pavementJunctionGeography, pavementJunctionContext);
		SpatialIndexManager.createIndex(pavementJunctionGeography, Junction.class);
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
		for(Junction j: roadNetwork.getNodes()) {
			
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
	
	@Test
	void testPavementJunctions() {
		try {
			setUpPedJunctions();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get a ped junction and check its attributes are not null
		for (Junction pJ: this.pavementJunctionGeography.getAllObjects()) {
			assert pJ.getp1pID() != null;
			assert pJ.getp1rlID() != null;
			assert pJ.getp2pID() != null;
			assert pJ.getp2rlID() != null;

			assert pJ.getv1pID() != null;
			assert pJ.getv1rlID() != null;
			assert pJ.getv2pID() != null;
			assert pJ.getv2rlID() != null;
		}
	}
	
	@Test
	void testPavementNetwork() {
		
		// Call methods that create the pavement junction geography and the pavement network
		try {
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Check all the junctions but this time via the network nodes.
		// Get a ped junction and check its attributes are not null
		for (Junction pJ: this.pavementNetwork.getNodes()) {
			assert pJ.getp1pID() != null;
			assert pJ.getp1rlID() != null;
			assert pJ.getp2pID() != null;
			assert pJ.getp2rlID() != null;

			assert pJ.getv1pID() != null;
			assert pJ.getv1rlID() != null;
			assert pJ.getv2pID() != null;
			assert pJ.getv2rlID() != null;
		}
		
	}

	/*
	 * T-junction straight over
	 */
	@Test
	void testTacticalHorizonEndJunctions1() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "9745D155-3C95-4CCD-BC65-0908D57FA83A_0";
		String rlOutHorzID = "A8675945-DE94-4E22-9905-B0623A326221_0";
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : this.roadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		// Now check the nodes as as expected
		String endJID1 = tacticalEndJunctions.get(0).getFID();
		String endJID2 = tacticalEndJunctions.get(1).getFID();
		
		boolean nodeCheck = false;
		if (endJID1.contentEquals("pave_node_89") & endJID2.contentEquals("pave_node_91")) {
			nodeCheck = true;
		}
		else if (endJID1.contentEquals("pave_node_91") & endJID2.contentEquals("pave_node_89")) {
			nodeCheck = true;
		}
		
		assert nodeCheck == true;
	}
	
	/*
	 * T-junction straight ahead, other direction
	 */
	@Test
	void testTacticalHorizonEndJunctions2() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "A8675945-DE94-4E22-9905-B0623A326221_0";
		String rlOutHorzID = "9745D155-3C95-4CCD-BC65-0908D57FA83A_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : this.roadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		// Now check the nodes as as expected
		String endJID1 = tacticalEndJunctions.get(0).getFID();
		String endJID2 = tacticalEndJunctions.get(1).getFID();
		
		boolean nodeCheck = false;
		if (endJID1.contentEquals("pave_node_90") & endJID2.contentEquals("pave_node_91")) {
			nodeCheck = true;
		}
		else if (endJID1.contentEquals("pave_node_91") & endJID2.contentEquals("pave_node_90")) {
			nodeCheck = true;
		}
		
		assert nodeCheck == true;
	}
	
	/*
	 * 4 way junction straight ahead
	 */
	@Test
	void testTacticalHorizonEndJunctions3() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0";
		String rlOutHorzID = "CF9F0CB7-1387-4C83-9D25-98F63CADBE26_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : this.roadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		// Now check the nodes as as expected
		String endJID1 = tacticalEndJunctions.get(0).getFID();
		String endJID2 = tacticalEndJunctions.get(1).getFID();
		
		boolean nodeCheck = false;
		if (endJID1.contentEquals("pave_node_119") & endJID2.contentEquals("pave_node_121")) {
			nodeCheck = true;
		}
		else if (endJID1.contentEquals("pave_node_121") & endJID2.contentEquals("pave_node_119")) {
			nodeCheck = true;
		}
		
		assert nodeCheck == true;
	}
	
	/*
	 * 4 way junction left turn
	 */
	@Test
	void testTacticalHorizonEndJunctions4() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0";
		String rlOutHorzID = "1DACEAB0-2BA5-4299-8D86-F854C2FAC565_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : this.roadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		// Now check the nodes as as expected
		String endJID1 = tacticalEndJunctions.get(0).getFID();
		String endJID2 = tacticalEndJunctions.get(1).getFID();
		
		boolean nodeCheck = false;
		if (endJID1.contentEquals("pave_node_119") & endJID2.contentEquals("pave_node_121")) {
			nodeCheck = true;
		}
		else if (endJID1.contentEquals("pave_node_121") & endJID2.contentEquals("pave_node_119")) {
			nodeCheck = true;
		}
		
		assert nodeCheck == true;
	}
	
	/*
	 * Complex covent garden junction straight ahead
	 */
	@Test
	void testTacticalHorizonEndJunctions5() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "BBABD5F1-74AC-4481-91C0-61D4C85ABD77_0";
		String rlOutHorzID = "3868DA68-A5D6-4B90-9E0C-4B117146CCFD_0";
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : this.roadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		// Now check the nodes as as expected
		String endJID1 = tacticalEndJunctions.get(0).getFID();
		String endJID2 = tacticalEndJunctions.get(1).getFID();
		
		boolean nodeCheck = false;
		if (endJID1.contentEquals("pave_node_44") & endJID2.contentEquals("pave_node_45")) {
			nodeCheck = true;
		}
		else if (endJID1.contentEquals("pave_node_45") & endJID2.contentEquals("pave_node_44")) {
			nodeCheck = true;
		}
		
		assert nodeCheck == true;
	}
	
	/*
	 * Complex covent garden junction missing pavement nodes
	 */
	@Test
	void testTacticalHorizonEndJunctions6() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "3868DA68-A5D6-4B90-9E0C-4B117146CCFD_0";		
		String rlOutHorzID = "A9B5D6A4-C673-4C4F-8DC4-98FB56A72974_0";
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : this.roadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		assert tacticalEndJunctions.size() > 0;
	}
	
	
	/*
	 * 4 way junction straight ahead
	 */
	@Test
	void testTacticalHorizonOutsideJunctions1() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0";
		String rlOutHorzID = "CF9F0CB7-1387-4C83-9D25-98F63CADBE26_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : this.roadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonOutsideJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		// Now check the nodes as as expected
		String endJID1 = tacticalEndJunctions.get(0).getFID();
		String endJID2 = tacticalEndJunctions.get(1).getFID();
		
		boolean nodeCheck = false;
		if (endJID1.contentEquals("pave_node_122") & endJID2.contentEquals("pave_node_120")) {
			nodeCheck = true;
		}
		else if (endJID1.contentEquals("pave_node_120") & endJID2.contentEquals("pave_node_122")) {
			nodeCheck = true;
		}
		
		assert nodeCheck == true;
	}
	
	/*
	 * 4 way junction left turn
	 */
	@Test
	void testTacticalHorizonOutsideJunctions2() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0";
		String rlOutHorzID = "1DACEAB0-2BA5-4299-8D86-F854C2FAC565_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : this.roadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonOutsideJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		// Now check the nodes as as expected
		String endJID1 = tacticalEndJunctions.get(0).getFID();
		String endJID2 = tacticalEndJunctions.get(1).getFID();
		
		boolean nodeCheck = false;
		if (endJID1.contentEquals("pave_node_122") & endJID2.contentEquals("pave_node_121")) {
			nodeCheck = true;
		}
		else if (endJID1.contentEquals("pave_node_121") & endJID2.contentEquals("pave_node_122")) {
			nodeCheck = true;
		}
		
		assert nodeCheck == true;
	}
	
	/*
	 * Complex covent garden junction straight ahead
	 */
	@Test
	void testTacticalHorizonOutsideJunctions3() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "BBABD5F1-74AC-4481-91C0-61D4C85ABD77_0";
		String rlOutHorzID = "3868DA68-A5D6-4B90-9E0C-4B117146CCFD_0";
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : this.roadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonOutsideJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		// Now check the nodes as as expected - only one junction returned in this case due to complicated junction
		String endJID1 = tacticalEndJunctions.get(0).getFID();
		
		boolean nodeCheck = false;
		if (endJID1.contentEquals("pave_node_35")) {
			nodeCheck = true;
		}
		assert nodeCheck == true;
	}
	
	/*
	 * Test the creation of a TacticalRoute object
	 * 
	 * This method tests that the tactical routes produces with a tactical horizon of one link are expected when
	 * using the minimise crossings and minimise distance heuristics
	 */
	@Test
	public void testSetupTacticalRoute1() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
			setUpODs("OD_pedestrian_nodes.shp");
			
			setUpCrossingAlternatives();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Initialise origin-destination pair to test
		OD o = null;
		OD d = null;
		for (OD od : this.odGeography.getAllObjects()) {
			if (od.getId() == 4) {
				o = od;
			}
			else if (od.getId() == 7) {
				d = od;
			}
		}
						
		boolean minimiseCrossings = true;
		Ped pMinCross = new Ped(geography, this.roadGeography, o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
		context.add(pMinCross);
        
        minimiseCrossings = false;
        Ped pMinDist = new Ped(geography, this.roadGeography, o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
        context.add(pMinCross);
        
        // Get the strategic path - will be the same for both pedestrians
        List<RoadLink> sP = pMinCross.getPathFinder().getStrategicPath();
		int horizonNLinks = 1;
		
		// Check the start and end pavement junctions of the route
		assert pMinCross.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_85");
		assert pMinDist.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_85");
		assert pMinCross.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_112");
		assert pMinDist.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_112");
		
        TacticalRoute tr = PedPathFinder.planTacticalPath(this.pavementNetwork, this.caGeography, this.roadGeography, horizonNLinks, pMinCross, sP, pMinCross.getPathFinder().getStartPavementJunction(), pMinCross.getPathFinder().getDestPavementJunction(), pMinCross.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                

		// Now validate the tactical route
		// Check that target junction is the no crossing junction as expected
		
		// Select which end junction to find tactical path to
		final String end1ID = "pave_node_87";
        assert tr.getCurrentJunction().getFID().contentEquals(end1ID);

		// Route path of size zero now means that path only ever had one link in
		assert tr.getRoutePath().size() == 0;
		
		// Since planning horizon is one link expect the remainder path to be empty also
		assert tr.getRemainderPath().size() == 0;
		
		// Produce tactical route for min distance ped 
        tr = PedPathFinder.planTacticalPath(this.pavementNetwork, this.caGeography, this.roadGeography, horizonNLinks, pMinDist, sP, pMinDist.getPathFinder().getStartPavementJunction(), pMinDist.getPathFinder().getDestPavementJunction(), pMinDist.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                

		
		final String end2ID = "pave_node_88";		
		
		// In this case, because current edge contains a primary crossing the current junction will be the default, no cross junction
		// The accumulator target junction will be the crossing junction
		assert tr.getCurrentJunction().getFID().contentEquals(end1ID);
		assert tr.getAccumulatorRoute().getTargetJunction().getFID().contentEquals(end2ID);
		assert tr.getRoutePath().size() == 0;
		assert tr.getRemainderPath().size() == 0;		
	}
	
	/*
	 * Test the combined effect of cost heuristics and planning horizon.
	 * 
	 * With short planning horizon both cost heuristics should give the same path.
	 * 
	 * With a longer planning horizon cost heuristics should give a different paths. 
	 */
	@Test
	public void testSetupTacticalRoute2() {
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
			setUpODs("OD_pedestrian_nodes.shp");
			
			setUpCrossingAlternatives();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Initialise origin-destination pair to test
		OD o = null;
		OD d = null;
		for (OD od : this.odGeography.getAllObjects()) {
			if (od.getId() == 8) {
				o = od;
			}
			else if (od.getId() == 7) {
				d = od;
			}
		}
		
		boolean minimiseCrossings = true;
		Ped pMinCross = new Ped(geography, this.roadGeography, o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
		context.add(pMinCross);
        
        minimiseCrossings = false;
        Ped pMinDist = new Ped(geography, this.roadGeography, o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
        context.add(pMinCross);
        
        // Get the strategic path - will be the same for both pedestrians
        List<RoadLink> sP = pMinCross.getPathFinder().getStrategicPath();
		int horizonNLinks = 1;
		
		// Check the start and end pavement junctions of the route
		assert pMinCross.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_88");
		assert pMinDist.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_88");
		assert pMinCross.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_112");
		assert pMinDist.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_112");
		
        TacticalRoute trMinCross = PedPathFinder.planTacticalPath(this.pavementNetwork, this.caGeography, this.roadGeography, horizonNLinks, pMinCross, sP, pMinCross.getPathFinder().getStartPavementJunction(), pMinCross.getPathFinder().getDestPavementJunction(), pMinCross.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                
        TacticalRoute trMinDist = PedPathFinder.planTacticalPath(this.pavementNetwork, this.caGeography, this.roadGeography, horizonNLinks, pMinDist, sP, pMinDist.getPathFinder().getStartPavementJunction(), pMinDist.getPathFinder().getDestPavementJunction(), pMinDist.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                

        // With tactical planning horizon of one link, min distance and min crossing tactical paths should be the same
		final String end1ID = "pave_node_79";
        assert trMinCross.getCurrentJunction().getFID().contentEquals(end1ID);
        assert trMinDist.getCurrentJunction().getFID().contentEquals(end1ID);

		// Route path of size zero now means that path only ever had one link in.
		// Since planning horizon is one link expect the remainder path to be empty also
		assert trMinCross.getRoutePath().size() == 0;
		assert trMinDist.getRoutePath().size() == 0;
		assert trMinCross.getRemainderPath().size() == 0;
		assert trMinDist.getRemainderPath().size() == 0;
		
		// Now re-plan with planning horizon equal to the whole strategic path
		horizonNLinks = sP.size();
        trMinCross = PedPathFinder.planTacticalPath(this.pavementNetwork, this.caGeography, this.roadGeography, horizonNLinks, pMinCross, sP, pMinCross.getPathFinder().getStartPavementJunction(), pMinCross.getPathFinder().getDestPavementJunction(), pMinCross.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                
        trMinDist = PedPathFinder.planTacticalPath(this.pavementNetwork, this.caGeography, this.roadGeography, horizonNLinks, pMinDist, sP, pMinDist.getPathFinder().getStartPavementJunction(), pMinDist.getPathFinder().getDestPavementJunction(), pMinDist.getPathFinder().getPrimaryCostHeuristic(), pMinDist.getPathFinder().getSecondaryCostHeuristic());                
        
        // Initial section of the path should be the same for both peds
        // Need to starting junction using the accumulator route since first link is a crossing link
		final String end2ID = "pave_node_81";
        assert trMinCross.getAccumulatorRoute().getTargetJunction().getFID().contentEquals(end2ID);
        assert trMinDist.getAccumulatorRoute().getTargetJunction().getFID().contentEquals(end2ID);
        
        // But rest of path will differ
        String [] expectedNodesMinCross = {end2ID, "pave_node_89", "pave_node_91", pMinCross.getPathFinder().getDestPavementJunction().getFID()};
        String [] expectedNodesMinDist = {end2ID, "pave_node_90", pMinCross.getPathFinder().getDestPavementJunction().getFID()};
        List<Junction> remainderPathNodesMinCross = trMinCross.getNetworkPathFinder().nodePathFromEdges(trMinCross.getRemainderPath(), trMinCross.getAccumulatorRoute().getTargetJunction());
        List<Junction> remainderPathNodesMinDist = trMinDist.getNetworkPathFinder().nodePathFromEdges(trMinDist.getRemainderPath(), trMinDist.getAccumulatorRoute().getTargetJunction());
		
		for (int i=0; i<Math.max(expectedNodesMinCross.length, remainderPathNodesMinCross.size()); i++) {
			assert remainderPathNodesMinCross.get(i).getFID().contentEquals(expectedNodesMinCross[i]);
		}
		
		for (int i=0; i<expectedNodesMinDist.length; i++) {
			assert remainderPathNodesMinDist.get(i).getFID().contentEquals(expectedNodesMinDist[i]);
		}
	}
	
	/*
	 * Testing the initialisation of a PedPathFinder object and initial update of tactical path via path finder. 
	 * 
	 * O Id = 4 D id = 1.
	 */
	@Test
	public void testPedPathFinder1() {
		
		// Setup environment
		try {
			setUpObjectGeography();

			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
			setUpCrossingAlternatives();
			
			setUpODs("OD_pedestrian_nodes.shp");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : this.odGeography.getAllObjects()) {
			if (i.getId() == 4) {
				o = i;
			}
			else if (i.getId() == 1) {
				d = i;
			}
		}
				
		// Initialise a pedestrian, this internally initialises a ped path finder
		boolean minimiseCrossings = false;
		Ped pedMinDist = new Ped(geography, this.roadGeography, o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
	 		
		// Check the start and end pavement junctions are as expected
		assert pedMinDist.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_85");
		assert pedMinDist.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_93");
		
		// Check that current junction is initially null
		assert pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction() == null;
				
		// Update tactical path and check the expected path is produced
		pedMinDist.getPathFinder().updateTacticalPath();
		
		assert pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_87");
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getTargetJunction().getFID().contentEquals("pave_node_88");
		
		// Repeat test for ped that minimises crossings
		minimiseCrossings = true;
		Ped pedMinCross = new Ped(geography, this.roadGeography, o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
		
		// Check the start and end pavement junctions are as expected
		assert pedMinCross.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_85");
		assert pedMinCross.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_93");
		
		// Check that current junction is initially null
		assert pedMinCross.getPathFinder().getTacticalPath().getCurrentJunction() == null;
				
		// Update tactical path and check the expected path is produced
		pedMinCross.getPathFinder().updateTacticalPath();
		
		assert pedMinCross.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_87");
		assert pedMinCross.getPathFinder().getTacticalPath().getAccumulatorRoute().getTargetJunction() == null;
	}
	
	/*
	 * Testing the initialisation of a PedPathFinder object. O Id = 2 D id = 1.
	 * 
	 * Test fails due to issues with strategic path finding model.
	 */
	@Test
	public void testPedPathFinder2() {
		
		// Setup environment
		try {
			setUpObjectGeography();

			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
			setUpCrossingAlternatives();
			
			setUpODs("OD_pedestrian_nodes.shp");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : this.odGeography.getAllObjects()) {
			if (i.getId() == 2) {
				o = i;
			}
			else if (i.getId() == 1) {
				d = i;
			}
		}
		
		
		// Set up ped path finder
		PedPathFinder ppf = new PedPathFinder(o, d, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
		
		// Check the strategic path is as expected
		String[] expectedRoadIDs = {"762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0", "56CF7BBA-28E4-4ACA-9F58-E096E88094FB_0", "B2B9D137-2BA4-4864-8350-2EDAA5910747_0"};
		for (int i=0;i<ppf.getStrategicPath().size(); i++) {
			assert ppf.getStrategicPath().get(i).getFID().contentEquals(expectedRoadIDs[i]);
		}
		
		// Check the start and end pavement junctions are as expected
		assert ppf.getStartPavementJunction().getFID().contentEquals("pave_node_121");
		assert ppf.getDestPavementJunction().getFID().contentEquals("pave_node_80");
		
		Ped p = new Ped(geography, this.roadGeography, o, d, 0.5, 1.0, 0.9, 3.0, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
		
		// Now test planning the first tactical path with this ped path finder object
		ppf.planTacticaAccumulatorPath(this.pavementNetwork, this.caGeography, this.roadGeography, p, ppf.getStrategicPath(), ppf.getStartPavementJunction(), ppf.getDestPavementJunction());
		
		// Check the current (default) and target tactical alternatives are as expected
		assert ppf.getTacticalPath().getCurrentTA().getEndJunction().getFID().contentEquals("pave_node_114");
		assert ppf.getTacticalPath().getTargetTA().getEndJunction().getFID().contentEquals("pave_node_114");
	}
	
	/*
	 * Testing the initialisation of a PedPathFinder object. O Id = 2 D id = 5.
	 * 
	 * Also tests updating the tactical path from a junction part way along the journey.
	 * 
	 */
	@Test
	public void testPedPathFinder3() {
		
		// Setup environment
		try {
			setUpObjectGeography();

			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
			setUpCrossingAlternatives();
			
			setUpODs("OD_pedestrian_nodes.shp");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : this.odGeography.getAllObjects()) {
			if (i.getId() == 2) {
				o = i;
			}
			else if (i.getId() == 5) {
				d = i;
			}
		}
		
		
		// Set up ped path finder
		PedPathFinder ppf = new PedPathFinder(o, d, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
		
		// Check the strategic path is as expected
		String[] expectedRoadIDs = {"762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0", "A8675945-DE94-4E22-9905-B0623A326221_0", "F4C0B1FB-762C-4492-BB0D-673CC4950CBE_0", "8A9E2D7B-3B48-4A19-B89A-0B4F4D516870_2"};
		for (int i=0;i<ppf.getStrategicPath().size(); i++) {
			assert ppf.getStrategicPath().get(i).getFID().contentEquals(expectedRoadIDs[i]);
		}
		
		// Check the start and end pavement junctions are as expected
		assert ppf.getStartPavementJunction().getFID().contentEquals("pave_node_121");
		assert ppf.getDestPavementJunction().getFID().contentEquals("pave_node_87");
		
		Ped p = new Ped(geography, this.roadGeography, o, d, 0.5, 1.0, 0.9, 3.0, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
		
		// Now test planning the first tactical path with this ped path finder object
		ppf.planTacticaAccumulatorPath(this.pavementNetwork, this.caGeography, this.roadGeography, p, ppf.getStrategicPath(), ppf.getStartPavementJunction(), ppf.getDestPavementJunction());
		
		// Check the current (default) and target tactical alternatives are as expected
		assert ppf.getTacticalPath().getCurrentTA().getEndJunction().getFID().contentEquals("pave_node_91");
		assert ppf.getTacticalPath().getTargetTA().getEndJunction().getFID().contentEquals("pave_node_91");
		
		// Test planning the second tactical path
		String[] expectedRouteJunctions = {"pave_node_114", "pave_node_112", "pave_node_91", "pave_node_89"};
		for (int i=0;i<expectedRouteJunctions.length;i++) {
			assert ppf.getTacticalPath().getTargetTA().getRouteJunctions().get(i).getFID().contentEquals(expectedRouteJunctions[i]);
		}
		
		Junction secondStart = ppf.getTacticalPath().getTargetTA().getOutsideJunction();
		List<RoadLink> updatedSP = ppf.getStrategicPath().subList(2, ppf.getStrategicPath().size());
		ppf.planTacticaAccumulatorPath(this.pavementNetwork, this.caGeography, this.roadGeography, p, updatedSP, secondStart, ppf.getDestPavementJunction());
		assert ppf.getTacticalPath().getCurrentTA().getEndJunction().getFID().contentEquals("pave_node_81");
		
	}
	
	/*
	 * Testing the initialisation of a PedPathFinder object. O Id = 3 D id = 5.
	 * 
	 * Tests updating the tactical path until the end of the journey.
	 * 
	 */
	@Test
	public void testPedPathFinder4() {
		
		// Setup environment
		try {
			setUpObjectGeography();

			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
			setUpCrossingAlternatives();
			
			setUpODs("OD_pedestrian_nodes.shp");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : this.odGeography.getAllObjects()) {
			if (i.getId() == 3) {
				o = i;
			}
			else if (i.getId() == 4) {
				d = i;
			}
		}
		
		
		// Set up ped path finder
		PedPathFinder ppf = new PedPathFinder(o, d, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
		
		// Check the strategic path is as expected
		String[] expectedRoadIDs = {"9745D155-3C95-4CCD-BC65-0908D57FA83A_0", "F4C0B1FB-762C-4492-BB0D-673CC4950CBE_0", "8A9E2D7B-3B48-4A19-B89A-0B4F4D516870_2", "8A9E2D7B-3B48-4A19-B89A-0B4F4D516870_1"};
		for (int i=0;i<ppf.getStrategicPath().size(); i++) {
			assert ppf.getStrategicPath().get(i).getFID().contentEquals(expectedRoadIDs[i]);
		}
		
		// Check the start and end pavement junctions are as expected
		assert ppf.getStartPavementJunction().getFID().contentEquals("pave_node_108");
		assert ppf.getDestPavementJunction().getFID().contentEquals("pave_node_85");
		
		Ped p = new Ped(geography, this.roadGeography, o, d, 0.5, 1.0, 0.9, 3.0, this.roadLinkGeography, this.roadNetwork, this.odGeography, this.pavementJunctionGeography, this.pavementNetwork);
		
		// Now test planning the first tactical path with this ped path finder object
		ppf.planTacticaAccumulatorPath(this.pavementNetwork, this.caGeography, this.roadGeography, p, ppf.getStrategicPath(), ppf.getStartPavementJunction(), ppf.getDestPavementJunction());
		
		// Check the current (default) and target tactical alternatives are as expected
		assert ppf.getTacticalPath().getCurrentTA().getEndJunction().getFID().contentEquals("pave_node_91");
		assert ppf.getTacticalPath().getTargetTA().getEndJunction().getFID().contentEquals("pave_node_91");
		
		// Test planning the second tactical path
		String[] expectedRouteJunctions = {"pave_node_91", "pave_node_90"};
		for (int i=0;i<expectedRouteJunctions.length;i++) {
			assert ppf.getTacticalPath().getTargetTA().getRouteJunctions().get(i).getFID().contentEquals(expectedRouteJunctions[i]);
		}
		
		Junction updatedStart = ppf.getTacticalPath().getTargetTA().getOutsideJunction();
		List<RoadLink> updatedSP = ppf.getStrategicPath().subList(1, ppf.getStrategicPath().size());
		ppf.planTacticaAccumulatorPath(this.pavementNetwork, this.caGeography, this.roadGeography, p, updatedSP, updatedStart, ppf.getDestPavementJunction());
		assert ppf.getTacticalPath().getCurrentTA().getEndJunction().getFID().contentEquals("pave_node_82");
		assert ppf.getTacticalPath().getCurrentTA().getOutsideJunction().getFID().contentEquals("pave_node_79");
		
		// Update again
		updatedStart = ppf.getTacticalPath().getTargetTA().getOutsideJunction();
		updatedSP = updatedSP.subList(1, updatedSP.size());
		ppf.planTacticaAccumulatorPath(this.pavementNetwork, this.caGeography, this.roadGeography, p, updatedSP, updatedStart, ppf.getDestPavementJunction());
		assert ppf.getTacticalPath().getCurrentTA().getEndJunction().getFID().contentEquals("pave_node_88");
		assert ppf.getTacticalPath().getCurrentTA().getOutsideJunction().getFID().contentEquals("pave_node_88");
		
		// Update a final time
		// Should get difference between target and current TA because now on final link
		updatedStart = ppf.getTacticalPath().getTargetTA().getOutsideJunction();
		updatedSP = updatedSP.subList(1, updatedSP.size());
		ppf.planTacticaAccumulatorPath(this.pavementNetwork, this.caGeography, this.roadGeography, p, updatedSP, updatedStart, ppf.getDestPavementJunction());
		assert ppf.getTacticalPath().getCurrentTA().getEndJunction().getFID().contentEquals("pave_node_83");
		assert ppf.getTacticalPath().getTargetTA().getEndJunction().getFID().contentEquals("pave_node_85");
		
		assert ppf.getTacticalPath().getCurrentTA().getOutsideJunction() == null;
		assert ppf.getTacticalPath().getTargetTA().getOutsideJunction() == null;
		
		// Finally test updating the targetTA as if it were chosen
		
		// First check the current junction
		assert ppf.getTacticalPath().getCurrentTA().getCurrentJunction().getFID().contentEquals("pave_node_83");
		
		// Then update the target TA
		ppf.getTacticalPath().getTargetTA().updatePathToEnd(ppf.getTacticalPath().getCurrentTA().getCurrentJunction());
		assert ppf.getTacticalPath().getTargetTA().getCurrentJunction().getFID().contentEquals("pave_node_85");
		
	}
}
