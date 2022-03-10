package repastInterSim.tests;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.engine.environment.RunState;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;
import repastInterSim.environment.OD;
import repastInterSim.environment.RoadLink;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.PedPathFinder;
import repastInterSim.pathfinding.TacticalRoute;

class PedPathFinderTest {
	
	Context<Object> context = new DefaultContext<Object>();
	
	String testGISDir = ".//data//test_gis_data//";
	String pedestrianRoadsPath = null;
	String vehicleRoadsPath = null;
	String roadLinkPath = null;
	String pavementLinkPath = null;
	String pedJPath = null;
	String serialisedLookupPath = null;
	

	
	public void addPedToWorld(Ped p, OD o) {
		Context<Object> c = RunState.getInstance().getMasterContext();
		Geography<Object> g = (Geography<Object>) c.getProjection(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY);
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(p.getRad());
		
		GISFunctions.moveAgentToGeometry(g, circle, p);
		p.setLoc();
	}
	

	
	@Test
	void testGetLinksWithinAngularDistance() throws Exception {
		
		EnvironmentSetup.setUpProperties();
		
		// Load links
		EnvironmentSetup.orRoadLinkGeography = EnvironmentSetup.setUpRoadLinks("test_strategic_path1.shp", GlobalVars.CONTEXT_NAMES.ROAD_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.ROAD_LINK_GEOGRAPHY);
		
		List<RoadLink> sP = new ArrayList<RoadLink>();
		EnvironmentSetup.orRoadLinkGeography.getAllObjects().forEach(sP::add);
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
		
		EnvironmentSetup.setUpProperties();
		
		// Load links
		EnvironmentSetup.orRoadLinkGeography = EnvironmentSetup.setUpRoadLinks("test_strategic_path2.shp", GlobalVars.CONTEXT_NAMES.ROAD_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.ROAD_LINK_GEOGRAPHY);
		
		List<RoadLink> sP = new ArrayList<RoadLink>();
		EnvironmentSetup.orRoadLinkGeography.getAllObjects().forEach(sP::add);
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
		EnvironmentSetup.setUpProperties();
		
		// Difference between "test_strategic_path3.shp" and "test_strategic_path2.shp"
		// is that "test_strategic_path3.shp" reverses the order of coords for one of the line strings 
		// compared to the others. This tests that angle is still correctly calculated.
		
		// Load links
		EnvironmentSetup.orRoadLinkGeography = EnvironmentSetup.setUpRoadLinks("test_strategic_path3.shp", GlobalVars.CONTEXT_NAMES.ROAD_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.ROAD_LINK_GEOGRAPHY);
		
		List<RoadLink> sP = new ArrayList<RoadLink>();
		EnvironmentSetup.orRoadLinkGeography.getAllObjects().forEach(sP::add);
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpPedJunctions();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get a ped junction and check its attributes are not null
		for (Junction pJ: EnvironmentSetup.pavementJunctionGeography.getAllObjects()) {
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Check all the junctions but this time via the network nodes.
		// Get a ped junction and check its attributes are not null
		for (Junction pJ: EnvironmentSetup.pavementNetwork.getNodes()) {
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "9745D155-3C95-4CCD-BC65-0908D57FA83A_0";
		String rlOutHorzID = "A8675945-DE94-4E22-9905-B0623A326221_0";
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(EnvironmentSetup.pavementNetwork, rlEndHorz, rlOutHorz);
		
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "A8675945-DE94-4E22-9905-B0623A326221_0";
		String rlOutHorzID = "9745D155-3C95-4CCD-BC65-0908D57FA83A_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(EnvironmentSetup.pavementNetwork, rlEndHorz, rlOutHorz);
		
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0";
		String rlOutHorzID = "CF9F0CB7-1387-4C83-9D25-98F63CADBE26_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(EnvironmentSetup.pavementNetwork, rlEndHorz, rlOutHorz);
		
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0";
		String rlOutHorzID = "1DACEAB0-2BA5-4299-8D86-F854C2FAC565_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(EnvironmentSetup.pavementNetwork, rlEndHorz, rlOutHorz);
		
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "BBABD5F1-74AC-4481-91C0-61D4C85ABD77_0";
		String rlOutHorzID = "3868DA68-A5D6-4B90-9E0C-4B117146CCFD_0";
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(EnvironmentSetup.pavementNetwork, rlEndHorz, rlOutHorz);
		
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "3868DA68-A5D6-4B90-9E0C-4B117146CCFD_0";		
		String rlOutHorzID = "A9B5D6A4-C673-4C4F-8DC4-98FB56A72974_0";
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonEndJunctions(EnvironmentSetup.pavementNetwork, rlEndHorz, rlOutHorz);
		
		assert tacticalEndJunctions.size() > 0;
	}
	
	
	/*
	 * 4 way junction straight ahead
	 */
	@Test
	void testTacticalHorizonOutsideJunctions1() {
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0";
		String rlOutHorzID = "CF9F0CB7-1387-4C83-9D25-98F63CADBE26_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonOutsideJunctions(EnvironmentSetup.pavementNetwork, rlEndHorz, rlOutHorz);
		
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0";
		String rlOutHorzID = "1DACEAB0-2BA5-4299-8D86-F854C2FAC565_0";		
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonOutsideJunctions(EnvironmentSetup.pavementNetwork, rlEndHorz, rlOutHorz);
		
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Manually set the two road links to get pavement junctions between
		String rlEndHorzID = "BBABD5F1-74AC-4481-91C0-61D4C85ABD77_0";
		String rlOutHorzID = "3868DA68-A5D6-4B90-9E0C-4B117146CCFD_0";
		
		RoadLink rlEndHorz = null;
		RoadLink rlOutHorz = null;
		for (RoadLink rl : EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(rlEndHorzID)) {
				rlEndHorz = rl;
				continue;
			}
			
			if (rl.getPedRLID().contentEquals(rlOutHorzID)) {
				rlOutHorz = rl;
				continue;
			}
		}
		
		List<Junction> tacticalEndJunctions = PedPathFinder.tacticalHorizonOutsideJunctions(EnvironmentSetup.pavementNetwork, rlEndHorz, rlOutHorz);
		
		// Now check the nodes as as expected - only one junction returned in this case due to complicated junction
		String endJID1 = tacticalEndJunctions.get(0).getFID();
		
		boolean nodeCheck = false;
		if (endJID1.contentEquals("pave_node_44")) {
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
			EnvironmentSetup.setUpProperties();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Initialise origin-destination pair to test
		OD o = null;
		OD d = null;
		for (OD od : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (od.getId() == 5) {
				o = od;
			}
			else if (od.getId() == 2) {
				d = od;
			}
		}
						
		boolean minimiseCrossings = true;
		Ped pMinCross = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
		EnvironmentSetup.context.add(pMinCross);
        
        minimiseCrossings = false;
        Ped pMinDist = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
        EnvironmentSetup.context.add(pMinCross);
        
        // Get the strategic path - will be the same for both pedestrians
        List<RoadLink> sP = pMinCross.getPathFinder().getStrategicPath();
		int horizonNLinks = 1;
		
		// Check the start and end pavement junctions of the route
		assert pMinCross.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_85");
		assert pMinDist.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_85");
		assert pMinCross.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_112");
		assert pMinDist.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_112");
		
        TacticalRoute tr = pMinCross.getPathFinder().planTacticalPath(EnvironmentSetup.pavementNetwork, EnvironmentSetup.caGeography, EnvironmentSetup.roadGeography, horizonNLinks, pMinCross, sP, pMinCross.getPathFinder().getStartPavementJunction(), pMinCross.getPathFinder().getDestPavementJunction(), pMinCross.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                

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
        tr = pMinDist.getPathFinder().planTacticalPath(EnvironmentSetup.pavementNetwork, EnvironmentSetup.caGeography, EnvironmentSetup.roadGeography, horizonNLinks, pMinDist, sP, pMinDist.getPathFinder().getStartPavementJunction(), pMinDist.getPathFinder().getDestPavementJunction(), pMinDist.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                

		
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
			EnvironmentSetup.setUpProperties();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Initialise origin-destination pair to test
		OD o = null;
		OD d = null;
		for (OD od : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (od.getId() == 6) {
				o = od;
			}
			else if (od.getId() == 2) {
				d = od;
			}
		}
		
		boolean minimiseCrossings = true;
		Ped pMinCross = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
		EnvironmentSetup.context.add(pMinCross);
        
        minimiseCrossings = false;
        Ped pMinDist = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
        EnvironmentSetup.context.add(pMinCross);
        
        // Get the strategic path - will be the same for both pedestrians
        List<RoadLink> sP = pMinCross.getPathFinder().getStrategicPath();
		int horizonNLinks = 1;
		
		// Check the start and end pavement junctions of the route
		assert pMinCross.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_88");
		assert pMinDist.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_88");
		assert pMinCross.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_112");
		assert pMinDist.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_112");
		
        TacticalRoute trMinCross = pMinCross.getPathFinder().planTacticalPath(EnvironmentSetup.pavementNetwork, EnvironmentSetup.caGeography, EnvironmentSetup.roadGeography, horizonNLinks, pMinCross, sP, pMinCross.getPathFinder().getStartPavementJunction(), pMinCross.getPathFinder().getDestPavementJunction(), pMinCross.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                
        TacticalRoute trMinDist = pMinDist.getPathFinder().planTacticalPath(EnvironmentSetup.pavementNetwork, EnvironmentSetup.caGeography, EnvironmentSetup.roadGeography, horizonNLinks, pMinDist, sP, pMinDist.getPathFinder().getStartPavementJunction(), pMinDist.getPathFinder().getDestPavementJunction(), pMinDist.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                

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
        trMinCross = pMinCross.getPathFinder().planTacticalPath(EnvironmentSetup.pavementNetwork, EnvironmentSetup.caGeography, EnvironmentSetup.roadGeography, horizonNLinks, pMinCross, sP, pMinCross.getPathFinder().getStartPavementJunction(), pMinCross.getPathFinder().getDestPavementJunction(), pMinCross.getPathFinder().getPrimaryCostHeuristic(), pMinCross.getPathFinder().getSecondaryCostHeuristic());                
        trMinDist = pMinDist.getPathFinder().planTacticalPath(EnvironmentSetup.pavementNetwork, EnvironmentSetup.caGeography, EnvironmentSetup.roadGeography, horizonNLinks, pMinDist, sP, pMinDist.getPathFinder().getStartPavementJunction(), pMinDist.getPathFinder().getDestPavementJunction(), pMinDist.getPathFinder().getPrimaryCostHeuristic(), pMinDist.getPathFinder().getSecondaryCostHeuristic());                
        
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
			EnvironmentSetup.setUpProperties();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == 5) {
				o = i;
			}
			else if (i.getId() == 7) {
				d = i;
			}
		}
				
		// Initialise a pedestrian, this internally initialises a ped path finder
		boolean minimiseCrossings = false;
		Ped pedMinDist = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
	 		
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
		Ped pedMinCross = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
		
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
	 * Testing the initialisation of a PedPathFinder object.
	 * 
	 * Test fails due to issue with how destination pavement network node is identified. Caused by road link route not including
	 * the link the nearest pavement junction to the destination is on.
	 * 
	 * Need to revise this process to ensure dest node is the node nearest the destination. Will require 
	 * distinguishing between cases where the final road link is included in the strategic path or not.
	 */
	@Test
	public void testPedPathFinder2() {
		
		// Setup environment
		try {
			EnvironmentSetup.setUpProperties();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == 8) {
				o = i;
			}
			else if (i.getId() == 7) {
				d = i;
			}
		}
		
		// Initialise a pedestrian, this internally initialises a ped path finder
		boolean minimiseCrossings = false;
		Ped pedMinDist = new Ped(o, d, GlobalVars.pedVavg, GlobalVars.pedMassAv, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, 180.0, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
				
		// Check the strategic path is as expected
		String[] expectedRoadIDs = {"762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0", "56CF7BBA-28E4-4ACA-9F58-E096E88094FB_0", "B2B9D137-2BA4-4864-8350-2EDAA5910747_0"};
		for (int i=0;i<pedMinDist.getPathFinder().getStrategicPath().size(); i++) {
			assert pedMinDist.getPathFinder().getStrategicPath().get(i).getFID().contentEquals(expectedRoadIDs[i]);
		}
		
		// Check the start and end pavement junctions are as expected
		assert pedMinDist.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_121");
		assert pedMinDist.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_80");
				
		// Now test planning the first tactical path with this ped path finder object
		pedMinDist.getPathFinder().updateTacticalPath();
		
		// Check the current (default) and target tactical alternatives are as expected
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getDefaultJunction().getFID().contentEquals("pave_node_114");
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getTargetJunction().getFID().contentEquals("pave_node_113");
	}
	
	/*
	 * Testing the initialisation of a PedPathFinder object.
	 * 
	 * Also tests updating the tactical path from a junction part way along the journey.
	 * 
	 */
	@Test
	public void testPedPathFinder3() {
		
		// Setup environment
		try {
			EnvironmentSetup.setUpProperties();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == 8) {
				o = i;
			}
			else if (i.getId() == 9) {
				d = i;
			}
		}
		
		// Initialise a pedestrian, this internally initialises a ped path finder
		boolean minimiseCrossings = true;
		Ped pedMinDist = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);		
				
		// Check the strategic path is as expected
		String[] expectedRoadIDs = {"762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0", "A8675945-DE94-4E22-9905-B0623A326221_0", "F4C0B1FB-762C-4492-BB0D-673CC4950CBE_0", "8A9E2D7B-3B48-4A19-B89A-0B4F4D516870_2"};
		for (int i=0;i<pedMinDist.getPathFinder().getStrategicPath().size(); i++) {
			assert pedMinDist.getPathFinder().getStrategicPath().get(i).getFID().contentEquals(expectedRoadIDs[i]);
		}
		
		// Check the start and end pavement junctions are as expected
		assert pedMinDist.getPathFinder().getStartPavementJunction().getFID().contentEquals("pave_node_121");
		assert pedMinDist.getPathFinder().getDestPavementJunction().getFID().contentEquals("pave_node_87");
				
		// Now test planning the first tactical path with this ped path finder object
		pedMinDist.getPathFinder().updateTacticalPath();
		
		// Check target junction is as expected
		assert pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_114");
		
		// Check remainder path is as expected
        List<RepastEdge<Junction>> remainderPath = pedMinDist.getPathFinder().getTacticalPath().getRemainderPath();
        List<Junction> remainderPathNodes = pedMinDist.getPathFinder().getTacticalPath().getNetworkPathFinder().nodePathFromEdges(remainderPath, pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction());
		
        String[] expectedRouteJunctions = {"pave_node_114", "pave_node_112", "pave_node_91"};
		for (int i=0;i< Math.max(expectedRouteJunctions.length, remainderPathNodes.size());i++) {
			assert remainderPathNodes.get(i).getFID().contentEquals(expectedRouteJunctions[i]);
		}
		
		// Test planning the second tactical path
		pedMinDist.getPathFinder().getTacticalPath().updateCurrentJunction();
		pedMinDist.getPathFinder().updateTacticalPath();
		
		// After this update, expect ped's current pavement link to be a direct crossing link. 
		// Check for this and manually call function to represent choice of crossing alternative being made.
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().isDirectCrossing();
		pedMinDist.getPathFinder().getTacticalPath().crossingStartedUpdateCurrentJunction();
		
		
		// This time no remainder path but two junctions in route path
		assert pedMinDist.getPathFinder().getTacticalPath().getRemainderPath().size() == 0;		
		assert pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_112");
		pedMinDist.getPathFinder().getTacticalPath().updateCurrentJunction();
		assert pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_91");
		
		// Update tactical path again
		pedMinDist.getPathFinder().getTacticalPath().updateCurrentJunction();
		assert pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction() == null;
		pedMinDist.getPathFinder().updateTacticalPath();
		
		// Again, direct crossing required, meaning that the current junction is unchanged while the ped chooses a crossing location
		assert pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_91");
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().isDirectCrossing();
		pedMinDist.getPathFinder().getTacticalPath().crossingStartedUpdateCurrentJunction();
		

		// Again expect zero remainder path but multiple nodes in route path
		// In this case there are two equal distance candidate paths and so the chosen path is random.
		// Account for this by checking for the different possible junctions
		assert pedMinDist.getPathFinder().getTacticalPath().getRemainderPath().size() == 0;
		assert pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_89") | pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_90");
		pedMinDist.getPathFinder().getTacticalPath().updateCurrentJunction();
		assert pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_81") | pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_82");
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
			EnvironmentSetup.setUpProperties();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == 3) {
				o = i;
			}
			else if (i.getId() == 5) {
				d = i;
			}
		}
		
		// Initialise a pedestrian, this internally initialises a ped path finder
		boolean minimiseCrossings = false;
		Ped pedMinDist = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);		
		
		// Need to give ped location in order to test updating tactical path following crossing choice
        EnvironmentSetup.context.add(pedMinDist);        
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(pedMinDist.getRad());		
		GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, circle, pedMinDist);
		pedMinDist.setLoc();
		
		// Set up ped path finder
		PedPathFinder ppf = pedMinDist.getPathFinder();
		
		// Check the strategic path is as expected
		String[] expectedRoadIDs = {"9745D155-3C95-4CCD-BC65-0908D57FA83A_0", "F4C0B1FB-762C-4492-BB0D-673CC4950CBE_0", "8A9E2D7B-3B48-4A19-B89A-0B4F4D516870_2", "8A9E2D7B-3B48-4A19-B89A-0B4F4D516870_1"};
		for (int i=0;i<ppf.getStrategicPath().size(); i++) {
			assert ppf.getStrategicPath().get(i).getFID().contentEquals(expectedRoadIDs[i]);
		}
		
		// Check the start and end pavement junctions are as expected
		assert ppf.getStartPavementJunction().getFID().contentEquals("pave_node_108");
		assert ppf.getDestPavementJunction().getFID().contentEquals("pave_node_85");
				
		// Now test planning the first tactical path with this ped path finder object
		ppf.updateTacticalPath();
		
		// Check the current (default) and target tactical alternatives are as expected
		assert ppf.getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_91");
		assert ppf.getTacticalPath().getAccumulatorRoute().getTargetJunction().getFID().contentEquals("pave_node_89");
		
		assert ppf.getTacticalPath().getRemainderPath().size() == 0;
		
		// Test updating route where a choice of crossing alternative is involved
		// Manually set which crossing is chosen to test that tactical route is updated as expected
		List<CrossingAlternative> cas = ppf.getTacticalPath().getAccumulatorRoute().getCAs();
		ppf.getTacticalPath().getAccumulatorRoute().setChosenCA(cas.get(0));
		ppf.getTacticalPath().crossingStartedUpdateCurrentJunction();
		
		// Test planning the second tactical path
		ppf.getTacticalPath().updateCurrentJunction();
		ppf.updateTacticalPath();
		
		// This time no remainder path but two junctions in route path
		assert ppf.getTacticalPath().getRemainderPath().size() == 0;		
		assert ppf.getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_81");
		ppf.getTacticalPath().updateCurrentJunction();
		assert ppf.getTacticalPath().getCurrentJunction() == null;
		
		// Update again
		ppf.updateTacticalPath();
		
		assert ppf.getTacticalPath().getRemainderPath().size() == 0;		
		assert ppf.getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_87");
		assert ppf.getTacticalPath().getAccumulatorRoute().getTargetJunction().getFID().contentEquals("pave_node_88");
		ppf.getTacticalPath().updateCurrentJunction();
		assert ppf.getTacticalPath().getCurrentJunction() == null;
		
		// Update tactical path for the last time (should not reach destination)
		ppf.updateTacticalPath();
		
		// This time, because a crossing wasn't chosen tactical route should be non-crossing
		// Had the ped chosen a crossing location previously would expect them to choose a crossing link for the last section of the journey 
		assert ppf.getTacticalPath().getRemainderPath().size() == 0;		
		assert ppf.getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_85");
		assert ppf.getTacticalPath().getAccumulatorRoute().isBlank();
		
		// Check the target coordiante is the destination
		assert ppf.getTacticalPath().getTargetCoordinate().equals2D(d.getGeom().getCoordinate());
		ppf.getTacticalPath().updateCurrentJunction();
		assert ppf.getTacticalPath().getCurrentJunction() == null;
		
	}
	
	/*
	 * Test updating the tactical route when at the final link of a strategic path in order to check that the destination coordinate is correctly returned
	 * 
	 * Test fails due to issues with identifying starting junction where starting coordinate is near the end of the current road link.
	 * 
	 */
	@Test
	public void testPedPathFinder5() {
		
		try {
			EnvironmentSetup.setUpProperties();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		OD o = null;
		OD d = null;
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId()==10) {
				o = i;
			}
			else if (i.getId()==8) {
				d = i;
			}
		}
		
		boolean minimiseCrossings = true;
		Ped p = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
		
		Junction startJ = p.getPathFinder().getStartPavementJunction();
		Junction destJ = p.getPathFinder().getDestPavementJunction();
		
		assert startJ.getFID().contentEquals("pave_node_91");
		assert destJ.getFID().contentEquals("pave_node_121");
		
		
		// Initialise tactical paths and test path junctions are correct as path gets updated.
		p.getPathFinder().updateTacticalPath();
		
		assert p.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_112");
		
		p.getPathFinder().getTacticalPath().updateCurrentJunction();
		assert p.getPathFinder().getTacticalPath().getCurrentJunction() == null;

		p.getPathFinder().updateTacticalPath();
		
		// Direct crossing required to get from pave_node_112 to pave_node_114. So check for this and that the junction is not updated until the crossing is chosen
		assert p.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_112");
		assert p.getPathFinder().getTacticalPath().getAccumulatorRoute().isDirectCrossing();
		p.getPathFinder().getTacticalPath().updateTargetCoordiante(); // Updating the target coordinate should not update the junction when the crossing hasn't been chosen.
		assert p.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_112");
		
		// Now manually update the junction as if a crossing had been chosen
		p.getPathFinder().getTacticalPath().crossingStartedUpdateCurrentJunction();
		
		assert p.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_114");
		p.getPathFinder().getTacticalPath().updateCurrentJunction();
		assert p.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals(destJ.getFID());
		assert p.getPathFinder().getTacticalPath().getTargetCoordinate().equals2D(d.getGeom().getCoordinate());	
	}
	
	/*
	 * Test that pedestrian gets added and removed from list of pedestrians on a road link as it progresses on its route.
	 */
	@Test
	void testPedsOnRoadLinkList() {
		// Setup the environment
		try {
			EnvironmentSetup.setUpProperties();
			
			EnvironmentSetup.setUpPedObstructions();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Choose origin and destination and initialise pedestrian
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == 3) {
				o = i;
			}
			else if (i.getId() == 2) {
				d = i;
			}
		}
		
		// Initialise a pedestrian, this internally initialises a ped path finder
		boolean minimiseCrossings = false;
		Ped pedMinDist = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);		
		
		// Need to give ped location in order to test updating tactical path following crossing choice
        EnvironmentSetup.context.add(pedMinDist);        
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(pedMinDist.getRad());		
		GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, circle, pedMinDist);
		pedMinDist.setLoc();
		
		
		// Before peds first step there should be no peds assigns to first road link in route
		assert pedMinDist.getPathFinder().getStrategicPath().get(0).getPeds().size()==0;
		
		// Then after steeping the ped the ped should be added to the road link
		try {
			pedMinDist.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		assert pedMinDist.getPathFinder().getStrategicPath().get(0).getPeds().size()==1;
		assert pedMinDist.getPathFinder().getStrategicPath().get(0).getPeds().get(0) == pedMinDist;
		
		// Now record first road link and keep stepping until ped has moved onto next road link. Check that record of which peds are on which links is updates
		RoadLink firstRL = pedMinDist.getPathFinder().getStrategicPath().get(0);
		while (pedMinDist.getPathFinder().getStrategicPath().get(0).getFID().contentEquals(firstRL.getFID())) {
			try {
				pedMinDist.step();
			} catch (Exception e) {
				e.printStackTrace();
			}			
		}
		
		assert firstRL.getPeds().size()==0;
		assert pedMinDist.getPathFinder().getStrategicPath().get(0).getPeds().size()==1;
		assert pedMinDist.getPathFinder().getStrategicPath().get(0).getPeds().get(0) == pedMinDist;
	}
	
	
	/*
	 * Test that flags that indicate whether pedestrian requires crossing, has chosen a crossing location
	 * and is crossing work as expected.
	 * 
	 * These are used by vehicle agents to identify which pedestrians to yield to.
	 */
	@Test
	void testPedCrossingFlags1() {
		// Setup the environment
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRandomDistributions(1);
			
			EnvironmentSetup.setUpPedObstructions();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();;
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Choose origin and destination and initialise pedestrian
		// Set the IDs of the road network junctions to travel to and get strategic path between these
		OD o = null;
		OD d = null;
		
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == 5) {
				o = i;
			}
			else if (i.getId() == 6) {
				d = i;
			}
		}
		
		// Create distance minimising and crossing minimising pedestrians
		Ped pedMinDist = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, false, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);		
		
		Ped pedMinCross = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, 30, 60, 1.0, 0.75, true, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
		
		// Step the peds once to initiase their routes
		try {
			addPedToWorld(pedMinCross, o);
			pedMinCross.step();
			
			addPedToWorld(pedMinDist, o);
			pedMinDist.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().crossingRequired() == true;
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().caChosen() == false;
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().isCrossing() == false;
		
		assert pedMinCross.getPathFinder().getTacticalPath().getAccumulatorRoute().crossingRequired() == false;
		assert pedMinCross.getPathFinder().getTacticalPath().getAccumulatorRoute().caChosen() == false;
		assert pedMinCross.getPathFinder().getTacticalPath().getAccumulatorRoute().isCrossing() == false;
		
		// Step the path finder until a crossing is chosen. This keeps the ped stationary so creates a situation where the ped is not next to the crossing coord
		while (pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().caChosen() == false) {
			try {
				pedMinDist.getPathFinder().step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// At this point crossing is chosen but ped hasn't started crossing yet
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().isCrossing() == false;
		
		// Once target coordinate is updated, that means ped has started crossing since next target coordinate is the crossing
		assert pedMinDist.getPathFinder().getTacticalPath().getTargetCoordinate().equals2D(pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getCrossingCoordinates().getLast());
		pedMinDist.getPathFinder().getTacticalPath().updateTargetCoordiante();
		
		// Ped should now be flagged as reaching crossing but will default to yielding and so will not start crossing
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().reachedCrossing() == true;
		assert pedMinDist.getYield() == true;
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().isCrossing() == false;
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().crossingRequired() == true;
		
		// step the ped along it should stop yielding, update it's position and be classed as starting to yield
		Coordinate prevLoc = pedMinDist.getLoc();
		try {
			pedMinDist.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		assert pedMinDist.getYield()==false;
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().isCrossing()==true;
		
		// Update target coordinate again and check that ped is no longer crossing
		pedMinDist.getPathFinder().getTacticalPath().updateTargetCoordiante();
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().isCrossing() == false;
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().crossingRequired() == false;
	}
	
	/*
	 * Test tactical path finding for a specific case where informal crossing is allowed. This test complements the test testOD60NoInformalCrossingPathTacticalPath
	 * which considers the case where informal crossing is not allowed.
	 */
	@Test
	public void testOD162InformalCrossingPathTacticalPath() {
		
		try {
			String origTestDir = EnvironmentSetup.testGISDir;
			EnvironmentSetup.testGISDir += "/clapham_common/";
			
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRandomDistributions(1);
			
			EnvironmentSetup.setUpPedObstructions();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternativesTollerance03.shp");
			
			EnvironmentSetup.setUpPedODs("OD_pedestrian_nodes.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
			EnvironmentSetup.testGISDir = origTestDir;
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}

		Ped p = EnvironmentSetup.createPedestrian(null, null, "od_162", "od_0", true);
		
		PedPathFinder ppf = p.getPathFinder();
		
		// Check the road network route
		String[] rnrIDs = {"or_link_313", "or_link_314", "or_link_266", "or_link_265", "or_link_233", "or_link_229", "or_link_205", "or_link_204", "or_link_191", "or_link_180", "or_link_179", "or_link_186"};
		for (int i=0; i<rnrIDs.length; i++) {
			assert rnrIDs[i].contentEquals(ppf.getStrategicPath().get(i).getFID());
		}
		
		// Two updates to get to point where tactical path is different under no informal crossing scenario
		ppf.updateTacticalPath();
		
		// Check tactical path is as expected
		// - since expect only one link in tacgtical path initially, after update tactical path should be empty and current link should match the expected first link in tactical path
		String expectedCurrentLink = "pave_link_595_597";
		assert ppf.getTacticalPath().getRoutePath().size()==0;
		assert ( (NetworkEdge<Junction>)ppf.getTacticalPath().getCurrentEdge()).getRoadLink().getFID().contentEquals(expectedCurrentLink);
		
		/// The reminder path should go up to the end of the planning horizon
		String[] expectedRemainderPath = {"pave_link_588_596", "pave_link_515_588", "pave_link_515_516","pave_link_451_516","pave_link_440_451","pave_link_393_440"};
		assert ppf.getTacticalPath().getRemainderPath().size() == expectedRemainderPath.length;
		for (int i=0; i<Math.max(expectedRemainderPath.length, ppf.getTacticalPath().getRemainderPath().size()); i++) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) ppf.getTacticalPath().getRemainderPath().get(i);
			assert ne.getRoadLink().getFID().contentEquals(expectedRemainderPath[i]);
		}
	}
	
	/*
	 * Test tactical path finding in cases where pedestrian is prevented from making informal crossings.
	 * 
	 * Specifically check that when a detour is required, the initial tactical path and path to end of ph are correctly identified and the road 
	 * network route is correctly updated.
	 */
	@Test
	public void testOD162NoInformalCrossingPathTacticalPath() {
		
		try {
			String origTestDir = EnvironmentSetup.testGISDir;
			EnvironmentSetup.testGISDir += "/clapham_common/";
			
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRandomDistributions(1);
			
			EnvironmentSetup.setUpPedObstructions();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternativesTollerance03.shp");
			
			EnvironmentSetup.setUpPedODs("OD_pedestrian_nodes.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
			EnvironmentSetup.testGISDir = origTestDir;
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}

		Network<Junction> pavementNetwork = SpaceBuilder.getNetwork(GlobalVars.CONTEXT_NAMES.PAVEMENT_NETWORK);
		Geography<CrossingAlternative> caG = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.CA_GEOGRAPHY);
		
		// Importantly, now edit the pavement network to remove road crossing links where no crossing infrastructure is available
		SpaceBuilder.setInformalCrossingStatus(false);
		SpaceBuilder.removeCrossingLinksFromPavementNetwork(pavementNetwork, caG);
		
		Ped p = EnvironmentSetup.createPedestrian(null, null, "od_162", "od_0", true);
		
		PedPathFinder ppf = p.getPathFinder();
		
		// Check the road network route
		String[] rnrIDs = {"or_link_313", "or_link_314", "or_link_266", "or_link_265", "or_link_233", "or_link_229", "or_link_205", "or_link_204", "or_link_191", "or_link_180", "or_link_179", "or_link_186"};
		for (int i=0; i<rnrIDs.length; i++) {
			assert rnrIDs[i].contentEquals(ppf.getStrategicPath().get(i).getFID());
		}
		
		ppf.updateTacticalPath();
		
		// Road network route should be different now that ped has to travel a different way to reach their destination
		String[] rnrIDsUpdated = {"or_link_319","or_link_321","or_link_323","or_link_332","or_link_330","or_link_331","or_link_328","or_link_327","or_link_326","or_link_267","or_link_265","or_link_233","or_link_229","or_link_205","or_link_204","or_link_191","or_link_180","or_link_179","or_link_186"};
		for (int i=0; i<rnrIDsUpdated.length; i++) {
			assert rnrIDsUpdated[i].contentEquals(ppf.getStrategicPath().get(i).getFID());
		}
		
		// Check tactical path is as expected
		String expectedCurrentEdge = "pave_link_597_606";
		assert ( (NetworkEdge<Junction>) ppf.getTacticalPath().getCurrentEdge()).getRoadLink().getFID().contentEquals(expectedCurrentEdge);
		
		String[] expectedTacticalPath = {};
		assert ppf.getTacticalPath().getRoutePath().size() == expectedTacticalPath.length;	
		
		String[] expectedRemainderPath = {"pave_link_606_608","pave_link_608_629","pave_link_623_629","pave_link_623_664","pave_link_627_664","pave_link_620_627","pave_link_617_620","pave_link_615_617","pave_link_516_615", "pave_link_451_516", "pave_link_440_451","pave_link_393_440"};
		assert ppf.getTacticalPath().getRemainderPath().size() == expectedRemainderPath.length;
		for (int i=0; i<Math.max(expectedRemainderPath.length, ppf.getTacticalPath().getRemainderPath().size()); i++) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) ppf.getTacticalPath().getRemainderPath().get(i);
			assert ne.getRoadLink().getFID().contentEquals(expectedRemainderPath[i]);
		}
		
		// This route path meets back up with the strategic path 3 links in, therefore expect strategic path to be reduced by 3 after next update
		while (ppf.getTacticalPath().getCurrentEdge()!=null) {
			ppf.getTacticalPath().updateCurrentJunction();
		}
		ppf.updateTacticalPath();
		assert ppf.getStrategicPath().size() == rnrIDsUpdated.length - 1;
	}

	
	/*
	 * Test No Informal crossing path finding from origin od_78. The tactical path extends all the way to the 
	 * destination junction without meeting back up with the strategic path.
	 */
	@Test
	public void testOD78NoInformalCrossingPathTacticalPath() {
		
		try {
			String origTestDir = EnvironmentSetup.testGISDir;
			EnvironmentSetup.testGISDir += "/clapham_common/";
			
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRandomDistributions(1);
			
			EnvironmentSetup.setUpPedObstructions();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternativesTollerance05.shp");
			
			EnvironmentSetup.setUpPedODs("OD_pedestrian_nodes.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
			EnvironmentSetup.testGISDir = origTestDir;
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}

		Network<Junction> pavementNetwork = SpaceBuilder.getNetwork(GlobalVars.CONTEXT_NAMES.PAVEMENT_NETWORK);
		Geography<CrossingAlternative> caG = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.CA_GEOGRAPHY);
		
		// Importantly, now edit the pavement network to remove road crossing links where no crossing infrastructure is available
		SpaceBuilder.setInformalCrossingStatus(false);
		SpaceBuilder.removeCrossingLinksFromPavementNetwork(pavementNetwork, caG);
		
		Ped p = EnvironmentSetup.createPedestrian("od_78", "od_0", GlobalVars.pedVavg, GlobalVars.pedMassAv, 0.7634, 1.728, 0.976, 7.22, 78, 86, 1.85, 0.75, true, 276.37);;
		
		PedPathFinder ppf = p.getPathFinder();
		
		// Check the road network route
		String[] rnrIDs = {"or_link_58","or_link_57","or_link_295","or_link_293","or_link_290","or_link_289","or_link_212","or_link_197","or_link_194","or_link_187","or_link_186"};
		for (int i=0; i<rnrIDs.length; i++) {
			assert rnrIDs[i].contentEquals(ppf.getStrategicPath().get(i).getFID());
		}
		
		ppf.updateTacticalPath();
		
		// Expect strategic path to have changed now
		String[] rnrIDs2 = {"or_link_297","or_link_299","or_link_305","or_link_270","or_link_259","or_link_257","or_link_258","or_link_255","or_link_253","or_link_254","or_link_206","or_link_204","or_link_191","or_link_180","or_link_179"};
		for (int i=0; i<rnrIDs2.length; i++) {
			assert rnrIDs2[i].contentEquals(ppf.getStrategicPath().get(i).getFID());
		}
		
		// Check tactical path is as expected
		String expectedCurrentEdge = "pave_link_113_566";
		assert ( (NetworkEdge<Junction>) ppf.getTacticalPath().getCurrentEdge()).getRoadLink().getFID().contentEquals(expectedCurrentEdge);
		
		String[] expectedTacticalPath = {};
		assert ppf.getTacticalPath().getRoutePath().size() == expectedTacticalPath.length;
		for (int i=0; i<Math.max(expectedTacticalPath.length, ppf.getTacticalPath().getRoutePath().size()); i++) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) ppf.getTacticalPath().getRoutePath().get(i);
			assert ne.getRoadLink().getFID().contentEquals(expectedTacticalPath[i]);
		}		
		
		String[] expectedRemainderPath = {"pave_link_566_578","pave_link_578_580","pave_link_521_580","pave_link_501_521","pave_link_501_504","pave_link_498_504","pave_link_498_500","pave_link_492_500","pave_link_488_492","pave_link_488_490","pave_link_398_490","pave_link_398_654","pave_link_395_654","pave_link_371_395","pave_link_371_372","pave_link_346_372","pave_link_346_363"};
		assert ppf.getTacticalPath().getRemainderPath().size() == expectedRemainderPath.length;
		for (int i=0; i<Math.max(expectedRemainderPath.length, ppf.getTacticalPath().getRemainderPath().size()); i++) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) ppf.getTacticalPath().getRemainderPath().get(i);
			assert ne.getRoadLink().getFID().contentEquals(expectedRemainderPath[i]);
		}
	}
	
	/*
	 * Test No Informal crossing path finding from origin od_78. The tactical path extends all the way to the 
	 * destination junction without meeting back up with the strategic path.
	 */
	//@Test
	public void testOD133NoInformalCrossingPathTacticalPath() {
		
		try {
			String origTestDir = EnvironmentSetup.testGISDir;
			EnvironmentSetup.testGISDir += "/clapham_common/";
			
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRandomDistributions(1);
			
			EnvironmentSetup.setUpPedObstructions();
			
			EnvironmentSetup.setUpRoads();
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternativesTollerance05.shp");
			
			EnvironmentSetup.setUpPedODs("OD_pedestrian_nodes.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
			EnvironmentSetup.testGISDir = origTestDir;
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}

		Network<Junction> pavementNetwork = SpaceBuilder.getNetwork(GlobalVars.CONTEXT_NAMES.PAVEMENT_NETWORK);
		Geography<CrossingAlternative> caG = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.CA_GEOGRAPHY);
		
		// Importantly, now edit the pavement network to remove road crossing links where no crossing infrastructure is available
		SpaceBuilder.setInformalCrossingStatus(false);
		SpaceBuilder.removeCrossingLinksFromPavementNetwork(pavementNetwork, caG);
		
		Ped p = EnvironmentSetup.createPedestrian("od_133", "od_0", GlobalVars.pedVavg, GlobalVars.pedMassAv, 0.7634, 1.728, 0.976, 7.22, 78, 86, 1.85, 0.75, true, 276.37);;
		
		PedPathFinder ppf = p.getPathFinder();
		
		ppf.updateTacticalPath();
		
		// Check tactical path is as expected
		String expectedCurrentEdge = "pave_link_113_566";
		assert ( (NetworkEdge<Junction>) ppf.getTacticalPath().getCurrentEdge()).getRoadLink().getFID().contentEquals(expectedCurrentEdge);
		
		String[] expectedTacticalPath = {"pave_link_566_578","pave_link_578_580","pave_link_521_580","pave_link_501_521","pave_link_501_504","pave_link_498_504","pave_link_498_500","pave_link_492_500","pave_link_488_492","pave_link_488_490","pave_link_398_490","pave_link_398_654","pave_link_395_654","pave_link_370_395","pave_link_370_373","pave_link_346_373","pave_link_346_363"};
		assert ppf.getTacticalPath().getRoutePath().size() == expectedTacticalPath.length;
		for (int i=0; i<Math.max(expectedTacticalPath.length, ppf.getTacticalPath().getRoutePath().size()); i++) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) ppf.getTacticalPath().getRoutePath().get(i);
			assert ne.getRoadLink().getFID().contentEquals(expectedTacticalPath[i]);
		}		
		
		String[] expectedRemainderPath = {};
		assert ppf.getTacticalPath().getRemainderPath().size() == expectedRemainderPath.length;
		for (int i=0; i<Math.max(expectedRemainderPath.length, ppf.getTacticalPath().getRemainderPath().size()); i++) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) ppf.getTacticalPath().getRemainderPath().get(i);
			assert ne.getRoadLink().getFID().contentEquals(expectedRemainderPath[i]);
		}
	}
	

}
