package repastInterSim.tests;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.OD;
import repastInterSim.environment.RoadLink;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.PedPathFinder;
import repastInterSim.pathfinding.TacticalRoute;

class TacticalRouteTest {

	Context<Object> context = new DefaultContext<Object>();;
	
	String testGISDir = ".//data//test_gis_data//";
	String pedestrianRoadsPath = null;
	String vehicleRoadsPath = null;
	String roadLinkPath = null;
	String pavementLinkPath = null;
	String pedJPath = null;
	String serialisedLookupPath = null;


	/*
	 * This tests a scenario where the strategic path is 2 links long.
	 * 
	 * The test checks that TacticalRoute object is initialised with the expected currentJunction and that the remainderPath leads to the end pavement junction
	 */
	@Test
	void testTacticalRouteSetup1(){
		
		
		// Setup environment
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			
			EnvironmentSetup.setUpORRoadNetwork(false);
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
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
			if (i.getId() == 9) {
				o = i;
			}
			else if (i.getId() == 7) {
				d = i;
			}
		}
		
		// Set up ped path finder
		boolean minimiseCrossings = false;		
		Ped p = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
		
        context.add(p);        
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(p.getRad());		
		GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, circle, p);
        p.setLoc();
        
        PedPathFinder ppf = p.getPathFinder();
        
		// Check the start and end pavement junctions are as expected
		assert ppf.getStartPavementJunction().getFID().contentEquals("pave_node_87");
		assert ppf.getDestPavementJunction().getFID().contentEquals("pave_node_93");
		
		// Now test planning the first tactical path with this ped path finder object
        int tacticalHorizonLinks = PedPathFinder.getNLinksWithinAngularDistance(ppf.getStrategicPath(), p.getpHorizon());
        TacticalRoute tr = PedPathFinder.planTacticalPath(EnvironmentSetup.pavementNetwork, EnvironmentSetup.caGeography, EnvironmentSetup.roadGeography, tacticalHorizonLinks, p, ppf.getStrategicPath(), ppf.getStartPavementJunction(), ppf.getDestPavementJunction(), ppf.getPrimaryCostHeuristic(), ppf.getSecondaryCostHeuristic());                
        
        
        // Test needs to get end junction of route path and use that rather than current junction, current junction could not connect to remainder path if first link is a crossing link
        Junction startRemainderPathJunction = null;
        if (tr.getAccumulatorRoute().getTargetJunction() != null) {
        	startRemainderPathJunction = tr.getAccumulatorRoute().getTargetJunction(); 
        }
        else {
        	startRemainderPathJunction = tr.getCurrentJunction();
        }
        
        List<RepastEdge<Junction>> remainderPath = tr.getRemainderPath();
        List<Junction> remainderPathNodes = tr.getNetworkPathFinder().nodePathFromEdges(remainderPath, startRemainderPathJunction);
        
		// Check the end junctions of the chosen remainder path - planning horizon extends to the desitination so end junction is the same as destination junction
		assert remainderPathNodes.get(remainderPathNodes.size()-1).getFID().contentEquals("pave_node_93");
		assert tr.getCurrentJunction().getFID().contentEquals("pave_node_79") | tr.getCurrentJunction().getFID().contentEquals("pave_node_81");
	}
	
	/*
	 * This tests a scenario where the strategic path is 1 link long and primary crossing is required o reach destination.
	 * 
	 * This means the setting of up an AccumulatorRoute is tested.
	 * 
	 */
	@Test
	void testTacticalRouteSetup2(){
		
		
		// Setup environment
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRandomDistributions();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			
			EnvironmentSetup.setUpORRoadNetwork(false);
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			EnvironmentSetup.setUpPedObstructions();
			
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
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
			if (i.getId() == 11) {
				o = i;
			}
			else if (i.getId() == 3) {
				d = i;
			}
		}
		
		// Set up ped path finder
		boolean minimiseCrossings = false;
		
		Ped p = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
		
        context.add(p);        
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(p.getRad());		
		GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, circle, p);
        p.setLoc();
        
        PedPathFinder ppf = p.getPathFinder();
		
		// Check the start and end pavement junctions are as expected
		assert ppf.getStartPavementJunction().getFID().contentEquals("pave_node_89");
		assert ppf.getDestPavementJunction().getFID().contentEquals("pave_node_108");
		
		// Now test planning the first tactical path with this ped path finder object
        int tacticalHorizonLinks = PedPathFinder.getNLinksWithinAngularDistance(ppf.getStrategicPath(), p.getpHorizon());
        TacticalRoute tr = PedPathFinder.planTacticalPath(EnvironmentSetup.pavementNetwork, EnvironmentSetup.caGeography, EnvironmentSetup.roadGeography, tacticalHorizonLinks, p, ppf.getStrategicPath(), ppf.getStartPavementJunction(), ppf.getDestPavementJunction(), ppf.getPrimaryCostHeuristic(), ppf.getSecondaryCostHeuristic());        
        
		assert tr.getCurrentJunction().getFID().contentEquals("pave_node_106"); // This is the default junction
		assert tr.getAccumulatorRoute().getTargetJunction().getFID().contentEquals("pave_node_108");
		
		// Now test accumulating activation
		tr.step();
	}
	
	@Test
	public void testGetCrossingAlternatives() {
		
		// Setup environment
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			
			EnvironmentSetup.setUpORRoadNetwork(false);
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
			
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// origin and destination don't matter, just choose arbitrarily
		OD o = null;
		OD d = null;
		for (OD i : EnvironmentSetup.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId()==3) {
				o = i;
			}
			else if (i.getId()==1) {
				d = i;
			}
		}
        boolean minimiseCrossings = false;
        Ped ped = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, EnvironmentSetup.pavementJunctionGeography, EnvironmentSetup.pavementNetwork);
        
		// Select which set of road links to get crossing alternatives for
		List<RoadLink> rls = new ArrayList<RoadLink>();
		String[] roadLinkIDs = {"A8675945-DE94-4E22-9905-B0623A326221_0", "9745D155-3C95-4CCD-BC65-0908D57FA83A_0"};
		for (RoadLink rl: EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			for (int i=0;i<roadLinkIDs.length; i++) {
				if(rl.getFID().contentEquals(roadLinkIDs[i])) {
					rls.add(rl);
				}
			}
		}
		
		// Get crossing alternatives
		List<CrossingAlternative> cas = TacticalRoute.getCrossingAlternatives(EnvironmentSetup.caGeography, rls, ped, EnvironmentSetup.roadGeography);
		
		// Validate crossing alternatives by checking types are as expected
		String[] expectedTypes = {"unsignalised", "unsignalised", "unmarked"};
		for (int i=0;i<expectedTypes.length; i++) {
			assert cas.get(i).getType().contentEquals(expectedTypes[i]);
		}
		
		// Repeat with just one road link
		rls = new ArrayList<RoadLink>();
		String[] roadLinkIDs2 = {"A8675945-DE94-4E22-9905-B0623A326221_0"};
		for (RoadLink rl: EnvironmentSetup.orRoadLinkGeography.getAllObjects()) {
			for (int i=0;i<roadLinkIDs2.length; i++) {
				if(rl.getFID().contentEquals(roadLinkIDs2[i])) {
					rls.add(rl);
				}
			}
		}
		
		// Get crossing alternatives
		cas = TacticalRoute.getCrossingAlternatives(EnvironmentSetup.caGeography, rls, ped, EnvironmentSetup.roadGeography);
		
		String[] expectedTypes2 = {"unsignalised", "unmarked"};
		for (int i=0;i<expectedTypes2.length; i++) {
			assert cas.get(i).getType().contentEquals(expectedTypes2[i]);
		}
	}
	

}
