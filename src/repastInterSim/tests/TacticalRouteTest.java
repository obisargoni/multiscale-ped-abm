package repastInterSim.tests;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.List;

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
	
	void setUpProperties() throws IOException {
		IO.readProperties();
	}
	
	void setUpObjectGeography() {
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		SpaceBuilder.geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, context, geoParams);
		SpaceBuilder.geography.setCRS(GlobalVars.geographyCRSString);
		context.add(SpaceBuilder.geography);
	}
	
	void setUpRoads() throws Exception {
		pedestrianRoadsPath = testGISDir + "topographicAreaPedestrian.shp";
		vehicleRoadsPath = testGISDir + "topographicAreaVehicle.shp";
		serialisedLookupPath = testGISDir + "road_link_roads_cache.serialised";
		
		// Get road geography
		Context<Road> testRoadContext = new RoadContext();
		GeographyParameters<Road> GeoParamsRoad = new GeographyParameters<Road>();
		SpaceBuilder.roadGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testRoadGeography", testRoadContext, GeoParamsRoad);
		SpaceBuilder.roadGeography.setCRS(GlobalVars.geographyCRSString);

		
		// Load vehicle origins and destinations
		try {
			GISFunctions.readShapefile(Road.class, vehicleRoadsPath, SpaceBuilder.roadGeography, testRoadContext);
			GISFunctions.readShapefile(Road.class, pedestrianRoadsPath, SpaceBuilder.roadGeography, testRoadContext);
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		SpatialIndexManager.createIndex(SpaceBuilder.roadGeography, Road.class);
		
		/**
		 * To model two way roads, two RoadLink objects are created with the same id. Therefore, when getting the RoadLink objects
		 * that this Road object is associated, allow for multiple RoadLink objects to be associated.
		 */
		
		for (Road r: SpaceBuilder.roadGeography.getAllObjects()) {
			List<RoadLink> roadLinks = new ArrayList<RoadLink>();
			for(RoadLink rl: SpaceBuilder.roadLinkGeography.getAllObjects()) {
				// Iterating over the vehicle road links (ITN) but using their corresponding ped road link (open road) id to check whether they belong to this vehicle polygon
				if (rl.getPedRLID().contentEquals(r.getRoadLinkID())) {
					roadLinks.add(rl);
				}
			}
			r.setRoadLinks(roadLinks);
		}

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
		SpaceBuilder.roadLinkGeography = setUpLinks(roadLinkFile);
	}
	
	void setUpPavementLinks(String linkFile) throws MalformedURLException, FileNotFoundException {
		SpaceBuilder.pavementLinkGeography = setUpLinks(linkFile);
	}
		
	void setUpODs(String odFile) throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> ODContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		SpaceBuilder.pedestrianDestinationGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", ODContext, GeoParamsOD);
		SpaceBuilder.pedestrianDestinationGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String testODFile = testGISDir + odFile;
		GISFunctions.readShapefile(OD.class, testODFile, SpaceBuilder.pedestrianDestinationGeography, ODContext);
		SpatialIndexManager.createIndex(SpaceBuilder.pedestrianDestinationGeography, OD.class);
	}
	
	void setUpPedObstructions() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<PedObstruction> pedObstructContext = new PedObstructionContext();
		GeographyParameters<PedObstruction> GeoParams = new GeographyParameters<PedObstruction>();
		SpaceBuilder.pedObstructGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedObstructGeography", pedObstructContext, GeoParams);
		SpaceBuilder.pedObstructGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testPedObstructFile = testGISDir + "boundaryPedestrianVehicleArea.shp";
		GISFunctions.readShapefile(PedObstruction.class, testPedObstructFile, SpaceBuilder.pedObstructGeography, pedObstructContext);
		SpatialIndexManager.createIndex(SpaceBuilder.pedObstructGeography, PedObstruction.class);
	}
	
	void setUpCrossingAlternatives() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<CrossingAlternative> caContext = new CAContext();
		GeographyParameters<CrossingAlternative> GeoParams = new GeographyParameters<CrossingAlternative>();
		SpaceBuilder.caGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("caGeography", caContext, GeoParams);
		SpaceBuilder.caGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testCAFile = testGISDir + "crossing_lines.shp";
		GISFunctions.readShapefile(CrossingAlternative.class, testCAFile, SpaceBuilder.caGeography, caContext);
		SpatialIndexManager.createIndex(SpaceBuilder.caGeography, CrossingAlternative.class);
	}
	
	void setUpRoadNetwork(boolean isDirected) {
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		Geography<Junction> junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", junctionContext, GeoParamsJunc);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);
		
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		SpaceBuilder.roadNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(SpaceBuilder.roadLinkGeography, junctionContext, junctionGeography, SpaceBuilder.roadNetwork);
	}
	
	void setUpPavementNetwork() {
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>("PAVEMENT_NETWORK", SpaceBuilder.pavementJunctionContext, false);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		SpaceBuilder.pavementNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(SpaceBuilder.pavementLinkGeography, SpaceBuilder.pavementJunctionContext, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementNetwork);
	}
	
	void setUpPedJunctions() throws Exception {
		setUpProperties();
		pedJPath = testGISDir + IO.getProperty("PavementJunctionsShapefile");
		
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
		for(Junction j: SpaceBuilder.roadNetwork.getNodes()) {
			
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
			shortestRoute = rnr.getShortestRoute(SpaceBuilder.roadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
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
	 * This tests a scenario where the strategic path is 2 links long.
	 * 
	 * The test checks that TacticalRoute object is initialised with the expected currentJunction and that the remainderPath leads to the end pavement junction
	 */
	@Test
	void testTacticalRouteSetup1(){
		
		
		// Setup environment
		try {
			setUpObjectGeography();
			
			setUpRoadLinks("mastermap-itn RoadLink Intersect Within with orientation.shp");

			setUpRoads();

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
		
		for (OD i : SpaceBuilder.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == 5) {
				o = i;
			}
			else if (i.getId() == 1) {
				d = i;
			}
		}
		
		// Set up ped path finder
		boolean minimiseCrossings = false;		
		Ped p = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementNetwork);
		
        context.add(p);        
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(p.getRad());		
		GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, circle, p);
        p.setLoc();
        
        PedPathFinder ppf = p.getPathFinder();
        
		// Check the start and end pavement junctions are as expected
		assert ppf.getStartPavementJunction().getFID().contentEquals("pave_node_87");
		assert ppf.getDestPavementJunction().getFID().contentEquals("pave_node_93");
		
		// Now test planning the first tactical path with this ped path finder object
        int tacticalHorizonLinks = PedPathFinder.getNLinksWithinAngularDistance(ppf.getStrategicPath(), p.getpHorizon());
        TacticalRoute tr = PedPathFinder.planTacticalPath(SpaceBuilder.pavementNetwork, SpaceBuilder.caGeography, SpaceBuilder.roadGeography, tacticalHorizonLinks, p, ppf.getStrategicPath(), ppf.getStartPavementJunction(), ppf.getDestPavementJunction(), ppf.getPrimaryCostHeuristic(), ppf.getSecondaryCostHeuristic());                
        
        
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
			setUpObjectGeography();
			
			setUpRoadLinks("mastermap-itn RoadLink Intersect Within with orientation.shp");

			setUpRoads();

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
		
		for (OD i : SpaceBuilder.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == 6) {
				o = i;
			}
			else if (i.getId() == 3) {
				d = i;
			}
		}
		
		// Set up ped path finder
		boolean minimiseCrossings = false;
		
		Ped p = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementNetwork);
		
        context.add(p);        
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(p.getRad());		
		GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, circle, p);
        p.setLoc();
        
        PedPathFinder ppf = p.getPathFinder();
		
		// Check the start and end pavement junctions are as expected
		assert ppf.getStartPavementJunction().getFID().contentEquals("pave_node_89");
		assert ppf.getDestPavementJunction().getFID().contentEquals("pave_node_108");
		
		// Now test planning the first tactical path with this ped path finder object
        int tacticalHorizonLinks = PedPathFinder.getNLinksWithinAngularDistance(ppf.getStrategicPath(), p.getpHorizon());
        TacticalRoute tr = PedPathFinder.planTacticalPath(SpaceBuilder.pavementNetwork, SpaceBuilder.caGeography, SpaceBuilder.roadGeography, tacticalHorizonLinks, p, ppf.getStrategicPath(), ppf.getStartPavementJunction(), ppf.getDestPavementJunction(), ppf.getPrimaryCostHeuristic(), ppf.getSecondaryCostHeuristic());        
        
		assert tr.getCurrentJunction().getFID().contentEquals("pave_node_106"); // This is the default junction
		assert tr.getAccumulatorRoute().getTargetJunction().getFID().contentEquals("pave_node_108");
		
		// Now test accumulating activation
		tr.step();
	}
	
	@Test
	public void testGetCrossingAlternatives() {
		
		// Setup environment
		try {
			setUpObjectGeography();
			
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPedJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
			
			setUpRoads();
			
			setUpODs("OD_pedestrian_nodes.shp");
			
			setUpCrossingAlternatives();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// origin and destination don't matter, just choose arbitrarily
		OD o = null;
		OD d = null;
		for (OD i : SpaceBuilder.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId()==2) {
				o = i;
			}
			else if (i.getId()==4) {
				d = i;
			}
		}
        boolean minimiseCrossings = false;
        Ped ped = new Ped(o, d, 0.5, 1.0, 0.9, 3.0, minimiseCrossings, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementNetwork);
        
		// Select which set of road links to get crossing alternatives for
		List<RoadLink> rls = new ArrayList<RoadLink>();
		String[] roadLinkIDs = {"A8675945-DE94-4E22-9905-B0623A326221_0", "9745D155-3C95-4CCD-BC65-0908D57FA83A_0"};
		for (RoadLink rl: SpaceBuilder.roadLinkGeography.getAllObjects()) {
			for (int i=0;i<roadLinkIDs.length; i++) {
				if(rl.getFID().contentEquals(roadLinkIDs[i])) {
					rls.add(rl);
				}
			}
		}
		
		// Get crossing alternatives
		List<CrossingAlternative> cas = TacticalRoute.getCrossingAlternatives(SpaceBuilder.caGeography, rls, ped, SpaceBuilder.roadGeography);
		
		// Validate crossing alternatives by checking types are as expected
		String[] expectedTypes = {"unsignalised", "unsignalised", "unmarked"};
		for (int i=0;i<expectedTypes.length; i++) {
			assert cas.get(i).getType().contentEquals(expectedTypes[i]);
		}
		
		// Repeat with just one road link
		rls = new ArrayList<RoadLink>();
		String[] roadLinkIDs2 = {"A8675945-DE94-4E22-9905-B0623A326221_0"};
		for (RoadLink rl: SpaceBuilder.roadLinkGeography.getAllObjects()) {
			for (int i=0;i<roadLinkIDs2.length; i++) {
				if(rl.getFID().contentEquals(roadLinkIDs2[i])) {
					rls.add(rl);
				}
			}
		}
		
		// Get crossing alternatives
		cas = TacticalRoute.getCrossingAlternatives(SpaceBuilder.caGeography, rls, ped, SpaceBuilder.roadGeography);
		
		String[] expectedTypes2 = {"unsignalised", "unmarked"};
		for (int i=0;i<expectedTypes2.length; i++) {
			assert cas.get(i).getType().contentEquals(expectedTypes2[i]);
		}
	}
	

}
