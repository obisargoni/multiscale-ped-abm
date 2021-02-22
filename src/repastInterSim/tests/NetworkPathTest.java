package repastInterSim.tests;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Stack;
import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.Transformer;
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
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.pathfinding.NetworkPath;
import repastInterSim.pathfinding.RoadNetworkRoute;
import repastInterSim.pathfinding.transformers.EdgeRoadLinkIDTransformer;
import repastInterSim.pathfinding.transformers.EdgeWeightTransformer;

class NetworkPathTest {
	
	Geography<RoadLink> roadLinkGeography = null;
	Network<Junction> roadNetwork = null;
	
	Context<Junction> pavementJunctionContext = null;
	Geography<RoadLink> pavementLinkGeography = null;
	Geography<Junction> pavementJunctionGeography = null;
	Network<Junction> pavementNetwork = null;
	
	String testGISDir = ".//data//test_gis_data//";
	String roadLinkPath;
	String pavementLinkPath = null;
	String pedJPath = null;

	void setUpProperties() throws IOException {
		IO.readProperties();
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
	
	void setUpPavementLinks(String linkFile) throws MalformedURLException, FileNotFoundException {
		pavementLinkGeography = setUpLinks(linkFile);
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
	}
	
	void setUpPavementJunctions() throws Exception {
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
	void testPathsBetweenNodePavementJunctions1() {
		try {
			setUpPavementJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Initialise the path finder
		NetworkPath<Junction> np = new NetworkPath<Junction>(this.pavementNetwork);
		
		// Get the start and end nodes to find paths between
		String sourceID = "pave_node_111";
		String targetID = "pave_node_112";
		Junction source = null;
		Junction target = null;
		for (Junction j : this.pavementJunctionGeography.getAllObjects()) {
			if (j.getFID().contentEquals(sourceID)) {
				source = j;
			}
			else if (j.getFID().contentEquals(targetID)) {
				target = j;
			}
		}
		
		String roadJuncID = source.getjuncNodeID();
		
		// Now get just the paths around the road node
		Predicate<Junction> filter = j -> j.getjuncNodeID().contentEquals(roadJuncID);
		List<Stack<RepastEdge<Junction>>> edgePaths = np.getSimplePaths(source, target, filter);
		
		assert edgePaths.size()==2;
		
		// Finally convert these paths to edge paths and check the edges
		String[] expectedLinks = {"pave_link_111_113", "pave_link_113_114", "pave_link_112_114"};
		for (Stack<RepastEdge<Junction>> edgePath:edgePaths) {
			if (edgePath.size()==1) {
				NetworkEdge<Junction> ne = (NetworkEdge<Junction>) edgePath.get(0);
				assert ne.getRoadLink().getFID().contentEquals("pave_link_111_112");
			}
			else {
				for (int i=0; i<edgePath.size(); i++) {
					NetworkEdge<Junction> ne = (NetworkEdge<Junction>)edgePath.get(i);
					assert ne.getRoadLink().getFID().contentEquals(expectedLinks[i]);
				}
			}
		}
	}
	
	@Test
	void testPathToSelf() {
		try {
			setUpPavementJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Initialise the path finder
		NetworkPath<Junction> np = new NetworkPath<Junction>(this.pavementNetwork);
		
		// Get the start and end nodes to find paths between
		String sourceID = "pave_node_114";
		Junction source = null;
		for (Junction j : this.pavementJunctionGeography.getAllObjects()) {
			if (j.getFID().contentEquals(sourceID)) {
				source = j;
			}
		}
		
		String roadJuncID = source.getjuncNodeID();
		
		// Now get just the paths around the road node
		Predicate<Junction> filter = j -> j.getjuncNodeID().contentEquals(roadJuncID);
		List<Stack<RepastEdge<Junction>>> selfPaths = np.getSimplePaths(source, source, filter);
		
		assert selfPaths.size()==1;
		
		// Check that the edge path for this single entry node path is empty
		Stack<RepastEdge<Junction>> edgePath = selfPaths.get(0);
		assert edgePath.size()==0;
	}
	
	@Test
	void testFilterNetwork() {
		try {
			setUpPavementJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		NetworkPath<Junction> np = new NetworkPath<Junction>(this.pavementNetwork);
		
		// Check initial number of nodes
		int nNodes = 0;
		for (Junction n : this.pavementNetwork.getNodes()) {
			nNodes++;
		}
		assert nNodes == 125;
		
		assert np.getGraph().getVertexCount() == 125;
		
		// Now filter to include all road nodes in the strategic path and check number of simple paths returned
		String[] tacticalPavementNodesArray = {"pave_node_79","pave_node_87","pave_node_88","pave_node_80","pave_node_81","pave_node_82","pave_node_83","pave_node_84","pave_node_85","pave_node_86","pave_node_89","pave_node_90","pave_node_91","pave_node_111","pave_node_112","pave_node_113","pave_node_114"};
		List<String> tacticalPavementNodeIDs = Arrays.asList(tacticalPavementNodesArray);
		
		//filter = j -> sPNodes.stream().anyMatch( n -> n.getFID().contentEquals(j.getjuncNodeID()));
		Predicate<Junction> filter = j -> tacticalPavementNodeIDs.stream().anyMatch(n -> n.contentEquals(j.getFID()));
		np.filterGraph(filter);
		assert np.getGraph().getVertexCount() == tacticalPavementNodeIDs.size();
		
		
		// Filter again this time using a collection
		List<Junction> tacticalPavementNodes = new ArrayList<Junction>();
		for (Junction j: this.pavementJunctionGeography.getAllObjects()) {
			if (tacticalPavementNodeIDs.stream().anyMatch(n -> n.contentEquals(j.getFID()))) {
				tacticalPavementNodes.add(j);
			}
		}
		
		assert tacticalPavementNodes.size() == tacticalPavementNodeIDs.size();
		
		np = new NetworkPath<Junction>(this.pavementNetwork);
		np.filterGraph(tacticalPavementNodes);
		assert np.getGraph().getVertexCount() == tacticalPavementNodes.size();
	}
	
	@Test
	void testSimplePathsAlongTacticalHorizon() {
		
		// setup road network and pavement network
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPavementJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String sourceID = "pave_node_85";
		String targetID = "pave_node_112";
		Junction source = null;
		Junction target = null;
		for (Junction j : this.pavementJunctionGeography.getAllObjects()) {
			if (j.getFID().contentEquals(sourceID)) {
				source = j;
			}
			else if (j.getFID().contentEquals(targetID)) {
				target = j;
			}
		}
		
		String startJunctionID = source.getjuncNodeID();
		String endJunctionID = target.getjuncNodeID();
		Coordinate o = null;
		Coordinate d = null;
		List<RoadLink> sP = planStrategicPath(o, d, startJunctionID, endJunctionID);
		
		List<Junction> sPNodes = new ArrayList<Junction>();
		for (RoadLink rl: sP) {
			rl.getJunctions().stream().forEach(sPNodes::add);
		}
		
		NetworkPath<Junction> np = new NetworkPath<Junction>(this.pavementNetwork);

		// First test that when filter excludes target node, zero simple paths returned
		String firstJunctionID = sPNodes.get(0).getFID();
		Predicate<Junction> filter = j -> j.getjuncNodeID().contentEquals(firstJunctionID);
		List<Stack<RepastEdge<Junction>>> paths = np.getSimplePaths(source, target, filter);
		assert paths.size() == 0;
		
		// Now change filter to include all road nodes in the strategic path and check number of simple paths returned
		filter = j -> sPNodes.stream().anyMatch( n -> n.getFID().contentEquals(j.getjuncNodeID()));
		paths = np.getSimplePaths(source, target, filter);
		assert paths.size() == 2076; // python +networkx gives me 2076
	}
	
	/*
	 * Test choosing the shortest of a collection of simple paths based on heruistic.
	 * 
	 * This test uses the edge length cost function and compares shortest path to the Dijkstras shortest path
	 */
	@Test
	public void testGetShortestOfMultiplePaths1() {
		// setup road network and pavement network
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPavementJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String sourceID = "pave_node_85";
		String targetID = "pave_node_112";
		Junction source = null;
		Junction target = null;
		for (Junction j : this.pavementJunctionGeography.getAllObjects()) {
			if (j.getFID().contentEquals(sourceID)) {
				source = j;
			}
			else if (j.getFID().contentEquals(targetID)) {
				target = j;
			}
		}
		
		String startJunctionID = source.getjuncNodeID();
		String endJunctionID = target.getjuncNodeID();
		Coordinate o = null;
		Coordinate d = null;
		List<RoadLink> sP = planStrategicPath(o, d, startJunctionID, endJunctionID);
		
		List<Junction> sPNodes = new ArrayList<Junction>();
		for (RoadLink rl: sP) {
			rl.getJunctions().stream().forEach(sPNodes::add);
		}
		
		NetworkPath<Junction> np = new NetworkPath<Junction>(this.pavementNetwork);
		
		Predicate<Junction> filter = j -> sPNodes.stream().anyMatch( n -> n.getFID().contentEquals(j.getjuncNodeID()));
		List<Stack<RepastEdge<Junction>>> paths = np.getSimplePaths(source, target, filter);
		
		Transformer<RepastEdge<Junction>, Integer> distanceTransformer = new EdgeWeightTransformer<Junction>();
		List<Stack<RepastEdge<Junction>>> shortestPaths = np.getShortestOfMultiplePaths(paths, distanceTransformer);
		
		assert shortestPaths.size() == 1;
		
		Stack<RepastEdge<Junction>> shortestPath = shortestPaths.get(0);
		
		// Check that this gives the same paths as dijkstras on the filtered network
		np.filterGraph(filter);
		List<RepastEdge<Junction>> shortestPathDijkstra = np.getShortestPath(source, target);
		
		for (int i=0; i<Math.max(shortestPath.size(), shortestPathDijkstra.size()); i++) {
			NetworkEdge<Junction> e1 = (NetworkEdge<Junction>) shortestPath.get(i);
			NetworkEdge<Junction> e2 = (NetworkEdge<Junction>) shortestPathDijkstra.get(i);
			assert e1.getRoadLink().getFID().contentEquals(e2.getRoadLink().getFID());
		}
	}
	
	/*
	 * Test choosing the shortest of a collection of simple paths based on heruistic.
	 * 
	 * This test compares shortest paths selected using distance based and crossing based heuristics.
	 */
	@Test
	public void testGetShortestOfMultiplePaths2() {
		// setup road network and pavement network
		try {
			setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
			setUpRoadNetwork(false);
			
			setUpPavementJunctions();
			setUpPavementLinks("pedNetworkLinks.shp");
			setUpPavementNetwork();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String sourceID = "pave_node_87";
		String targetID = "pave_node_112";
		Junction source = null;
		Junction target = null;
		for (Junction j : this.pavementJunctionGeography.getAllObjects()) {
			if (j.getFID().contentEquals(sourceID)) {
				source = j;
			}
			else if (j.getFID().contentEquals(targetID)) {
				target = j;
			}
		}
		
		String startJunctionID = source.getjuncNodeID();
		String endJunctionID = target.getjuncNodeID();
		Coordinate o = null;
		Coordinate d = null;
		List<RoadLink> sP = planStrategicPath(o, d, startJunctionID, endJunctionID);
		
		List<Junction> sPNodes = new ArrayList<Junction>();
		for (RoadLink rl: sP) {
			rl.getJunctions().stream().forEach(sPNodes::add);
		}
		
		NetworkPath<Junction> np = new NetworkPath<Junction>(this.pavementNetwork);
		
		Predicate<Junction> filter = j -> sPNodes.stream().anyMatch( n -> n.getFID().contentEquals(j.getjuncNodeID()));
		List<Stack<RepastEdge<Junction>>> paths = np.getSimplePaths(source, target, filter);
		
		Transformer<RepastEdge<Junction>, Integer> distanceTransformer = new EdgeWeightTransformer<Junction>();
		Transformer<RepastEdge<Junction>, Integer> crossingTransformer = new EdgeRoadLinkIDTransformer<Junction>();
		
		
		List<Stack<RepastEdge<Junction>>> shortestDistancePaths = np.getShortestOfMultiplePaths(paths, distanceTransformer);
		List<Stack<RepastEdge<Junction>>> leastCrossingsPaths = np.getShortestOfMultiplePaths(paths, crossingTransformer);
		
		assert shortestDistancePaths.size() == 1;
		assert leastCrossingsPaths.size() == 1;
		
		Stack<RepastEdge<Junction>> shortestDistancePath = shortestDistancePaths.get(0);
		Stack<RepastEdge<Junction>> leastCrossingPath = leastCrossingsPaths.get(0);
		
		// First, get edges of expected path and check that the expected path is included in the collection of simple paths
		Stack<RepastEdge<Junction>> expectedMinCrossPath = new Stack<RepastEdge<Junction>>();
		String [] expectedEdgesMinCross = {"pave_link_87_81", "pave_link_81_89","pave_link_89_91","pave_link_91_112"};
		for (int i=0; i<expectedEdgesMinCross.length; i++) {
			String rlID = expectedEdgesMinCross[i];
			for (RepastEdge<Junction> re: this.pavementNetwork.getEdges()) {
				NetworkEdge<Junction> e = (NetworkEdge<Junction>) re;
				if (e.getRoadLink().getFID().contentEquals(rlID)) {
					expectedMinCrossPath.add(e);
				}
			}
		}
		
		boolean containsExpected = paths.stream().anyMatch(p -> p.containsAll(expectedMinCrossPath));
		assert containsExpected;
		
		// Now compare path length of the expected path to path length of returned shortest path
		int expectedLength = NetworkPath.getIntPathLength(expectedMinCrossPath, crossingTransformer);
		int returnedLength = NetworkPath.getIntPathLength(leastCrossingPath, crossingTransformer);
		assert expectedLength==returnedLength;
		
        // But rest of path will differ
        String [] expectedNodesMinDist = {"pave_node_87", "pave_node_81", "pave_node_90", "pave_node_112"};
        String [] expectedNodesMinCross = {"pave_node_87", "pave_node_81", "pave_node_89", "pave_node_91", "pave_node_112"};

        List<Junction> pathNodesMinDist = np.nodePathFromEdges(shortestDistancePath, source);
        List<Junction> pathNodesMinCross = np.nodePathFromEdges(leastCrossingPath, source);
		
		for (int i=0; i<Math.max(expectedNodesMinDist.length, pathNodesMinDist.size()); i++) {
			assert pathNodesMinDist.get(i).getFID().contentEquals(expectedNodesMinDist[i]);
		}
		
		for (int i=0; i<Math.max(expectedNodesMinCross.length, pathNodesMinCross.size()); i++) {
			assert pathNodesMinCross.get(i).getFID().contentEquals(expectedNodesMinCross[i]);
		}
		
	}
}
