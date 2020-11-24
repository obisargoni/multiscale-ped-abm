package repastInterSim.tests;

import static org.junit.jupiter.api.Assertions.*;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.List;
import java.util.Stack;
import java.util.function.Predicate;

import org.junit.jupiter.api.Test;

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

class NetworkPathTest {
	
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
	
	void setUpPavementLinks(String linkFile) throws MalformedURLException, FileNotFoundException {
		pavementLinkGeography = setUpLinks(linkFile);
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

}
