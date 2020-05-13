package tests;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;

import pathfinding.RoadNetworkRoute;
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
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.environment.contexts.VehicleDestinationContext;
import repastInterSim.main.GlobalVars;

public class RoadNetworkRouteTest {
	
	private RoadNetworkRoute rnr;
	
	private Boolean isDirected = true;
	
	private Context<RoadLink> roadLinkContext;
	
	private Context<OD> testODContext;
	
	private Context<Junction> junctionContext;
	private Network<Junction> roadNetwork;
	
	@BeforeEach
    public void setUp() throws Exception {
		
	    // Load agents from shapefiles
		String TestDataDir = ".//data//test_gis_data//";
		
		roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
		
		junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		Geography<Junction> junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", junctionContext, GeoParamsJunc);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);
		
		testODContext = new VehicleDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> testODGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", testODContext, GeoParamsOD);
		testODGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// 1. Load the road links
		String roadLinkFile = TestDataDir + "mastermap-itn RoadLink Intersect Within with orientation.shp";
		GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
		String vehicleDestinationsFile = TestDataDir + "OD_vehicle_nodes_intersect_within.shp";
		GISFunctions.readShapefile(OD.class, vehicleDestinationsFile, testODGeography, testODContext);
		
		// 2. roadNetwork
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		roadNetwork = builder.buildNetwork();
		GISFunctions.buildGISRoadNetwork(roadLinkGeography, junctionContext,junctionGeography, roadNetwork);
		
		
		// Initialise the road network route
		Coordinate o = testODContext.getObjects(OD.class).get(10).getGeom().getCoordinate();
		Coordinate d = testODContext.getObjects(OD.class).get(6).getGeom().getCoordinate();

		rnr = new RoadNetworkRoute(o , d);
    }

	
	@Test
	public void testGetShortestRoute() {
		
		// Define my starting and ending junctions to test
		List<Junction> currentJunctions = new ArrayList<Junction>();
		List<Junction> destJunctions = new ArrayList<Junction>();
		for(Junction j: junctionContext.getObjects(Junction.class)) {
			
			// Set the test current junctions 
			if (j.getFID().contentEquals("osgb4000000029971605")) {
				currentJunctions.add(j);
			}
			
			// Set the test destination junctions
			if (j.getFID().contentEquals("osgb4000000029970678")) {
				destJunctions.add(j);
			}
		}
		
		Junction[] routeEndpoints = new Junction[2];
		
		// Get shortest Route according to the Route class
		List<RepastEdge<Junction>> shortestRoute = null;
		try {
			shortestRoute = rnr.getShortestRoute(roadNetwork, currentJunctions, destJunctions, routeEndpoints);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
		// Assert this matches what I expect - what I expect
		String[] expectedRoadLinkRouteFIDs = {"osgb4000000030419069", "osgb4000000030344373", "osgb5000005123276401",
												"osgb5000005123276405", "osgb4000000030343876", "osgb5000005156082205", 
												"osgb4000000030343920"};
		
		for(int i=0; i< shortestRoute.size(); i++) {
			NetworkEdge<Junction> e = (NetworkEdge<Junction>)shortestRoute.get(i);
			String expectedFID = expectedRoadLinkRouteFIDs[i];
			assert e.getRoadLink().getFID().contentEquals(expectedFID);
		}
		
	}

}
