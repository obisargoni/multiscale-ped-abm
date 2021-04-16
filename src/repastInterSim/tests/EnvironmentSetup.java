package repastInterSim.tests;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
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
import repastInterSim.agent.Vehicle;
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
import repastInterSim.environment.contexts.PedObstructionPointsContext;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.RoadNetworkRoute;

public class EnvironmentSetup {
	
	static String testGISDir = ".//data//test_gis_data//";
	static String pedestrianRoadsPath = testGISDir + "topographicAreaPedestrian.shp";
	static String vehicleRoadsPath = testGISDir + "topographicAreaVehicle.shp";
	static String serialisedLookupPath = testGISDir + "road_link_roads_cache.serialised";
	
	public EnvironmentSetup() {
		// TODO Auto-generated constructor stub
	}
	
	static void setUpProperties() throws IOException {
		IO.readProperties();
		EnvironmentSetup.clearCaches();
		SpaceBuilder.fac = new GeometryFactory();
	}
	
	static void setUpObjectGeography() {
		SpaceBuilder.context = new DefaultContext<Object>();
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		SpaceBuilder.geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, SpaceBuilder.context, geoParams);
		SpaceBuilder.geography.setCRS(GlobalVars.geographyCRSString);
		SpaceBuilder.context.add(SpaceBuilder.geography);
	}
	
	static void setUpRoads() throws Exception {
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
	}
	
	static Geography<RoadLink> setUpRoadLinks(String roadLinkFile) throws MalformedURLException, FileNotFoundException {
		String roadLinkPath = testGISDir + roadLinkFile;
		
		// Initialise test road link geography and context
		Context<RoadLink> roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> rlG = GeographyFactoryFinder.createGeographyFactory(null).createGeography("orRoadLinkGeography", roadLinkContext, GeoParams);
		rlG.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(RoadLink.class, roadLinkPath, rlG, roadLinkContext);
		SpatialIndexManager.createIndex(rlG, RoadLink.class);
		
		return rlG;
	}
	
	static void setUpORRoadLinks() throws Exception {
		SpaceBuilder.orRoadLinkGeography = setUpRoadLinks("open-roads RoadLink Intersect Within simplify angles.shp");
	}
	
	static void setUpITNRoadLinks() throws Exception {
		SpaceBuilder.roadLinkGeography = setUpRoadLinks("mastermap-itn RoadLink Intersect Within with orientation.shp");
	}
	
	static void setUpPavementLinks(String linkFile) throws MalformedURLException, FileNotFoundException {
		SpaceBuilder.pavementLinkGeography = setUpRoadLinks(linkFile);
	}
		
	static Geography<OD> setUpODs(String odFile, String geogName) throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> ODContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> odG = GeographyFactoryFinder.createGeographyFactory(null).createGeography(geogName, ODContext, GeoParamsOD);
		odG.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String testODFile = testGISDir + odFile;
		GISFunctions.readShapefile(OD.class, testODFile, odG, ODContext);
		SpatialIndexManager.createIndex(odG, OD.class);
		return odG;		
	}
	
	static void setUpPedODs() throws MalformedURLException, FileNotFoundException {
		SpaceBuilder.pedestrianDestinationGeography = setUpODs("OD_pedestrian_nodes_test.shp", "pedODGeography");
	}
	
	static void setUpVehicleODs() throws MalformedURLException, FileNotFoundException {
		SpaceBuilder.vehicleDestinationGeography = setUpODs("OD_vehicle_nodes_intersect_within.shp", "vehicleODGeography");		
	}
	
	static void setUpVehicleODs(String file) throws MalformedURLException, FileNotFoundException {
		SpaceBuilder.vehicleDestinationGeography = setUpODs(file, "vehicleODGeography");		
	}
	
	static void setUpPedObstructions() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<PedObstruction> pedObstructContext = new PedObstructionContext();
		GeographyParameters<PedObstruction> GeoParams = new GeographyParameters<PedObstruction>();
		SpaceBuilder.pedObstructGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedObstructGeography", pedObstructContext, GeoParams);
		SpaceBuilder.pedObstructGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testPedObstructFile = testGISDir + IO.getProperty("PedestrianObstructionShapefile");
		GISFunctions.readShapefile(PedObstruction.class, testPedObstructFile, SpaceBuilder.pedObstructGeography, pedObstructContext);
		SpatialIndexManager.createIndex(SpaceBuilder.pedObstructGeography, PedObstruction.class);
	}
	
	static void setUpPedObstructionPoints() throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<PedObstruction> pedObstructionPointsContext = new PedObstructionPointsContext();
		GeographyParameters<PedObstruction> GeoParams = new GeographyParameters<PedObstruction>();
		SpaceBuilder.pedObstructionPointsGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedObstructionPointsGeography", pedObstructionPointsContext, GeoParams);
		SpaceBuilder.pedObstructionPointsGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testPedObstructFile = testGISDir + IO.getProperty("PedestrianObstructionPointsShapefile");
		GISFunctions.readShapefile(PedObstruction.class, testPedObstructFile, SpaceBuilder.pedObstructionPointsGeography, pedObstructionPointsContext);
		SpatialIndexManager.createIndex(SpaceBuilder.pedObstructionPointsGeography, PedObstruction.class);
	}
	
	static void setUpCrossingAlternatives(String caFile) throws MalformedURLException, FileNotFoundException {
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		Context<CrossingAlternative> caContext = new CAContext();
		GeographyParameters<CrossingAlternative> GeoParams = new GeographyParameters<CrossingAlternative>();
		SpaceBuilder.caGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("caGeography", caContext, GeoParams);
		SpaceBuilder.caGeography.setCRS(GlobalVars.geographyCRSString);
		
		
		// Load ped obstructions data
		String testCAFile = testGISDir + caFile;
		GISFunctions.readShapefile(CrossingAlternative.class, testCAFile, SpaceBuilder.caGeography, caContext);
		SpatialIndexManager.createIndex(SpaceBuilder.caGeography, CrossingAlternative.class);
	}
	
	static Network<Junction> setUpRoadNetwork(boolean isDirected, Geography<RoadLink> rlG, String name) {
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		Geography<Junction> junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", junctionContext, GeoParamsJunc);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);
		
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(name,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		Network<Junction> rN = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(rlG, junctionContext, junctionGeography, rN);
		
		return rN;
	}
	
	static void setUpORRoadNetwork(boolean isDirected) {
		
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		SpaceBuilder.orJunctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.OR_JUNCTION_GEOGRAPHY, junctionContext, GeoParamsJunc);
		SpaceBuilder.orJunctionGeography.setCRS(GlobalVars.geographyCRSString);
		
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.OR_ROAD_NETWORK,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		SpaceBuilder.orRoadNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(SpaceBuilder.orRoadLinkGeography, junctionContext, SpaceBuilder.orJunctionGeography, SpaceBuilder.orRoadNetwork);
	}
	
	static void setUpITNRoadNetwork(boolean isDirected) {
		SpaceBuilder.roadNetwork = setUpRoadNetwork(isDirected, SpaceBuilder.roadLinkGeography, GlobalVars.CONTEXT_NAMES.ROAD_NETWORK);
	}
	
	static void setUpPavementNetwork() {
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>("PAVEMENT_NETWORK", SpaceBuilder.pavementJunctionContext, false);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		SpaceBuilder.pavementNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(SpaceBuilder.pavementLinkGeography, SpaceBuilder.pavementJunctionContext, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementNetwork);
	}
	
	static void setUpPedJunctions() throws Exception {
		String pedJPath = testGISDir + IO.getProperty("PavementJunctionsShapefile");
		
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
		for(Junction j: SpaceBuilder.orRoadNetwork.getNodes()) {
			
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
			shortestRoute = rnr.getShortestRoute(SpaceBuilder.orRoadNetwork, currentJunctions, destJunctions, routeEndpoints, false);
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
	 * Code taken from SpaceBuilder that connects Roads with Road Links and visa versa.
	 */
	static void assocaiteRoadsWithRoadLinks() {
		// Link road with itn and OR road links
		// Also assigns the Road objects to the road links. This enables lookups between OR and ITN road links, through the road objects.
		for (Road r: SpaceBuilder.roadGeography.getAllObjects()) {
			List<RoadLink> roadLinks = new ArrayList<RoadLink>();
			for(RoadLink rl: SpaceBuilder.roadLinkGeography.getAllObjects()) {
				// Iterating over the vehicle road links (ITN) but using their corresponding ped road link (open road) id to check whether they belong to this vehicle polygon
				if (rl.getPedRLID().contentEquals(r.getRoadLinkID())) {
					roadLinks.add(rl);
					rl.getRoads().add(r);
				}
			}
			
			RoadLink orLink = null;
			for(RoadLink rl: SpaceBuilder.orRoadLinkGeography.getAllObjects()) {
				// Iterating over the vehicle road links (ITN) but using their corresponding ped road link (open road) id to check whether they belong to this vehicle polygon
				if (rl.getFID().contentEquals(r.getRoadLinkID())) {
					orLink = rl;
					orLink.getRoads().add(r);
					break;
				}
			}
			
			r.setRoadLinks(roadLinks);
			r.setORRoadLink(orLink);
		}
	}
	
	static Vehicle createVehicle(OD o, OD d) {
		Vehicle V = new Vehicle(GlobalVars.maxVehicleSpeed, GlobalVars.defaultVehicleAcceleration, GlobalVars.initialVehicleSpeed, o, d);
		SpaceBuilder.context.add(V);
		Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		GeometryFactory fac = new GeometryFactory();
		Point pt = fac.createPoint(oCoord);
		Geometry vehicleCircle = pt.buffer(2);
		GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, vehicleCircle, V);
		V.setLoc();
		return V;
	}
	
	static Vehicle createVehicle(String oID, String dID) {
		OD o = null;
		OD d = null;
		for (OD od: SpaceBuilder.vehicleDestinationGeography.getAllObjects()) {
			if (od.getFID().contentEquals(oID)) {
				o = od;
			}
			else if (od.getFID().contentEquals(dID)) {
				d = od;
			}
		}
		return createVehicle(o,d);
	}
	
	static Ped createPedestrian(int oID, int dID, boolean minimisesCrossing) {
		double s = 0.0;
		double m = 0.0;
		double alpha = 0.5;
		double lambda = 1.0;
		double gamma = 0.9;
		double epsilon = 3.0;
		double pH = 20.0;
		
		return createPedestrian(oID, dID, s, m, alpha, lambda, gamma, epsilon, minimisesCrossing, pH);
	}
	
	static Ped createPedestrian(int oID, int dID, Double s, Double m, Double alpha, Double lambda, Double gamma, Double epsilon, boolean minimiseCrossings, Double pH) {
		
		OD o = null;
		OD d = null;
		
		for (OD i : SpaceBuilder.pedestrianDestinationGeography.getAllObjects()) {
			if (i.getId() == oID) {
				o = i;
			}
			else if (i.getId() == dID) {
				d = i;
			}
		}
		
		Ped p = new Ped(o, d, s, m, alpha, lambda, gamma, epsilon, minimiseCrossings, pH, SpaceBuilder.pavementJunctionGeography, SpaceBuilder.pavementNetwork);

        SpaceBuilder.context.add(p);        
        Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		Point pt = GISFunctions.pointGeometryFromCoordinate(oCoord);
		Geometry circle = pt.buffer(p.getRad());		
		GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, circle, p);
		p.setLoc();
		return p;
	}
	
	public static void clearCaches() {
		RoadNetworkRoute.clearCaches();
		SpatialIndexManager.clearCaches();
	}

}
