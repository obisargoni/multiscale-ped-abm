package repastInterSim.tests;

import static org.junit.jupiter.api.Assertions.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.context.Context;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.OD;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
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
		String testODFile = testGISDir + "OD_pedestrian_nodes.shp";
		GISFunctions.readShapefile(OD.class, testODFile, odGeography, ODContext);
		
	}
}
