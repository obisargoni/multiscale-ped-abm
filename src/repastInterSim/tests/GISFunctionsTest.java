package repastInterSim.tests;


import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;

import repast.simphony.context.Context;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.OD;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;

class GISFunctionsTest {
	
	private String TestDataDir = ".//data//test_gis_data//";
	
	private Context<RoadLink> roadLinkContext;
	
	private Context<OD> pedestrianDestinationContext;

	void setUp(String lineDataFile) throws Exception {
		
	    // Initialise contexts and geographies used by all tests	
		roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("roadLinkGeography", roadLinkContext, GeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);
		
		pedestrianDestinationContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		Geography<OD> pedestrianDestinationGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("pedestrianDestinationGeography", pedestrianDestinationContext, GeoParamsOD);
		pedestrianDestinationGeography.setCRS(GlobalVars.geographyCRSString);
		
		// 1. Load road network data
		String roadLinkFile = TestDataDir + lineDataFile;
		GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, roadLinkContext);
		SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);
		
		String testODFile = TestDataDir + "parity_test_OD.shp";
		GISFunctions.readShapefile(OD.class, testODFile, pedestrianDestinationGeography, pedestrianDestinationContext);
		SpatialIndexManager.createIndex(pedestrianDestinationGeography, OD.class);
	}
	
	int getNumberIntersectingCoords(String lineDataFile) {
		try {
			setUp(lineDataFile);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get number of intersecting points and compare to expected value
		Coordinate[] odCoords = new Coordinate[pedestrianDestinationContext.getObjects(OD.class).size()];
		int i = 0;
		for (OD od : pedestrianDestinationContext.getObjects(OD.class)) {
			odCoords[i] = od.getGeom().getCoordinate();
			i++;
		}
		LineString ODLine = GISFunctions.lineStringGeometryFromCoordinates(odCoords);
		Geometry[] ODLineGeom = {ODLine};

		
		// Loop through road links in the rout and count number of times the ODLine intersects
		Geometry[] rlGeoms = new Geometry[roadLinkContext.getObjects(RoadLink.class).size()];
		i = 0;
		for (RoadLink rl: roadLinkContext.getObjects(RoadLink.class)) {
			rlGeoms[i] = rl.getGeom();
			i++;
		}
		
		// Now calculate number of intersecting coordinates
		int nIntersections = GISFunctions.calculateNIntersectionCoords(ODLineGeom, rlGeoms);
		
		return nIntersections;
	}


	@Test
	void testCalculateNIntersectionCoords() {
		
		// Now calculate number of intersecting coordinates
		int nIntersections1 = getNumberIntersectingCoords("parity_test_lines1.shp");
		int nIntersections2 = getNumberIntersectingCoords("parity_test_lines2.shp");
		int nIntersections3 = getNumberIntersectingCoords("parity_test_lines3.shp");
		
		assert nIntersections1 == 1;
		assert nIntersections2 == 2;
		assert nIntersections3 == 3;
	}

}
