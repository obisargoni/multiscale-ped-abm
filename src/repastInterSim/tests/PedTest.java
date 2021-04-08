package repastInterSim.tests;

import static org.junit.jupiter.api.Assertions.*;

import java.io.IOException;
import java.util.List;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Polygon;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;

class PedTest {
	
	Context<Object> context = new DefaultContext<Object>();
	
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

	@Test
	void testDistanceToObject() {
		// Setup the environment
		try {
			IO.readProperties();
			SpaceBuilder.fac = new GeometryFactory();
			
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpPedObstructions();

			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp");
			EnvironmentSetup.setUpPavementNetwork();
						
			EnvironmentSetup.setUpPedODs();
			EnvironmentSetup.setUpVehicleODs("mastermap-itn RoadNode Intersect Within.shp");
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Create pedestrian and point it towards it's destination (which in this case is just across the road)
		Ped pedMinDist = EnvironmentSetup.createPedestrian(3,4,false);
		
		// Move ped to position and bearing that has caused an error in the simulation
		Coordinate newLoc = new Coordinate(530512.2083577294, 180906.40573744595);
		pedMinDist.setLoc(newLoc);
		pedMinDist.setBearing(2.8797932657906435);		
		
		// Now that ped is pointiing in a direction can compare distance to nearest object using the full search method and method that searches 
		// just objects in field of vision
		
		// Get filters groups of objects
        Polygon fieldOfVisionApprox = pedMinDist.getPedestrianFieldOfVisionPolygon();
        double start = System.currentTimeMillis();
        Iterable<Object> mobileAgentsInArea = SpaceBuilder.geography.getObjectsWithin(fieldOfVisionApprox.getEnvelopeInternal());
        double durMAAgent = System.currentTimeMillis() - start;
        Iterable<PedObstruction> obstructionsInArea = SpaceBuilder.pedObstructGeography.queryInexact(fieldOfVisionApprox.getEnvelopeInternal());
        double durObstGeog = System.currentTimeMillis() - durMAAgent;
        List<PedObstruction> pedObsInSearchArea = SpatialIndexManager.findIntersectingObjects(SpaceBuilder.pedObstructGeography, fieldOfVisionApprox);
		double durObstSI = System.currentTimeMillis() - durObstGeog;
		
		System.out.print("Duration get mobile agents:" + durMAAgent + "\n");
		System.out.print("Duration get ped obst by geography:" + durObstGeog + "\n");
		System.out.print("Duration get ped obst by spatial index:" + durObstSI + "\n");
		
		// Get distance to nearest object and compare
		List<Double> fovAngles = pedMinDist.sampleFoV();
		
		for(int i=0; i<fovAngles.size(); i++) {
			double d1 = pedMinDist.distanceToObject(fovAngles.get(i));
			double d2 = pedMinDist.distanceToObject(fovAngles.get(i), mobileAgentsInArea, obstructionsInArea);
			double d3 = pedMinDist.distanceToObject(fovAngles.get(i), mobileAgentsInArea, pedObsInSearchArea);
			assert Double.compare(d1, d2) == 0;
			assert Double.compare(d2, d3) == 0;
			assert Double.compare(d3, d1) == 0;
		}
	}

}
