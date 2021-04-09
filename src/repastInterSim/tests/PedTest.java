package repastInterSim.tests;

import java.io.IOException;
import java.util.List;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Polygon;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repastInterSim.agent.Ped;
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
		double b = 2.8797932657906435;
		pedMinDist.setBearing(b);	
		
		// Get distance to nearest object and compare
		List<Double> fovAngles = pedMinDist.sampleFoV();
		
		// Get filters groups of objects
        Polygon fieldOfVisionApprox = pedMinDist.getPedestrianFieldOfVisionPolygon(b);
        List<Geometry> obstGeoms = pedMinDist.getObstacleGeometries(fieldOfVisionApprox);
		
		// Now that ped is pointiing in a direction can compare distance to nearest object using the full search method and method that searches 
		// just objects in field of vision
		
		// Speed test difference between original distanceToObject method and methods that uses pre-filtered obstacle geometries.
		
        double start = System.currentTimeMillis();
        double[] d1s = new double[fovAngles.size()];
		for(int i=0; i<fovAngles.size(); i++) {
			d1s[i] = pedMinDist.distanceToObject(fovAngles.get(i));
		}
        double durObst = System.currentTimeMillis() - start;
        
        double[] d2s = new double[fovAngles.size()];
		for(int i=0; i<fovAngles.size(); i++) {
			d2s[i] = pedMinDist.distanceToObject(fovAngles.get(i), obstGeoms);
		}
        double durObstGeog = System.currentTimeMillis() - (start + durObst);
        
		System.out.print("Duration get ped obst all:" + durObst + "\n");
		System.out.print("Duration get ped obst by geography:" + durObstGeog + "\n");
		
		// Now check that all distances returned match
		for(int i=0; i<fovAngles.size(); i++) {
			assert Double.compare(d1s[i], d2s[i]) == 0;
		}

	}

}
