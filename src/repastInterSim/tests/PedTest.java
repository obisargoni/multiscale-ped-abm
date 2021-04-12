package repastInterSim.tests;

import java.io.IOException;
import java.util.List;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repastInterSim.agent.Ped;
import repastInterSim.environment.GISFunctions;
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

	
	void testDistanceToObjectFromCoordBreaing(Coordinate c, double b) {
		// Setup the environment
		try {
			IO.readProperties();
			SpaceBuilder.fac = new GeometryFactory();
			
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpPedObstructions();
			EnvironmentSetup.setUpPedObstructionPoints();

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
        Point pt = GISFunctions.pointGeometryFromCoordinate(c);
		Geometry pGeomNew = pt.buffer(pedMinDist.getRad());
        GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, pGeomNew, pedMinDist);
		pedMinDist.setLoc();
		pedMinDist.setBearing(b);	
		
		// Get distance to nearest object and compare
		List<Double> fovAngles = pedMinDist.sampleFoV();
		
		// Get filters groups of objects
        Polygon fieldOfVisionApprox = pedMinDist.getPedestrianFieldOfVisionPolygon(b);
        List<Geometry> obstGeomsPoints = pedMinDist.getObstacleGeometries(fieldOfVisionApprox, SpaceBuilder.pedObstructionPointsGeography);
        
		// Now that ped is pointiing in a direction can compare distance to nearest object using the full search method and method that searches 
		// just objects in field of vision
		
		// Speed test difference between original distanceToObject method and methods that uses pre-filtered obstacle geometries.
		
        double start = System.currentTimeMillis();
        double[] d1s = new double[fovAngles.size()];
		for(int i=0; i<fovAngles.size(); i++) {
			d1s[i] = pedMinDist.distanceToObject(fovAngles.get(i));
		}
        double durObst = System.currentTimeMillis() - start;
        
    	// Then sample field of vision and initialise arrays of distances that correspond to each angle
        double[] d2s = new double[fovAngles.size()-1];
    	double[] displacementDistances = new double[fovAngles.size()-1];
    	
    	// Initialise values as -1 so can identify default values
    	for (int i=0; i<d2s.length; i++) {
    		d2s[i] = -1;
    		displacementDistances[i] = -1;
    	}
    	
    	// First find the minimum displacement distance and associated angle for obstruction geometries
    	pedMinDist.dispalcementDistancesToGeometries(obstGeomsPoints, fovAngles, d2s, displacementDistances);
    	
    	double durObstGeog = System.currentTimeMillis() - (start + durObst);
        
		System.out.print("Duration get ped obst all:" + durObst + "\n");
		System.out.print("Duration get ped obst by geography:" + durObstGeog + "\n");
		
		// Now check that all distances returned match
		int count = 0;
		for(int i=0; i<fovAngles.size()-1; i++) {
			// Only compare for angles where there is an obstacle
			if ( (d2s[i] < 0) | (d1s[i] == 10.0)) {
				continue;
			}
			double diff = Math.abs(d1s[i] - d2s[i]);
			//assert diff < 1;
			System.out.print(i+", "+diff + "\n");
			count++;
		}
		System.out.print("Count: "+count);

	}
	
	void testDistanceToObject1() {
		Coordinate c = new Coordinate(530512.2083577294, 180906.40573744595);
		double b = 2.8797932657906435;
		
		testDistanceToObjectFromCoordBreaing(c, b);
	}
	
	@Test
	void testDistanceToObject2() {
		Coordinate c = new Coordinate(530510.0, 180909.0);
		double b = (3*Math.PI)/2;
		
		testDistanceToObjectFromCoordBreaing(c, b);
	}
	
	@Test
	void testDistanceToObject3() {
		Coordinate c = new Coordinate(530509.6389832983, 180908.11179611267);
		double b = 3.9209401504483683;
		
		testDistanceToObjectFromCoordBreaing(c, b);
	}
	
	// Actually need to compare chosen directions
	@Test
	void testChosenDirection1() {
		
	}

}
