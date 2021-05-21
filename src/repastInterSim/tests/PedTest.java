package repastInterSim.tests;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.operation.distance.DistanceOp;

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
	
	Ped createPedAtLocation(boolean minimisesDistance, Coordinate c, double b) {
		// Create pedestrian and point it towards it's destination (which in this case is just across the road)
		Ped ped = EnvironmentSetup.createPedestrian(3,4,minimisesDistance);
		
		// Move ped to position and bearing that has caused an error in the simulation
        Point pt = GISFunctions.pointGeometryFromCoordinate(c);
		Geometry pGeomNew = pt.buffer(ped.getRad());
        GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, pGeomNew, ped);
		ped.setLoc();
		ped.setBearing(b);
		return ped;
	}

	
	public HashMap<String, double[]> wrapperDispalcementDistancesToGeometries(Coordinate c, double b, boolean intersects) {
		// Setup the environment
		try {
			EnvironmentSetup.setUpProperties();
			
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
        List<Geometry> obstGeomsPoints = null;		
        
    	// Then sample field of vision and initialise arrays of distances that correspond to each angle
        double[] d2s = null;
    	double[] displacementDistances = null;
    	
    	// First find the minimum displacement distance and associated angle for obstruction geometries
    	double start = System.currentTimeMillis();
    	if(intersects) {
    		obstGeomsPoints = pedMinDist.getObstacleGeometries(fieldOfVisionApprox, SpaceBuilder.pedObstructGeography);
            d2s = new double[fovAngles.size()];
        	displacementDistances = new double[fovAngles.size()];
    		pedMinDist.displacementDistancesToObstacleGeometries(obstGeomsPoints, fovAngles, d2s, displacementDistances);
    	}
    	else {
    		obstGeomsPoints = pedMinDist.getObstacleGeometries(fieldOfVisionApprox, SpaceBuilder.pedObstructionPointsGeography);
            d2s = new double[fovAngles.size()-1];
        	displacementDistances = new double[fovAngles.size()-1];
        	
        	// Initialise values as -1 so can identify default values
        	for (int i=0; i<d2s.length; i++) {
        		d2s[i] = -1;
        		displacementDistances[i] = -1;
        	}
    		pedMinDist.dispalcementDistancesToPointGeometries(obstGeomsPoints, fovAngles, d2s, displacementDistances);
    	}
    	
    	double durObstGeog = System.currentTimeMillis() - start;
        System.out.print("Duration get ped obst by geography:" + durObstGeog + "\n");
        
        HashMap<String, double[]> output = new HashMap<String, double[]>();
        double[] angles = new double[fovAngles.size()];
        for (int i=0; i<fovAngles.size(); i++) {
        	angles[i] = fovAngles.get(i);
        }
        output.put("angles", angles);
        output.put("distances", d2s);
        output.put("dispDistances", displacementDistances);
        
        return output;
	}
	
	void validateOutput(HashMap<String, double[]> output, double[] expectedAngles, double[] expectedDistances) {
		for (int i=0; i< Math.max(expectedAngles.length, output.get("angles").length); i++) {
			assert Double.compare(output.get("angles")[i], expectedAngles[i]) == 0;
		}
		
		for (int i=0; i< Math.max(expectedDistances.length, output.get("distances").length); i++) {
			assert Math.abs(output.get("distances")[i] - expectedDistances[i]) < 0.000001;
		}
	}
	
	/*
	 * Test that ped far from any walls doesn't identify obstacle geometries
	 */
	@Test
	void testDistanceToObject1() {
		Coordinate c = new Coordinate(530522.0, 180918);
		double b = 3.9209401504483683;
		
		HashMap<String, double[]> output = wrapperDispalcementDistancesToGeometries(c, b, false);
		
		// Set up expected angles
		List<Double> eA = new ArrayList<Double>();
		double theta = (2*Math.PI*75) / 360; // Copied from ped init function
		double angRes = (2*Math.PI) / (36 / 3); // Also copied from init function
		double a = b - theta;
		while (a <= b+theta) {
			eA.add(a);
			a+=angRes;
		}
		
		double[] expectedAngles = new double[eA.size()];
		for(int i=0; i<expectedAngles.length; i++) {
			expectedAngles[i] = eA.get(i);
		}
		
		double[] expectedDistances = new double[eA.size()-1];
		for(int i=0; i<expectedDistances.length; i++) {
			expectedDistances[i] = -1.0;
		}
		
		validateOutput(output, expectedAngles, expectedDistances);
	}
	
	/*
	 * Test that ped far from any walls doesn't identify obstacle geometries.
	 * 
	 * Test the intersect method.
	 */
	@Test
	void testDistanceToObject1b() {
		Coordinate c = new Coordinate(530522.0, 180918);
		double b = 3.9209401504483683;
		
		HashMap<String, double[]> output = wrapperDispalcementDistancesToGeometries(c, b, true);
		
		// Set up expected angles
		List<Double> eA = new ArrayList<Double>();
		double theta = (2*Math.PI*75) / 360; // Copied from ped init function
		double angRes = (2*Math.PI) / (36 / 3); // Also copied from init function
		double a = b - theta;
		while (a <= b+theta) {
			eA.add(a);
			a+=angRes;
		}
		
		double[] expectedAngles = new double[eA.size()];
		for(int i=0; i<expectedAngles.length; i++) {
			expectedAngles[i] = eA.get(i);
		}
		
		double[] expectedDistances = new double[eA.size()];
		for(int i=0; i<expectedDistances.length; i++) {
			expectedDistances[i] = 10.0;
		}
		
		validateOutput(output, expectedAngles, expectedDistances);
	}
	
	/*
	 * Test ped next to wall identifies obstacles points at expected angles and distances.
	 */
	@Test
	void testDistanceToObject2() {
		Coordinate c = new Coordinate(530509.6389832983, 180908.11179611267);
		double b = 3.9209401504483683;
		
		HashMap<String, double[]> output = wrapperDispalcementDistancesToGeometries(c, b, false);
		
		double[] expectedAngles = {2.611943211452621, 3.6429872589076853, 4.143583854519489, 4.520772370555282, 4.765214246949521, 5.229937089444116};
		double[] expectedDistances = {-1.0, 0.5471421433667476, 0.2539724955402363, 0.20371170329094668, 0.1932528727861567};
		
		for (int i=0; i< Math.max(expectedAngles.length, output.get("angles").length); i++) {
			assert Double.compare(output.get("angles")[i], expectedAngles[i]) == 0;
		}
		
		for (int i=0; i< Math.max(expectedDistances.length, output.get("distances").length); i++) {
			assert Double.compare(output.get("distances")[i], expectedDistances[i]) == 0;
		}
	}
	
	/*
	 * Test ped next to wall identifies obstacles points at expected angles and distances.
	 * 
	 * Use the intersect method.
	 */
	@Test
	void testDistanceToObject2b() {
		Coordinate c = new Coordinate(530509.6389832983, 180908.11179611267);
		double b = 3.9209401504483683;
		
		HashMap<String, double[]> output = wrapperDispalcementDistancesToGeometries(c, b, true);
		
		// Set up expected angles
		List<Double> eA = new ArrayList<Double>();
		double theta = (2*Math.PI*75) / 360; // Copied from ped init function
		double angRes = (2*Math.PI) / (36 / 3); // Also copied from init function
		double a = b - theta;
		while (a <= b+theta) {
			eA.add(a);
			a+=angRes;
		}
		
		double[] expectedAngles = new double[eA.size()];
		for(int i=0; i<expectedAngles.length; i++) {
			expectedAngles[i] = eA.get(i);
		}
		
		double[] expectedDistances = {10.0, 10.0, 0.33716113674871473, 0.0583762406829704, 0.007112777681278104, 0.019415176460662875};
		
		validateOutput(output, expectedAngles, expectedDistances);
	}
	
	@Test
	void testDistanceOpDetectsContactWithPed() {
		
		// Create three pedestrians,two that intersect and one that doesnt.
		Coordinate c1 = new Coordinate(530509.6389832983, 180908.11179611267);
		Point p1 = GISFunctions.pointGeometryFromCoordinate(c1);
		double r1 = 0.2;
		Geometry g1 = p1.buffer(r1);
		
		Coordinate c2 = new Coordinate(530509.846, 180908.255);
		Point p2 = GISFunctions.pointGeometryFromCoordinate(c2);
		double r2 = 0.5;
		Geometry g2 = p2.buffer(r2);

		
		Coordinate c3 = new Coordinate(530512, 180907);
		Point p3 = GISFunctions.pointGeometryFromCoordinate(c3);
		double r3 = 0.5;
		Geometry g3 = p3.buffer(r3);
		
		Coordinate c4 = new Coordinate(530513, 180907);
		Point p4 = GISFunctions.pointGeometryFromCoordinate(c4);
		double r4 = 0.5;
		Geometry g4 = p4.buffer(r4);
		
		// Now check distOp and intersects
		DistanceOp dist12 = new DistanceOp(g1, g2);
		assert (dist12.distance()==0) & (g1.intersects(g2));
		
		DistanceOp dist13 = new DistanceOp(g1, g3);
		assert (dist13.distance()>0) & (g1.intersects(g3)==false);
		
		// g3 and g4 should just touch, does identify as touching
		DistanceOp dist34 = new DistanceOp(g3, g4);
		assert (dist34.distance()==0) & (g3.intersects(g4)==true);
		
	}
	
	@Test
	void testcontactDistancePedGeomWithPedObstruction() {
		
		// Create linestring object that represents a wall
		Coordinate lc1 = new Coordinate(530511, 180900);
		Coordinate lc2 = new Coordinate(530511, 180910);
		Coordinate[] lcs = {lc1, lc2};
		Geometry l = GISFunctions.lineStringGeometryFromCoordinates(lcs);

		
		Coordinate c3 = new Coordinate(530512, 180907);
		Point p3 = GISFunctions.pointGeometryFromCoordinate(c3);
		double r3 = 0.5;
		Geometry g3 = p3.buffer(r3);
		
		Coordinate c4 = new Coordinate(530511.3, 180907);
		Point p4 = GISFunctions.pointGeometryFromCoordinate(c4);
		double r4 = 0.5;
		Geometry g4 = p4.buffer(r4);
		
		// Now check distOp and intersects
		DistanceOp distl3 = new DistanceOp(l, g3);
		assert (distl3.distance()==0.5) & (l.intersects(g3)==false);
		
		// Check dist op identifies the same coordiante that intersection does
		Coordinate expectedC = new Coordinate(530511.0, 180906.60249647763);
		DistanceOp distl4 = new DistanceOp(l, g4);
		assert (distl4.distance()==0) & (l.intersects(g4)==true);
		Coordinate i = distl4.nearestPoints()[0];
		assert  (Math.abs(i.x - expectedC.x) < 0.0001) & (Math.abs(i.y - expectedC.y) < 0.0001);
		
		// Also check that intersecting coords is as expected
		Coordinate[] intCoords = l.intersection(g4).getCoordinates();
		assert intCoords.length==2;
		assert (intCoords[0].equals2D(i)) | (intCoords[1].equals2D(i));
		
		// And mid point between intersecting coords is as expected
		Coordinate expectedMid = new Coordinate(530511.0, 180907.0);
		Coordinate mid = GISFunctions.midwayBetweenTwoCoordinates(intCoords[0], intCoords[1]);
		assert mid.equals2D(expectedMid);
		
	}
	
	
	@Test
	void testPedsPassingContactAccel() {
		
		// Setup the environment
		try {
			EnvironmentSetup.setUpProperties();
			
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
		
		Ped p17 = EnvironmentSetup.createPedestrian(1, 3, 0.8098254410639371, 56.15750059331304, 0.5, 0.1, 0.9, 4.0, true, 20.0);
		assert (p17.getRad() - 0.17549218935410324) < 0.00001;
		// Step to initialise
		try {
			p17.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// step p17 until it is on last section
		while(p17.getCurrentRoad().getRoadLinkID().contentEquals("9745D155-3C95-4CCD-BC65-0908D57FA83A_0")==false) {
			try {
				p17.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		assert p17.getPathFinder().getTacticalPath().getCurrentJunction().getFID().contentEquals("pave_node_108");
		
		// Move ped 17 to location it should be at when ped 33 is added to simulation
		Coordinate c = new Coordinate(530428.678340819, 180841.2802190967);
        Point pt = GISFunctions.pointGeometryFromCoordinate(c);
		Geometry pGeomNew = pt.buffer(p17.getRad());
        GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, pGeomNew, p17);
		p17.setLoc();
		p17.setBearing(4.016138701200504);
		double[] v17 = {-0.6213405465259059, -0.5193776759134808};
		p17.setV(v17);
		
		
		// Live OD data id 6 is test OD data id 13
		Ped p33 = EnvironmentSetup.createPedestrian(3,13, 0.698355745646177, 61.108214657651644, 0.5, 0.1, 0.9, 4.0, true, 20.0);
		assert (p33.getRad() - 0.1909631708051614) < 0.00001;
		try {
			p33.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Now step both peds along and try to recreate bug
		int count = 0; // After 30 steps ped 33 get's push through a wall. when count =29 things go wrong
		while((p33.getCurrentRoad()!=null) & (p17.getCurrentRoad()!=null)) {
			try {
				// Break if p17 gets to destination
				if (p17.getPathFinder().getTacticalPath().getCurrentJunction()==null) {
					break;
				}
				p17.step();
				p33.step();
				count++;
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		System.out.print(count);
		assert p33.getCurrentRoad()!=null;

	}

}
