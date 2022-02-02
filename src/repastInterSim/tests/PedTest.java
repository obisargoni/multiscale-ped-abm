package repastInterSim.tests;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.operation.distance.DistanceOp;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.space.gis.Geography;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.OD;
import repastInterSim.environment.Road;
import repastInterSim.environment.UnmarkedCrossingAlternative;
import repastInterSim.environment.Vector;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.TacticalRoute;

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
		Ped ped = EnvironmentSetup.createPedestrian(3,4, null, null, minimisesDistance);
		
		// Move ped to position and bearing that has caused an error in the simulation
        Point pt = GISFunctions.pointGeometryFromCoordinate(c);
		Geometry pGeomNew = pt.buffer(ped.getRad());
        GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, pGeomNew, ped);
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
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp", GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_GEOGRAPHY);
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
		Ped pedMinDist = EnvironmentSetup.createPedestrian(3,4, null, null, false);
		
		// Move ped to position and bearing that has caused an error in the simulation
        Point pt = GISFunctions.pointGeometryFromCoordinate(c);
		Geometry pGeomNew = pt.buffer(pedMinDist.getRad());
        GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, pGeomNew, pedMinDist);
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
    		obstGeomsPoints = pedMinDist.getObstacleGeometries(fieldOfVisionApprox, EnvironmentSetup.pedObstructGeography);
            d2s = new double[fovAngles.size()];
        	displacementDistances = new double[fovAngles.size()];
    		pedMinDist.displacementDistancesToObstacleGeometries(obstGeomsPoints, fovAngles, d2s, displacementDistances);
    	}
    	else {
    		obstGeomsPoints = pedMinDist.getObstacleGeometries(fieldOfVisionApprox, EnvironmentSetup.pedObstructionPointsGeography);
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
		
		
		// Check distance op with point and intersects
		DistanceOp dist12p = new DistanceOp(p1, g2);
		assert (dist12p.distance()<=(r1+r2)) & (g1.intersects(g2));
		
		DistanceOp dist13p = new DistanceOp(p1, g3);
		assert (dist13.distance()>(r1+r3)) & (g1.intersects(g3)==false);
		
		// g3 and g4 should just touch, does identify as touching
		DistanceOp dist34p = new DistanceOp(p3, g4);
		assert (dist34.distance()<=(r3+r4)) & (g3.intersects(g4)==true);
		
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
	
	/*
	 * Test whether using DistanceOp with the ped's centroid can identify the midpoint 
	 * of the peds circle geometry with line obstructions
	 */
	@Test
	void testcontactDistancePedCentroidWithPedObstruction1() {
		
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
		DistanceOp distl3 = new DistanceOp(l, p3);
		assert (distl3.distance()>r3) & (l.intersects(g3)==false);
		
		// Check dist op identifies the same coordiante that intersection does
		DistanceOp distl4 = new DistanceOp(l, p4);
		assert (distl4.distance()<=r4) & (l.intersects(g4)==true);
		
		// Also check that intersecting coords is as expected
		Coordinate i = distl4.nearestPoints()[0];
		Coordinate[] intCoords = l.intersection(g4).getCoordinates();
		assert intCoords.length==2;
		Coordinate expectedMid = new Coordinate(530511.0, 180907.0);
		Coordinate mid = GISFunctions.midwayBetweenTwoCoordinates(intCoords[0], intCoords[1]);
		assert mid.equals2D(expectedMid);
		
		// Expected nearest point is the mid point between the intersecting geometries
		assert  (Math.abs(i.x - expectedMid.x) < 0.0001) & (Math.abs(i.y - expectedMid.y) < 0.0001);
	}
	
	
	/*
	 * Further tests of DistanceOp using peds centroid.
	 * 
	 * Test that it identifies when ped just touches barrier.
	 * 
	 * Test that distance given by DistanceOp is the distance from the ped centroid to the nearest
	 * point on the line.
	 */
	@Test
	void testcontactDistancePedCentroidWithPedObstruction2() {
		
		// Create linestring object that represents a wall
		Coordinate lc1 = new Coordinate(530511, 180910);
		Coordinate lc2 = new Coordinate(530521, 180910);
		Coordinate[] lcs = {lc1, lc2};
		Geometry l = GISFunctions.lineStringGeometryFromCoordinates(lcs);

		
		Coordinate c3 = new Coordinate(530515, 180910.5);
		Point p3 = GISFunctions.pointGeometryFromCoordinate(c3);
		double r3 = 0.5;
		Geometry g3 = p3.buffer(r3);
		
		Coordinate c4 = new Coordinate(530515, 180910.2);
		Point p4 = GISFunctions.pointGeometryFromCoordinate(c4);
		double r4 = 0.5;
		Geometry g4 = p4.buffer(r4);

		
		// Now check distOp and intersects
		DistanceOp distl3 = new DistanceOp(l, p3);
		assert (distl3.distance()<=r3) & (l.intersects(g3)==true);
		
		// Check dist op identifies the same coordiante that intersection does
		DistanceOp distl4 = new DistanceOp(l, p4);
		assert (distl4.distance()<=r4) & (l.intersects(g4)==true);
		
		// Also check that intersecting coords is as expected
		Coordinate i = distl4.nearestPoints()[0];
		Coordinate[] intCoords = l.intersection(g4).getCoordinates();
		assert intCoords.length==2;
		Coordinate expectedMid = new Coordinate(530515.0, 180910.0);
		Coordinate mid = GISFunctions.midwayBetweenTwoCoordinates(intCoords[0], intCoords[1]);
		assert mid.equals2D(expectedMid);
		
		// Expected nearest point is the mid point between the intersecting geometries
		assert  (Math.abs(i.x - expectedMid.x) < 0.0001) & (Math.abs(i.y - expectedMid.y) < 0.0001);
		assert Double.compare(distl4.distance(), mid.distance(c4))==0;
	}
	
	
	@Test
	void testPedsPassingContactAccel() {
		
		// Setup the environment
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRandomDistributions(1);
			
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpPedObstructions();
			EnvironmentSetup.setUpPedObstructionPoints();

			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp", GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_GEOGRAPHY);
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
		
		Ped p17 = EnvironmentSetup.createPedestrian(1, 3, 0.8098254410639371, 56.15750059331304, 0.5, 0.1, 0.9, 4.0, 30, true, 20.0);
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
        GISFunctions.moveAgentToGeometry(EnvironmentSetup.geography, pGeomNew, p17);
		p17.setLoc();
		p17.setBearing(4.016138701200504);
		double[] v17 = {-0.6213405465259059, -0.5193776759134808};
		p17.setV(v17);
		
		
		// Live OD data id 6 is test OD data id 13
		Ped p33 = EnvironmentSetup.createPedestrian(3,13, 0.698355745646177, 61.108214657651644, 0.5, 0.1, 0.9, 4.0, 30, true, 20.0);
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
	
	/*
	 * Method to test that ped correctly updates it route when making a direct crossing.
	 * 
	 * Ped should be prevented from updating whilst choosing a crossing alternative, then should update once crossing alternative is chosen.
	 */
	@Test
	void testPedestrianDirectCrossing1() {
		
		// Setup the environment
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRandomDistributions(1);
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpPedObstructions();

			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp", GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_GEOGRAPHY);
			EnvironmentSetup.setUpPavementNetwork();
						
			EnvironmentSetup.setUpPedODs();
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Initialise a pedestrian, this internally initialises a ped path finder
		boolean minimiseCrossings = true;
		Ped pedMinDist = EnvironmentSetup.createPedestrian(8,9, null, null, minimiseCrossings);
						
		// Check the strategic path is as expected
		String[] expectedRoadIDs = {"762DB27A-3B61-4EAA-B63E-6F1B0BD80D98_0", "A8675945-DE94-4E22-9905-B0623A326221_0", "F4C0B1FB-762C-4492-BB0D-673CC4950CBE_0", "8A9E2D7B-3B48-4A19-B89A-0B4F4D516870_2"};
		for (int i=0;i<pedMinDist.getPathFinder().getStrategicPath().size(); i++) {
			assert pedMinDist.getPathFinder().getStrategicPath().get(i).getFID().contentEquals(expectedRoadIDs[i]);
		}
		
		// Check crossing type.
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getChosenCA()==null;
		
		// Step the ped until a crossing is required
		while (pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().isBlank()) {
			try {
				pedMinDist.step();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		// Expect a direct crossing, which means that the default accumulator junction is the same as the peds current junction
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().isDirectCrossing();
		
		// Now check that ped doesn't pdate its current junction until a crossin alternative is chosen
		Junction currentJ = pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction();
		Junction targetJ = pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getTargetJunction();
		Junction defaultJ = pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getDefaultJunction();
		
		assert defaultJ.getGeom().equals(currentJ.getGeom());
		
		// step ped until crossing chosen
		while( (pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().caChosen()==false) & (pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().crossingRequired()) ) {
			try {
				pedMinDist.step();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		// Check crossing type is now not none
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getChosenCA()!=null;
				
		// Now currentJunction should be the targetJunction 
		assert targetJ.getGeom().equals(pedMinDist.getPathFinder().getTacticalPath().getCurrentJunction().getGeom());
		
		// the accumulator route should have coordinates to go to
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getCrossingCoordinates().size()>0;
		
		// Check that crossing type stays not none until crossing is over
		while ( pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getCrossingCoordinates().size()>0 ) {
			try {
				pedMinDist.step();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getChosenCA()!=null;
		}
	}
	
	/*
	 * Test that pedestrian agent initially yields when reaching the start of a marked crossing be crosses at the next tick since peds have right on way at marked crossings.
	 */
	@Test
	void testPedestrianYieldMarkedCrossing1() {
		// Setup the environment
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpPedObstructions();

			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp", GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_GEOGRAPHY);
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
		
		// Create pedestrian that will cross first road link
		Ped ped = EnvironmentSetup.createPedestrian(3,4, null, null, false);
		ped.getPathFinder().updateTacticalPath(); // initialise peds route
		
		// Create a vehicle that moves along same link as pedestrian
		Vehicle v = EnvironmentSetup.createVehicle("osgb4000000029970447", "osgb4000000029970446");
		try {
			v.getRoute().setRoute();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		v.setCurrentRoadLinkAndQueuePos(v.getRoute().getRoadsX().get(0)); //  Add vehicle to first road link in its route.
				
		// Identify crossing alternatives
		List<CrossingAlternative> cas = TacticalRoute.getCrossingAlternatives(EnvironmentSetup.caGeography, ped.getPathFinder().getStrategicPath().subList(0, 1), ped, EnvironmentSetup.roadGeography);
		
		// Identify marked crossing
		CrossingAlternative mCA = cas.stream().filter(ca -> ca.getType().contentEquals("unsignalised")).collect(Collectors.toList()).get(0);
		
		ped.getPathFinder().getTacticalPath().getAccumulatorRoute().setChosenCA(mCA);
		ped.getPathFinder().step(); // Updates the tactical route to incorporate the crossing
		
		// ped initially isn't set to yield
		assert ped.getYield()==false;
		
		// Step ped until it reaches it's crossing
		while (ped.getPathFinder().getTacticalPath().getAccumulatorRoute().getCrossingCoordinates().size()==2) {
			try {
				ped.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// ped should initially yield. isCrossing is et to true currently when ped reaches start of crossing
		assert ped.getYield()==true;
		assert ped.isCrossing()==true;
		
		// Now set ped speed such that vehicle and ped should collide at the middle of the crossing
		Coordinate crossingMid = GISFunctions.midwayBetweenTwoCoordinates(mCA.getC1(), mCA.getC2());
		double b = GISFunctions.bearingBetweenCoordinates(v.getLoc(), crossingMid);
		v.setBearing(b);
		
		// Set vehicle velocity to ensure a conflict with ped if ped were to cross
		double vDist = v.getLoc().distance(crossingMid)- GlobalVars.vehicleLength/2;
		double pDist = ped.getLoc().distance(crossingMid);
		double pedSpeed = ped.getV0(); // Ped's desired walking speed
		double tPed = pDist / pedSpeed;
		double vSpeed = vDist / tPed ;
		v.setSpeed(vSpeed);
		
		// but after an additional step should continue to cross bc at marked crossing peds trusts vehicle will yield
		try {
			ped.step();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		assert ped.getYield()==false;
		assert ped.isCrossing()==true;
	}
	
	/*
	 * Test that pedestrian agent initially yields when reaching the start of an unmarked crossing and will continue to yield until a sufficient gap in vehicle traffic appears.
	 */
	@Test
	void testPedestrianYieldUnMarkedCrossing1() {
		// Setup the environment
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpPedObstructions();

			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpORRoadNetwork(false);
			
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.setUpPedJunctions();
			EnvironmentSetup.setUpPavementLinks("pedNetworkLinks.shp", GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_CONTEXT, GlobalVars.CONTEXT_NAMES.PAVEMENT_LINK_GEOGRAPHY);
			EnvironmentSetup.setUpPavementNetwork();
						
			EnvironmentSetup.setUpPedODs();
			EnvironmentSetup.setUpVehicleODs("mastermap-itn RoadNode Intersect Within.shp");
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Create pedestrian that will cross first road link
		Ped ped = EnvironmentSetup.createPedestrian(3,4, null, null, false);
		ped.getPathFinder().updateTacticalPath(); // initialise peds route
		
		// Create a vehicle that moves along same link as pedestrian
		Vehicle v = EnvironmentSetup.createVehicle("osgb4000000029970447", "osgb4000000029970446");
		try {
			v.getRoute().setRoute();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		v.setCurrentRoadLinkAndQueuePos(v.getRoute().getRoadsX().get(0)); //  Add vehicle to first road link in its route.
				
		// Identify crossing alternatives
		List<CrossingAlternative> cas = TacticalRoute.getCrossingAlternatives(EnvironmentSetup.caGeography, ped.getPathFinder().getStrategicPath().subList(0, 1), ped, EnvironmentSetup.roadGeography);
		
		// Identify marked crossing
		CrossingAlternative umCA = cas.stream().filter(ca -> ca.getType().contentEquals("unmarked")).collect(Collectors.toList()).get(0);
		
		ped.getPathFinder().getTacticalPath().getAccumulatorRoute().setChosenCA(umCA);
		ped.getPathFinder().step(); // Updates the tactical route to incorporate the crossing
		
		// Initially ped is set to not yield
		assert ped.getYield()==false;
		
		// Step ped until it reaches it's crossing
		while (ped.getPathFinder().getTacticalPath().getAccumulatorRoute().getCrossingCoordinates().size()==2) {
			try {
				ped.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Now that ped has reached crossing (therefore is crossing is true) ped yields
		assert ped.getYield() == true;
		assert ped.isCrossing()==true;
		
		
		// Now set ped speed such that vehicle and ped should collide at the middle of the crossing
		Coordinate crossingMid = GISFunctions.midwayBetweenTwoCoordinates(umCA.getC1(), umCA.getC2());
		double b = GISFunctions.bearingBetweenCoordinates(v.getLoc(), crossingMid);
		v.setBearing(b);
		
		// Set vehicle velocity to ensure a conflict with ped if ped were to cross
		double vDist = v.getLoc().distance(crossingMid)- GlobalVars.vehicleLength/2;
		double pDist = ped.getLoc().distance(crossingMid);
		double pedSpeed = ped.getV0(); // Ped's desired walking speed
		double tPed = pDist / pedSpeed;
		double vSpeed = vDist / tPed ;
		v.setSpeed(vSpeed);
		
		// after an addition step ped still yields due to non null ttc with vehicle
		try {
			ped.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// since vehicle ttc currently not null ped will continue to yield
		double[] pLoc = {ped.getLoc().x, ped.getLoc().y};
    	double[] pV = {ped.getV0()*Math.sin(umCA.getCrossingBearing()), ped.getV0()*Math.cos(umCA.getCrossingBearing())}; 
		HashMap<Vehicle, Double> ttcs = umCA.vehicleTTCs(pLoc, pV);
		assert ttcs.get(v)!=null;
		assert ped.getYield()==true;
		assert ped.isCrossing()==true;
		
		// If vehicle comes to a stop, ttcs become null and ped should stop yielding
		v.setSpeed(0.0);
		try {
			ped.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		ttcs = umCA.vehicleTTCs(pLoc, pV);
		assert ttcs.get(v)==null;
		assert ped.getYield()==false;
		assert ped.isCrossing()==true;
	}
}
