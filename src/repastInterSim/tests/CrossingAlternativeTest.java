package repastInterSim.tests;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.UnmarkedCrossingAlternative;
import repastInterSim.environment.Vector;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.TacticalRoute;

class CrossingAlternativeTest {

	@Test
	void testLoadCrossingAlternatives1() {
		try {
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		assert SpaceBuilder.caGeography != null;
		assert SpaceBuilder.caGeography.size() == 5;
		
		CrossingAlternative ca1 = null;
		for (CrossingAlternative ca_: SpaceBuilder.caGeography.getAllObjects()) {
			if(ca_.getID() == 3) {
				ca1 = ca_;
				break;
			}
		}
		
		// Test ca1 setup as expected
		assert ca1.getRoadLinkID().contentEquals("CCBB6F91-5D00-46AA-8835-CC06A8D97B91_0");
		
		String[] expectedITNLinks = {"osgb4000000030238817","osgb4000000030238818"};
		for (int i=0; i<expectedITNLinks.length; i++) {
			assert ca1.getITNRoadLinkIDs()[i].contentEquals(expectedITNLinks[i]);
		}
		
		char[] phase1 = {'r','r'};
		char[] phase2 = {'g','g'};
		char[][] expectedPhases = {phase1, phase2};
		for(int i=0; i< expectedPhases.length; i++) {
			for (int j=0; j<expectedPhases[i].length; j++) {
				assert ca1.getPhases()[i][j] == expectedPhases[i][j];
			}
		}
		
		int[] expectedDurations = {10,30};
		for (int i=0; i<expectedDurations.length; i++) {
			assert ca1.getPhaseDurations()[i] == expectedDurations[i];
		}
		
		// Finally check initial state of signal
		assert ca1.getState("osgb4000000030238817") == 'r';
		assert ca1.getState("osgb4000000030238818") == 'r';
		assert ca1.getState("wrong_link_id") == 'u';
	}
	
	@Test
	void testSignalPhaseChange() {
		try {
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		assert SpaceBuilder.caGeography != null;
		assert SpaceBuilder.caGeography.size() == 5;
		
		CrossingAlternative ca1030 = null;
		CrossingAlternative ca530a = null;
		CrossingAlternative ca530b = null;
		for (CrossingAlternative ca_: SpaceBuilder.caGeography.getAllObjects()) {
			if(ca_.getID() == 3) {
				ca1030 = ca_;
			}
			else if (ca_.getID() == 4) {
				ca530a = ca_;
			}
			else if (ca_.getID() == 5) {
				ca530b = ca_;
			}
		}
		
		// Check initial states of the signals then check that state changes as expected
		assert ca1030.getState("osgb4000000030238817") == 'r';
		assert ca1030.getState("osgb4000000030238818") == 'r';
		
		assert ca530a.getState("osgb4000000030238839") == 'r';
		assert ca530b.getState("osgb4000000030418774") == 'r';
		
		// After 5 ie on the 6th tick) ticks two of the signals should change to green
		for (int i=0; i<6; i++) {
			ca1030.step();
			ca530a.step();
			ca530b.step();
		}
		
		assert ca1030.getState("osgb4000000030238817") == 'r';
		assert ca1030.getState("osgb4000000030238818") == 'r';
		
		assert ca530a.getState("osgb4000000030238839") == 'g';
		assert ca530b.getState("osgb4000000030418774") == 'g';
		
		// After 5 more ticks the other signal should also be green
		for (int i=0; i<6; i++) {
			ca1030.step();
			ca530a.step();
			ca530b.step();
		}

		assert ca1030.getState("osgb4000000030238817") == 'g';
		assert ca1030.getState("osgb4000000030238818") == 'g';
		
		assert ca530a.getState("osgb4000000030238839") == 'g';
		assert ca530b.getState("osgb4000000030418774") == 'g';
		
		// After a further 25 steps 530 signals should revert to red
		for (int i=0; i<26; i++) {
			ca1030.step();
			ca530a.step();
			ca530b.step();
		}
		
		assert ca1030.getState("osgb4000000030238817") == 'g';
		assert ca1030.getState("osgb4000000030238818") == 'g';
		
		assert ca530a.getState("osgb4000000030238839") == 'r';
		assert ca530b.getState("osgb4000000030418774") == 'r';
		
		// After 5 more ticks 1030 signal should be on red whilst the others are back to green
		for (int i=0; i<6; i++) {
			ca1030.step();
			ca530a.step();
			ca530b.step();
		}
		assert ca1030.getState("osgb4000000030238817") == 'r';
		assert ca1030.getState("osgb4000000030238818") == 'r';
		
		assert ca530a.getState("osgb4000000030238839") == 'g';
		assert ca530b.getState("osgb4000000030418774") == 'g';
		
	}
	
	@Test
	void testUnsignalisedCrossings() {
		try {
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		assert SpaceBuilder.caGeography != null;
		assert SpaceBuilder.caGeography.size() == 5;
		
		CrossingAlternative us1 = null;
		CrossingAlternative us2 = null;
		for (CrossingAlternative ca_: SpaceBuilder.caGeography.getAllObjects()) {
			if(ca_.getID() == 1) {
				us1 = ca_;
			}
			else if (ca_.getID() == 2) {
				us2 = ca_;
			}
		}
		
		// Test that these crossing alternatives return 'u' when trying to get their state
		assert us1.getState("osgb4000000030238839") == 'u';
		assert us2.getState("osgb4000000030343781") == 'u';
		
		assert us1.getPhaseDurations() == null;
		assert us1.getPhases() == null;
		
	}
	
	@Test
	void testUnmarkedCrossingAlternativeCoordinates1() {
		// Setup environment
		try {
			EnvironmentSetup.setUpProperties();
			
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
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Coordinate pedLoc = new Coordinate(530420.5, 180821.69);
		double bearing = 2*Math.PI * (55.0/360); // 55 degrees
		Ped p = EnvironmentSetup.createPedAtLocation(4, 2, false, pedLoc, bearing);
		
		UnmarkedCrossingAlternative caU = new UnmarkedCrossingAlternative();
		caU.setPed(p);
		
		Coordinate c1 = caU.getC1();
		Coordinate c2 = caU.getC2();
		
		assert (Math.abs(c1.x-pedLoc.x) < 0.000001) & (Math.abs(c1.y - pedLoc.y) < 0.0000001 );
		assert c2.equals2D(new Coordinate(530416.1087234715,180827.5097640739));
		
	}
	
	@Test
	void testUnmarkedCrossingAlternativeCoordinates2() {
		// Setup environment
		try {
			EnvironmentSetup.setUpProperties();
			
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
			
			EnvironmentSetup.setUpCrossingAlternatives("crossing_lines.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Coordinate pedLoc = new Coordinate(530669.1, 180776.85);
		double bearing = 2*Math.PI * (294.0/360); // 294 degrees
		Ped p = EnvironmentSetup.createPedAtLocation(14, 13, false, pedLoc, bearing);
		
		UnmarkedCrossingAlternative caU = new UnmarkedCrossingAlternative();
		caU.setPed(p);
		
		Coordinate c1 = caU.getC1();
		Coordinate c2 = caU.getC2();
		
		assert (Math.abs(c1.x-pedLoc.x) < 0.000001) & (Math.abs(c1.y - pedLoc.y) < 0.0000001 );
		assert c2.equals2D(new Coordinate(530661.25, 180772.85));
		
	}
	
	@Test
	void testUnmarkedCrossingVehicleFlow() {
		// Setup the environment
		try {
			IO.readProperties();
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
		
		// Create pedestrian that will cross first road link
		Ped ped = EnvironmentSetup.createPedestrian(3,4,false);
		
		// Create a vehicle that moves along same link as pedestrian
		Vehicle v = EnvironmentSetup.createVehicle("osgb4000000029970447", "osgb4000000029970446");
		
		// Step the vehicle and pedestrian once to initialise them
		try {
			v.step();
			ped.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Identify crossing alternatives
		List<CrossingAlternative> cas = TacticalRoute.getCrossingAlternatives(SpaceBuilder.caGeography, ped.getPathFinder().getStrategicPath().subList(0, 1), ped, SpaceBuilder.roadGeography);
		
		// Get unmarked crossing alternative
		UnmarkedCrossingAlternative uc = (UnmarkedCrossingAlternative) cas.stream().filter(ca -> ca.getType().contentEquals("unmarked")).collect(Collectors.toList()).get(0);
		
		assert uc.getC1().equals2D(ped.getLoc());
		
		// Now check vehicle flow of crossing as vehicle approaches
		double crossingTime = (uc.getC1().distance(uc.getC2())) / ped.getSpeed();
		
		while (v.getLoc().distance(uc.getC1()) - v.getSpeed()*crossingTime > 1) {
			assert uc.getvFlow() == 0;
			try {
				v.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Now after last step vehicle is within reach of the crossing
		assert uc.getvFlow() == 1;
		
		// Vehicle flow remains 1 until vehicle passes crossing
		// Test this by moving vehicle along to the end of its route
		while (GISFunctions.coordInFront(v.getLoc(), v.getBearing(), uc.getC1())) {
			assert uc.getvFlow()==1;
			try {
				v.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Check that vehicle is still on road but that vehicle flow passing through crossing has fallen to zero
		int count = 0;
		for (RoadLink orRL: SpaceBuilder.orRoadLinkGeography.getAllObjects()) {
			if(orRL.getFID().contentEquals(uc.getRoadLinkID())) {
				List<RoadLink> itnLinks = SpaceBuilder.orToITN.get(orRL);
				for (int i=0; i<itnLinks.size(); i++) {
					count += itnLinks.get(i).getVehicleCount();
				}
			}
		}
		assert count>0;
		assert uc.getvFlow()==0;
		
	}
	
	@Test
	void testCATTC() {
		// Setup the environment
		try {
			IO.readProperties();
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
		
		// Create pedestrian that will cross first road link
		Ped ped = EnvironmentSetup.createPedestrian(3,4,false);
		
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
		List<CrossingAlternative> cas = TacticalRoute.getCrossingAlternatives(SpaceBuilder.caGeography, ped.getPathFinder().getStrategicPath().subList(0, 1), ped, SpaceBuilder.roadGeography);
		
		// Place ped at edge of a crossing, wlaking in the direction of the other side of the crossing.
		CrossingAlternative ca = cas.stream().filter(ca_ -> ca_.getType().contentEquals("unsignalised")).collect(Collectors.toList()).get(0);
		ped.setLoc(ca.getC1());
		double[] vCross = {ca.getC2().x-ca.getC1().x, ca.getC2().y-ca.getC1().y};
		ped.setBearing(Vector.angleBetweenNorthAndUnitVector(Vector.unitV(vCross)));
		
		// Find mid point
		Coordinate crossingMid = GISFunctions.midwayBetweenTwoCoordinates(ca.getC1(), ca.getC2());
		
		// Point vehicle towards crossing midpoint
		double[] vectorToMidPoint = {crossingMid.x-v.getLoc().x, crossingMid.y-v.getLoc().y};
		v.setBearing(Vector.angleBetweenNorthAndUnitVector(Vector.unitV(vectorToMidPoint)));
		
		// Get TTC for this ped and this vehicle when ped is stationary
		HashMap<Vehicle, Double> ttcs = ca.vehicleTTCs(ped);
		
		// Because ped is stationary expect null ttc value
		assert Math.abs(0.0-ped.getSpeed())<0.0000001;
		assert (ttcs.size() == 1);
		assert ttcs.values().stream().allMatch(ttc -> ttc==null);
		
		// Now set ped speed such that vehicle and ped should collide at the middle of the crossing
		double vDist = v.getLoc().distance(crossingMid)- GlobalVars.vehicleLength/2;
		double pDist = ped.getLoc().distance(crossingMid);
		v.setSpeed(3.0);
		double tVeh = vDist / v.getSpeed(); // Time it will take front of vehicle to reach crossing.
		double pedSpeed = pDist / tVeh; // Speed ped should walk at to reach centre of crossing at about the same time as the vehicle.
		double[] pV = {pedSpeed*Math.sin(ped.getBearing()), pedSpeed*Math.cos(ped.getBearing())}; 
		ped.setV(pV);
		
		ttcs = ca.vehicleTTCs(ped);
		assert ttcs.size() == 1;
		assert ttcs.values().stream().allMatch(ttc -> ttc!=null);
		assert Math.abs(ttcs.get(v) - tVeh) < 0.0000001; 
		
		
		// Test again by placing centre of the vehicle at the crossing midpoint. Make vehicle stationary
		v.setLoc(crossingMid);
		v.setSpeed(0.0);
		v.setBearing(ped.getBearing()+Math.PI/2);
		tVeh = (pDist - GlobalVars.vehicleWidth/2) / ped.getSpeed();
		
		ttcs = ca.vehicleTTCs(ped);
		assert Math.abs(ttcs.get(v) - tVeh) < 0.0000001;
		
		// Test that if vehicle just moves out of the way ttc is null again
		double vehSpeed = GlobalVars.vehicleLength/2 / tVeh + 0.000001;
		v.setSpeed(vehSpeed);
		ttcs = ca.vehicleTTCs(ped);
		assert ttcs.values().stream().allMatch(ttc -> ttc==null);
		
		// Add another vehicle to the road and check that two ttc values are returned
		Vehicle v2 = EnvironmentSetup.createVehicle("osgb4000000029970447", "osgb4000000029970446");
		try {
			v2.getRoute().setRoute();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		v2.setCurrentRoadLinkAndQueuePos(v2.getRoute().getRoadsX().get(0)); //  Add vehicle to first road link in its route.
		
		ttcs = ca.vehicleTTCs(ped);
		assert ttcs.size()==2;
		assert ttcs.values().stream().allMatch(ttc -> ttc==null);
		
		// Move ped and test getting ttc via unmarked crossing alternative
		UnmarkedCrossingAlternative uca = (UnmarkedCrossingAlternative) cas.stream().filter(ca_ -> ca_.getType().contentEquals("unmarked")).collect(Collectors.toList()).get(0);
		
		ped.setLoc(new Coordinate(530451, 180855));
		double[] vCross2 = {uca.getC2().x-uca.getC1().x, uca.getC2().y-uca.getC1().y};
		ped.setBearing(Vector.angleBetweenNorthAndUnitVector(Vector.unitV(vCross2)));
		crossingMid = GISFunctions.midwayBetweenTwoCoordinates(uca.getC1(), uca.getC2());
		
		// Point vehicle towards crossing midpoint
		double[] vectorToMidPoint2 = {crossingMid.x-v2.getLoc().x, crossingMid.y-v2.getLoc().y};
		v2.setBearing(Vector.angleBetweenNorthAndUnitVector(Vector.unitV(vectorToMidPoint2)));
		
		vDist = v2.getLoc().distance(crossingMid)- GlobalVars.vehicleLength/2;
		pDist = ped.getLoc().distance(crossingMid);
		v2.setSpeed(3.0);
		tVeh = vDist / v2.getSpeed(); // Time it will take front of vehicle to reach crossing.
		pedSpeed = pDist / tVeh; // Speed ped should walk at to reach centre of crossing at about the same time as the vehicle.
		double[] pV2 = {pedSpeed*Math.sin(ped.getBearing()), pedSpeed*Math.cos(ped.getBearing())}; 
		ped.setV(pV2);
		
		ttcs = uca.vehicleTTCs(ped);
		assert ttcs.size()==2;
		assert ttcs.get(v) == null;
		assert Math.abs(ttcs.get(v2)-tVeh)<0.00000001;
		
	}
}
