package repastInterSim.tests;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;

import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.UnmarkedCrossingAlternative;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.TacticalRoute;

class TestCrossingAlternative {
	
	
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
		assert uc.getRoad().getRoadLinksVehicleCount()>0;
		assert uc.getvFlow()==0;
		
	}
}
