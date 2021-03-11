package repastInterSim.tests;

import java.io.IOException;
import java.util.List;

import org.junit.jupiter.api.Test;


import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;

class VehicleTest {
	
	/*
	 * Test that get crossing pedestrians returns an empty list when there are no pedestrians
	 */
	@Test
	void testGetCrossingPedestrians1() {
		// Set up environment
		try {
			IO.readProperties();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpVehicleODs("OD_vehicle_nodes_intersect_within.shp");
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Choose origin and destination
		String originID = "osgb4000000029970681";
		String destID = "osgb4000000029971717";
		
		// Create vehicle
		Vehicle v = EnvironmentSetup.createVehicle(originID, destID);
		
		// Initialise route and current link of vehicle		
		try {
			v.getRoute().setRoute();
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		v.setCurrentRoadLinkAndQueuePos(v.getRoute().getRoadsX().get(0));
		
		assert v.getCrossingPedestrians().size()==0; 
	}
	
	@Test
	void testVehicleDriveLeader() {
		
		// Set up environment
		try {
			IO.readProperties();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpVehicleODs("OD_vehicle_nodes_intersect_within.shp");
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Choose origin and destination
		String originID = "osgb4000000029970681";
		String destID = "osgb4000000029971717";
		
		// Create vehicle
		Vehicle v = EnvironmentSetup.createVehicle(originID, destID);
		
		// Initialise route and current link of vehicle		
		try {
			v.getRoute().setRoute();
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		v.setCurrentRoadLinkAndQueuePos(v.getRoute().getRoadsX().get(0));
		
		// Check that no vehicle is in front and that acceleration is equal to the default acceleration value
		Vehicle vinfront = v.getVehicleInFront(); 
		assert vinfront == null;
		assert v.carFollowingAcceleration(vinfront) == GlobalVars.defaultVehicleAcceleration;
		
		// Check speed after one step
		try {
			v.step();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert v.getSpeed() == GlobalVars.initialVehicleSpeed + GlobalVars.defaultVehicleAcceleration;
	}
	
	/*
	 * Test that leader vehicle is identified when on the same link as follower vehicle. Test that acceleration is set so that follower vehicle
	 * matches speed of leader after one tick - this is an overly simplistic car following model.
	 * 
	 * Check that follower vehicle decelerates when leader vehicle has stopped.
	 */
	@Test
	void testVehicleDriveFollower1() {
		
		// Set up environment
		try {
			IO.readProperties();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpVehicleODs("OD_vehicle_nodes_intersect_within.shp");
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Choose origin and destination
		String originID = "osgb4000000029970681";
		String destID = "osgb4000000029971717";
		
		// Create leader vehicle and step forward for some ticks
		Vehicle vLeader = EnvironmentSetup.createVehicle(originID, destID);
		
		// Step the leader vehicle ahead a few ticks
		for (int i=0; i<5; i++) {
			try {
				vLeader.step();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		// Now add follower vehicle
		Vehicle vFollower = EnvironmentSetup.createVehicle(originID, destID);
		
		try {
			vFollower.getRoute().setRoute();
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		vFollower.setCurrentRoadLinkAndQueuePos(vFollower.getRoute().getRoadsX().get(0));
		
		// Check vehicle in front is the leader vehicle.
		Vehicle vinfront = vFollower.getVehicleInFront();
		assert vinfront.getID() == vLeader.getID();
		
		// Acceleration is greater than default bc driver speeds up to catch up with vehicle ahead
		assert vFollower.carFollowingAcceleration(vinfront) > GlobalVars.defaultVehicleAcceleration;
		
		// Check speed after one step
		try {
			vFollower.step();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert vFollower.getSpeed() > (GlobalVars.initialVehicleSpeed + GlobalVars.defaultVehicleAcceleration);
		
		
		// Now reduce speed of leader vehicle and check that follower vehicle also responds
		vLeader.setSpeed(0.0);
		
		// Check vehicle in front is the leader vehicle.
		vinfront = vFollower.getVehicleInFront();
		assert vinfront.getID() == vLeader.getID();
		
		// Acceleration should still be default because car in front is moving at a higher speed. No need for this vehicle to slow down.
		assert vFollower.carFollowingAcceleration(vinfront) < 0.0;
		
		// Check speed after one step
		try {
			vFollower.step();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert vFollower.getSpeed() - vLeader.getSpeed() < 0.0000001;
	}
	
	
	/*
	 * Test vehicle following across road links.
	 */
	@Test
	void testVehicleDriveFollower2() {
		
		// Set up environment
		try {
			IO.readProperties();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpVehicleODs("OD_vehicle_nodes_intersect_within.shp");
			EnvironmentSetup.setUpITNRoadNetwork(true);
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Choose origin and destination
		String originID = "osgb4000000029970681";
		String destID = "osgb4000000029971717";
		
		// Create leader vehicle and step forward for some ticks
		Vehicle vLeader = EnvironmentSetup.createVehicle(originID, destID);
		
		// Step once to initialise route
		try {
			vLeader.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Step the leader vehicle ahead until it moves onto the next road link
		String currentRLID = vLeader.getRoute().getRoadsX().get(0).getFID();
		while (vLeader.getRoute().getRoadsX().get(0).getFID().contentEquals(currentRLID)) {
			try {
				vLeader.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		assert vLeader.getRoute().getRoadsX().get(0).getFID().contentEquals(currentRLID) == false;
		
		// Now add follower vehicle
		Vehicle vFollower = EnvironmentSetup.createVehicle(originID, destID);
		
		try {
			vFollower.getRoute().setRoute();
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		vFollower.setCurrentRoadLinkAndQueuePos(vFollower.getRoute().getRoadsX().get(0));
		
		// Check vehicle at end of next link is the leader, but that vehicle in front is null bc leader vehicle is too far away
		Vehicle vEnd = vFollower.getVehicleAtEndOfNextRoadLink();
		assert vEnd.getID() == vLeader.getID();
		
		Vehicle vinfront = vFollower.getVehicleInFront();
		assert vinfront == null;
		
		// Acceleration should be the default value bc leader vehicle is far enough away not to follow
		assert vFollower.carFollowingAcceleration(vinfront) == GlobalVars.defaultVehicleAcceleration;
		
		// Check speed after one step
		try {
			vFollower.step();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert vFollower.getSpeed() == (GlobalVars.initialVehicleSpeed + GlobalVars.defaultVehicleAcceleration);
		
		
		// Now step follower vehicle until it is about to enter next road link and check that it adjusts speed to match leader vehicle
		currentRLID = vFollower.getRoute().getRoadsX().get(0).getFID();
		boolean keepStepping = true;
		while (keepStepping) {
			try {
				vFollower.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			// Stop stepping if about to enter next road link
			if (vFollower.getRoute().getRoadsX().get(1).getFID().contentEquals(currentRLID)==false) {
				keepStepping = false;
			}
		}
		
		// Set leader vehicle speed to that follower has to slow down
		vLeader.setSpeed(0);
		
		// Check vehicle in front is the leader vehicle.
		vinfront = vFollower.getVehicleInFront();
		assert vinfront.getID() == vLeader.getID();
		
		// Acceleration should be -ve bc follower vehicle needs to slow down.
		assert vFollower.carFollowingAcceleration(vinfront) < 0;
		
	}
	
	/*
	 * Test that a vehicle agents identifies a crossing pedestrian and adjusts its acceleration in response.
	 */
	@Test
	void testVehicleYieldsCrossingPedestrian1() {
		
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
		Ped pedMinDist = EnvironmentSetup.createPedestrian(3,4,false);
		
		// Create a vehicle that moves along same link as pedestrian
		Vehicle v = EnvironmentSetup.createVehicle("osgb4000000029970447", "osgb4000000029970446");
		
		// Step the vehicle and pedestrian once to initialise them
		try {
			v.step();
			pedMinDist.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Check that the vehicle doesn't yet perceive any peds crossing
		assert v.getCrossingPedestrians().size() == 0;
		
		// Step the path finder until a crossing is chosen. This keeps the ped stationary so creates a situation where the ped is not next to the crossing coord
		while (pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().caChosen() == false) {
			try {
				pedMinDist.getPathFinder().step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Vehicle still shouldn't perceive any crossing pedestrians
		assert v.getCrossingPedestrians().size() == 0;
		
		// Now update ped's current coordinate, which means that ped has reach first crossing coordinate and therefore has started crossing
		pedMinDist.getPathFinder().getTacticalPath().updateTargetCoordiante();
		assert v.getCrossingPedestrians().size() == 1;
		
		// Vehicle not yet close enough to pedestrian to yield, so ped yielding acceleration is set to max value
		double pedAcc = v.pedYieldingAcceleration(v.getCrossingPedestrians());
		assert pedAcc == Double.MAX_VALUE;
		
		// Steping vehicle along brings it closer to the pedestrian
		while (v.getLoc().distance(pedMinDist.getLoc()) > v.getDMax()) {
			try {
				v.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		assert v.getLoc().distance(pedMinDist.getLoc())< v.getDMax();
		
		// Now vehicle should yield to ped
		pedAcc = v.pedYieldingAcceleration(v.getCrossingPedestrians());
		assert pedAcc < GlobalVars.defaultVehicleAcceleration;		
	}
	
	/*
	 * Test that a vehicle agents identifies a crossing pedestrian but does not yield when pedestrian is behind the vehicle
	 */
	@Test
	void testVehicleYieldsCrossingPedestrian2() {
		
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
		Ped pedMinDist = EnvironmentSetup.createPedestrian(3,4,false);
		
		// Create a vehicle that moves along same link as pedestrian
		Vehicle v = EnvironmentSetup.createVehicle("osgb4000000029970447", "osgb4000000029970446");
		
		// Step the vehicle and pedestrian once to initialise them
		try {
			v.step();
			pedMinDist.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Check that the vehicle doesn't yet perceive any peds crossing
		assert v.getCrossingPedestrians().size() == 0;
		
		// Step the path finder until a crossing is chosen. This keeps the ped stationary so creates a situation where the ped is not next to the crossing coord
		while (pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().caChosen() == false) {
			try {
				pedMinDist.getPathFinder().step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Now walk pedestrian along until it is about to start crossing
		while (pedMinDist.getLoc().distance(pedMinDist.getPathFinder().getTacticalPath().getTargetCoordinate()) > 2.5) {
			try {
				pedMinDist.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Vehicle still shouldn't perceive any crossing pedestrians
		assert v.getCrossingPedestrians().size() == 0;
		
		// Now step vehicle along until it passes the pedestrian
		double prevDist = v.getLoc().distance(pedMinDist.getLoc());
		double distDiff = 0.0; 
		while (distDiff != 1.0) {
			try {
				v.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			distDiff = Math.signum(v.getLoc().distance(pedMinDist.getLoc()) - prevDist);
			prevDist = v.getLoc().distance(pedMinDist.getLoc());
		}
		
		// Now update ped's current coordinate, which means that ped has reach first crossing coordinate and therefore has started crossing
		pedMinDist.getPathFinder().getTacticalPath().updateTargetCoordiante();
		assert v.getCrossingPedestrians().size() == 1;
		
		// Now vehicle should not yield to ped despite being close enough to detect ped
		assert v.getLoc().distance(pedMinDist.getLoc()) < v.getDMax();
		double pedAcc = v.pedYieldingAcceleration(v.getCrossingPedestrians());
		assert pedAcc == Double.MAX_VALUE;		
	}
	
	/*
	 * Test that a vehicle agent identifies and yields to a crossing alternative when the crossing signals red.
	 */
	@Test
	void testVehicleYieldsCrossingAlternative1() {
		
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
			
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
			
			EnvironmentSetup.assocaiteRoadsWithRoadLinks();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Create a vehicle that travels along a link with a crossing on
		Vehicle v = EnvironmentSetup.createVehicle("osgb4000000029970431", "osgb4000000029971717");
		
		// Step the vehicle once to initialise
		try {
			v.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// Check that vehicle identifies crossing alternative
		List<CrossingAlternative> cas = v.getRoadLinkCrossingAlterantives(v.getRoute().getRoadsX().get(0).getFID());
		assert cas.size()==1;
		CrossingAlternative ca = cas.get(0);
		
		// Check that while vehicle is further than perception distance away acceleration crossing alternative does not affect acceleration
		while (v.getLoc().distance(ca.getSignalLoc())>v.getDMax()) {
			assert v.crossingAlternativeAcceleration(cas) == Double.MAX_VALUE;
			try {
				v.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Now that vehicle is close enough to perceive crossing alternative check that it adjusts its acceleration accordingly
		double caAcc = v.crossingAlternativeAcceleration(cas); 
		assert caAcc < GlobalVars.defaultVehicleAcceleration;
		
		// Step the vehicle and check that its acceleration is the acceleration due to crossing alternative
		try {
			v.step();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		assert v.getAcc() == caAcc; 
		
		// Now step the crossing alternative so that it changes from red to green and check that vehicle acceleration reverts to default value
		for(int i=0; i<ca.getPhaseDurations()[0]+1; i++){
			ca.step();
		}
		
		caAcc = v.crossingAlternativeAcceleration(cas); 
		assert caAcc == Double.MAX_VALUE;
	}

}
