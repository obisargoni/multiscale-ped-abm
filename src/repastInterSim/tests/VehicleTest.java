package repastInterSim.tests;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;

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
			EnvironmentSetup.setUpProperties();
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpVehicleODs("OD_vehicle_nodes_intersect_within.shp");
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
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
		assert v.carFollowingSafeSpeed(vinfront) == Double.MAX_VALUE;
		assert v.getDesiredSpeed(vinfront, new ArrayList<Ped>(), new ArrayList<CrossingAlternative>()) == v.getSpeed() + GlobalVars.defaultVehicleAcceleration*GlobalVars.stepToTimeRatio;
		
		// Check speed after one step
		try {
			v.step();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert v.getSpeed() == GlobalVars.initialVehicleSpeed + GlobalVars.defaultVehicleAcceleration*GlobalVars.stepToTimeRatio;
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
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpVehicleODs("OD_vehicle_nodes_intersect_within.shp");
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
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
		
		// Speed is greater than initial value bc driver speeds up to catch up with vehicle ahead
		double sFollow = vFollower.carFollowingSafeSpeed(vinfront);
		assert sFollow > vFollower.getSpeed();
		
		// However, safe car following speed is greater than the speed vehicle can reach in one time step, so desired speed will be current speed + default acceleration
		double vDes = vFollower.getDesiredSpeed(vinfront, new ArrayList<Ped>(), new ArrayList<CrossingAlternative>()); 
		assert vDes == vFollower.getSpeed() + GlobalVars.defaultVehicleAcceleration*GlobalVars.stepToTimeRatio;
		
		
		// Check speed after one step speed is the desired speed
		try {
			vFollower.step();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert vFollower.getSpeed() ==  vDes;
		
		
		// Now reduce speed of leader vehicle and check that follower vehicle also responds
		vLeader.setSpeed(0.0);
		
		// Check vehicle in front is the leader vehicle.
		vinfront = vFollower.getVehicleInFront();
		assert vinfront.getID() == vLeader.getID();
		
		// Safe following speed should now be less than what it was before
		assert vFollower.carFollowingSafeSpeed(vinfront) < sFollow;
		
		// However, since safe following speed is still greater than can be achieve in a single time step actual speed increases by default acceleration only.
		try {
			vFollower.step();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert vFollower.getSpeed() == vDes + GlobalVars.defaultVehicleAcceleration*GlobalVars.stepToTimeRatio;
	}
	
	
	/*
	 * Test vehicle following across road links.
	 */
	@Test
	void testVehicleDriveFollower2() {
		
		// Set up environment
		try {
			EnvironmentSetup.setUpProperties();
			EnvironmentSetup.setUpObjectGeography();
			EnvironmentSetup.setUpRoads();
			EnvironmentSetup.setUpITNRoadLinks();
			EnvironmentSetup.setUpORRoadLinks();
			EnvironmentSetup.setUpVehicleODs("OD_vehicle_nodes_intersect_within.shp");
			EnvironmentSetup.setUpITNRoadNetwork(true);
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
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
		
		// Need to step once more to move vehicle onto next link
		try {
			vLeader.step();
		} catch (Exception e) {
			e.printStackTrace();
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
		
		// Safe following speed should be max value bc leader vehicle is far enough away not to follow
		assert vFollower.carFollowingSafeSpeed(vinfront) == Double.MAX_VALUE;
		
		// Check speed after one step
		try {
			vFollower.step();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		assert vFollower.getSpeed() == (GlobalVars.initialVehicleSpeed + GlobalVars.defaultVehicleAcceleration*GlobalVars.stepToTimeRatio);
		
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
		
		// Speed should be lower than current speed since follower vehicle needs to slow down.
		assert vFollower.carFollowingSafeSpeed(vinfront) < vFollower.getSpeed();
		
	}
	
	/*
	 * Test that a vehicle agents identifies a crossing pedestrian and adjusts its acceleration in response.
	 */
	@Test
	void testVehicleYieldsCrossingPedestrian1() {
		
		// Setup the environment
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
		
		
		// In ped's step function after stepping the PathFinder ped checks if target coordinate needs updating. 
		// Do that once here in case ped chooses unmarked crossing, which means that target coordinate is current position, and results in null bearing
    	// Finally update the target coordinate if current target coordinate has been reached
    	if (pedMinDist.getLoc().distance(pedMinDist.getPathFinder().getTacticalPath().getTargetCoordinate()) < 0.5) {
    		pedMinDist.getPathFinder().getTacticalPath().updateTargetCoordiante();
    	}
		
		// Now move ped along until it reaches its first crossing point, as indicated by it's distance from the second crossing point
    	Coordinate firstCrossCoord = pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getCrossingCoordinates().getLast(); 
		while (pedMinDist.getLoc().distance(firstCrossCoord)>6.0) {
			try {
				pedMinDist.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}		
		assert pedMinDist.getPathFinder().getTacticalPath().getAccumulatorRoute().getCrossingCoordinates().size() == 1;
		assert v.getCrossingPedestrians().size() == 1;
		
		// Vehicle not yet close enough to pedestrian to yield, so safe ped speed  is set to max value
		double pedSpd = v.pedYieldingSafeSpeed(v.getCrossingPedestrians());
		assert pedSpd == Double.MAX_VALUE;
		
		// Steping vehicle along brings it closer to the pedestrian
		while (v.getLoc().distance(pedMinDist.getLoc()) > v.getDMax()) {
			try {
				v.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		assert v.getLoc().distance(pedMinDist.getLoc()) < v.getDMax();
		
		// Now safe following speed should be lower
		assert v.pedYieldingSafeSpeed(v.getCrossingPedestrians()) < pedSpd;
		
		// How ever safe following speed still not low enough to make vehicle slow down so remains at current (max) speed
		assert v.getDesiredSpeed(null, v.getCrossingPedestrians(), new ArrayList<CrossingAlternative>()) == v.getSpeed();		
		
		// Step vehicle ahead some more to check that it eventually slows down
		// Steping vehicle along brings it closer to the pedestrian
		while (v.getLoc().distance(pedMinDist.getLoc()) > 5) {
			try {
				v.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Desired speed should now be lower than accelerated speed safe following speed still not low enough to prevent vehicle from accelerating
		assert v.pedYieldingSafeSpeed(v.getCrossingPedestrians()) < v.getSpeed();
		assert v.getDesiredSpeed(null, v.getCrossingPedestrians(), new ArrayList<CrossingAlternative>()) == v.pedYieldingSafeSpeed(v.getCrossingPedestrians());
		
	}
	
	/*
	 * Test that a vehicle agents identifies a crossing pedestrian but does not yield when pedestrian is behind the vehicle
	 */
	@Test
	void testVehicleYieldsCrossingPedestrian2() {
		
		// Setup the environment
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
		
		// Now walk pedestrian along until it is about to start crossing. If ped chooses unmarked crossing it will be at the location of first crossing point
		// so while loop won't be entered.
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
		
		// Now update ped's current coordinate, which means that ped has reached first crossing coordinate and therefore has started crossing
		pedMinDist.getPathFinder().getTacticalPath().updateTargetCoordiante();
		assert v.getCrossingPedestrians().size() == 1;
		
		// Now vehicle should not yield to ped despite being close enough to detect ped
		assert v.getLoc().distance(pedMinDist.getLoc()) < v.getDMax();
		double pedSpd = v.pedYieldingSafeSpeed(v.getCrossingPedestrians());
		assert pedSpd == Double.MAX_VALUE;		
	}
	
	/*
	 * Test that a vehicle agent identifies and yields to a crossing alternative when the crossing signals red.
	 */
	@Test
	void testVehicleYieldsCrossingAlternative1() {
		
		// Setup the environment
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
			assert v.crossingAlternativeSafeSpeed(cas) == Double.MAX_VALUE;
			try {
				v.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Now that vehicle is close enough to perceive crossing alternative
		double caSpd = v.crossingAlternativeSafeSpeed(cas); 
		assert caSpd < Double.MAX_VALUE;
		
		// Safe speed still greater than current speed so vehicle doesn't slow down
		assert v.getDesiredSpeed(null, v.getCrossingPedestrians(), cas) == v.getSpeed();
		
		// Move vehicle along so that it is closer to the crossing
		while (v.getLoc().distance(ca.getSignalLoc())>5) {
			try {
				v.step();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		// Now safe speed is lower and is the speed the vehicle must follow
		assert v.crossingAlternativeSafeSpeed(cas) < caSpd;
		assert v.getDesiredSpeed(null, v.getCrossingPedestrians(), cas) == v.crossingAlternativeSafeSpeed(cas);
		 
		
		// Now step the crossing alternative so that it changes from red to green and check that vehicle acceleration reverts to default value
		for(int i=0; i<ca.getPhaseDurations()[0]+1; i++){
			ca.step();
		}
		
		caSpd = v.crossingAlternativeSafeSpeed(cas); 
		assert caSpd == Double.MAX_VALUE;
	}

}
