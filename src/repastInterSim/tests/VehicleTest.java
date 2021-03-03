package repastInterSim.tests;


import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.OD;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadLinkContext;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;

class VehicleTest {
		
	Context<Object> context = new DefaultContext<Object>();;
	
	String testGISDir = ".//data//test_gis_data//";
	String roadLinkPath = null;
	String pedestrianRoadsPath = null;
	String vehicleRoadsPath = null;
	String pavementLinkPath = null;
	String pedJPath = null;
	String serialisedLookupPath = null;
	
	void setUpProperties() throws IOException {
		IO.readProperties();
	}
	
	void setUpObjectGeography() {
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		SpaceBuilder.geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, context, geoParams);
		SpaceBuilder.geography.setCRS(GlobalVars.geographyCRSString);
		context.add(SpaceBuilder.geography);
	}
	
	void setUpRoads() throws Exception {
		pedestrianRoadsPath = testGISDir + "topographicAreaPedestrian.shp";
		vehicleRoadsPath = testGISDir + "topographicAreaVehicle.shp";
		serialisedLookupPath = testGISDir + "road_link_roads_cache.serialised";
		
		// Get road geography
		Context<Road> testRoadContext = new RoadContext();
		GeographyParameters<Road> GeoParamsRoad = new GeographyParameters<Road>();
		SpaceBuilder.roadGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testRoadGeography", testRoadContext, GeoParamsRoad);
		SpaceBuilder.roadGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		try {
			GISFunctions.readShapefile(Road.class, vehicleRoadsPath, SpaceBuilder.roadGeography, testRoadContext);
			GISFunctions.readShapefile(Road.class, pedestrianRoadsPath, SpaceBuilder.roadGeography, testRoadContext);
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		SpatialIndexManager.createIndex(SpaceBuilder.roadGeography, Road.class);
	}
	
	Geography<RoadLink> setUpRoadLinks(String roadLinkFile) throws MalformedURLException, FileNotFoundException {
		roadLinkPath = testGISDir + roadLinkFile;
		
		// Initialise test road link geography and context
		SpaceBuilder.roadLinkContext = new RoadLinkContext();
		GeographyParameters<RoadLink> GeoParams = new GeographyParameters<RoadLink>();
		Geography<RoadLink> rlG = GeographyFactoryFinder.createGeographyFactory(null).createGeography("orRoadLinkGeography", SpaceBuilder.roadLinkContext, GeoParams);
		rlG.setCRS(GlobalVars.geographyCRSString);
				
		GISFunctions.readShapefile(RoadLink.class, roadLinkPath, rlG, SpaceBuilder.roadLinkContext);
		SpatialIndexManager.createIndex(rlG, RoadLink.class);
		
		return rlG;
	}
	
	void setUpITNRoadLinkGeography() throws Exception {
		SpaceBuilder.roadLinkGeography = setUpRoadLinks("mastermap-itn RoadLink Intersect Within with orientation.shp");
	}

	void setUpODs(String odFile) throws MalformedURLException, FileNotFoundException {
		
		// Initialise OD context and geography
		Context<OD> ODContext = new PedestrianDestinationContext();
		GeographyParameters<OD> GeoParamsOD = new GeographyParameters<OD>();
		SpaceBuilder.vehicleDestinationGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("testODGeography", ODContext, GeoParamsOD);
		SpaceBuilder.vehicleDestinationGeography.setCRS(GlobalVars.geographyCRSString);
		
		// Load vehicle origins and destinations
		String testODFile = testGISDir + odFile;
		GISFunctions.readShapefile(OD.class, testODFile, SpaceBuilder.vehicleDestinationGeography, ODContext);
		SpatialIndexManager.createIndex(SpaceBuilder.vehicleDestinationGeography, OD.class);
	}
	
	void setUpRoadNetwork(boolean isDirected) {
		Context<Junction> junctionContext = new JunctionContext();
		GeographyParameters<Junction> GeoParamsJunc = new GeographyParameters<Junction>();
		SpaceBuilder.junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("junctionGeography", junctionContext, GeoParamsJunc);
		SpaceBuilder.junctionGeography.setCRS(GlobalVars.geographyCRSString);
		
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, isDirected);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		SpaceBuilder.roadNetwork = builder.buildNetwork();
		
		GISFunctions.buildGISRoadNetwork(SpaceBuilder.roadLinkGeography, junctionContext, SpaceBuilder.junctionGeography, SpaceBuilder.roadNetwork);
	}
	
	private Vehicle createVehicle(OD o, OD d) {
		Vehicle V = new Vehicle(GlobalVars.maxVehicleSpeed, GlobalVars.defaultVehicleAcceleration, GlobalVars.initialVehicleSpeed, o, d);
		context.add(V);
		Coordinate oCoord = o.getGeom().getCentroid().getCoordinate();
		GeometryFactory fac = new GeometryFactory();
		Point pt = fac.createPoint(oCoord);
		Geometry vehicleCircle = pt.buffer(2);
		GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, vehicleCircle, V);
		V.setLoc();
		return V;
	}
	
	private Vehicle createVehicle(String oID, String dID) {
		OD o = null;
		OD d = null;
		for (OD od: SpaceBuilder.vehicleDestinationGeography.getAllObjects()) {
			if (od.getFID().contentEquals(oID)) {
				o = od;
			}
			else if (od.getFID().contentEquals(dID)) {
				d = od;
			}
		}
		return createVehicle(o,d);
	}
	
	@Test
	void testVehicleDriveLeader() {
		
		// Set up environment
		try {
			IO.readProperties();
			setUpObjectGeography();
			setUpITNRoadLinkGeography();
			setUpODs("OD_vehicle_nodes_intersect_within.shp");
			setUpRoadNetwork(true);
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
		Vehicle v = createVehicle(originID, destID);
		
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
		assert v.setAccFollowing(vinfront) == GlobalVars.defaultVehicleAcceleration;
		
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
			setUpObjectGeography();
			setUpITNRoadLinkGeography();
			setUpODs("OD_vehicle_nodes_intersect_within.shp");
			setUpRoadNetwork(true);
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
		Vehicle vLeader = createVehicle(originID, destID);
		
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
		Vehicle vFollower = createVehicle(originID, destID);
		
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
		assert vFollower.setAccFollowing(vinfront) > GlobalVars.defaultVehicleAcceleration;
		
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
		assert vFollower.setAccFollowing(vinfront) < 0.0;
		
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
			setUpObjectGeography();
			setUpITNRoadLinkGeography();
			setUpODs("OD_vehicle_nodes_intersect_within.shp");
			setUpRoadNetwork(true);
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
		Vehicle vLeader = createVehicle(originID, destID);
		
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
		Vehicle vFollower = createVehicle(originID, destID);
		
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
		assert vFollower.setAccFollowing(vinfront) == GlobalVars.defaultVehicleAcceleration;
		
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
		assert vFollower.setAccFollowing(vinfront) < 0;
		
	}

}
