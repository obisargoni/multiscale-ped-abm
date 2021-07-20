/*
©Copyright 2012 Nick Malleson
This file is part of RepastCity.

RepastCity is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RepastCity is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RepastCity.  If not, see <http://www.gnu.org/licenses/>.
*/

package repastInterSim.pathfinding;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.operation.distance.DistanceOp;

import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repast.simphony.space.graph.ShortestPath;
import repastInterSim.environment.Cacheable;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.OD;
import repastInterSim.environment.Road;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;
/**
 * Create routes around a GIS road network. The <code>setRoute</code> function actually finds the route and can be
 * overridden by subclasses to create different types of Route. See the method documentation for details of how routes
 * are calculated.
 * 
 * <p>
 * A "unit of travel" is the distance that an agent can cover in one iteration (one square on a grid environment or the
 * distance covered at walking speed in an iteration on a GIS environment). This will change depending on the type of
 * transport the agent is using. E.g. if they are in a car they will be able to travel faster, similarly if they are
 * travelling along a transport route they will cover more ground.
 * </p>
 * 
 * @author Nick Malleson
 */
public class RoadNetworkRoute implements Cacheable {

	private static Logger LOGGER = Logger.getLogger(RoadNetworkRoute.class.getName());

	static {
		// Route.routeCache = new Hashtable<CachedRoute, CachedRoute>();
	}
	
	protected Coordinate origin;
	protected Coordinate destination;
	
	protected Junction[] routeEndpoints;

	/*
	 * The route consists of a list of roads which describe how to get to the destination from the origin.
	 */
	private int currentPosition;
	protected List<Double> routeSpeedsX;
	protected List<RoadLink> roadsX;

	// Record which function has added each coord, useful for debugging
	protected List<String> routeDescriptionX;
	
	// Records whether origin and destination lie on opposite sides of the road network route
	private int routeParity;

	/*
	 * Cache every coordinate which forms a road so that Route.onRoad() is quicker. Also save the Road(s) they are part
	 * of, useful for the agent's awareness space (see getRoadFromCoordCache()).
	 */
	private static volatile Map<Coordinate, List<RoadLink>> coordCache;
	
	
	/*
	 * Store which pavement junction each OD is closest to
	 */
	private static volatile NearestPavementJunctionCoordCache odPaveJuncCache;
	// To stop threads competing for the cache:
	private static Object odPaveJuncCacheLock = new Object();
	/*
	 * Cache of road link ID to the road objects that have these IDs. Use to find the pedestrian and vehicle polygons links to a road link
	 */
	private static volatile RoadLinkRoadsCache roadLinkRoadsCache = null;
	// To stop threads competing for the cache:
	private static Object roadLinkRoadsCacheLock = new Object();
	
	// File names used for cache
	private static String gisDir;
	private static File odFile;
	private static File vehcileRoadsFile;
	private static File pedestrianRoadsFile;
	private static File pavementJunctionFile;	

	private static File serialisedLoc;
	private static File paveJuncSeriealizedLoc;
	
	private static GeometryFactory geomFac = new GeometryFactory();



	/**
	 * Create a new route object
	 * 
	 * @param origin
	 * 		The origin coordinate of the route
	 * @param destination
	 * 		The destination coordinate of the route
	 */
	public RoadNetworkRoute(Coordinate origin, Coordinate destination) {
		
		// Assert that origin and destination coordinates are not null
		try {
			checkNotNull(origin, destination);
		} catch (RoutingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}		
		
		this.origin = origin;
		this.destination = destination;
		
		// File names used for caches
		gisDir = IO.getProperty(GlobalVars.GISDataDirectory);
		odFile = new File(gisDir + IO.getProperty(GlobalVars.PedestrianDestinationsFile));
		vehcileRoadsFile = new File(gisDir + IO.getProperty(GlobalVars.VehicleRoadShapefile));
		pedestrianRoadsFile = new File(gisDir + IO.getProperty(GlobalVars.PedestrianRoadShapefile));
		pavementJunctionFile = new File(gisDir + IO.getProperty(GlobalVars.PavementJunctionShapeFile));
		
		serialisedLoc = new File(gisDir + IO.getProperty(GlobalVars.ODORRoadLinkCoordsCache));
		paveJuncSeriealizedLoc = new File(gisDir + IO.getProperty(GlobalVars.ODPavementJunctionCache));
	}
	
	/**
	 * Find the shortest route from the input origin junctions to the input destination junctions. If multiple origin or destination junctions are given
	 * the shortest of the routes between OD pairs is set are the road link route.
	 * 
	 * @throws Exception
	 */
	public void setRoadLinkRoute(List<Junction> currentJunctions, List<Junction> destJunctions) throws Exception {
		long time = System.nanoTime();

		this.roadsX = new Vector<RoadLink>();
		this.routeDescriptionX = new Vector<String>();
		this.routeSpeedsX = new Vector<Double>();

		// No route cached, have to create a new one (and cache it at the end).
		try {			
			/*
			 * Now have possible routes (2 origin junctions, 2 destination junctions max, less if directed roads don't allow junctions to be  accessed) need to pick which junctions
			 * form shortest route
			 */
			this.routeEndpoints = new Junction[2];
			List<RepastEdge<Junction>> shortestPath = getShortestRoute(SpaceBuilder.orRoadNetwork, currentJunctions, destJunctions, routeEndpoints);

			/*
			 * Add the road links that make up the shortest path to the class attribute lists
			 */
			this.addPathToRoute(shortestPath);

		} catch (RoutingException e) {
			/*
			LOGGER.log(Level.SEVERE, "Route.setRoute(): Problem creating route for " + " going from " + currentCoord.toString() + " to " + this.destination.toString());
					*/
			throw e;
		}
		// Cache the route and route speeds
		// List<Coordinate> routeClone = Cloning.copy(theRoute);
		// LinkedHashMap<Coordinate, Double> routeSpeedsClone = Cloning.copy(this.routeSpeeds);
		// cachedRoute.setRoute(routeClone);
		// cachedRoute.setRouteSpeeds(routeSpeedsClone);

		// cachedRoute.setRoute(this.routeX, this.roadsX, this.routeSpeedsX, this.routeDescriptionX);
		// synchronized (Route.routeCache) {
		// // Same cached route is both value and key
		// Route.routeCache.put(cachedRoute, cachedRoute);
		// }
		// TempLogger.out("...Route cacheing new route with unique id " + cachedRoute.hashCode());
		/*
		LOGGER.log(Level.FINER, "Route Finished planning route for " + "with "
				+ this.routeX.size() + " coords in " + (0.000001 * (System.nanoTime() - time)) + "ms.");
				*/
	}
	

	/**
	 * Find a route from the origin to the destination. A route is a list of Coordinates which describe the route to a
	 * destination restricted to a road network.
	 * 
	 * @throws Exception
	 */
	public void setRoadLinkRoute() throws Exception {
		Coordinate currentCoord = this.origin;		
		Coordinate destCoord = this.destination;


		// No route cached, have to create a new one (and cache it at the end).
		try {
			
			/*
			 * Find the nearest junctions to our current position (road endpoints)
			 */

			// Start by Finding the road that this coordinate is on
			RoadLink currentRoad = findNearestObject(currentCoord, SpaceBuilder.orRoadLinkGeography, null,
					GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE);
			// Find which Junction is closest to us on the road.
			List<Junction> currentJunctions = currentRoad.getJunctions();

			/* Find the nearest Junctions to our destination (road endpoints) */

			// Find the road that this coordinate is on
			RoadLink destRoad = findNearestObject(destCoord, SpaceBuilder.orRoadLinkGeography, null,
					GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE);
			// Find which Junction connected to the edge is closest to the coordinate.
			List<Junction> destJunctions = destRoad.getJunctions();
			
			setRoadLinkRoute(currentJunctions, destJunctions);

		} catch (RoutingException e) {
			/*
			LOGGER.log(Level.SEVERE, "Route.setRoute(): Problem creating route for " + " going from " + currentCoord.toString() + " to " + this.destination.toString());
					*/
			throw e;
		}
		// Cache the route and route speeds
		// List<Coordinate> routeClone = Cloning.copy(theRoute);
		// LinkedHashMap<Coordinate, Double> routeSpeedsClone = Cloning.copy(this.routeSpeeds);
		// cachedRoute.setRoute(routeClone);
		// cachedRoute.setRouteSpeeds(routeSpeedsClone);

		// cachedRoute.setRoute(this.routeX, this.roadsX, this.routeSpeedsX, this.routeDescriptionX);
		// synchronized (Route.routeCache) {
		// // Same cached route is both value and key
		// Route.routeCache.put(cachedRoute, cachedRoute);
		// }
		// TempLogger.out("...Route cacheing new route with unique id " + cachedRoute.hashCode());
		/*
		LOGGER.log(Level.FINER, "Route Finished planning route for " + "with "
				+ this.routeX.size() + " coords in " + (0.000001 * (System.nanoTime() - time)) + "ms.");
				*/
	}
	
	/**
	 * Find a route from the origin to the destination by first finding the nearest junctions to the pedestrian then finding the shortest road network
	 * path between those junctions. 
	 * 
	 * A route is a list of Coordinates which describe the route to a
	 * destination restricted to a road network.
	 * 
	 * @throws Exception
	 */
	public Junction[] setRoadLinkRoute(Geography<Junction> pavementJunctionsGeography, Network<Junction> pavementNetwork) throws Exception {
		Coordinate currentCoord = this.origin;		
		Coordinate destCoord = this.destination;


		// No route cached, have to create a new one (and cache it at the end).
		Junction currentPaveJ;
		Junction destPaveJ;
		RoadLink currentRoad = null;
		RoadLink destRoad = null;
			
		// Find the nearest pavement junctions to start and end of journey. Use these to find start and end road junctions
		currentPaveJ = getNearestPavementJunctionToOD(currentCoord);
		destPaveJ = getNearestPavementJunctionToOD(destCoord);
		try {
			// Find the road link corresponding to the pavement polygon the agents start and end pavement junctions are on
			// I think this is quite slow and needs rethinking, eg using a cache
			String currentORLinkFID = GISFunctions.getCoordinateRoad(currentCoord, currentPaveJ.getGeom().getCoordinate()).getRoadLinkID();
			String destORLinkFID = GISFunctions.getCoordinateRoad(destCoord, destPaveJ.getGeom().getCoordinate()).getRoadLinkID();
			
			for (RoadLink rl: SpaceBuilder.orRoadLinkGeography.getAllObjects()) {
				if (rl.getFID().contentEquals(currentORLinkFID)) {
					currentRoad = rl;
				}
				
				if (rl.getFID().contentEquals(destORLinkFID)) {
					destRoad = rl;
				}
			}
		} catch (RoutingException e) {
			throw e;
		}
		
		
		List<Junction> currentORJ = new ArrayList<Junction>();
		List<Junction> destORJ = new ArrayList<Junction>();
		
		// Loop through road junction geography to select start and end road nodes
		for (Junction rJ: SpaceBuilder.orJunctionGeography.getAllObjects()) {
			if(rJ.getFID().contentEquals(currentPaveJ.getjuncNodeID())) {
				currentORJ.add(rJ);
			}
			if(rJ.getFID().contentEquals(destPaveJ.getjuncNodeID())) {
				destORJ.add(rJ);
			}
		}
		
		assert currentORJ.size() == 1;
		assert destORJ.size() == 1;
			
		setRoadLinkRoute(currentORJ, destORJ);
		
		// If no roads in path this means start and end pavement junction are on the same link. Add current road to route manually
		Junction startORJ = null;
		Junction endORJ = null;
		
		if(this.roadsX.size()==0) {
			this.addToRoute(currentRoad, currentRoad.getEdge().getSpeed(), "current road manually added");
		}
		
		if (currentRoad.getFID().contentEquals(destRoad.getFID())) {			
			// Identify the road junctions that bookend the route by checking distance to start coord
			double dSource = currentCoord.distance(currentRoad.getEdge().getSource().getGeom().getCoordinate());
			double dTarget = currentCoord.distance(currentRoad.getEdge().getTarget().getGeom().getCoordinate());
			
			if (dSource < dTarget) {
				startORJ = currentRoad.getEdge().getSource();
				endORJ = currentRoad.getEdge().getTarget();
			}
			else {
				startORJ = currentRoad.getEdge().getTarget();
				endORJ = currentRoad.getEdge().getSource();
			}			
		}
		else {
			// In this case identify junctions the bookend route by finding which junction is not in the route end points

			// Find the OR junctions the pedestrian agent walks away from to start its route, and the junction that it walks to at the end of the route (the other junction connected to its current road link)
			startORJ = currentORJ.stream().filter(j -> !j.getFID().contentEquals(this.routeEndpoints[0].getFID())).collect(Collectors.toList()).get(0);
			
			// Repeat with the end OR junction to find the end pavement junction
			endORJ = destORJ.stream().filter(j -> !j.getFID().contentEquals(this.routeEndpoints[1].getFID())).collect(Collectors.toList()).get(0);
		}
		
		// Now check whether the links the origin and destination lie on are included in the route and if not add them in
		Junction startPaveJ = null;
		if (!this.roadsX.get(0).getFID().contentEquals(currentRoad.getFID())) {
			
			// Check whether the road the Origin belongs to, curentRoad, connects with the start of the road link route
			// If is does prefix the route with the currentRoad. This check is required to esure the links in the road link route form a path.
			for ( RepastEdge<Junction> re: SpaceBuilder.orRoadNetwork.getEdges(currentORJ.get(0))) {
				if (re.equals(currentRoad.getEdge())) {
					this.prefixToRoute(currentRoad, currentRoad.getEdge().getSpeed(), "Prefixing road link origin on manually");
					
					// Set the current pavement junction to be that connected to the nearest pavement junction but at the other end of the starting road link
					
					// Get road junction that is not near to the pavement junction
					Junction roadNode = currentRoad.getJunctions().stream().filter(n -> !n.getFID().contentEquals(currentPaveJ.getjuncNodeID())).collect(Collectors.toList()).get(0);
					
					startPaveJ = sameSideOppEndPavementJunction(pavementNetwork, currentPaveJ, roadNode);
				}
			}
		}
		
		if (startPaveJ==null){
			startPaveJ = currentPaveJ;
		}
		
		// Do the same for the end of the route
		Junction endPaveJ = null;
		if (!this.roadsX.get(this.roadsX.size()-1).getFID().contentEquals(destRoad.getFID())) {
			
			// Again, add the destination road to the route if it connects to end junction
			for ( RepastEdge<Junction> re: SpaceBuilder.orRoadNetwork.getEdges(destORJ.get(0))) {
				if (re.equals(destRoad.getEdge())) {
					this.addToRoute(destRoad, destRoad.getEdge().getSpeed(), "Adding road link destination on manually");
					
					// Update the dest pavement junction to be that connected to the nearest pavement junction but at the other end of the ending road link
					
					// Get road junction that is not near to the pavement junction
					Junction roadNode = destRoad.getJunctions().stream().filter(n -> !n.getFID().contentEquals(destPaveJ.getjuncNodeID())).collect(Collectors.toList()).get(0);
					
					// Find connecting pavement node that is associated with this road node
					endPaveJ = sameSideOppEndPavementJunction(pavementNetwork, destPaveJ, roadNode);
				}
			}
		}
		
		if (endPaveJ==null) {
			endPaveJ = destPaveJ;
		}
		
		// Check that a route has actually been created
		checkListSizes();
		
		setRouteParity();

		
		// Cache the route and route speeds
		// List<Coordinate> routeClone = Cloning.copy(theRoute);
		// LinkedHashMap<Coordinate, Double> routeSpeedsClone = Cloning.copy(this.routeSpeeds);
		// cachedRoute.setRoute(routeClone);
		// cachedRoute.setRouteSpeeds(routeSpeedsClone);

		// cachedRoute.setRoute(this.routeX, this.roadsX, this.routeSpeedsX, this.routeDescriptionX);
		// synchronized (Route.routeCache) {
		// // Same cached route is both value and key
		// Route.routeCache.put(cachedRoute, cachedRoute);
		// }
		// TempLogger.out("...Route cacheing new route with unique id " + cachedRoute.hashCode());
		/*
		LOGGER.log(Level.FINER, "Route Finished planning route for " + "with "
				+ this.routeX.size() + " coords in " + (0.000001 * (System.nanoTime() - time)) + "ms.");
				*/
		Junction[] routeEnds = {startPaveJ, endPaveJ};
		return routeEnds;
		
	}
	
	/*
	 * Find connecting pavement node that is associated with this road node and that is on the same side of the road
	 * 
	 * @param Network<Junction> pavementNetwork
	 * 		The pavement network
	 * @param Junction pJ
	 * 		The pavement junction to find an adjacent pavement junction to.
	 * @param Junction rJ
	 * 		The road junction indicating the opposite end of the road.
	 */
	public Junction sameSideOppEndPavementJunction(Network<Junction> pavementNetwork, Junction pJ, Junction rJ) {
		Junction adjacent=null;
		for(RepastEdge<Junction> e: pavementNetwork.getEdges(pJ)) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) e;
			if (ne.getRoadLink().getPedRLID().contentEquals("")) {
				
				Junction candidate = null;
				if (ne.getSource().getFID().contentEquals(pJ.getFID())) {
					candidate = ne.getTarget();
				}
				else {
					candidate = ne.getSource();
				}
				
				if (candidate.getjuncNodeID().contentEquals(rJ.getFID())) {
					adjacent = candidate;
				}
			}
		}
		return adjacent;
	}
	

	private void checkListSizes() {
		assert this.roadsX.size() > 0 && this.routeDescriptionX.size() == this.routeSpeedsX.size()
				&& this.roadsX.size() == this.routeDescriptionX.size() : this.roadsX.size()
				+ "," + this.routeDescriptionX.size() + "," + this.routeSpeedsX.size();
	}
	
	/**
	 * Loops through the edges in the shortest paths and adds their corresponding road links to the route list.
	 * 
	 * @param shortestPath
	 * @param startingJunction
	 *            The junction the path starts from, this is required so that the algorithm knows which road coordinate
	 *            to add first (could be first or last depending on the order that the road coordinates are stored
	 *            internally).
	 * @throws RoutingException
	 */
	private void addPathToRoute(List<RepastEdge<Junction>> shortestPath)
			throws RoutingException {
		double time = System.nanoTime();
		if (shortestPath.size() < 1) {
			// This could happen if the agent's destination is on the same road
			// as the origin
			return;
		}
		// Lock the currentAgent so that NetworkEdge obejcts know what speed to use (depends on transport available to
		// the specific agent).
		synchronized (GlobalVars.TRANSPORT_PARAMS.currentBurglarLock) {

			// Iterate over all edges in the route adding coords and weights as appropriate
			NetworkEdge<Junction> e;
			RoadLink r;
			// Loop through the edges in the shortest path and add the corresponding road links for each to the RoadNetworkRoute's list of road links
			for (int i = 0; i < shortestPath.size(); i++) {
				e = (NetworkEdge<Junction>) shortestPath.get(i);
				r = e.getRoadLink();
				
				double speed = e.getSpeed();
				if (speed < 1)
					speed = 1;
				
				addToRoute(r, speed, "getRouteBetweenJuctions - on road");
				
			// Finished!
			/*
			LOGGER.log(Level.FINER, "getRouteBetweenJunctions (" + (0.000001 * (System.nanoTime() - time)) + "ms");
			*/
			}
		}
		return;
	} // getRouteBetweenJunctions


	/**
	 * Used to add details to the route. This should be used rather than updating
	 * individual lists because it makes sure that all lists stay in sync
	 * 
	 * @param r
	 *            The road that the coordinate is part of
	 * @param speed
	 *            The speed that the road can be travelled along
	 * @param description
	 *            A description of why the coordinate has been added
	 */
	protected void addToRoute(RoadLink r, double speed, String description) {
		this.roadsX.add(r);
		this.routeSpeedsX.add(speed);
		this.routeDescriptionX.add(description);
	}
	
	/**
	 * Used to details to the start of the route route. This should be used rather than updating
	 * individual lists because it makes sure that all lists stay in sync
	 * 
	 * @param r
	 *            The road that the coordinate is part of
	 * @param speed
	 *            The speed that the road can be travelled along
	 * @param description
	 *            A description of why the coordinate has been added
	 */
	protected void prefixToRoute(RoadLink r, double speed, String description) {
		this.roadsX.add(0, r);
		this.routeSpeedsX.add(0, speed);
		this.routeDescriptionX.add(0, description);
	}
	
	protected void removeNextFromRoute() {
		removeFromRoute(0);
	}
	
	protected void removeFromRoute(int index) {
		this.roadsX.remove(index);
		this.routeSpeedsX.remove(index);
		this.routeDescriptionX.remove(index);
	}
	
	protected void updateRouteRoadDescription(RoadLink r, String newDesc) {
		int i = this.roadsX.indexOf(r);
		this.routeDescriptionX.set(i, newDesc);
	}
	
	
	
	/*
	 * Get the Road objects linked to a given Road Link ID. These are the vehicle and pedestrian polygons that for the street segment
	 * corresponding to the road link ID.
	 * 
	 * @param roadLinkID
	 * 			The ID of the road link to get the Road objects for.
	 * 
	 */
	public static synchronized List<Road> getRoadLinkRoads(String roadLinkID, Geography<Road> roadGeography) throws Exception {
		// double time = System.nanoTime();
		
		// Don't bother with the cache for now
		synchronized (roadLinkRoadsCacheLock) {
			if (roadLinkRoadsCache == null) {
				/*
				LOGGER.log(Level.FINE, "Route.getRoadLinkRoads called for first time, "
						+ "creating cache of all roads and the buildings which are on them ...");
						*/
				// Create a new cache object, this will be read from disk if
				// possible (which is why the getInstance() method is used
				// instead of the constructor.
				//String gisDir = IO.getProperty(GlobalVars.GISDataDirectory);
				//File vehcileRoadsFile = new File(gisDir + IO.getProperty(GlobalVars.VehicleRoadShapefile));
				//File pedestrianRoadsFile = new File(gisDir + IO.getProperty(GlobalVars.PedestrianRoadShapefile));
				//File serialisedLoc = new File(gisDir + IO.getProperty(GlobalVars.RoadLinkRoadsCache));

				roadLinkRoadsCache = RoadLinkRoadsCache.getInstance(roadGeography, vehcileRoadsFile, pedestrianRoadsFile, serialisedLoc);
			} // if not cached
		} // synchronized
		return roadLinkRoadsCache.get(roadLinkID);
	}
	
	/*
	 * Wrapper method that returns only the pedestrian priority roads associated to a particular road link.
	 * 
	 * @param roadLinkID
	 * 			The String ID of the road link to get pedestrian priority roads for.
	 *
	 * @return
	 * 		List of Roads
	 */
	public static synchronized List<Road> getRoadLinkPedestrianRoads(String roadLinkID, Geography<Road> roadGeography) throws Exception {
		List<Road> roads = getRoadLinkRoads(roadLinkID, roadGeography);
		
		// Remove non pedestrian roads
		for (int i = 0; i< roads.size(); i++) {
			if (!roads.get(i).getPriority().contentEquals("pedestrian")) {
				roads.remove(i);
			}
		}
		return roads;
	}
	
	/*
	 * Wrapper function for initialising road link coord and getting roads. This allow the initialisation to be tested in test cass.
	 */
	public static synchronized List<Road> getRoadLinkRoads(Geography<Road> g, File vRF, File pRF, File sL, String roadLinkID) throws Exception {
		
		RoadLinkRoadsCache roadLinkRoadsCache = RoadLinkRoadsCache.getInstance(g, vRF, pRF, sL);

		List<Road> roads = roadLinkRoadsCache.get(roadLinkID);
		
		return roads;
	}
	
	/*
	 * Wrapper method that returns only the pedestrian priority roads associated to a particular road link. 
	 * Different set of arguments allow the method to be tested because references to SpaceBuilder are not required.
	 * 
	 */
	public static synchronized List<Road> getRoadLinkPedestrianRoads(Geography<Road> g, File vRF, File pRF, File sL, String roadLinkID) throws Exception {
		List<Road> roads = getRoadLinkRoads(g, vRF, pRF, sL, roadLinkID);
		
		// Remove non pedestrian roads
		for (int i = 0; i< roads.size(); i++) {
			if (!roads.get(i).getPriority().contentEquals("pedestrian")) {
				roads.remove(i);
			}
		}
		return roads;
	}
	
	/*
	 * Use cache to find the nearest pavement junction to the input coordinate.
	 * 
	 * @param Coordinate c
	 * 			The coordinate to find nearest junction to.
	 * 
	 */
	public static synchronized Junction getNearestPavementJunctionToOD(Coordinate c) throws Exception {
		// double time = System.nanoTime();
		
		// Don't bother with the cache for now
		synchronized (odPaveJuncCacheLock) {
			if (odPaveJuncCache == null) {

				// Create a new cache object, this will be read from disk if
				// possible (which is why the getInstance() method is used
				// instead of the constructor.
				odPaveJuncCache = NearestPavementJunctionCoordCache.getInstance(SpaceBuilder.pedestrianDestinationGeography, odFile, SpaceBuilder.pavementJunctionGeography, pavementJunctionFile, paveJuncSeriealizedLoc, geomFac);
			} // if not cached
		} // synchronized
		return odPaveJuncCache.get(c, SpaceBuilder.pavementJunctionGeography);
	}
	
	/*
	 * A version of the methods for getting the nearest junction to a coordinate that is used when testing since can control which data files to use.
	 */
	public static synchronized Junction getNearestpavementJunctionToOD(Coordinate c, File odFile, File pavementJunctionFile, File paveJuncSeriealizedLoc) throws Exception {
		// Don't bother with the cache for now
		synchronized (odPaveJuncCacheLock) {
			if (odPaveJuncCache == null) {
				odPaveJuncCache = NearestPavementJunctionCoordCache.getInstance(SpaceBuilder.pedestrianDestinationGeography, odFile, SpaceBuilder.pavementJunctionGeography, pavementJunctionFile, paveJuncSeriealizedLoc, geomFac);
			} // if not cached
		} // synchronized
		return odPaveJuncCache.get(c, SpaceBuilder.pavementJunctionGeography);
	}

	/**
	 * Finds the shortest route between multiple origin and destination junctions. Will return the shortest path and
	 * also, via two parameters, can return the origin and destination junctions which make up the shortest route.
	 * 
	 * @param currentJunctions
	 *            An array of origin junctions
	 * @param destJunctions
	 *            An array of destination junctions
	 * @param routeEndpoints
	 *            An array of size 2 which can be used to store the origin (index 0) and destination (index 1) Junctions
	 *            which form the end points of the shortest route.
	 * @param longestShortest
	 * 				Boolean. Select whether or not to return the longest of the shortest routes between pairs of origin and destination junctions.
	 * 				Longest path is useful because it includes the road link of the origin and destination.
	 * @return the shortest route between the origin and destination junctions
	 * @throws Exception
	 */
	public List<RepastEdge<Junction>> getShortestRoute(Network<Junction> net, Iterable<Junction> currentJunctions, Iterable<Junction> destJunctions,
			Junction[] routeEndpoints, boolean longestShortest) throws Exception {
		double time = System.nanoTime();
		synchronized (GlobalVars.TRANSPORT_PARAMS.currentBurglarLock) {
			
			Double comparisonPathLength;
			int compVal = 0;
			if (longestShortest) {
				comparisonPathLength = Double.MIN_VALUE;
				compVal = +1;
			}
			else {
				comparisonPathLength = Double.MAX_VALUE;
				compVal = -1;
			}
			Double pathLength = 0.0;
			ShortestPath<Junction> p;
			List<RepastEdge<Junction>> shortestPath = null;
			for (Junction o : currentJunctions) {
				for (Junction d : destJunctions) {
					if (o == null || d == null) {
						/*
						LOGGER.log(Level.WARNING, "Route.getShortestRoute() error: either the destination or origin "
								+ "junction is null. This can be caused by disconnected roads. It's probably OK"
								+ "to ignore this as a route should still be created anyway.");
								*/
					} else {
						p = new ShortestPath<Junction>(net);
						pathLength = p.getPathLength(o,d);
						int comp = Integer.signum(pathLength.compareTo(comparisonPathLength));
						if (comp == compVal) {
							comparisonPathLength = pathLength;
							shortestPath = p.getPath(o,d);
							routeEndpoints[0] = o;
							routeEndpoints[1] = d;
						}
						// TODO See if the shortestpath bug has been fixed, would make this unnecessary
						// this removes the projection listener from the network. Projection listener has something to do with listening for events.
						p.finalize();
						p = null;
					} // if junc null
				} // for dest junctions
			} // for origin junctions
			if (shortestPath == null) {
				String debugString = "Route.getShortestRoute() could not find a route. Looking for the shortest route between :\n";
				for (Junction j : currentJunctions)
					debugString += "\t" + j.toString() + ", roads: " + j.getRoads().toString() + "\n";
				for (Junction j : destJunctions)
					debugString += "\t" + j.toString() + ", roads: " + j.getRoads().toString() + "\n";
				throw new RoutingException(debugString);
			}
			/*
			LOGGER.log(Level.FINER, "Route.getShortestRoute (" + (0.000001 * (System.nanoTime() - time))
					+ "ms) found shortest path " + "(length: " + shortestPathLength + ") from "
					+ routeEndpoints[0].toString() + " to " + routeEndpoints[1].toString());
					*/
			return shortestPath;
		} // synchronized
	}
	
	/**
	 * Finds the shortest route between multiple origin and destination junctions. Will return the shortest path and
	 * also, via two parameters, can return the origin and destination junctions which make up the shortest route.
	 * 
	 * @param currentJunctions
	 *            An array of origin junctions
	 * @param destJunctions
	 *            An array of destination junctions
	 * @param routeEndpoints
	 *            An array of size 2 which can be used to store the origin (index 0) and destination (index 1) Junctions
	 *            which form the end points of the shortest route.
	 * @throws Exception
	 */
	public List<RepastEdge<Junction>> getShortestRoute(Network<Junction> net, Iterable<Junction> currentJunctions, Iterable<Junction> destJunctions,
			Junction[] routeEndpoints) throws Exception {
		return getShortestRoute(net, currentJunctions, destJunctions, routeEndpoints, false);
	}
	
	/*
	 * Method to calculate the number of times a line between two coordinates intersects road links.
	 */
	public static int calculateNRouteIntersections(Coordinate o, Coordinate d, List<RoadLink> roadLinkRoute) {		
		// Create linestring connecting origin to destination
		Coordinate[] odCoords = {o, d};
		LineString ODLine = GISFunctions.lineStringGeometryFromCoordinates(odCoords);
		Geometry[] ODLineGeom = {ODLine};

		
		// Loop through road links in the rout and count number of times the ODLine intersects
		Geometry[] rlGeoms = new Geometry[roadLinkRoute.size()];
		for (int i = 0; i< roadLinkRoute.size(); i++) {
			RoadLink rl = roadLinkRoute.get(i);
			rlGeoms[i] = rl.getGeom();
		}
				
		// Count number of intersections
		int nIntersections = GISFunctions.calculateNIntersectionCoords(ODLineGeom, rlGeoms);

		return nIntersections;
	}
	
	/*
	 * Method to calculate the route parity. The parity relates to the topology of the road network route in
	 * relation to the start and end points. Considering the road network route as a path that divides the space,
	 * even parity means the origin and destination are in the same region, odd parity means they are in different regions.
	 * 
	 * This is calculated by counting the number of times a straight line from the origin to the destination intersects the road
	 * network path.
	 */
	public static int calculateRouteParity(Coordinate o, Coordinate d, List<RoadLink> roadLinkRoute) {
		int p;
		
		// Count number of intersections
		int nIntersections = calculateNRouteIntersections(o, d, roadLinkRoute);

		p = nIntersections % 2;
		return p;
	}
	
	/*
	 * Identify whether the final destination lies within the agent's planning horizon. If it does return the 
	 * parity of the straight line connecting current location to destination, indicating whether a primary crossing is required.
	 * Otherwise return null
	 * 
	 * @param Coordinate o
	 * 		The start coordinate for parity calculation
	 * @param Coordinate d
	 * 		The end coordinate for parity calculation
	 * @param List<RoadLink> sP
	 * 		The strategic path
	 * @param int pH
	 * 		The planning horizon of the agent. Indicates how many links of the strategic path the agent is able to plan into the future.
	 */
	public static Integer boundedRouteParity(Coordinate o, Coordinate d, List<RoadLink> sP, int pH) {
		Integer parity = null;
		
		// Destination lies within planning horizon
		if (sP.size() <= pH) {
			parity = calculateRouteParity(o, d, sP);
		}
		else {
			parity = null;
		}
		
		return parity;
	}
	
	/*
	 * Calculate the interior angle between the geometries of two road links
	 * 
	 * 
	 */

	private static void checkNotNull(Object... args) throws RoutingException {
		for (Object o : args) {
			if (o == null) {
				throw new RoutingException("An input argument is null");
			}
		}
		return;
	}


	private void printRoute() {
		StringBuilder out = new StringBuilder();
		out.append("Printing route. Current position in list is "
				+ this.currentPosition + " ('" + this.routeDescriptionX.get(this.currentPosition) + "')");
		for (int i = 0; i < this.roadsX.size(); i++) {
			out.append("\t"+ this.routeSpeedsX.get(i).toString() + "\t" 
			+ this.roadsX.get(i) + "\t" + this.routeDescriptionX.get(i));
		}
		/*
		LOGGER.info(out.toString());
		*/
	}

	
	/**
	 * Find the nearest object in the given geography to the coordinate.
	 * 
	 * @param <T>
	 * @param x
	 *            The coordinate to search from
	 * @param geography
	 *            The given geography to look through
	 * @param closestPoints
	 *            An optional List that will be populated with the closest points to x (i.e. the results of
	 *            <code>distanceOp.closestPoints()</code>.
	 * @param searchDist
	 *            The maximum distance to search for objects in. Small distances are more efficient but larger ones are
	 *            less likely to find no objects.
	 * @return The nearest object.
	 * @throws RoutingException
	 *             If an object cannot be found.
	 */
	public static synchronized <T> T findNearestObject(Coordinate x, Geography<T> geography,
			List<Coordinate> closestPoints, GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE searchDist)
			throws RoutingException {
		if (x == null) {
			throw new RoutingException("The input coordinate is null, cannot find the nearest object");
		}

		T nearestObject = SpatialIndexManager.findNearestObject(geography, x, closestPoints, searchDist);

		if (nearestObject == null) {
			throw new RoutingException("Couldn't find an object close to these coordinates:\n\t" + x.toString());
		} else {
			return nearestObject;
		}
	}
	
	/**
	 * Find the object in the given geography that intersect the coordinate.
	 * 
	 * @param <T>
	 * @param x
	 *            The coordinate to search from
	 * @param geography
	 *            The given geography to look through
	 * @return List of intersecting objects
	 * @throws RoutingException
	 *             If an object cannot be found.
	 */
	public static synchronized <T> List<T> findIntersectingObjects(Coordinate x, Geography<T> geography)
			throws RoutingException {
		if (x == null) {
			throw new RoutingException("The input coordinate is null, cannot find the nearest object");
		}

		List<T> intersectingObjects = SpatialIndexManager.findIntersectingObjects(geography, x);

		if (intersectingObjects.size() == 0) {
			throw new RoutingException("Couldn't find an object close to these coordinates:\n\t" + x.toString());
		} else {
			return intersectingObjects;
		}
	}

	/**
	 * Returns the angle of the vector from p0 to p1 relative to the x axis
	 * <p>
	 * The angle will be between -Pi and Pi. I got this directly from the JUMP program source.
	 * 
	 * @return the angle (in radians) that p0p1 makes with the positive x-axis.
	 */
	public static synchronized double angle(Coordinate p0, Coordinate p1) {
		double dx = p1.x - p0.x;
		double dy = p1.y - p0.y;

		return Math.atan2(dy, dx);
	}


	/**
	 * The coordinate the route is targeting
	 * 
	 * @return the destination
	 */
	public Coordinate getDestination() {
		return this.destination;
	}

	/**
	 * Maintain a cache of all coordinates which are part of a road segment. Store the coords and all the road(s) they
	 * are part of.
	 * 
	 * @param coord
	 *            The coordinate which should be part of a road geometry
	 * @return The road(s) which the coordinate is part of or null if the coordinate is not part of any road
	 */
	private List<RoadLink> getRoadFromCoordCache(Coordinate coord) {

		populateCoordCache(SpaceBuilder.orRoadLinkGeography); // Check the cache has been populated
		return coordCache.get(coord);
	}

	/**
	 * Test if a coordinate is part of a road segment.
	 * 
	 * @param coord
	 *            The coordinate which we want to test
	 * @return True if the coordinate is part of a road segment
	 */
	private boolean coordOnRoad(Coordinate coord) {
		populateCoordCache(SpaceBuilder.orRoadLinkGeography); // check the cache has been populated
		return coordCache.containsKey(coord);
	}

	private synchronized static void populateCoordCache(Geography<RoadLink> rlG) {

		double time = System.nanoTime();
		if (coordCache == null) { // Fist check cache has been created
			coordCache = new HashMap<Coordinate, List<RoadLink>>();
			/*
			LOGGER.log(Level.FINER,
					"Route.populateCoordCache called for first time, creating new cache of all Road coordinates.");
					*/
		}
		if (coordCache.size() == 0) { // Now popualte it if it hasn't already
										// been populated
			/*
			LOGGER.log(Level.FINER, "Route.populateCoordCache: is empty, creating new cache of all Road coordinates.");
			*/

			for (RoadLink r : rlG.getAllObjects()) {
				for (Coordinate c : r.getGeom().getCoordinates()) {
					if (coordCache.containsKey(c)) {
						coordCache.get(c).add(r);
					} else {
						List<RoadLink> l = new ArrayList<RoadLink>();
						l.add(r);
						// TODO Need to put *new* coordinate here? Not use
						// existing one in memory?
						coordCache.put(new Coordinate(c), l);
					}
				}
			}
			/*
			LOGGER.log(Level.FINER, "... finished caching all road coordinates (in " + 0.000001
					* (System.nanoTime() - time) + "ms)");
					*/
		}
	}

	public static void clearCaches() {
		if (coordCache != null)
			coordCache.clear();
			coordCache = null;
		if (odPaveJuncCache != null) {
			odPaveJuncCache.clear();
			odPaveJuncCache = null;
		}
		
	}
	
	public List<RoadLink> getRoadsX(){
		return (List<RoadLink>)this.roadsX;
	}
	
	public List<String> getRouteDescriptionX(){
		return this.routeDescriptionX;
	}
	
	public List<Double> getRouteSpeedsX(){
		return this.routeSpeedsX;
	}


	public int getRouteParity() {
		return routeParity;
	}


	public void setRouteParity() {
		int p = calculateRouteParity(this.origin, this.destination, this.roadsX);
		this.routeParity = p;
	}

	/**
	 * Will add the given buildings to the awareness space of the Burglar who is being controlled by this Route.
	 * 
	 * @param buildings
	 *            A list of buildings
	 */
	/*
	protected <T> void passedObject(T object, Class<T> clazz) {
		List<T> list = new ArrayList<T>(1);
		list.add(object);
	}
	*/

}


/* ************************************************************************ */

/**
 * Caches the nearest pavement junctions to each coordiante in the pedestrian ODs file. This is used by pedestrian agents when first identifying the nearest pavement 
 * and road junctions.
 * 
 *
 * This class can be serialised so that if the GIS data doesn't change it doesn't have to be re-calculated each time.
 * 
 * @author Obi Sargoni (adapted from Nick Malleon's code)
 */
class NearestPavementJunctionCoordCache implements Serializable {

	private static Logger LOGGER = Logger.getLogger(NearestPavementJunctionCoordCache.class.getName());

	private static final long serialVersionUID = 1L;
	private Hashtable<Coordinate, Junction> theCache; // The actual cache
	// Check that the road/building data hasn't been changed since the cache was
	// last created
	private File odFile;
	private File pavementJunctionFile;
	// The location that the serialised object might be found.
	private File serialisedLoc;
	// The time that this cache was created, can be used to check data hasn't
	// changed since
	private long createdTime;
	
	private GeometryFactory geomFac;

	private NearestPavementJunctionCoordCache(Geography<OD> odGeography, File odFile, Geography<Junction> pavementJunctionGeography, 
			File pavementJunctionFile, File serialisedLoc, GeometryFactory geomFac) throws Exception {

		this.odFile = odFile;
		this.pavementJunctionFile = pavementJunctionFile;
		this.serialisedLoc = serialisedLoc;
		this.theCache = new Hashtable<Coordinate, Junction>();
		this.geomFac = geomFac;

		LOGGER.log(Level.FINE, "NearestPavementJunctionCoordCache() creating new cache with data (and modification date):\n\t"
				+ this.odFile.getAbsolutePath() + " (" + new Date(this.odFile.lastModified()) + ") \n\t"
				+ this.pavementJunctionFile.getAbsolutePath() + " (" + new Date(this.pavementJunctionFile.lastModified()) + "):\n\t"
				+ this.serialisedLoc.getAbsolutePath());
		
		populateCache(odGeography, pavementJunctionGeography);
		this.createdTime = new Date().getTime();
		//serialise();
	}

	public void clear() {
		this.theCache.clear();
	}

	
	private void populateCache(Geography<OD> odGeography, Geography<Junction> pavementJunctionGeography)
			throws Exception {
		double time = System.nanoTime();
		theCache = new Hashtable<Coordinate, Junction>();
		// Iterate over every building and find the nearest road point
		for (OD od : odGeography.getAllObjects()) {
			Geometry odGeom = od.getGeom();
			Junction nearestJ = RoadNetworkRoute.findNearestObject(odGeom.getCoordinate(), pavementJunctionGeography, null,	GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE);
			theCache.put(odGeom.getCoordinate(), nearestJ);
		}// for Buildings
		LOGGER.log(Level.FINER, "Finished caching nearest roads (" + (0.000001 * (System.nanoTime() - time)) + "ms)");
	} // if nearestRoadCoordCache = null;

	/**
	 * 
	 * @param c
	 * @return
	 * @throws Exception
	 */
	public Junction get(Coordinate c, Geography<Junction> pavementJunctionGeography) throws Exception {
		if (c == null) {
			throw new Exception("NearestPavementJunctionCoordCache.get() error: the given coordinate is null.");
		}
		double time = System.nanoTime();
		Junction nearestJ = this.theCache.get(c);
		if (nearestJ != null) {
			LOGGER.log(Level.FINER, "NearestRoadCoordCache.get() (using cache) - ("
					+ (0.000001 * (System.nanoTime() - time)) + "ms)");
			return nearestJ;
		}
		// If get here then the coord is not in the cache, agent not starting their journey from an OD already cached.
		// Search for it manually. Search all junctions in the vicinity, looking for the junction which is nearest the agent.
		double minDist = Double.MAX_VALUE;
		Junction nearestJunc = null;
		Point coordGeom = this.geomFac.createPoint(c);

		// Note: could use an expanding envelope that starts small and gets bigger
		double bufferDist = GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE.dist;
		double bufferMultiplier = 1.0;
		Envelope searchEnvelope = coordGeom.buffer(bufferDist * bufferMultiplier).getEnvelopeInternal();
		StringBuilder debug = new StringBuilder(); // in case the operation fails

		for (Junction j : pavementJunctionGeography.getObjectsWithin(searchEnvelope)) {

			DistanceOp distOp = new DistanceOp(coordGeom, pavementJunctionGeography.getGeometry(j));
			double thisDist = distOp.distance();
			// BUG?: if an agent is on a really long road, the long road will not be found by getObjectsWithin because
			// it is not within the buffer
			debug.append("\tpavement junction ").append(j.toString()).append(" is ").append(thisDist).append(
					" distance away (at closest point). ");

			if (thisDist < minDist) {
				minDist = thisDist;
				nearestJunc = j;
			} // if thisDist < minDist
			debug.append("\n");

		} // for near junctions

		if (nearestJunc != null) {
			LOGGER.log(Level.FINER, "NearestPavementJunctionCoordCache.get() (not using cache) - ("
					+ (0.000001 * (System.nanoTime() - time)) + "ms)");
			return nearestJunc;
		}
		/* IF HERE THEN ERROR, PRINT DEBUGGING INFO */
		StringBuilder debugIntro = new StringBuilder(); // Some extra info for debugging
		debugIntro.append("NearestPavementJunctionCoordCache.get() error: couldn't find a junction to return.\n");
		Iterable<Junction> juncs = pavementJunctionGeography.getObjectsWithin(searchEnvelope);
		debugIntro.append("Looking for nearest road coordinate around ").append(c.toString()).append(".\n");
		debugIntro.append("pavementJunctionGeography.getObjectsWithin() returned ").append(
				SpaceBuilder.sizeOfIterable(juncs) + " roads, printing debugging info:\n");
		debugIntro.append(debug);
		throw new Exception(debugIntro.toString());

	}

	private void serialise() throws IOException {
		double time = System.nanoTime();
		FileOutputStream fos = null;
		ObjectOutputStream out = null;
		try {
			if (!this.serialisedLoc.exists())
				this.serialisedLoc.createNewFile();
			fos = new FileOutputStream(this.serialisedLoc);
			out = new ObjectOutputStream(fos);
			out.writeObject(this);
			out.close();
		} catch (IOException ex) {
			if (serialisedLoc.exists()) {
				// delete to stop problems loading incomplete file next time
				serialisedLoc.delete();
			}
			throw ex;
		}
		LOGGER.log(Level.FINE, "... serialised NearestPavementJunctionCoordCache to " + this.serialisedLoc.getAbsolutePath()
				+ " in (" + 0.000001 * (System.nanoTime() - time) + "ms)");
	}

	/**
	 * Used to create a new NearestPavementJunctionCoordCache object. This function is used instead of the constructor directly so
	 * that the class can check if there is a serialised version on disk already. If not then a new one is created and
	 * returned.
	 * 
	 * @param odGeography
	 * @param odFile
	 * @param pavementJunctionGeography
	 * @param pavemetnJunctionFile
	 * @param serialisedLoc
	 * @param geomFac
	 * @return
	 * @throws Exception
	 */
	public synchronized static NearestPavementJunctionCoordCache getInstance(Geography<OD> odGeography, File odFile,
			Geography<Junction> pavementJunctionGeography, File pavementJunctionFile, File serialisedLoc, GeometryFactory geomFac) throws Exception {
		double time = System.nanoTime();
		// See if there is a cache object on disk.
		if (serialisedLoc.exists()) {
			FileInputStream fis = null;
			ObjectInputStream in = null;
			NearestPavementJunctionCoordCache ncc = null;
			try {
				fis = new FileInputStream(serialisedLoc);
				in = new ObjectInputStream(fis);
				ncc = (NearestPavementJunctionCoordCache) in.readObject();
				in.close();

				// Check that the cache is representing the correct data and the
				// modification dates are ok
				if (!odFile.getAbsolutePath().equals(ncc.odFile.getAbsolutePath())
						|| !pavementJunctionFile.getAbsolutePath().equals(ncc.pavementJunctionFile.getAbsolutePath())
						|| odFile.lastModified() > ncc.createdTime || pavementJunctionFile.lastModified() > ncc.createdTime) {
					LOGGER.log(Level.FINE, "NearestPavementJunctionCoordCache, found serialised object but it doesn't match the "
							+ "data (or could have different modification dates), will create a new cache.");
				} else {
					LOGGER.log(Level.FINER, "NearestPavementJunctionCoordCache, found serialised cache, returning it (in "
							+ 0.000001 * (System.nanoTime() - time) + "ms)");
					return ncc;
				}
			} catch (IOException ex) {
				if (serialisedLoc.exists())
					serialisedLoc.delete(); // delete to stop problems loading incomplete file next tinme
				throw ex;
			} catch (ClassNotFoundException ex) {
				if (serialisedLoc.exists())
					serialisedLoc.delete();
				throw ex;
			}

		}

		// No serialised object, or got an error when opening it, just create a new one
		return new NearestPavementJunctionCoordCache(odGeography, odFile, pavementJunctionGeography, pavementJunctionFile, serialisedLoc, geomFac);
	}

}


/**
 * Caches the road objects associated to a road link ID. Usful when wanting the pedestrian or vehicle spaces connected to a particular road link.
 * <p>
 * This class can be serialised so that if the GIS data doesn't change it doesn't have to be re-calculated each time.
 * 
 * @author Obi Thompson Sargoni (original code by Nick Malleson)
 */
class RoadLinkRoadsCache implements Serializable {

	private static Logger LOGGER = Logger.getLogger(RoadLinkRoadsCache.class.getName());

	private static final long serialVersionUID = 1L;
	private Hashtable<String, List<Road>> theCache; // The actual cache
	// Check that the road/building data hasn't been changed since the cache was
	// last created
	private File vehicleRoadsFile;
	private File pedestrianRoadsFile;
	// The location that the serialised object might be found.
	private File serialisedLoc;
	// The time that this cache was created, can be used to check data hasn't
	// changed since
	private long createdTime;

	private RoadLinkRoadsCache(Geography<Road> roadGeography, File vehicleRoadsFile, File pedestrianRoadsFile, File serialisedLoc)
			throws Exception {

		this.vehicleRoadsFile = vehicleRoadsFile;
		this.pedestrianRoadsFile = pedestrianRoadsFile;
		this.serialisedLoc = serialisedLoc;
		this.theCache = new Hashtable<String, List<Road>>();

		LOGGER.log(Level.FINE, "RoadLinkRoadsCache() creating new cache with data (and modification date):\n\t"
				+ this.vehicleRoadsFile.getAbsolutePath() + " (" + new Date(this.vehicleRoadsFile.lastModified()) + "):\n\t"
				+ this.pedestrianRoadsFile.getAbsolutePath() + " (" + new Date(this.pedestrianRoadsFile.lastModified()) + "):\n\t"
				+ this.serialisedLoc.getAbsolutePath());
		
		// Don't populate the cache for now
		populateCache(roadGeography);
		this.createdTime = new Date().getTime();
		serialise();
	}

	public void clear() {
		this.theCache.clear();
	}

	
	private void populateCache(Geography<Road> roadGeography)
			throws Exception {
		double time = System.nanoTime();
		theCache = new Hashtable<String, List<Road>>();
		// Iterate over every road object and group them by their road link IDs
		for (Road r : roadGeography.getAllObjects()) {
			String rlID = r.getRoadLinkID();
			
			// Check if cache doesn't contain this id, create new list and add to the cache
			if(!theCache.containsKey(rlID)) {
				List<Road> roads = new ArrayList<Road>();
				theCache.put(rlID, roads);
			}
			
			// Add the road to the list of roads under this ID
			theCache.get(rlID).add(r);
		}// for Buildings
		LOGGER.log(Level.FINER, "Finished caching nearest roads (" + (0.000001 * (System.nanoTime() - time)) + "ms)");
	} // if nearestRoadCoordCache = null;

	/**
	 * 
	 * @param c
	 * @return
	 * @throws Exception
	 */
	public List<Road> get(String rlID) throws Exception {
		if (rlID == null) {
			throw new Exception("RoadNetworkRoute.RoadLinkRoadsCache.get() error: the given coordinate is null.");
		}
		double time = System.nanoTime();
		List<Road> roads = this.theCache.get(rlID);
		if (roads != null) {
			LOGGER.log(Level.FINER, "RoadLinkRoadsCache.get() (using cache) - ("
					+ (0.000001 * (System.nanoTime() - time)) + "ms)");
			return roads;
		}
		
		/* IF HERE THEN ERROR, PRINT DEBUGGING INFO */
		StringBuilder debugIntro = new StringBuilder(); // Some extra info for debugging
		debugIntro.append("Route.RoadLinkRoadsCache.get() error: couldn't find a list of roads to return.\n");
		throw new Exception(debugIntro.toString());

	}

	private void serialise() throws IOException {
		double time = System.nanoTime();
		FileOutputStream fos = null;
		ObjectOutputStream out = null;
		try {
			if (!this.serialisedLoc.exists())
				this.serialisedLoc.createNewFile();
			fos = new FileOutputStream(this.serialisedLoc);
			out = new ObjectOutputStream(fos);
			out.writeObject(this);
			out.close();
		} catch (IOException ex) {
			if (serialisedLoc.exists()) {
				// delete to stop problems loading incomplete file next time
				serialisedLoc.delete();
			}
			throw ex;
		}
		LOGGER.log(Level.FINE, "... serialised NearestRoadCoordCache to " + this.serialisedLoc.getAbsolutePath()
				+ " in (" + 0.000001 * (System.nanoTime() - time) + "ms)");
	}

	/**
	 * Used to create a new BuildingsOnRoadCache object. This function is used instead of the constructor directly so
	 * that the class can check if there is a serialised version on disk already. If not then a new one is created and
	 * returned.
	 * 
	 * @param roadGeography
	 * @param vehicleRoadsFile
	 * @param pedestrianRoadsFile
	 * @param serialisedLoc
	 * @return
	 * @throws Exception
	 */
	public synchronized static RoadLinkRoadsCache getInstance(Geography<Road> roadGeography, File vehicleRoadsFile, File pedestrianRoadsFile, File serialisedLoc) throws Exception {
		double time = System.nanoTime();
		// See if there is a cache object on disk.
		if (serialisedLoc.exists()) {
			FileInputStream fis = null;
			ObjectInputStream in = null;
			RoadLinkRoadsCache rlr = null;
			try {

				fis = new FileInputStream(serialisedLoc);
				in = new ObjectInputStream(fis);
				rlr = (RoadLinkRoadsCache) in.readObject();
				in.close();

				// Check that the cache is representing the correct data and the
				// modification dates are ok
				if (!vehicleRoadsFile.getAbsolutePath().equals(rlr.vehicleRoadsFile.getAbsolutePath())
						|| !pedestrianRoadsFile.getAbsolutePath().equals(rlr.pedestrianRoadsFile.getAbsolutePath())
						|| vehicleRoadsFile.lastModified() > rlr.createdTime || pedestrianRoadsFile.lastModified() > rlr.createdTime) {
					LOGGER.log(Level.FINE, "BuildingsOnRoadCache, found serialised object but it doesn't match the "
							+ "data (or could have different modification dates), will create a new cache.");
				} else {
					LOGGER.log(Level.FINER, "RoadLinkRoadsCache, found serialised cache, returning it (in "
							+ 0.000001 * (System.nanoTime() - time) + "ms)");
					return rlr;
				}
			} catch (IOException ex) {
				if (serialisedLoc.exists())
					serialisedLoc.delete(); // delete to stop problems loading incomplete file next tinme
				throw ex;
			} catch (ClassNotFoundException ex) {
				if (serialisedLoc.exists())
					serialisedLoc.delete();
				throw ex;
			}

		}

		// No serialised object, or got an error when opening it, just create a new one
		return new RoadLinkRoadsCache(roadGeography, vehicleRoadsFile, pedestrianRoadsFile, serialisedLoc);
	}

}


/**
 * Convenience class for creating deep copies of lists/maps (copies the values stored as well). Haven't made this
 * generic because need access to constructors to create new objects (e.g. new Coord(c))
 */
final class Cloning {

	public static List<Coordinate> copy(List<Coordinate> in) {

		List<Coordinate> out = new ArrayList<Coordinate>(in.size());
		for (Coordinate c : in) {
			// TODO Check this Coordinate constructor does what I expect it to
			out.add(new Coordinate(c));
		}
		return out;
	}

	// Not used now that route speeds are a list, not a map
	// public static LinkedHashMap<Coordinate, Double>
	// copy(LinkedHashMap<Coordinate, Double> in) {
	//
	// LinkedHashMap<Coordinate, Double> out = new LinkedHashMap<Coordinate,
	// Double>(in.size());
	// for (Coordinate c : in.keySet()) {
	// out.put(c, in.get(c));
	// }
	// return out;
	// }

}
