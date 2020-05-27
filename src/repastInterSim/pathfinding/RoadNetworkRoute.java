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

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.operation.distance.DistanceOp;

import cern.colt.Arrays;
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
	 * Cache the nearest road Coordinate to every building for efficiency (agents usually/always need to get from the
	 * centroids of houses to/from the nearest road).
	 */
	private static volatile NearestRoadCoordCache nearestRoadCoordCache;
	/*
	 * Store which road every building is closest to. This is used to efficiently add buildings to the agent's awareness
	 * space
	 */
	//private static volatile BuildingsOnRoadCache buildingsOnRoadCache;
	// To stop threads competing for the cache:
	private static Object buildingsOnRoadCacheLock = new Object();
	/*
	 * Cache of road link ID to the road objects that have these IDs. Use to find the pedestrian and vehicle polygons links to a road link
	 */
	private static volatile RoadLinkRoadsCache roadLinkRoadsCache;
	// To stop threads competing for the cache:
	private static Object roadLinkRoadsCacheLock = new Object();


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
	}
	

	/**
	 * Find a route from the origin to the destination. A route is a list of Coordinates which describe the route to a
	 * destination restricted to a road network. The algorithm consists of three major parts:
	 * <ol>
	 * <li>Find out if the agent is on a road already, if not then move to the nearest road segment</li>
	 * <li>Get from the current location (probably mid-point on a road) to the nearest junction</li>
	 * <li>Travel to the junction which is closest to our destination (using Dijkstra's shortest path)</li>
	 * <li>Get from the final junction to the road which is nearest to the destination
	 * <li>
	 * <li>Move from the road to the destination</li>
	 * </ol>
	 * 
	 * @throws Exception
	 */
	public void setRoadLinkRoute() throws Exception {
		long time = System.nanoTime();

		this.roadsX = new Vector<RoadLink>();
		this.routeDescriptionX = new Vector<String>();
		this.routeSpeedsX = new Vector<Double>();
		
		Coordinate currentCoord = this.origin;		
		Coordinate destCoord = this.destination;


		// No route cached, have to create a new one (and cache it at the end).
		try {
			/*
			 * See if the current position and the destination are on road segments. If the destination is not on a road
			 * segment we have to move to the closest road segment, then onto the destination.
			 */
			boolean destinationOnRoad = true;
			Coordinate finalDestination = null;
			if (!coordOnRoad(currentCoord)) {
				/*
				 * Not on a road so the first coordinate to add to the route is the point on the closest road segment.
				 */
				currentCoord = getNearestRoadCoord(currentCoord);
			}
			if (!coordOnRoad(destCoord)) {
				/*
				 * Not on a road, so need to set the destination to be the closest point on a road, and set the
				 * destinationOnRoad boolean to false so we know to add the final dest coord at the end of the route
				 */
				destinationOnRoad = false;
				finalDestination = destCoord; // Added to route at end of alg.
				destCoord = getNearestRoadCoord(destCoord);
			}


			// Find the road link the agent is currently on
			List<RoadLink> currItersectingRoads = findIntersectingObjects(currentCoord, SpaceBuilder.roadLinkGeography);

			// For each of the road links the agent intersects, find which junctions can be accessed. ie is the target junction.
			// These are the junctions to input into the shortest path function
			Map<Junction, RoadLink> currentAccessMap = new HashMap<Junction, RoadLink>();
			for (RoadLink rl: currItersectingRoads) {
				NetworkEdge<Junction> e = rl.getEdge();
				currentAccessMap.put(e.getTarget(), rl);
			}
			
			// Find the road link the agent is currently on
			List<RoadLink> destIntersectingRoads = findIntersectingObjects(destCoord, SpaceBuilder.roadLinkGeography);

			// For each of the road links the agent intersects, find which junctions the destCoord can be accessed from. ie is the source junction.
			// These are the junctions to input into the shortest path function
			Map<Junction, RoadLink> destAccessMap = new HashMap<Junction, RoadLink>();
			for (RoadLink rl: destIntersectingRoads) {
				NetworkEdge<Junction> e = rl.getEdge();
				destAccessMap.put(e.getSource(), rl);
			}
			
			
			/*
			 * Now have possible routes (2 origin junctions, 2 destination junctions max, less if directed roads don't allow junctions to be  accessed) need to pick which junctions
			 * form shortest route
			 */
			Junction[] routeEndpoints = new Junction[2];
			List<RepastEdge<Junction>> shortestPath = getShortestRoute(SpaceBuilder.roadNetwork, currentAccessMap.keySet(), destAccessMap.keySet(), routeEndpoints);

			Junction currentJunction = routeEndpoints[0];
			Junction destJunction = routeEndpoints[1];

			/*
			 * Add the road links that make up the shortest path to the class attribute lists
			 */
			this.addPathToRoute(shortestPath, currentJunction);


			// Check that a route has actually been created
			checkListSizes();
			
			setRouteParity();

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
	private void addPathToRoute(List<RepastEdge<Junction>> shortestPath, Junction startingJunction)
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
	
	/**
	 * Find the nearest coordinate which is part of a Road. Returns the coordinate which is actually the closest to the
	 * given coord, not just the corner of the segment which is closest. Uses the DistanceOp class which finds the
	 * closest points between two geometries.
	 * <p>
	 * When first called, the function will populate the 'nearestRoadCoordCache' which calculates where the closest road
	 * coordinate is to each building. The agents will commonly start journeys from within buildings so this will
	 * improve efficiency.
	 * </p>
	 * 
	 * @param inCoord
	 *            The coordinate from which to find the nearest road coordinate
	 * @return the nearest road coordinate
	 * @throws Exception
	 */
	private synchronized Coordinate getNearestRoadCoord(Coordinate inCoord) throws Exception {
		// double time = System.nanoTime();
		
		// Don't bother with the cache for now
		synchronized (buildingsOnRoadCacheLock) {
			if (nearestRoadCoordCache == null) {
				/*
				LOGGER.log(Level.FINE, "Route.getNearestRoadCoord called for first time, "
						+ "creating cache of all roads and the buildings which are on them ...");
						*/
				// Create a new cache object, this will be read from disk if
				// possible (which is why the getInstance() method is used
				// instead of the constructor.
				String gisDir = IO.getProperty(GlobalVars.GISDataDirectory);
				File buildingsFile = new File(gisDir + IO.getProperty(GlobalVars.BuildingShapefile));
				File roadsFile = new File(gisDir + IO.getProperty(GlobalVars.RoadShapefile));
				File serialisedLoc = new File(gisDir + IO.getProperty(GlobalVars.BuildingsRoadsCoordsCache));

				nearestRoadCoordCache = NearestRoadCoordCache.getInstance(SpaceBuilder.vehicleDestinationGeography,
						buildingsFile, SpaceBuilder.roadLinkGeography, roadsFile, serialisedLoc, new GeometryFactory());
			} // if not cached
		} // synchronized
		return nearestRoadCoordCache.get(inCoord);
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
	 *            which form the endpoints of the shortest route.
	 * @return the shortest route between the origin and destination junctions
	 * @throws Exception
	 */
	public List<RepastEdge<Junction>> getShortestRoute(Network<Junction> net, Iterable<Junction> currentJunctions, Iterable<Junction> destJunctions,
			Junction[] routeEndpoints) throws Exception {
		double time = System.nanoTime();
		synchronized (GlobalVars.TRANSPORT_PARAMS.currentBurglarLock) {
			
			double shortestPathLength = Double.MAX_VALUE;
			double pathLength = 0;
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
						if (pathLength < shortestPathLength) {
							shortestPathLength = pathLength;
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
	
	/*
	 * Method to calculate the route parity. The parity relates to the topology of the road network route in
	 * relation to the start and end points. Considering the road network route as a path that divides the space,
	 * even parity means the origin and destination are in the same region, odd parity means they are in different regions.
	 * 
	 * This is calculated by counting the number of times a straight line from the origin to the destination intersects the road
	 * network path.
	 */
	public int calculateRouteParity(Coordinate o, Coordinate d, List<RoadLink> roadLinkRoute) {
		int p;
		
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

		p = nIntersections % 2;
		return p;
	}

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

		populateCoordCache(); // Check the cache has been populated
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
		populateCoordCache(); // check the cache has been populated
		return coordCache.containsKey(coord);
	}

	private synchronized static void populateCoordCache() {

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

			for (RoadLink r : SpaceBuilder.roadLinkContext.getObjects(RoadLink.class)) {
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



	public void clearCaches() {
		if (coordCache != null)
			coordCache.clear();
		if (nearestRoadCoordCache != null) {
			nearestRoadCoordCache.clear();
			nearestRoadCoordCache = null;
		}
		
		// Ignoring building cache for now
		/*
		if (buildingsOnRoadCache != null) {
			buildingsOnRoadCache.clear();
			buildingsOnRoadCache = null;
		}
		*/
	}
	
	public List<RoadLink> getRoadsX(){
		return this.roadsX;
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
 * Caches the nearest road Coordinate to every building for efficiency (agents usually/always need to get from the
 * centroids of houses to/from the nearest road).
 * <p>
 * This class can be serialised so that if the GIS data doesn't change it doesn't have to be re-calculated each time.
 * 
 * @author Nick Malleson
 */
class NearestRoadCoordCache implements Serializable {

	private static Logger LOGGER = Logger.getLogger(NearestRoadCoordCache.class.getName());

	private static final long serialVersionUID = 1L;
	private Hashtable<Coordinate, Coordinate> theCache; // The actual cache
	// Check that the road/building data hasn't been changed since the cache was
	// last created
	private File buildingsFile;
	private File roadsFile;
	// The location that the serialised object might be found.
	private File serialisedLoc;
	// The time that this cache was created, can be used to check data hasn't
	// changed since
	private long createdTime;

	private GeometryFactory geomFac;

	private NearestRoadCoordCache(Geography<OD> buildingEnvironment, File buildingsFile,
			Geography<RoadLink> roadLinkEnvironment, File roadsFile, File serialisedLoc, GeometryFactory geomFac)
			throws Exception {

		this.buildingsFile = buildingsFile;
		this.roadsFile = roadsFile;
		this.serialisedLoc = serialisedLoc;
		this.theCache = new Hashtable<Coordinate, Coordinate>();
		this.geomFac = geomFac;

		LOGGER.log(Level.FINE, "NearestRoadCoordCache() creating new cache with data (and modification date):\n\t"
				+ this.buildingsFile.getAbsolutePath() + " (" + new Date(this.buildingsFile.lastModified()) + ") \n\t"
				+ this.roadsFile.getAbsolutePath() + " (" + new Date(this.roadsFile.lastModified()) + "):\n\t"
				+ this.serialisedLoc.getAbsolutePath());
		
		// Don't populate the cache for now
		populateCache(buildingEnvironment, roadLinkEnvironment);
		this.createdTime = new Date().getTime();
		serialise();
	}

	public void clear() {
		this.theCache.clear();
	}

	
	private void populateCache(Geography<OD> buildingEnvironment, Geography<RoadLink> roadLinkEnvironment)
			throws Exception {
		double time = System.nanoTime();
		theCache = new Hashtable<Coordinate, Coordinate>();
		// Iterate over every building and find the nearest road point
		for (OD b : buildingEnvironment.getAllObjects()) {
			List<Coordinate> nearestCoords = new ArrayList<Coordinate>();
			Geometry bGeom = b.getGeom();
			RoadNetworkRoute.findNearestObject(bGeom.getCoordinate(), roadLinkEnvironment, nearestCoords,
					GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE);
			// Two coordinates returned by closestPoints(), need to find the one
			// which isn't the building coord
			Coordinate nearestPoint = null;
			for (Coordinate c : nearestCoords) {
				if (!c.equals(bGeom.getCoordinate())) {
					nearestPoint = c;
					break;
				}
			} // for nearestCoords
			if (nearestPoint == null) {
				throw new Exception("Route.getNearestRoadCoord() error: couldn't find a road coordinate which "
						+ "is close to building " + b.toString());
			}
			theCache.put(bGeom.getCoordinate(), nearestPoint);
		}// for Buildings
		LOGGER.log(Level.FINER, "Finished caching nearest roads (" + (0.000001 * (System.nanoTime() - time)) + "ms)");
	} // if nearestRoadCoordCache = null;

	/**
	 * 
	 * @param c
	 * @return
	 * @throws Exception
	 */
	public Coordinate get(Coordinate c) throws Exception {
		if (c == null) {
			throw new Exception("Route.NearestRoadCoordCache.get() error: the given coordinate is null.");
		}
		double time = System.nanoTime();
		Coordinate nearestCoord = this.theCache.get(c);
		if (nearestCoord != null) {
			LOGGER.log(Level.FINER, "NearestRoadCoordCache.get() (using cache) - ("
					+ (0.000001 * (System.nanoTime() - time)) + "ms)");
			return nearestCoord;
		}
		// If get here then the coord is not in the cache, agent not starting their journey from a house, search for
		// it manually. Search all roads in the vicinity, looking for the point which is nearest the person
		double minDist = Double.MAX_VALUE;
		Coordinate nearestPoint = null;
		Point coordGeom = this.geomFac.createPoint(c);

		// Note: could use an expanding envelope that starts small and gets bigger
		double bufferDist = GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE.dist;
		double bufferMultiplier = 1.0;
		Envelope searchEnvelope = coordGeom.buffer(bufferDist * bufferMultiplier).getEnvelopeInternal();
		StringBuilder debug = new StringBuilder(); // incase the operation fails

		for (RoadLink r : SpaceBuilder.roadLinkGeography.getObjectsWithin(searchEnvelope)) {

			DistanceOp distOp = new DistanceOp(coordGeom, SpaceBuilder.roadLinkGeography.getGeometry(r));
			double thisDist = distOp.distance();
			// BUG?: if an agent is on a really long road, the long road will not be found by getObjectsWithin because
			// it is not within the buffer
			debug.append("\troad ").append(r.toString()).append(" is ").append(thisDist).append(
					" distance away (at closest point). ");

			if (thisDist < minDist) {
				minDist = thisDist;
				Coordinate[] closestPoints = distOp.nearestPoints();
				// Two coordinates returned by closestPoints(), need to find the
				// one which isn''t the coord parameter
				debug.append("Closest points (").append(closestPoints.length).append(") are: ").append(
						Arrays.toString(closestPoints));
				nearestPoint = (c.equals(closestPoints[0])) ? closestPoints[1] : closestPoints[0];
				debug.append("Nearest point is ").append(nearestPoint.toString());
				nearestPoint = (c.equals(closestPoints[0])) ? closestPoints[1] : closestPoints[0];
			} // if thisDist < minDist
			debug.append("\n");

		} // for nearRoads

		if (nearestPoint != null) {
			LOGGER.log(Level.FINER, "NearestRoadCoordCache.get() (not using cache) - ("
					+ (0.000001 * (System.nanoTime() - time)) + "ms)");
			return nearestPoint;
		}
		/* IF HERE THEN ERROR, PRINT DEBUGGING INFO */
		StringBuilder debugIntro = new StringBuilder(); // Some extra info for debugging
		debugIntro.append("Route.NearestRoadCoordCache.get() error: couldn't find a coordinate to return.\n");
		Iterable<RoadLink> roads = SpaceBuilder.roadLinkGeography.getObjectsWithin(searchEnvelope);
		debugIntro.append("Looking for nearest road coordinate around ").append(c.toString()).append(".\n");
		debugIntro.append("roadLinkEnvironment.getObjectsWithin() returned ").append(
				SpaceBuilder.sizeOfIterable(roads) + " roads, printing debugging info:\n");
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
		LOGGER.log(Level.FINE, "... serialised NearestRoadCoordCache to " + this.serialisedLoc.getAbsolutePath()
				+ " in (" + 0.000001 * (System.nanoTime() - time) + "ms)");
	}

	/**
	 * Used to create a new BuildingsOnRoadCache object. This function is used instead of the constructor directly so
	 * that the class can check if there is a serialised version on disk already. If not then a new one is created and
	 * returned.
	 * 
	 * @param buildingEnv
	 * @param buildingsFile
	 * @param roadEnv
	 * @param roadsFile
	 * @param serialisedLoc
	 * @param geomFac
	 * @return
	 * @throws Exception
	 */
	public synchronized static NearestRoadCoordCache getInstance(Geography<OD> buildingEnv, File buildingsFile,
			Geography<RoadLink> roadEnv, File roadsFile, File serialisedLoc, GeometryFactory geomFac) throws Exception {
		double time = System.nanoTime();
		// See if there is a cache object on disk.
		if (serialisedLoc.exists()) {
			FileInputStream fis = null;
			ObjectInputStream in = null;
			NearestRoadCoordCache ncc = null;
			try {

				fis = new FileInputStream(serialisedLoc);
				in = new ObjectInputStream(fis);
				ncc = (NearestRoadCoordCache) in.readObject();
				in.close();

				// Check that the cache is representing the correct data and the
				// modification dates are ok
				if (!buildingsFile.getAbsolutePath().equals(ncc.buildingsFile.getAbsolutePath())
						|| !roadsFile.getAbsolutePath().equals(ncc.roadsFile.getAbsolutePath())
						|| buildingsFile.lastModified() > ncc.createdTime || roadsFile.lastModified() > ncc.createdTime) {
					LOGGER.log(Level.FINE, "BuildingsOnRoadCache, found serialised object but it doesn't match the "
							+ "data (or could have different modification dates), will create a new cache.");
				} else {
					LOGGER.log(Level.FINER, "NearestRoadCoordCache, found serialised cache, returning it (in "
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
		return new NearestRoadCoordCache(buildingEnv, buildingsFile, roadEnv, roadsFile, serialisedLoc, geomFac);
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
