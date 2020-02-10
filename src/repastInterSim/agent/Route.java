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

package repastInterSim.agent;

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

import org.apache.commons.lang3.ArrayUtils;
import org.geotools.referencing.GeodeticCalculator;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.operation.distance.DistanceOp;

import cern.colt.Arrays;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.RepastEdge;
import repast.simphony.space.graph.ShortestPath;
import repastInterSim.environment.Cacheable;
import repastInterSim.environment.Destination;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
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
 * travelling along a transort route they will cover more ground.
 * </p>
 * 
 * @author Nick Malleson
 */
public class Route implements Cacheable {

	private static Logger LOGGER = Logger.getLogger(Route.class.getName());

	static {
		// Route.routeCache = new Hashtable<CachedRoute, CachedRoute>();
	}
	protected Geography<Object> geography;
	
	protected MobileAgent mA;

	protected Coordinate destination;

	/*
	 * The route consists of a list of coordinates which describe how to get to the destination. Each coordinate might
	 * have an attached 'speed' which acts as a multiplier and is used to indicate whether or not the agent is
	 * travelling along a transport route (i.e. if a coordinate has an attached speed of '2' the agent will be able to
	 * get to the next coordinate twice as fast as they would do if they were walking). The current position incicate
	 * where in the lists of coords the agent is up to. Other attribute information about the route can be included as
	 * separate arrays with indices that match those of the 'route' array below.
	 */
	private int currentPosition;
	protected List<Coordinate> routeX;
	protected List<Double> routeSpeedsX;
	/*
	 * This maps route coordinates to their containing Road, used so that when travelling we know which road/community
	 * the agent is on. private
	 */
	protected List<RoadLink> roadsX;

	// Record which function has added each coord, useful for debugging
	protected List<String> routeDescriptionX;
	

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
	 * Store a route once it has been created, might be used later (note that the same object acts as key and value).
	 */
	// TODO Re-think route caching, would be better to cache the whole Route object
	// private static volatile Map<CachedRoute, CachedRoute> routeCache;
	// /** Store a route distance once it has been created */
	// private static volatile Map<CachedRouteDistance, Double> routeDistanceCache;


	/**
	 * Create a new route object
	 * 
	 * @param geography
	 * 		The geography projection that the mobile agent this route belongs to is in
	 * @param mA
	 * 		The mobile agent this route belongs to
	 * @param destination
	 * 		The destination coordinate of the route
	 */
	public Route(Geography<Object> geography, MobileAgent mA, Coordinate destination) {
		this.geography = geography;
		this.mA = mA;
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
	public void setRoute() throws Exception {
		long time = System.nanoTime();

		this.routeX = new Vector<Coordinate>();
		this.roadsX = new Vector<RoadLink>();
		this.routeDescriptionX = new Vector<String>();
		this.routeSpeedsX = new Vector<Double>();
		
		// Don't create a route if at destination
		if (atDestination()) {
			return;
		}

		Coordinate currentCoord = mA.getLoc();		
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
				addToRoute(currentCoord, RoadLink.nullRoad, 1, "setRoute() initial");
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
			List<RoadLink> currItersectingRoads = Route.findIntersectingObjects(currentCoord, SpaceBuilder.roadLinkGeography);

			// For each of the road links the agent intersects, find which junctions can be accessed. ie is the target junction.
			// These are the junctions to input into the shortest path function
			Map<Junction, RoadLink> currentAccessMap = new HashMap<Junction, RoadLink>();
			for (RoadLink rl: currItersectingRoads) {
				NetworkEdge<Junction> e = rl.getEdge();
				currentAccessMap.put(e.getTarget(), rl);
			}
			
			// Find the road link the agent is currently on
			List<RoadLink> destIntersectingRoads = Route.findIntersectingObjects(destCoord, SpaceBuilder.roadLinkGeography);

			// For each of the road links the agent intersects, find which junctions the destCoord can be accessed from. ie is the source junction.
			// These are the junctions to input into the shortest path function
			Map<Junction, RoadLink> destAccessMap = new HashMap<Junction, RoadLink>();
			for (RoadLink rl: destIntersectingRoads) {
				NetworkEdge<Junction> e = rl.getEdge();
				destAccessMap.put(e.getSource(), rl);
			}
			
			/*
			
			Old method that found nearest road rather than intersecting road links. 
			Might still be useful if vehicle agents are ever not on roads to begin with.
			
			GeometryFactory geomFac = new GeometryFactory();
			Point currentPoint = geomFac.createPoint(currentCoord);
			Point destPoint = geomFac.createPoint(destCoord);
			 
			RoadLink currentRoad = Route.findNearestObject(currentCoord, SpaceBuilder.roadLinkGeography, null,
					GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE);
			
			// Find which junctions this road is connected to
			List<Junction> currentJunctions = currentRoad.getJunctions();
			
			for (int i = 0; i < currentJunctions.size(); i++) {
				Junction j = currentJunctions.get(i);
				List<RoadLink> jRoads = j.getRoads();
				
				for (int k = 0; k < jRoads.size(); k++) {
					RoadLink jR = jRoads.get(k);
					NetworkEdge<Junction> e = jR.getEdge();
					if ((jR.getGeom().intersects(currentPoint)) & (e.getTarget().equals(j))) {
						currentAccessMap.put(j, jR);
					}
				}
			}
			

			// Find the road that this coordinate is on
			RoadLink destRoad = Route.findNearestObject(destCoord, SpaceBuilder.roadLinkGeography, null,
					GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE);
			
			// Find which Junctions are connected to this road
			List<Junction> destJunctions = destRoad.getJunctions();
			
			// Find which junctions can be accessed by RoadLinks the agent is currently intersection
			// This time junction needs to be the source of the link so agent can travel from last junction in path to its final destination
			Map<Junction, RoadLink> destAccessMap = new HashMap<Junction, RoadLink>();
			for (int i = 0; i < destJunctions.size(); i++) {
				Junction j = destJunctions.get(i);
				List<RoadLink> jRoads = j.getRoads();
				
				for (int k = 0; k < jRoads.size(); k++) {
					RoadLink jR = jRoads.get(k);
					if ((jR.getGeom().intersects(destPoint)) & (jR.getEdge().getSource().equals(j))) {
						destAccessMap.put(j, jR);
					}
				}
			}
			*/
			
			
			/*
			 * Now have possible routes (2 origin junctions, 2 destination junctions max, less if directed roads don't allow junctions to be  accessed) need to pick which junctions
			 * form shortest route
			 */
			Junction[] routeEndpoints = new Junction[2];
			List<RepastEdge<Junction>> shortestPath = getShortestRoute(currentAccessMap.keySet(), destAccessMap.keySet(), routeEndpoints);
			// NetworkEdge<Junction> temp = (NetworkEdge<Junction>)
			// shortestPath.get(0);
			Junction currentJunction = routeEndpoints[0];
			Junction destJunction = routeEndpoints[1];

			/* Add the coordinates describing how to get to the nearest junction */
			List<Coordinate> tempCoordList = new Vector<Coordinate>();
			this.getCoordsAlongRoad(currentCoord, currentJunction.getGeom().getCoordinate(), currentAccessMap.get(currentJunction), true, tempCoordList);
			addToRoute(tempCoordList, currentAccessMap.get(currentJunction), 1, "getCoordsAlongRoad (toJunction)");

			/*
			 * Add the coordinates and speeds etc which describe how to move along the chosen path
			 */
			this.getRouteBetweenJunctions(shortestPath, currentJunction, false);

			/*
			 * Add the coordinates describing how to get from the final junction to the destination.
			 */

			tempCoordList.clear();
			this.getCoordsAlongRoad(destJunction.getGeom().getCoordinate(),destCoord, destAccessMap.get(destJunction), false, tempCoordList);
			addToRoute(tempCoordList, destAccessMap.get(destJunction), 1, "getCoordsAlongRoad (fromJunction)");

			if (!destinationOnRoad) {
				addToRoute(finalDestination, RoadLink.nullRoad, 1, "setRoute final");
			}

			// Check that a route has actually been created
			checkListSizes();

			// If the algorithm was better no coordinates would have been duplicated
			// removePairs();

			// Check lists are still the same size.
			checkListSizes();

		} catch (RoutingException e) {
			LOGGER.log(Level.SEVERE, "Route.setRoute(): Problem creating route for " + this.mA.toString()
					+ " going from " + currentCoord.toString() + " to " + this.destination.toString());
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

		LOGGER.log(Level.FINER, "Route Finished planning route for " + this.mA.toString() + "with "
				+ this.routeX.size() + " coords in " + (0.000001 * (System.nanoTime() - time)) + "ms.");

		// Finished, just check that the route arrays are all in sync
		assert this.roadsX.size() == this.routeX.size() && this.routeDescriptionX.size() == this.routeSpeedsX.size()
				&& this.roadsX.size() == this.routeDescriptionX.size();
	}
	
	/**
	 * Find a route from the origin to the destination. A route is a list of Coordinates which describe the route to a
	 * destination restricted to a road network. This algorithm differs from the setRoute() algorithm by only adding
	 * coordinates of junctions along the route to the route coordinates. This allows pedestrian agent freedom to navigate
	 * between junctions following a separate model of movement. The algorithm consists of three major parts:
	 * <ol>
	 * <li>Get from the current location (probably mid-point on a road) to the nearest junction</li>
	 * <li>Travel to the junction which is closest to our destination (using Dijkstra's shortest path)</li>
	 * <li>Move from the road to the destination</li>
	 * </ol>
	 * 
	 * @throws Exception
	 */
	public void setPedestrianRoute() throws Exception {
		long time = System.nanoTime();

		this.routeX = new Vector<Coordinate>();
		this.roadsX = new Vector<RoadLink>();
		this.routeDescriptionX = new Vector<String>();
		this.routeSpeedsX = new Vector<Double>();
		
		// Don't create a route if at destination
		if (atDestination()) {
			return;
		}

		Coordinate currentCoord = mA.getLoc();
		Coordinate destCoord = this.destination;



		// No route cached, have to create a new one (and cache it at the end).
		try {
			/*
			 * See if the current position and the destination are on road segments. If the destination is not on a road
			 * segment we have to move to the closest road segment, then onto the destination.
			 */
			boolean destinationOnRoad = true;
			Coordinate finalDestination = null;
			
			// Add the starting coordinate to the route
			//addToRoute(currentCoord, RoadLink.nullRoad, 1, "setRoute() initial");

			/*
			 * Find the nearest junctions to our current position (road endpoints)
			 */

			// Start by Finding the road that this coordinate is on
			/*
			 * TODO EFFICIENCY: often the agent will be creating a new route from a building so will always find the
			 * same road, could use a cache. Even better, could implement a cache in FindNearestObject() method!
			 */
			RoadLink currentRoad = Route.findNearestObject(currentCoord, SpaceBuilder.roadLinkGeography, null,
					GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE);
			// Find which Junction is closest to us on the road.
			List<Junction> currentJunctions = currentRoad.getJunctions();

			/* Find the nearest Junctions to our destination (road endpoints) */

			// Find the road that this coordinate is on
			RoadLink destRoad = Route.findNearestObject(destCoord, SpaceBuilder.roadLinkGeography, null,
					GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.LARGE);
			// Find which Junction connected to the edge is closest to the coordinate.
			List<Junction> destJunctions = destRoad.getJunctions();
			/*
			 * Now have four possible routes (2 origin junctions, 2 destination junctions) need to pick which junctions
			 * form shortest route
			 */
			Junction[] routeEndpoints = new Junction[2];
			List<RepastEdge<Junction>> shortestPath = getShortestRoute(currentJunctions, destJunctions, routeEndpoints);
			// NetworkEdge<Junction> temp = (NetworkEdge<Junction>)
			// shortestPath.get(0);
			Junction currentJunction = routeEndpoints[0];
			Junction destJunction = routeEndpoints[1];

			/*
			 * Add the coordinates and speeds etc which describe how to move along the chosen path
			 * 
			 * Replace this with a simple loop to get the junctions?
			 */
			this.getRouteBetweenJunctions(shortestPath, currentJunction, true);

			addToRoute(destCoord, RoadLink.nullRoad, 1, "setRoute final");

			// Check that a route has actually been created
			checkListSizes();

			// If the algorithm was better no coordinates would have been duplicated
			// removePairs();

			// Check lists are still the same size.
			checkListSizes();

		} catch (RoutingException e) {
			LOGGER.log(Level.SEVERE, "Route.setRoute(): Problem creating route for " + this.mA.toString()
					+ " going from " + currentCoord.toString() + " to " + this.destination.toString());
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

		LOGGER.log(Level.FINER, "Route Finished planning route for " + this.mA.toString() + "with "
				+ this.routeX.size() + " coords in " + (0.000001 * (System.nanoTime() - time)) + "ms.");

		// Finished, just check that the route arrays are all in sync
		assert this.roadsX.size() == this.routeX.size() && this.routeDescriptionX.size() == this.routeSpeedsX.size()
				&& this.roadsX.size() == this.routeDescriptionX.size();
	}
	
	

	private void checkListSizes() {
		assert this.roadsX.size() > 0 && this.roadsX.size() == this.routeX.size()
				&& this.routeDescriptionX.size() == this.routeSpeedsX.size()
				&& this.roadsX.size() == this.routeDescriptionX.size() : this.routeX.size() + "," + this.roadsX.size()
				+ "," + this.routeDescriptionX.size() + "," + this.routeSpeedsX.size();

	}

	/**
	 * Convenience function that can be used to add details to the route. This should be used rather than updating
	 * individual lists because it makes sure that all lists stay in sync
	 * 
	 * @param coord
	 *            The coordinate to add to the route
	 * @param nullRoad
	 *            The road that the coordinate is part of
	 * @param speed
	 *            The speed that the road can be travelled along
	 * @param description
	 *            A description of why the coordinate has been added
	 */
	protected void addToRoute(Coordinate coord, RoadLink nullRoad, double speed, String description) {
		this.routeX.add(coord);
		this.roadsX.add(nullRoad);
		this.routeSpeedsX.add(speed);
		this.routeDescriptionX.add(description);
	}

	/**
	 * A convenience for adding to the route that will add a number of coordinates with the same description, road and
	 * speed.
	 * 
	 * @param coord
	 *            A list of coordinates to add to the route
	 * @param road
	 *            The road that the coordinates are part of
	 * @param speed
	 *            The speed that the road can be travelled along
	 * @param description
	 *            A description of why the coordinates have been added
	 */
	protected void addToRoute(List<Coordinate> coords, RoadLink road, double speed, String description) {
		for (Coordinate c : coords) {
			this.routeX.add(c);
			this.roadsX.add(road);
			this.routeSpeedsX.add(speed);
			this.routeDescriptionX.add(description);
		}
	}
	
	protected void removeNextFromRoute() {
		removeFromRoute(0);
	}
	
	protected void removeFromRoute(int index) {
		this.routeX.remove(index);
		this.roadsX.remove(index);
		this.routeSpeedsX.remove(index);
		this.routeDescriptionX.remove(index);
	}
	
	protected void updateRouteCoordDescription(Coordinate c, String newDesc) {
		int i = this.routeX.indexOf(c);
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
				LOGGER.log(Level.FINE, "Route.getNearestRoadCoord called for first time, "
						+ "creating cache of all roads and the buildings which are on them ...");
				// Create a new cache object, this will be read from disk if
				// possible (which is why the getInstance() method is used
				// instead of the constructor.
				String gisDir = SpaceBuilder.getProperty(GlobalVars.GISDataDirectory);
				File buildingsFile = new File(gisDir + SpaceBuilder.getProperty(GlobalVars.BuildingShapefile));
				File roadsFile = new File(gisDir + SpaceBuilder.getProperty(GlobalVars.RoadShapefile));
				File serialisedLoc = new File(gisDir + SpaceBuilder.getProperty(GlobalVars.BuildingsRoadsCoordsCache));

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
	private List<RepastEdge<Junction>> getShortestRoute(Iterable<Junction> currentJunctions, Iterable<Junction> destJunctions,
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
						LOGGER.log(Level.WARNING, "Route.getShortestRoute() error: either the destination or origin "
								+ "junction is null. This can be caused by disconnected roads. It's probably OK"
								+ "to ignore this as a route should still be created anyway.");
					} else {
						p = new ShortestPath<Junction>(SpaceBuilder.roadNetwork);
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
			LOGGER.log(Level.FINER, "Route.getShortestRoute (" + (0.000001 * (System.nanoTime() - time))
					+ "ms) found shortest path " + "(length: " + shortestPathLength + ") from "
					+ routeEndpoints[0].toString() + " to " + routeEndpoints[1].toString());
			return shortestPath;
		} // synchronized
	}

	/**
	 * Calculates the coordinates required to move an agent from their current position to the destination along a given
	 * road. The algorithm to do this is as follows:
	 * <ol>
	 * <li>Starting from the destination coordinate, record each vertex and check inside the booundary of each line
	 * segment until the destination point is found.</li>
	 * <li>Return all but the last vertex, this is the route to the destination.</li>
	 * </ol>
	 * A boolean allows for two cases: heading towards a junction (the endpoint of the line) or heading away from the
	 * endpoint of the line (this function can't be used to go to two midpoints on a line)
	 * 
	 * @param currentCoord
	 * @param destinationCoord
	 * @param currentRoad
	 * @param toJunction
	 *            whether or not we're travelling towards or away from a Junction
	 * @param coordList
	 *            A list which will be populated with the coordinates that the agent should follow to move along the
	 *            road.
	 * @param roadList
	 *            A list of roads associated with each coordinate.
	 * @throws Exception
	 */
	private void getCoordsAlongRoad(Coordinate currentCoord, Coordinate destinationCoord, RoadLink currentRoad,
			boolean toJunction, List<Coordinate> coordList) throws RoutingException {

		Route.checkNotNull(currentCoord, destinationCoord, currentRoad, coordList);

		double time = System.nanoTime();
		Coordinate[] roadCoords = currentRoad.getGeom().getCoordinates();

		// Check that the either the destination or current coordinate are actually part of the road
		boolean currentCorrect = false, destinationCorrect = false;
		for (int i = 0; i < roadCoords.length; i++) {
			if (toJunction && destinationCoord.equals(roadCoords[i])) {
				destinationCorrect = true;
				break;
			} else if (!toJunction && currentCoord.equals(roadCoords[i])) {
				currentCorrect = true;
				break;
			}
		} // for

		if (!(destinationCorrect || currentCorrect)) {
			String roadCoordsString = "";
			for (Coordinate c : roadCoords)
				roadCoordsString += c.toString() + " - ";
			throw new RoutingException("Neigher the origin or destination nor the current"
					+ "coordinate are part of the road '" + currentRoad.toString() + "' (person '" + this.mA.toString()
					+ "').\n" + "Road coords: " + roadCoordsString + "\n" + "\tOrigin: " + currentCoord.toString()
					+ "\n" + "\tDestination: " + destinationCoord.toString()+ " )\n " + "Heading " + (toJunction ? "to" : "away from")
					+ " a junction, so " + (toJunction ? "destination" : "origin")
					+ " should be part of a road segment");
		}

		// Might need to reverse the order of the road coordinates
		if (toJunction && !destinationCoord.equals(roadCoords[roadCoords.length - 1])) {
			// If heading towards a junction, destination coordinate must be at end of road segment
			ArrayUtils.reverse(roadCoords);
		} else if (!toJunction && !currentCoord.equals(roadCoords[0])) {
			// If heading away form junction current coord must be at beginning of road segment
			ArrayUtils.reverse(roadCoords);
		}
		
		GeometryFactory geomFac = new GeometryFactory();
		Point destinationPointGeom = geomFac.createPoint(destinationCoord);
		Point currentPointGeom = geomFac.createPoint(currentCoord);
		// If still false at end then algorithm hasn't worked
		boolean foundAllCoords = false;
		search: for (int i = 0; i < roadCoords.length - 1; i++) {
			Coordinate[] segmentCoords = new Coordinate[] { roadCoords[i], roadCoords[i + 1] };
			// Draw a small buffer around the line segment and look for the coordinate within the buffer
			Geometry buffer = geomFac.createLineString(segmentCoords).buffer(GlobalVars.GEOGRAPHY_PARAMS.BUFFER_DISTANCE.SMALL.dist);
			if (!toJunction) {
				/* If heading away from a junction, keep adding road coords until we find the destination */
				coordList.add(roadCoords[i]);
				if (destinationPointGeom.within(buffer)) {
					coordList.add(destinationCoord);
					foundAllCoords = true;
					break search;
				}
			} else if (toJunction) {
				/*
				 * If heading towards a junction: find the curent coord, add it to the route, then add all the remaining
				 * coords which make up the road segment
				 */
				if (currentPointGeom.within(buffer)) {
					for (int j = i + 1; j < roadCoords.length; j++) {
						coordList.add(roadCoords[j]);
					}
					coordList.add(destinationCoord);
					foundAllCoords = true;
					break search;
				}
			}
		} // for
		if (foundAllCoords) {
			LOGGER.log(Level.FINER, "getCoordsAlongRoad (" + (0.000001 * (System.nanoTime() - time)) + "ms)");
			return;
		} else { // If we get here then the route hasn't been created
			// A load of debugging info
			String error = "Route: getCoordsAlongRoad: could not find destination coordinates "
					+ "along the road.\n\tHeading *" + (toJunction ? "towards" : "away from")
					+ "* a junction.\n\t Person: " + this.mA.toString() + ")\n\tRoad causing problems: " + currentRoad.toString()
					+ "\n\tRoad vertex coordinates: " + Arrays.toString(roadCoords);
			throw new RoutingException(error);

		}
	}

	private static void checkNotNull(Object... args) throws RoutingException {
		for (Object o : args) {
			if (o == null) {
				throw new RoutingException("An input argument is null");
			}
		}
		return;
	}

	/**
	 * Returns all the coordinates that describe how to travel along a path, restricted to road coordinates. In some
	 * cases the route wont have an associated road, this occurs if the route is part of a transport network. In this
	 * case just the origin and destination coordinates are added to the route.
	 * 
	 * @param shortestPath
	 * @param startingJunction
	 *            The junction the path starts from, this is required so that the algorithm knows which road coordinate
	 *            to add first (could be first or last depending on the order that the road coordinates are stored
	 *            internally).
	 * @return the coordinates as a mapping between the coord and its associated speed (i.e. how fast the agent can
	 *         travel to the next coord) which is dependent on the type of edge and the agent (e.g.
	 *         driving/walking/bus). LinkedHashMap is used to guarantee the insertion order of the coords is maintained.
	 * @throws RoutingException
	 */
	private void getRouteBetweenJunctions(List<RepastEdge<Junction>> shortestPath, Junction startingJunction, boolean junctionCoordsOnly)
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
			// Use sourceFirst to represent whether or not the edge's source does actually represent the start of the
			// edge (agent could be going 'forwards' or 'backwards' over edge
			boolean sourceFirst;
			for (int i = 0; i < shortestPath.size(); i++) {
				e = (NetworkEdge<Junction>) shortestPath.get(i);
				r = e.getRoadLink();
				
				// No coords in route yet, compare the source to the starting junction
				sourceFirst = (e.getSource().getGeom().getCoordinate().equals(r.getGeom().getCoordinates()[0])) ? true : false;

				/*
				 * Now add the coordinates describing how to move along the road. If there is no road associated with
				 * the edge (i.e. it is a transport route) then just add the source/dest coords. Note that the shared
				 * coordinates between two edges will be added twice, these must be removed later
								
				 * Get the speed that the agent will be able to travel along this edge (depends on the transport
				 * available to the agent and the edge). Some speeds will be < 1 if the agent shouldn't be using this
				 * edge but doesn't have any other way of getting to the destination. in these cases set speed to 1
				 * (equivalent to walking).
				 */
				double speed = e.getSpeed();
				if (speed < 1)
					speed = 1;

				if (junctionCoordsOnly == true) { // Or we only want to add the coordinates of junctions to the route
					if (sourceFirst) {
						this.addToRoute(e.getSource().getGeom().getCoordinate(), r, speed, "getRouteBetweenJunctions - junctions only");
						this.addToRoute(e.getTarget().getGeom().getCoordinate(), r, -1, "getRouteBetweenJunctions - junctions only");
						// (Note speed = -1 used because we don't know the weight to the next
						// coordinate - this can be removed later)
					} else {
						this.addToRoute(e.getTarget().getGeom().getCoordinate(), r, speed, "getRouteBetweenJunctions - junctions only");
						this.addToRoute(e.getSource().getGeom().getCoordinate(), r, -1, "getRouteBetweenJunctions - junctions only");
					}
				}
				else {
					// This edge is a road, add all the coords which make up its geometry
					Coordinate[] roadCoords = r.getGeom().getCoordinates();
					if (roadCoords.length < 2)
						throw new RoutingException("Route.getRouteBetweenJunctions: for some reason road " + "'"
								+ r.toString() + "' doesn't have at least two coords as part of its geometry ("
								+ roadCoords.length + ")");
					
					// Make sure the coordinates of the road are added in the correct order
					if (!sourceFirst) {
						ArrayUtils.reverse(roadCoords);
					}

					// Add all the road geometry's coords
					for (int j = 0; j < roadCoords.length; j++) {
						this.addToRoute(roadCoords[j], r, speed, "getRouteBetweenJuctions - on road");
						// (Note that last coord will have wrong weight)
					} // for roadCoords.length
				} // if road!=null
			}
			// Check all lists are still the same size.
			assert this.roadsX.size() == this.routeX.size()
					&& this.routeDescriptionX.size() == this.routeSpeedsX.size()
					&& this.roadsX.size() == this.routeDescriptionX.size();

			// Check all lists are still the same size.
			assert this.roadsX.size() == this.routeX.size()
					&& this.routeDescriptionX.size() == this.routeSpeedsX.size()
					&& this.roadsX.size() == this.routeDescriptionX.size();

			// Finished!
			LOGGER.log(Level.FINER, "getRouteBetweenJunctions (" + (0.000001 * (System.nanoTime() - time)) + "ms");
			return;
		} // synchronized
	} // getRouteBetweenJunctions

	/**
	 * This might not work for pedestrians since they might not move to the exact location of the destination.
	 * 
	 * Determine whether or not the person associated with this Route is at their destination. Compares their current
	 * coordinates to the destination coordinates (must be an exact match).
	 * 
	 * @return True if the person is at their destination
	 * @throws TransformException 
	 * @throws MismatchedDimensionException 
	 */
	public boolean atDestination() throws MismatchedDimensionException, TransformException {
		return SpaceBuilder.getAgentGeometry(geography, mA).getCoordinate().equals(this.destination);
	}


	private void printRoute() {
		StringBuilder out = new StringBuilder();
		out.append("Printing route (" + this.mA.toString() + "). Current position in list is "
				+ this.currentPosition + " ('" + this.routeDescriptionX.get(this.currentPosition) + "')");
		for (int i = 0; i < this.routeX.size(); i++) {
			out.append("\t(" + this.mA.toString() + ") " + this.routeX.get(i).toString() + "\t"
					+ this.routeSpeedsX.get(i).toString() + "\t" + this.roadsX.get(i) + "\t"
					+ this.routeDescriptionX.get(i));
		}
		LOGGER.info(out.toString());
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
			LOGGER.log(Level.FINER,
					"Route.populateCoordCache called for first time, creating new cache of all Road coordinates.");
		}
		if (coordCache.size() == 0) { // Now popualte it if it hasn't already
										// been populated
			LOGGER.log(Level.FINER, "Route.populateCoordCache: is empty, creating new cache of all Road coordinates.");

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

			LOGGER.log(Level.FINER, "... finished caching all road coordinates (in " + 0.000001
					* (System.nanoTime() - time) + "ms)");
		}
	}


	/**
	 * Calculate the distance (in meters) between two Coordinates, using the coordinate reference system that the
	 * roadGeography is using. For efficiency it can return the angle as well (in the range -0 to 2PI) if returnVals
	 * passed in as a double[2] (the distance is stored in index 0 and angle stored in index 1).
	 * 
	 * @param c1
	 * @param c2
	 * @param returnVals
	 *            Used to return both the distance and the angle between the two Coordinates. If null then the distance
	 *            is just returned, otherwise this array is populated with the distance at index 0 and the angle at
	 *            index 1.
	 * @return The distance between Coordinates c1 and c2.
	 */
	public static synchronized double distance(Coordinate c1, Coordinate c2, double[] returnVals) {
		// TODO check this now, might be different way of getting distance in new Simphony
		GeodeticCalculator calculator = new GeodeticCalculator(SpaceBuilder.roadLinkGeography.getCRS());
		calculator.setStartingGeographicPoint(c1.x, c1.y);
		calculator.setDestinationGeographicPoint(c2.x, c2.y);
		double distance = calculator.getOrthodromicDistance();
		if (returnVals != null && returnVals.length == 2) {
			returnVals[0] = distance;
			double angle = Math.toRadians(calculator.getAzimuth()); // Angle in range -PI to PI
			// Need to transform azimuth (in range -180 -> 180 and where 0 points north)
			// to standard mathematical (range 0 -> 360 and 90 points north)
			if (angle > 0 && angle < 0.5 * Math.PI) { // NE Quadrant
				angle = 0.5 * Math.PI - angle;
			} else if (angle >= 0.5 * Math.PI) { // SE Quadrant
				angle = (-angle) + 2.5 * Math.PI;
			} else if (angle < 0 && angle > -0.5 * Math.PI) { // NW Quadrant
				angle = (-1 * angle) + 0.5 * Math.PI;
			} else { // SW Quadrant
				angle = -angle + 0.5 * Math.PI;
			}
			returnVals[1] = angle;
		}
		return distance;
	}

	/**
	 * Converts a distance lat/long distance (e.g. returned by DistanceOp) to meters. The calculation isn't very
	 * accurate because (probably) it assumes that the distance is between two points that lie exactly on a line of
	 * longitude (i.e. one is exactly due north of the other). For this reason the value shouldn't be used in any
	 * calculations which is why it's returned as a String.
	 * 
	 * @param dist
	 *            The distance (as returned by DistanceOp) to convert to meters
	 * @return The approximate distance in meters as a String (to discourage using this approximate value in
	 *         calculations).
	 * @throws Exception
	 * @see com.vividsolutions.jts.operation.distance.DistanceOp
	 */
	public static synchronized String distanceToMeters(double dist) throws Exception {
		// Works by creating two coords (close to a randomly chosen object) which are a certain distance apart
		// then using similar method as other distance() function
		GeodeticCalculator calculator = new GeodeticCalculator(SpaceBuilder.roadLinkGeography.getCRS());
		Coordinate c1 = SpaceBuilder.junctionContext.getRandomObject().getGeom().getCoordinate();
		calculator.setStartingGeographicPoint(c1.x, c1.y);
		calculator.setDestinationGeographicPoint(c1.x, c1.y + dist);
		return String.valueOf(calculator.getOrthodromicDistance());
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
	
	public List<Coordinate> getRouteX(){
		return this.routeX;
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
	
	public Coordinate getRouteXCoordinate(Integer i) {
		return this.routeX.get(i);
	}
	
	public void removeRouteXCoordinate(Coordinate c) {
		routeX.remove(c);
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
		this.mA.addToMemory(list, clazz);
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

	private NearestRoadCoordCache(Geography<Destination> buildingEnvironment, File buildingsFile,
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

	
	private void populateCache(Geography<Destination> buildingEnvironment, Geography<RoadLink> roadLinkEnvironment)
			throws Exception {
		double time = System.nanoTime();
		theCache = new Hashtable<Coordinate, Coordinate>();
		// Iterate over every building and find the nearest road point
		for (Destination b : buildingEnvironment.getAllObjects()) {
			List<Coordinate> nearestCoords = new ArrayList<Coordinate>();
			Geometry bGeom = b.getGeom();
			Route.findNearestObject(bGeom.getCoordinate(), roadLinkEnvironment, nearestCoords,
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
	public synchronized static NearestRoadCoordCache getInstance(Geography<Destination> buildingEnv, File buildingsFile,
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
