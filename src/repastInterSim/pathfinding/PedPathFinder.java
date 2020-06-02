package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.geotools.coverage.grid.GridCoverage2D;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;

import repast.simphony.space.gis.Geography;
import repastInterSim.agent.GridRoute;
import repastInterSim.agent.Ped;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.OD;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class PedPathFinder {
	
	private OD origin;
	private OD destination;
	private Geography<Object> geography;
	
	private List<RoadLink> strategicPath;
	private GridRoute tacticalPath = new GridRoute();
	
	private Coordinate nextTacticalPathCoord;

	private Coordinate nextCrossingCoord;
	private RoadLink currentRoadLink;

	public PedPathFinder(Geography<Object> g, OD o, OD d) {
		this.geography = g;
		this.origin = o;
		this.destination = d;
		
		planStrategicPath();
	}
	
	/**
	 * Initialise a new road link routing object that can be used to find a path along a topological road network.
	 * Use this to identify the shortest path through the network and assign this path to this classes' strategic path attribute.
	 * 
	 */
	public void planStrategicPath() {
		// Initialise road network route - needs to ne non-directed for pedestrians! fix this
		RoadNetworkRoute rnr = new RoadNetworkRoute(origin.getGeom().getCoordinate(), destination.getGeom().getCoordinate());
		
		// Find shortest path using road network route
		try {
			rnr.setRoadLinkRoute();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get path of road links and set this as the strategic path
		this.strategicPath = rnr.getRoadsX();
		
		// The strategic path can be empty when origin and destination are very close
		if(this.strategicPath.isEmpty() == false) {
			this.currentRoadLink = this.strategicPath.get(0);
		}

	}
	
	/*
	 * Set the next tactical coordinate. If there are no coordinates in the tactical path, plan a new tactical path. Otherwise get the next coordinate
	 * in the tactical path.
	 * 
	 * @param gSPM
	 * 			The map from integers used to indicate the road user priority of grid cells to the perceived cost of moving through those grid cells. Used for routing on a grid.
	 * @param tacticalOriginCoord
	 * 			The start coordinate of the tactical route.
	 */
	public void updateTacticalPathCoordinate(HashMap<Integer, Double> gSPM, Coordinate tacticalOriginCoord) {
		
		// If reached the end of one section of the route, or if route has just been created, need to produce next set of route coordinates.
		if(this.tacticalPath.getRouteX().size() == 0) {
			Coordinate tacticalDestCoord = tacticalDestinationCoordinate(tacticalOriginCoord);

			// Both use the roadLinkCoordX[0] to set, consider passing in as parameter?
			planTacticaPath(gSPM, tacticalOriginCoord, tacticalDestCoord);
		}
		
    	// If crossing coord is null, check the route for upcoming crossing coord
		// Importantly, do this before any coords are removed from the route
    	if(this.nextCrossingCoord == null) {
    		this.nextCrossingCoord = this.tacticalPath.getNextRouteCrossingCoord();
    	}
    	
		this.nextTacticalPathCoord = this.tacticalPath.getRouteX().get(0);
		this.tacticalPath.removeNextFromRoute();
    }
	
	/*
	 * Plan a tactical level path from an origin coordinate to a destination coordinate. The tactical path is planned using the
	 * GridRoute class.
	 * 
	 * @param gSPM
	 * 			The map from integers used to indicate the road user priority of grid cells to the perceived cost of moving through those grid cells. Used for routing on a grid.
	 * @param tacticalOriginCoord
	 * 			The start coordinate of the tactical route.
	 * @param tacticalDestCoord
	 * 			The end coordinate of the tactical route.
	 */
	public void planTacticaPath(HashMap<Integer, Double> gSPM, Coordinate tacticalOriginCoord, Coordinate tacticalDestCoord) {
		
		GridCoverage2D grid = this.geography.getCoverage(GlobalVars.CONTEXT_NAMES.BASE_COVERAGE);
		GridRoute tP = new GridRoute(grid, gSPM, tacticalOriginCoord, tacticalDestCoord, true);
		
    	// Get updated set of route coords to follow to next road link coordinate
		tP.setGroupedGridPath();
    	
    	// Adds coordinates to route from next section of grid path
		tP.setNextRouteSection();
		
		this.tacticalPath= tP;
		
		// Remove this road link from the strategic path, no longer the next road link
		if (this.strategicPath.isEmpty() == false) {
			this.strategicPath.remove(0);
		}
	}
	
	/*
	 * Identify the destination coordinate for the tactical path. This is currently defined as the coordinate of the end junction of the
	 * current road link in the road link path, if the agent has reached the final road link in the path, the final destination.
	 */
	private Coordinate tacticalDestinationCoordinate() {
		Coordinate tacticalDestCoord;
		if (this.strategicPath.isEmpty() | this.strategicPath.size() == 1) {
			tacticalDestCoord = this.destination.getGeom().getCoordinate();
		}
		else {
			this.currentRoadLink = this.strategicPath.get(0);
			
			// Need to use the network edge associated to the road link to make sure the end of the road link is used as the destination
			tacticalDestCoord = this.currentRoadLink.getEdge().getTarget().getGeom().getCoordinate();
		}
		
		return tacticalDestCoord;
	}

	public static HashMap<String, List<Coordinate>> getTacticalDestinationCoodinateOptions(Coordinate oC, List<Road> cPR, List<RoadLink> sPS, Geography<PedObstruction> obstructGeography){
		
		// Find which road involved crossing the road link
		HashMap<String, List<Coordinate>> tacticalDestinationOptions = new HashMap<String, List<Coordinate>>();
		tacticalDestinationOptions.put("cross", new ArrayList<Coordinate>());
		tacticalDestinationOptions.put("nocross", new ArrayList<Coordinate>());
		
		for(Road r: cPR) {
			// Get candidate destination coordiante from pedestrian road
			Coordinate c = farthestUnobstructedRoadCoordinate(oC, r.getGeom(), obstructGeography);
			
			// Null coordinate returned when it is not possible to see a coordinate on a ped road without obstruction. Skip these
			if (c==null) {
				continue;
			}
			// Check parity - number of times road link is intersected travelling from origin to centroid of polygon.
			// Tells us whether primary crossing is needed to reach polygon
			int parity = RoadNetworkRoute.calculateRouteParity(oC, c, sPS);
			
			if (parity ==1) {
				tacticalDestinationOptions.get("cross").add(c);
			}
			else {
				tacticalDestinationOptions.get("nocross").add(c);
			}
		}
		
		return tacticalDestinationOptions;
	}
	
	/*
	 *Method for finding either the nearest or farthest coordinate of a geometry from a starting point, that doesn't obstruct a geometry in the obstruction
	 *geography
	 * 
	 * @param Coordinate originCoord
	 * @param Geometry rGeom
	 * @param Geography<T> obstructionGeography
	 * @param String nearOrFar
	 */
	public static <T> Coordinate xestUnobstructedRoadCoordinate(Coordinate originCoord, Geometry rGeom, Geography<T> obstructionGeography, String nearOrFar) {
		Coordinate destCoord = null;
		Double destDist = null;
		int compVal = 0;
		if (nearOrFar.contentEquals("near")) {
			destDist = Double.MAX_VALUE;
			compVal = -1;
		}
		else if (nearOrFar.contentEquals("far")){
			destDist = 0.0;
			compVal = +1;
		}
		else {
			return destCoord;
		}
		
		
		// Loop through coordinates of road geometry
		// Find the farthest coordinate not blocked by a pedestrian obstruction
		Coordinate[] rGeomCoords = rGeom.getCoordinates();
		for(Coordinate c: rGeomCoords) {
			
			// Check for obstructions
			Coordinate[] lineCoords = {originCoord, c};
			LineString pathLine = GISFunctions.lineStringGeometryFromCoordinates(lineCoords);

			// Check if line passes through a ped obstruction
			Boolean isObstructingObjects = GISFunctions.doesIntersectGeographyObjects(pathLine, obstructionGeography);
			if(!isObstructingObjects) {
				Double cDist = c.distance(originCoord);
				int comp = Integer.signum(cDist.compareTo(destDist));
				if (comp == compVal) {
					destDist = cDist;
					destCoord = c;
				}
			}
		}
		return destCoord;
	}
	
	/*
	 *Method for finding the farthest coordinate of a geometry from a starting point, that doesn't obstruct a geometry in the obstruction
	 *geography
	 * 
	 * @param Coordinate originCoord
	 * @param Geometry rGeom
	 * @param Geography<T> obstructionGeography
	 */
	public static <T> Coordinate farthestUnobstructedRoadCoordinate(Coordinate originCoord, Geometry rGeom, Geography<T> obstructionGeography) {		
		Coordinate destCoord = xestUnobstructedRoadCoordinate(originCoord, rGeom, obstructionGeography, "far");
		return destCoord;
	}
	
	/*
	 *Method for finding the nearest coordinate of a geometry from a starting point, that doesn't obstruct a geometry in the obstruction
	 *geography
	 * 
	 * @param Coordinate originCoord
	 * @param Geometry rGeom
	 * @param Geography<T> obstructionGeography
	 */
	public static <T> Coordinate nearestUnobstructedRoadCoordinate(Coordinate originCoord, Geometry rGeom, Geography<T> obstructionGeography) {		
		Coordinate destCoord = xestUnobstructedRoadCoordinate(originCoord, rGeom, obstructionGeography, "near");
		return destCoord;
	}
	
	
	public List<RoadLink> getStrategicPath() {
		return this.strategicPath;
	}
	
	public GridRoute getTacticalPath() {
		return this.tacticalPath;
	}

	public Coordinate getNextTacticalPathCoord() {
		return nextTacticalPathCoord;
	}
	
	public Coordinate getNextCrossingCoord() {
		return nextCrossingCoord;
	}
	
	public void setNextCrossingCoord(Coordinate c) {
		this.nextCrossingCoord = c;
	}
	
	public OD getOrigin() {
		return origin;
	}

	public void setOrigin(OD origin) {
		this.origin = origin;
	}
	
	public RoadLink getCurrentRoadLink() {
		return this.currentRoadLink;
	}

}
