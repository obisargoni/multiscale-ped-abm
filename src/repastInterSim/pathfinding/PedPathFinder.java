package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

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
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.TacticalAlternative;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class PedPathFinder {
	
	private Ped ped;
	
	private OD origin;
	private OD destination;
	private Geography<Object> geography;
	
	private List<RoadLink> strategicPath;
	private GridRoute tacticalPath = new GridRoute();
	
	private Coordinate nextTacticalPathCoord = null;

	private Coordinate nextCrossingCoord;	
	private String prevCoordType = "";

	public PedPathFinder(Geography<Object> g, OD o, OD d) {
		this.geography = g;
		this.origin = o;
		this.destination = d;
		
		planStrategicPath();
	}
	
	public PedPathFinder(Ped p) {
		this.ped = p;
		this.geography = p.getGeography();
		this.origin = p.getOrigin();
		this.destination = p.getDestination();
		
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
	public void updateTacticalPathCoordinate() {
				
		// If reached the end of one section of the route, or if route has just been created, need to produce next set of route coordinates.
		if(this.tacticalPath.getRouteX().size() == 0) {
			
			Coordinate tacticalOriginCoord = this.ped.getLoc();
			
			// if previous destination coordinate reached the end of a link, remove the old link from the strategic path
			if (prevCoordType.contentEquals("not_intersects_next") | prevCoordType.contentEquals("intersects_next")) {
				// Remove this road link from the strategic path, no longer the next road link
				if (this.strategicPath.isEmpty() == false) {
					this.strategicPath.remove(0);
				}
			}
			
			// Update pedestrians perceptions of road link
    		String currentRoadLinkID = getCurrentRoadLinkID(this.ped.getLoc());
    		String nextRoadLinkID = getNextRoadLinkID();
    		
        	// Update perception of vehicle priority space cost
        	// triggered by reaching end of tactical path so that new tactical path always uses newly calculated costs of vehicle priority space
    		HashMap<Integer, Double> gSPM = this.ped.calculateDynamicGridSummandPriorityMap(currentRoadLinkID);
			
			// Identify new tactical destination coordinate
			Coordinate tacticalDestCoord = null;
			
			// "not_intersects_next" indicates that the previous tactical destination coord reached the end of one road link but doesn't touch the start of the next road link 
			// in the strategic path. In this case find the nearest coordinate that touches a ped road on the next strategic path link
			if (prevCoordType.contentEquals("not_intersects_next")) {
				tacticalDestCoord = chooseTacticalDestinationCoordinate(tacticalOriginCoord, this.destination.getGeom().getCoordinate(), SpaceBuilder.roadGeography, SpaceBuilder.pedObstructGeography, this.ped.getpHorizon(), this.strategicPath, true);
				prevCoordType = "start_next";
			}
			// Otherwise either haven't reached end of road link or previous dest coord touches next road link. Find farthest coord along road link
			else {
				tacticalDestCoord = chooseTacticalDestinationCoordinate(tacticalOriginCoord, this.destination.getGeom().getCoordinate(), SpaceBuilder.roadGeography, SpaceBuilder.pedObstructGeography, this.ped.getpHorizon(), this.strategicPath, false);
				prevCoordType = checkCoordinateIntersectingRoads(tacticalDestCoord, SpaceBuilder.roadGeography, currentRoadLinkID, nextRoadLinkID);
			}
			
			// Get path to that coordinate
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
	}
	
	
	/*
	 * When agent reaches a tactical destination coordinate they must identify the next tactical destination coordinate. This involves interaction between
	 * the strategic and tactical pathfinding levels as an agents needs to identify whether it has reached the end of one road link.
	 * 
	 * To do this, the tactical destination coordinate needs to be categorised in terms of whether it is at the end of a road link and at the stat of the next road link,
	 * at the end of a road link and not yet at the start of the next road link, or not yet at the end of a road link.
	 * 
	 * Ideally also want to be able to tell if at the start of a road link but not sure how to do this now.
	 * 
	 * @param Coordinate c
	 * 		The coordinate to get information on
	 * @param Geography<Road>
	 * 		The projection containing road objects
	 * @param String cRL
	 * 		The ID of the current road link
	 * @param String nRL
	 * 		The ID of the next road link in the strategic path.
	 *  
	 * @return String
	 * 		Indicates the location of coordinate in reference to the strategic path
	 */
	public static String checkCoordinateIntersectingRoads(Coordinate c, Geography<Road> roadGeography, String cRL, String nRL) {
				
		// Get the intersecting road objects
		List<Road> roads = SpatialIndexManager.findIntersectingObjects(roadGeography, c);
		
		// If all roads have the ID of the current road link then coordinate is not at the end of the road link
		Boolean notAtEnd = roads.stream().allMatch(r -> r.getRoadLinkID().contentEquals(cRL));
		if(notAtEnd) {
			return "not_at_end";
		}
		
		// Check if intersecting roads include the next link in the strategic path
		Boolean intersectsNext = roads.stream().anyMatch(r -> r.getRoadLinkID().contentEquals(nRL));
		
		if (intersectsNext) {
			return "intersects_next";
		}
		else {
			return "not_intersects_next";
		}
	}
	
	/*
	 * Identify the destination coordinate for the tactical path. This is currently defined as the coordinate of the end junction of the
	 * current road link in the road link path, if the agent has reached the final road link in the path, the final destination.
	 * 
	 * @param Coordinate originCoord
	 * 		The starting coordinate
	 * @param String prevDestType
	 * 		A string indicating the type of destination coordinate the last tactical destination was. This is used to determine how next
	 * destination options are identified.
	 */
	public static Coordinate chooseTacticalDestinationCoordinate(Coordinate originCoord, Coordinate destCoord, Geography<Road> roadGeography, Geography<PedObstruction> pedObstructGeography, int pH, List<RoadLink> sP, Boolean secondaryCrossing) {
		
		Coordinate tacticalDestCoord = null;
		
		// If reached end of strategic path tactical destination is final route destination
		if (sP.isEmpty() | sP.size() == 1) {
			tacticalDestCoord = destCoord;
		}
		else {			
			// Get the ID of the current road link
			String roadLinkID = sP.get(0).getFID();
			
			// Get pedestrian roads linked to this road link
			List<Road> pedRoads = null;
			try {
				pedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(roadLinkID, roadGeography);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			// If not at a secondary crossing, destination coordinate options are the farthest coordinates on the ped roads associated to the road link
			if (secondaryCrossing==false) {
				// Get tactical destination options
				ArrayList<TacticalAlternative> alternatives = getTacticalDestinationAlternatives(originCoord, pedRoads, sP, destCoord, pH, pedObstructGeography, true);
				
				// Write a comparator that sorts on multiple attributes of alternative
				Comparator<TacticalAlternative> comparator = Comparator.comparing(TacticalAlternative::getParityS,Comparator.nullsLast(Comparator.naturalOrder()))
																		.thenComparing(TacticalAlternative::getParityT, Comparator.nullsLast(Comparator.naturalOrder()))
																		.thenComparing(TacticalAlternative::getCostT, Comparator.reverseOrder());
				
				// Now sort the tactical alternatives and choose appropriate one
				List<TacticalAlternative> sortedAlternatives = alternatives.stream().sorted(comparator).collect(Collectors.toList());
				tacticalDestCoord = sortedAlternatives.get(0).getC();
			}
			
			// If about to perform a secondary crossing destination options are nearest coordinates. Always select the no crossing option
			else if (secondaryCrossing==true) {
				// Get tactical destination options
				ArrayList<TacticalAlternative> alternatives = getTacticalDestinationAlternatives(originCoord, pedRoads, sP, destCoord, pH, pedObstructGeography, false);
				
				// Always select only the nearest no cross option
				tacticalDestCoord = alternatives.stream().filter(ta -> ta.getParityT() == 0).sorted((ta1,ta2) -> ta1.getCostT().compareTo(ta2.getCostT())).collect(Collectors.toList()).get(0).getC();
			}
		}
		
		return tacticalDestCoord;
	}
	
	
	/*
	 * Given a starting coordinate, ie the current position of the agent, the pedestrian roads the agent is considering as destinations, the road link associated to these
	 * pedestrian roads and the geography projection of pedestrian obstructions.
	 * 
	 * @param Coordinate oC
	 * 		The origin coordinate
	 * @param List<Road> pR
	 * 		The pedestrian roads
	 * @param List<RoadLink> sPS
	 * 		A segment of the strategic path that consists of the link associated to the pedestrian roads
	 * @param Geography<PedObstruction> obstructGeography
	 * 		The geography projection containing obstructions to pedestrian movement
	 * 
	 * @return HashMap<String, List<Coordinate>>
	 */
	public static ArrayList<TacticalAlternative> getTacticalDestinationAlternatives(Coordinate oC, List<Road> pR, List<RoadLink> sP, Coordinate d, int pH, Geography<PedObstruction> obstructGeography, Boolean far){
		
		ArrayList<TacticalAlternative> alternatives = new ArrayList<TacticalAlternative>();
		
		for(Road r: pR) {
			// Get candidate destination coordiante from pedestrian road
			Coordinate c = xestUnobstructedRoadCoordinate(oC, r.getGeom(), obstructGeography, far);
			
			// Null coordinate returned when it is not possible to see a coordinate on a ped road without obstruction. Skip these
			if (c==null) {
				continue;
			}
			
			// Create new alternative
			TacticalAlternative tc = new TacticalAlternative(c);
			
			// Find parity of coordinate at tactical scale
			// Tells us whether primary crossing is needed to reach coordinate
			tc.setParityT(RoadNetworkRoute.calculateRouteParity(oC, c, sP.subList(0, 1)));
			
			// Find distance to tactical coordinate
			tc.setCostT(c.distance(oC));
			
			// Find parity at strategic level - identifies whether primary crossing required to get from tactical alternative to destination
			// Takes into account bounded rationality
			tc.setParityS(RoadNetworkRoute.boundedRouteParity(c, d, sP, pH));
			
			// Add to collection of choice alternatives
			alternatives.add(tc);
		}
		return alternatives;
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
	public static <T> Coordinate xestUnobstructedRoadCoordinate(Coordinate originCoord, Geometry rGeom, Geography<T> obstructionGeography, Boolean far) {
		Coordinate destCoord = null;
		Double destDist = null;
		int compVal = 0;
		if (far == false) {
			destDist = Double.MAX_VALUE;
			compVal = -1;
		}
		else if (far == true){
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
		Coordinate destCoord = xestUnobstructedRoadCoordinate(originCoord, rGeom, obstructionGeography, true);
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
		Coordinate destCoord = xestUnobstructedRoadCoordinate(originCoord, rGeom, obstructionGeography, false);
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
	
	public String getCurrentRoadLinkID(Coordinate c) {
		
		String currentRoadLinkID = null;
		
		// For pedestrian ODs that are close to one another on the same road link the currentRoadLink is null
		// In this case pedestrian identifies the road link to estimate vehicle priority costs by using their current location
		if (this.strategicPath.isEmpty()) {
			Road currentRoad;
			try {
				currentRoad = GISFunctions.getCoordinateRoad(c);
	    		currentRoadLinkID = currentRoad.getRoadLinkID();
			} catch (RoutingException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		else {
			currentRoadLinkID = this.strategicPath.get(0).getFID();
		}
		
		return currentRoadLinkID;
	}
	
	private String getNextRoadLinkID() {
		if(this.strategicPath.isEmpty() | this.strategicPath.size()<2) {
			return "";
		}
		else {
			return this.strategicPath.get(1).getFID();
		}
	}

}
