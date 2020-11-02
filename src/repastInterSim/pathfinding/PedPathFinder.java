package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Stack;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.geotools.coverage.grid.GridCoverage2D;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.LineString;

import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;
import repastInterSim.environment.OD;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class PedPathFinder {
	
	private Ped ped;
	
	private OD origin;
	private OD destination;
	private Geography<Object> geography;
	
	private Geography<PedObstruction> obstructGeography;
	private Geography<Road> rGeography;
	
	
	private List<RoadLink> strategicPath;
	private AccumulatorRoute tacticalPath = new AccumulatorRoute();
	
	private Coordinate nextCrossingCoord;	
	private String prevCoordType = "";

	public PedPathFinder(Geography<Object> g, OD o, OD d) {
		init(g, o, d);
	}
	
	public PedPathFinder(Ped p) {
		this.ped = p;
		init(p.getGeography(), p.getOrigin(), p.getDestination());
	}
	
	private void init(Geography<Object> g, OD o, OD d) {
		this.geography = g;
		this.origin = o;
		this.destination = d;
		
		this.obstructGeography = SpaceBuilder.pedObstructGeography;
		this.rGeography = SpaceBuilder.roadGeography;
		
		planStrategicPath(this.origin, this.destination, SpaceBuilder.orRoadLinkGeography, SpaceBuilder.orRoadNetwork, SpaceBuilder.pedestrianDestinationGeography);
	}
	
	/**
	 * Initialise a new road link routing object that can be used to find a path along a topological road network.
	 * Use this to identify the shortest path through the network and assign this path to this classes' strategic path attribute.
	 * 
	 */
	public void planStrategicPath(Coordinate oC, Coordinate dC, Geography<RoadLink> rlG, Network<Junction> pedNet, Geography<OD> odG) {
		// Initialise road network route - needs to ne non-directed for pedestrians! fix this
		RoadNetworkRoute rnr = new RoadNetworkRoute(oC, dC, rlG, pedNet, odG);
		
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
	 * Set the start and destination pavement network junctions from the road network route
	 */
	public static Junction[] strategicPathPavementJunctions(Network<Junction> pavementNetwork, RoadNetworkRoute rnr, Coordinate oC, Coordinate dC) {
		
		Junction[] rnrEndpoints = rnr.routeEndpoints;
		
		// Choose the candidate associated to the first road link in the strategic path and that is closest to the origin
		Junction originPavementJunction = null;
		double d = Double.MAX_VALUE;
		for (Junction j: pavementNetwork.getNodes()) {
			if (j.getjuncNodeID().contentEquals(rnrEndpoints[0].getjuncNodeID())) {
				RoadLink startLink = rnr.getRoadsX().get(0);
				if (j.getv1rlID().contentEquals(startLink.getFID()) | j.getv2rlID().contentEquals(startLink.getFID())) {
					Double dj = oC.distance(j.getGeom().getCoordinate());
					if (dj < d) {
						d = dj;
						originPavementJunction = j;
					}
				}
			}
		}
		
		Junction destPavementJunction = null;
		d = Double.MAX_VALUE;
		for (Junction j: pavementNetwork.getNodes()) {
			if (j.getjuncNodeID().contentEquals(rnrEndpoints[1].getjuncNodeID())) {
				RoadLink endLink = rnr.getRoadsX().get(rnr.getRoadsX().size()-1);
				if (j.getv1rlID().contentEquals(endLink.getFID()) | j.getv2rlID().contentEquals(endLink.getFID())) {
					Double dj = oC.distance(j.getGeom().getCoordinate());
					if (dj < d) {
						d = dj;
						destPavementJunction = j;
					}
				}
			}
		}
		
		Junction[] odPavementJunctions = {originPavementJunction, destPavementJunction}; 
		return odPavementJunctions;
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
		
		// Calculate number of links in planning horizon
		int nLinks = getNLinksWithinAngularDistance(this.strategicPath, this.ped.getpHorizon());
		
		// "not_intersects_next" indicates that the previous tactical destination coord reached the end of one road link but doesn't touch the start of the next road link 
		// in the strategic path. In this case find the nearest coordinate that touches a ped road on the next strategic path link
		if (prevCoordType.contentEquals("not_intersects_next")) {
			tacticalDestCoord = chooseTacticalDestinationCoordinate(tacticalOriginCoord, this.destination.getGeom().getCoordinate(), this.rGeography, this.obstructGeography, nLinks, this.strategicPath, true);
			prevCoordType = "start_next";
		}
		// Otherwise either haven't reached end of road link or previous dest coord touches next road link. Find farthest coord along road link
		else {
			tacticalDestCoord = chooseTacticalDestinationCoordinate(tacticalOriginCoord, this.destination.getGeom().getCoordinate(), this.rGeography, this.obstructGeography, nLinks, this.strategicPath, false);
			prevCoordType = checkCoordinateIntersectingRoads(tacticalDestCoord, this.rGeography, currentRoadLinkID, nextRoadLinkID);
		}
		
		Coordinate defaultDest = defaultDestinationCoordinate(SpaceBuilder.pedJunctionGeography, this.strategicPath.subList(0, nLinks), this.ped.getLoc());
		
		// Set the pedestrian's default destination
		this.ped.setDefaultDestination(defaultDest);
		
		// Initialise Accumulator Route that agent will use to navigate along the planning horizon
		planTacticaAccumulatorPath(SpaceBuilder.caGeography, this.strategicPath.subList(0, nLinks), this.ped, SpaceBuilder.roadGeography, tacticalDestCoord, defaultDest);
    }
	
	/*
	 * Returns the ID of the road network node connecting two road links. Returns null if there isn't a shared node.
	 */
	public static String connectingNodeID(RoadLink rl1, RoadLink rl2) {
		
		// Find the common node
		if (rl1.getMNodeFID().contentEquals(rl2.getMNodeFID()) | rl1.getMNodeFID().contentEquals(rl2.getPNodeFID())) {
			return rl1.getMNodeFID();
		}
		else if (rl1.getPNodeFID().contentEquals(rl2.getPNodeFID()) | rl1.getPNodeFID().contentEquals(rl2.getMNodeFID())) {
			return rl1.getPNodeFID();
		}
		else {
			return null;
		}
	}
	
	/*
	 * Given a road node ID return the pavement network junctions associated to this node. These are the pavement junctions
	 * around a road network node.
	 */
	public static List<Junction> roadNodePavementJunctions(Network<Junction> pavementNetwork, String roadNodeID) {
		 List<Junction> nodeJunctions = new ArrayList<Junction>();
		 for(Junction j:pavementNetwork.getNodes()) {
			 if (j.getjuncNodeID().contentEquals(roadNodeID)) {
				 nodeJunctions.add(j);
			 }
		 }
		 return nodeJunctions;
	}
	
	/*
	 * Get the pavement network junctions that are at the end of the tactical planning horizon.
	 * 
	 */
	public static List<Junction> tacticalHorizonEndJunctions(Network<Junction> pavementNetwork, RoadLink rlEndHorz, RoadLink rlOutHorz) {
		
		HashMap<String, List<Junction>> tacticalEndJunctions = tacticalHorizonJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		return tacticalEndJunctions.get("end");
	}
	
	/*
	 * Get the pavement network junctions that are at the start of the first link beyond the tactical planning horizon.
	 * 
	 */
	public static List<Junction> tacticalHorizonOutsideJunctions(Network<Junction> pavementNetwork, RoadLink rlEndHorz, RoadLink rlOutHorz) {
		
		HashMap<String, List<Junction>> tacticalEndJunctions = tacticalHorizonJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		return tacticalEndJunctions.get("outside");
	}
	
	/*
	 * Get the pavement network junctions that at the intersection between the final road link in the tactical planning horizin and
	 * the next link in the strategic path that lies outside the tactical planning horizon.
	 * 
	 * Choose whether to return the junctions connected to the link at the end of the horizon or connected to the next link beyong the horizon.
	 * 
	 * @param Network<Junction> pavementNetwork
	 * 			The network of pavement junctions used to connect pedestrian pavements
	 * @param RoadLink rlEndHorz
	 * 			The road link at the end of the tactical planning horizon
	 * @param RoadLink rlOutHorz
	 * 			The road link on the other side of the tactical planning horizon
	 * @param boolean outside
	 * 			Selects whether to return the junctions at the end of the planning horizon or outside the planning horizon at the start of the next link
	 */
	public static HashMap<String, List<Junction>> tacticalHorizonJunctions(Network<Junction> pavementNetwork, RoadLink rlEndHorz, RoadLink rlOutHorz) {
		
		// First get the road network node id connecting these links
		String nodeID = connectingNodeID(rlEndHorz, rlOutHorz);
		
		// Get the pavement junctions linked to this road node
		List<Junction> pavementJunctions =  roadNodePavementJunctions(pavementNetwork, nodeID);

		// Loop over pavement junctions and select those touching the road link at the end of the planning horizon
		HashMap<String, List<Junction>> tacticalEndJunctions = new HashMap<String, List<Junction>>();
		tacticalEndJunctions.put("end", new ArrayList<Junction>());
		tacticalEndJunctions.put("outside", new ArrayList<Junction>());
		for (Junction j: pavementJunctions) {
			if (j.getv1rlID().contentEquals(rlEndHorz.getPedRLID()) | j.getv2rlID().contentEquals(rlEndHorz.getPedRLID())) {
				tacticalEndJunctions.get("end").add(j);
			}
			if (j.getv1rlID().contentEquals(rlOutHorz.getPedRLID()) | j.getv2rlID().contentEquals(rlOutHorz.getPedRLID())) {
				tacticalEndJunctions.get("outside").add(j);
			}
		}
		
		return tacticalEndJunctions;
	}
	
	public static TacticalRoute setupTacticalRoute(NetworkPath<Junction> nP, List<RoadLink> sP, Junction eJ, List<Junction> outsideJunctions, Junction currentJ, Junction destJ) {
		
		TacticalRoute tr = new TacticalRoute();
		
		// Get path to end junctions
		List<RepastEdge<Junction>> path = nP.getShortestPath(currentJ, eJ);
		
		// Get path from end junction to the junction at the start of the first link outside the tactical planning horizon
		List<RepastEdge<Junction>> pathToOutside = new ArrayList<RepastEdge<Junction>>();
		Stack<Junction> nodesToOutside = new Stack<Junction>();
		double minLength = Double.MAX_VALUE;
				
		// Calculate path from the end junction to each possible outside junction. Select the shortest path that does not involve a primary crossing
		// This is the path that stays on the same side of the road and moves to the start of the next route section
		for (Junction j: outsideJunctions) {
			
			Predicate<Junction> nodeFilter = n -> n.getjuncNodeID().contentEquals(eJ.getjuncNodeID());
			List<Stack<Junction>> nodePaths = nP.getSimplePaths(eJ, j, nodeFilter);
			
			// Loop through these paths and find the shortest one that does not make a primary crossing
			for (Stack<Junction> nodePath : nodePaths) {
				List<RepastEdge<Junction>> edgePath = nP.edgePathFromNodes(nodePath);
				if (containsPrimaryCrossing(edgePath, sP) == false) {
					double pathLength = nP.getPathLength(edgePath);
					if (pathLength < minLength) {
						nodesToOutside = nodePath;
						pathToOutside = edgePath;
						minLength = pathLength;
					}
				}
			}
		}
		
		// Get junctions to go via in orget to get from end junction to outside junction, including the end junction itself
		for (Junction j : nodesToOutside) {
			tr.addJunction(j);
		}
		
		// Combine the two paths and add to the tactical route
		for (RepastEdge<Junction> e: pathToOutside) {
			path.add(e);
		}
		
		tr.setRoutePath(path);
		
		// Finally, if the destination junction is known, calculate the path from the last junction added to the tactical route to the destination junction
		// This is recorded separately as the path required to complete the journey
		if (destJ != null) {
			tr.setRouteRemainderPath(nP.getShortestPath(nodesToOutside.get(nodesToOutside.size()-1), destJ));
		}
		
		return tr;
	}
	private static boolean containsPrimaryCrossing(List<RepastEdge<Junction>> path, List<RoadLink> sP) {
		boolean cPC = false;
		for (RepastEdge<Junction> e: path) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) e;
			for (RoadLink rl : sP) {
				if (rl.getPedRLID().contentEquals(ne.getRoadLink().getPedRLID())) {
					cPC = true;
				}
			}
		}		
		return cPC;
	}
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
		
		//this.tacticalPath= tP;
	}
	
	/*
	 * Plan a tactical level path using the accumulator crossing choice path finding model.
	 */
	public void planTacticaAccumulatorPath(Geography<CrossingAlternative> caG, List<RoadLink> sP, Ped p, Geography<Road> rG, Coordinate tacticalDestCoord, Coordinate defaultDestCoord) {
		
		// Get crossing alternatives within planning horizon
		List<CrossingAlternative> cas = getCrossingAlternatives(caG, sP, p, rG, tacticalDestCoord);
		
		// Get length of the planning horizon, this is used in the accumulator route
		double pHLength = 0;
		for (RoadLink rl: sP) {
			pHLength += rl.getGeom().getLength();
		}
		
		// Initialse AccumulatorRoute with strategic path planning horizon and origin and destination coordiantes, crossing alternatives
		AccumulatorRoute accRoute = new AccumulatorRoute(p, defaultDestCoord, cas, pHLength);
		
		this.tacticalPath = accRoute;
	}
	
	public void step() {
		
		// Only update the accumulator path if the ped agent is walking along a road linked to the strategic path,
		// ie not crossing over another road
		if (pedOnStrategicPathRoadLink()) {
			this.tacticalPath.step();
		}
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
		
		// If number of links in planning horizon equal to or greater than length of strategic path then destination is within planning horizon.
		// Set the destination as the tactical destination coordinate
		if (pH >= sP.size()) {
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
				ArrayList<TacticalRoute> alternatives = getTacticalDestinationAlternatives(originCoord, pedRoads, sP, destCoord, pH, pedObstructGeography, true);
				
				// Write a comparator that sorts on multiple attributes of alternative
				Comparator<TacticalRoute> comparator = Comparator.comparing(TacticalRoute::getParityS,Comparator.nullsLast(Comparator.naturalOrder()))
																		.thenComparing(TacticalRoute::getParityT, Comparator.nullsLast(Comparator.naturalOrder()))
																		.thenComparing(TacticalRoute::getCostT, Comparator.reverseOrder());
				
				// Now sort the tactical alternatives and choose appropriate one
				List<TacticalRoute> sortedAlternatives = alternatives.stream().sorted(comparator).collect(Collectors.toList());
				tacticalDestCoord = sortedAlternatives.get(0).getC();
			}
			
			// If about to perform a secondary crossing destination options are nearest coordinates. Always select the no crossing option
			else if (secondaryCrossing==true) {
				// Get tactical destination options
				ArrayList<TacticalRoute> alternatives = getTacticalDestinationAlternatives(originCoord, pedRoads, sP, destCoord, pH, pedObstructGeography, false);
				
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
	public static ArrayList<TacticalRoute> getTacticalDestinationAlternatives(Coordinate oC, List<Road> pR, List<RoadLink> sP, Coordinate d, int pH, Geography<PedObstruction> obstructGeography, Boolean far){
		
		ArrayList<TacticalRoute> alternatives = new ArrayList<TacticalRoute>();
		
		for(Road r: pR) {
			// Get candidate destination coordiante from pedestrian road
			Coordinate c = GISFunctions.xestUnobstructedGeomCoordinate(oC, r.getGeom(), obstructGeography, far);
			
			// Null coordinate returned when it is not possible to see a coordinate on a ped road without obstruction. Skip these
			if (c==null) {
				continue;
			}
			
			// Create new alternative
			TacticalRoute tc = new TacticalRoute(c);
			
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
	 * Identify the crossing alternatives that lie on the given road links. Prepare these crossing alternatives
	 * for use in the accumulator choice model.
	 * 
	 * @param Geography<CrossingAlternative> caG
	 * 		The geography containing all crossing alternative objects
	 * @param List<RoadLink> sP
	 * 		The road links to identify crossing alternatives along
	 * @param Ped p
	 * 		The pedestrian agent perceiving the crossing alternatives
	 * @param Geography<Road> rG
	 * 		The Geography containing Road objects
	 */
	public static List<CrossingAlternative> getCrossingAlternatives(Geography<CrossingAlternative> caG, List<RoadLink> sP, Ped p, Geography<Road> rG, Coordinate dest){
		
		List<CrossingAlternative> cas = new ArrayList<CrossingAlternative>();
		
		// Loop through open road sp section road links and crossing alternatives and add any crossing alternatives that match the road link id to the list
		for (RoadLink rl: sP) {
			for (CrossingAlternative ca: caG.getAllObjects()) {
				if (ca.getRoadLinkID().contentEquals(rl.getFID())) {
					
					// Set up this crossing alternative
					ca.setPed(p);
					ca.setRoadGeography(rG);
					ca.setDestination(dest);
					ca.setStrategicPathSection(sP);
					ca.setRoad();
					
					// Add to list
					cas.add(ca);
				}
			}
		}
		
		// Add unmarked crossing alternative to lsit
		// Don't set road for unmarked crossings
		CrossingAlternative caU = new CrossingAlternative();
		caU.setType("unmarked");
		caU.setPed(p);
		caU.setRoadGeography(rG);
		caU.setDestination(dest);
		caU.setStrategicPathSection(sP);
		cas.add(caU);
		
		return cas;		
	}
	
	public static int getNLinksWithinAngularDistance(List<RoadLink> sP, double angDistThreshold) {
		// Initialise number of links in planning horizon as one (current link always in planning horizon)
		int nLinks = 1;
		
		// Initialise ang dist
		double angDist = 0.0;
		
		// Loop through road links in strategic path adding up the angular distance between each one
		for (int i = 0; i < sP.size() - 1; i++) {
			RoadLink rl1 = sP.get(i);
			RoadLink rl2 = sP.get(i+1);
			
			// Calculate angular distance to next road link
			angDist += GISFunctions.angleBetweenConnectedLineStrings((LineString) rl1.getGeom(), (LineString) rl2.getGeom());
			
			// If distance within threshold, add next link to horizon
			if (angDist > angDistThreshold) {
				break;
			}
			
			// If angular distance threshold not exceeded add 1 to number of links within horizon
			nLinks++;

		}
		return nLinks;
	}
	
	public static List<RoadLink> getLinksWithinAngularDistance(List<RoadLink> sP, double angDistThreshold){		
		// Calculate number of links that are within planning horizon
		int nLinks = getNLinksWithinAngularDistance(sP, angDistThreshold);
		
		// use this to get get planning horizon as list of links
		return sP.subList(0, nLinks);
	}
	
	/*
	 * Used to check whether the pedestrian agents is currently walking on a road
	 * polygon alongside a link in the strategic path or not.
	 * 
	 * When a pedestrian agent is making a secondary crossing it is not on a road beside
	 * a strategic path link.
	 */
	public Boolean pedOnStrategicPathRoadLink() {
		Boolean pedOnSLink = false;
		Road r = this.ped.getCurrentRoad();
		for (RoadLink rl: this.strategicPath) {
			if (rl.getFID().contentEquals(r.getRoadLinkID())) {
				pedOnSLink = true;
			}
		}
		return pedOnSLink;
	}
	
	
	public List<RoadLink> getStrategicPath() {
		return this.strategicPath;
	}
	
	public AccumulatorRoute getTacticalPath() {
		return this.tacticalPath;
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
