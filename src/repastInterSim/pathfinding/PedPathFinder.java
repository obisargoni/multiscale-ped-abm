package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Collections;
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
	
	private List<RoadLink> strategicPath;
	private Junction[] spPavementJunctionEndpoints;
	private AccumulatorRoute tacticalPath = new AccumulatorRoute();
	
	private Coordinate nextCrossingCoord;	

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
		
		planStrategicPath(this.origin.getGeom().getCoordinate(), this.destination.getGeom().getCoordinate(), SpaceBuilder.orRoadLinkGeography, SpaceBuilder.orRoadNetwork, SpaceBuilder.pedestrianDestinationGeography, SpaceBuilder.pavementNetwork);
	}
	
	public void step() {
		
		// Only update the accumulator path if the ped agent is walking along a road linked to the strategic path,
		// ie not crossing over another road
		if (pedOnStrategicPathRoadLink()) {
			this.tacticalPath.step();
		}
	}
	
	/**
	 * Initialise a new road link routing object that can be used to find a path along a topological road network.
	 * Use this to identify the shortest path through the network and assign this path to this classes' strategic path attribute.
	 * 
	 */
	public void planStrategicPath(Coordinate oC, Coordinate dC, Geography<RoadLink> rlG, Network<Junction> orNetwork, Geography<OD> odG, Network<Junction> paveNetwork) {
		// Initialise road network route - needs to ne non-directed for pedestrians! fix this
		RoadNetworkRoute rnr = new RoadNetworkRoute(oC, dC, rlG, orNetwork, odG);
		
		// Find shortest path using road network route
		try {
			rnr.setRoadLinkRoute();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get path of road links and set this as the strategic path
		this.strategicPath = rnr.getRoadsX();
		
		this.spPavementJunctionEndpoints = strategicPathPavementJunctions(paveNetwork, rnr, oC, dC);
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
		// Initialise Accumulator Route that agent will use to navigate along the planning horizon
		planTacticaAccumulatorPath(SpaceBuilder.pavementNetwork, SpaceBuilder.caGeography, SpaceBuilder.roadGeography, this.ped, this.strategicPath, this.spPavementJunctionEndpoints[0], this.spPavementJunctionEndpoints[1]);
    }
	
	/*
	 * Plan a tactical level path using the accumulator crossing choice path finding model.
	 */
	public void planTacticaAccumulatorPath(Network<Junction> pavementNetwork, Geography<CrossingAlternative> caG, Geography<Road> rG, Ped p, List<RoadLink> sP, Junction currentJ, Junction destJ) {
		
		// Calculate number of links in planning horizon
		int nLinks = getNLinksWithinAngularDistance(sP, p.getpHorizon());
		List<RoadLink> tacticalHorizon = this.strategicPath.subList(0, nLinks);
		
		// First identify tactical route alternatives
		List<TacticalRoute> trs = tacticalRoutes(pavementNetwork, sP, p.getpHorizon(), currentJ, destJ);
		
		// Sort routes based on the length of the path to the end of tactical horizon
		List<String> strategiRoadLinkIDs = sP.stream().map(rl->rl.getPedRLID()).collect(Collectors.toList());
		PavementRoadLinkTransformer<Junction> transformer = new PavementRoadLinkTransformer<Junction>(strategiRoadLinkIDs, Double.MAX_VALUE);
		List<Double> pathLengths = trs.stream().map(tr -> NetworkPath.getPathLength(tr.getRoutePath(), transformer)).collect(Collectors.toList());
		
		// Choose tactical route alternative
		// Default to choosing alternative with fewest primary crossings required to complete tactical horizon
		// By definition this is also the default TacticalRoute the agent walks towards
		TacticalRoute chosenTR = trs.get(pathLengths.indexOf(Collections.min(pathLengths)));
		
		// If chosen to cross road, identify crossing options and initialise accumulator route
		List<CrossingAlternative> cas = new ArrayList<CrossingAlternative>();
		if (NetworkPath.getPathLength(chosenTR.getRoutePath(), transformer)<Double.MAX_VALUE) {
			// No primary crossings required
		}
		else {
			// Get crossing alternatives within planning horizon
			cas = getCrossingAlternatives(caG, tacticalHorizon, p, rG, chosenTR.getRouteJunctions().get(0).getGeom().getCoordinate());
		}
		
		// Get length of the planning horizon, this is used in the accumulator route
		double pHLength = 0;
		for (RoadLink rl: sP) {
			pHLength += rl.getGeom().getLength();
		}
		
		//Initialise Accumulator route with the chosen tactical route also set as the default.
		AccumulatorRoute accRoute = new AccumulatorRoute(p, cas, pHLength, chosenTR, chosenTR);
		
		this.tacticalPath = accRoute;
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
	
	public static List<TacticalRoute> tacticalRoutes(Network<Junction> pavementNetwork, List<RoadLink> sP, Double pH, Junction currentJ, Junction destJ) {
		// Get road link ID of link at end of planning horizon and first strategic path road link outside of planning horizon
		List<RoadLink> tacticalPlanHorz = PedPathFinder.getLinksWithinAngularDistance(sP, pH);
		RoadLink rlEndHorz = tacticalPlanHorz.get(tacticalPlanHorz.size()-1);
		RoadLink rlOutHorz = sP.get(tacticalPlanHorz.size());
		
		// Get horizon junctions
		HashMap<String, List<Junction>> horizonJunctions = tacticalHorizonJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
		
		// get path to each of the junctions at the end of the planning horizon
		List<Junction> endJunctions = horizonJunctions.get("end");
		
		NetworkPath<Junction> nP = new NetworkPath<Junction>(pavementNetwork);
		
		// For each of the end junctions create a TacticalRoute object representing the route via this junction
		List<TacticalRoute> tacticalRoutes = new ArrayList<TacticalRoute>();
		for (Junction eJ: endJunctions) {
			TacticalRoute tr = setupTacticalRoute(nP, sP, eJ, horizonJunctions.get("outside"), currentJ, destJ);
			tacticalRoutes.add(tr);
		}
		return tacticalRoutes;
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
	
	/*
	 * Identify the crossing alternatives that lie on the given road links. Prepare these crossing alternatives
	 * for use in the accumulator choice model.
	 * 
	 * @param Geography<CrossingAlternative> caG
	 * 		The geography containing all crossing alternative objects
	 * @param List<RoadLink> rls
	 * 		The road links to identify crossing alternatives along, typically the tactical planning horizon section of the strategic path.
	 * @param Ped p
	 * 		The pedestrian agent perceiving the crossing alternatives
	 * @param Geography<Road> rG
	 * 		The Geography containing Road objects
	 */
	public static List<CrossingAlternative> getCrossingAlternatives(Geography<CrossingAlternative> caG, List<RoadLink> rls, Ped p, Geography<Road> rG, Coordinate dest){
		
		List<CrossingAlternative> cas = new ArrayList<CrossingAlternative>();
		
		// Loop through open road sp section road links and crossing alternatives and add any crossing alternatives that match the road link id to the list
		for (RoadLink rl: rls) {
			for (CrossingAlternative ca: caG.getAllObjects()) {
				if (ca.getRoadLinkID().contentEquals(rl.getFID())) {
					
					// Set up this crossing alternative
					ca.setPed(p);
					ca.setRoadGeography(rG);
					ca.setDestination(dest);
					ca.setStrategicPathSection(rls);
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
		caU.setStrategicPathSection(rls);
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
}
