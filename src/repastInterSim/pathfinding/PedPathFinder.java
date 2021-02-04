package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Stack;
import java.util.stream.Collectors;

import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.Transformer;

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
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.UnmarkedCrossingAlternative;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.SpaceBuilder;

public class PedPathFinder {
	
	private Ped ped;
	
	private OD origin;
	private OD destination;
		
	private List<RoadLink> strategicPath;
	int tacticalHorizonLinks = 0;
	private Junction startPavementJunction;
	private Junction destPavementJunction;
	private TacticalAlternative tacticalPath = new TacticalAlternative();
	
	private Coordinate nextCrossingCoord;	

	public PedPathFinder(OD o, OD d, Geography<RoadLink> rlG, Network<Junction> orNetwork, Geography<OD> odG, Geography<Junction> paveG, Network<Junction> paveNetwork) {
		init(o, d, rlG, orNetwork, odG, paveG, paveNetwork);
	}
	
	public PedPathFinder(Ped p, Geography<RoadLink> rlG, Network<Junction> orNetwork, Geography<OD> odG, Geography<Junction> paveG, Network<Junction> paveNetwork) {
		this.ped = p;
		init(p.getOrigin(), p.getDestination(), rlG, orNetwork, odG, paveG, paveNetwork);
	}
	
	private void init(OD o, OD d, Geography<RoadLink> rlG, Network<Junction> orNetwork, Geography<OD> odG, Geography<Junction> paveG, Network<Junction> paveNetwork) {
		this.origin = o;
		this.destination = d;
				
		planStrategicPath(this.origin.getGeom().getCoordinate(), this.destination.getGeom().getCoordinate(), rlG, orNetwork, odG, paveG, paveNetwork);
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
	public void planStrategicPath(Coordinate oC, Coordinate dC, Geography<RoadLink> rlG, Network<Junction> orNetwork, Geography<OD> odG, Geography<Junction> paveG, Network<Junction> paveNetwork) {
		RoadNetworkRoute rnr = new RoadNetworkRoute(oC, dC, rlG, orNetwork, odG);
		
		// initialise object to record the start and end pavement junctions of the route
		Junction[] routeEnds = null;
		
		// Find shortest path using road network route
		try {
			routeEnds = rnr.setRoadLinkRoute(paveG, paveNetwork);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get path of road links and set this as the strategic path
		this.strategicPath = rnr.getRoadsX();
		
		
		this.startPavementJunction = routeEnds[0];
		this.destPavementJunction = routeEnds[1];
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
	public void updateTacticalPath() {
		
		// First update the strategic path by removing the links that formed part of the previous tactical planning horizon
		for (int i = 0; i < tacticalHorizonLinks; i++) {
			this.strategicPath.remove(0);
		}
		
		// If no tactical path has been set use the strategic path start junction, otherwise set the start junction as the end junction of previous tactical path
		Junction startJunction = null;
		if (this.tacticalPath.isBlank()) {
			startJunction = this.startPavementJunction;
		}
		else {
			startJunction = this.tacticalPath.getFinalJunction();
		}
		
		// Initialise Accumulator Route that agent will use to navigate along the planning horizon, and update the number of links in the tactical planning horizon
		tacticalHorizonLinks = planTacticalPath(SpaceBuilder.pavementNetwork, SpaceBuilder.caGeography, SpaceBuilder.roadGeography, this.ped, this.strategicPath, startJunction, this.destPavementJunction);
    }
	
	/*
	 * Plan a tactical level path using the accumulator crossing choice path finding model.
	 */
	public int planTacticalPath(Network<Junction> pavementNetwork, Geography<CrossingAlternative> caG, Geography<Road> rG, Ped p, List<RoadLink> sP, Junction currentJ, Junction destJ) {
		
		// Calculate number of links in planning horizon
		int tacticalNLinks = getNLinksWithinAngularDistance(sP, p.getpHorizon());
		
		boolean destInPlanningHorizon = false;
		if (tacticalNLinks == sP.size()) {
			destInPlanningHorizon = true;
		}
		
		// Get road link ID of link at end of planning horizon and first strategic path road link outside of planning horizon
		List<Junction> endJunctions = null;
		if (destInPlanningHorizon==false) {
			RoadLink rlEndHorz = sP.get(tacticalNLinks-1);
			RoadLink rlOutHorz = sP.get(tacticalNLinks);
			HashMap<String, List<Junction>> horizonJunctions = tacticalHorizonJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
			endJunctions = horizonJunctions.get("end");
		}
		else {
			endJunctions = new ArrayList<Junction>();
			endJunctions.add(destJ);
		}
		// Choose path to end of tactical horizon
		// Need to get two junctions at the end of the first link in strategic path
		List<Junction> endFirstLinkJunctions = null;
		if (sP.size()>1) {
			HashMap<String, List<Junction>> firstLinkJunctions = tacticalHorizonJunctions(pavementNetwork,  sP.get(0), sP.get(1));
			endFirstLinkJunctions = firstLinkJunctions.get("end");
		}
		else {
			endFirstLinkJunctions = new ArrayList<Junction>();
			endFirstLinkJunctions.add(destJ);
		}
		List<RepastEdge<Junction>> tacticalPath = chooseTacticalPath(pavementNetwork, currentJ, endJunctions, heuristic1, heuristic2);
		
		// Create tactical alternative from this path		
		TacticalAlternative tr = setupChosenTacticalAlternative(sP, tacticalNLinks, tacticalPath, currentJ, endFirstLinkJunctions, caG, rG, p);
		
		this.tacticalPath = tr;
		
		return tacticalNLinks;
	}
	
	/*
	 * Method that chooses the tactical path. Does by finding all simple paths to the junction(s) at the end of the
	 * planning horizon using a filtered version of the pavement graph. 
	 * 
	 * These paths are then ranked, first using heurisitc 1 then heuristic 2. The joint shortest paths by these measures are added to a list
	 * and the chosen path is selected at random from that list.
	 * 
	 * @param networkpath<Junction> nP
	 * 		NetworkPath instance used for finding simple paths. Graph paths found on is stored internally.
	 * @param Junction currentJunction
	 * 		The junction to start paths from.
	 * @param colelction<Junction> targetJunctions
	 * 		The destination junctions to find paths to.
	 * @param Transformer<RepastEdge<Junction>> heuristic1
	 * 		Primary metric to measure and rank path length by.
	 * @param Transformer<RepastEdge<Junction>> heuristic2
	 * 		Secondary metric to measure and rank path length by.
	 * 
	 * @returns List<RepastEdge<Junction>>
	 */
	public List<RepastEdge<Junction>> chooseTacticalPath(Network<Junction> pavementNetwork, Junction currentJ, Collection<Junction> targetJunctions, Transformer<RepastEdge<Junction>,Double> heuristic1, Transformer<RepastEdge<Junction>,Double> heuristic2) {
		
		NetworkPath<Junction> nP = new NetworkPath<Junction>(pavementNetwork);

		List<Stack<RepastEdge<Junction>>> candidatePaths = new ArrayList<Stack<RepastEdge<Junction>>>();
		
		// Loop through simple paths to target junctions. Identify those with shortest path length and add to the licat of candicate paths
		Integer minPathLength1 = Integer.MAX_VALUE;
		Integer minPathLength2 = Integer.MAX_VALUE;
		for (Junction tJ: targetJunctions) {
			
			List<Stack<RepastEdge<Junction>>> simplePaths = nP.getSimplePaths(currentJ, tJ);
			
			for (Stack<RepastEdge<Junction>> path : simplePaths) {
				Integer pathLength1 = (int) NetworkPath.getPathLength(path, heuristic1);
				
				if (pathLength1 < minPathLength1) {
					candidatePaths.clear();
					minPathLength1 = pathLength1;
					candidatePaths.add(path);
					
					minPathLength2 = (int) NetworkPath.getPathLength(path, heuristic2);
				}
				// If paths are tied on 1st distance heuristic, use 2nd to discriminate
				else if (pathLength1 == minPathLength1) {
					Integer pathLength2 = (int) NetworkPath.getPathLength(path, heuristic2);
					if (pathLength2 < minPathLength2) {
						candidatePaths.clear();
						minPathLength2 = pathLength2;
						candidatePaths.add(path);
					}
					
					else if (pathLength2 == minPathLength2) {
						candidatePaths.add(path);
					}
				}
			}
		}
		
		// Any paths in candidatePaths have equally low path length when measured using both heuristic 1 and heuristic 2.
		// To choose between these we choose at random
	    Random rand = new Random();
	    List<RepastEdge<Junction>> chosenPath = candidatePaths.get(rand.nextInt(candidatePaths.size()));
	    return chosenPath;
	}
	
	/*
	 * New method for setting up a tactical alternative. The method takes the chosen tactical path along with the strategic path and uses this to 
	 * set up the tactical alternative, which requires identifying at which points in the tactical path crossing locations need to be chosen and how to choose 
	 * crossing locations at those points
	 */
	public static TacticalAlternative setupChosenTacticalAlternative(List<RoadLink> sP, int tacticalNLinks, List<RepastEdge<Junction>> tacticalPath, Junction currentJ, List<Junction> endFirstLinkHorizonJunctions, Geography<CrossingAlternative> caG, Geography<Road> rG, Ped p) {
				
		// Need to split the chosen tactical path into three section sections
		List<RepastEdge<Junction>> initTacticalPath = new ArrayList<RepastEdge<Junction>>(); 
		List<RepastEdge<Junction>> firstLinkTacticalPath = new ArrayList<RepastEdge<Junction>>();
		List<RepastEdge<Junction>> remainderTacticalPath = new ArrayList<RepastEdge<Junction>>();
		
		
		// Split firstLinkTacticalPath from remainderTacticalPath using endFirstLinkHorizonJunctions
		Junction prev = currentJ;
		boolean reachedEndJunction = false;
		int i = 0;
		while (reachedEndJunction == false) {
			RepastEdge<Junction> e = tacticalPath.get(i);
			
			Junction next = null;
			if (e.getSource().equals(prev)) {
				next = e.getTarget();
			}
			else {
				next = e.getSource();
			}
			
			// Check whether previous node / next node is one of the end junctions
			boolean matchPrev = endFirstLinkHorizonJunctions.stream().anyMatch(j -> j.getFID().contentEquals(prev.getFID()));
			boolean matchNext = endFirstLinkHorizonJunctions.stream().anyMatch(j -> j.getFID().contentEquals(next.getFID()));
			
			// Based on these matches decide whether to include the current edge in the path section that goes to the end of the first link
			if (matchPrev & !matchNext) {
				reachedEndJunction = true;
			}
			else if (matchPrev & matchNext) {
				firstLinkTacticalPath.add(e);
				i++;
				reachedEndJunction = true;
			}
			else {
				firstLinkTacticalPath.add(e);
				i++;
			}
		}
		
		remainderTacticalPath = tacticalPath.subList(i, tacticalPath.size()-1);

		
		// Initialise the tactical alternative - sets the path
		TacticalAlternative tr = new TacticalAlternative(p, sP, tacticalNLinks, initTacticalPath, firstLinkTacticalPath, remainderTacticalPath, currentJ, caG, rG);
		tr.updateCurrentJunction();
		
		return tr;
		
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
	 * Get the pavement network junctions around the input node id. Group these junctions according to whether they lie on the input road link or not.
	 * 
	 * @param Network<Junction> pavementNetwork
	 * 			The network of pavement junctions used to connect pedestrian pavements
	 * @param RoadLink rlEndHorz
	 * 			The road link at the end of the tactical planning horizon
	 * @param String endNodeID
	 * 			The road node to get pavement junctions for
	 */
	public static HashMap<String, List<Junction>> tacticalHorizonJunctions(Network<Junction> pavementNetwork, RoadLink rl, String endNodeID) {
		
		// Get the pavement junctions linked to this road node
		List<Junction> pavementJunctions =  roadNodePavementJunctions(pavementNetwork, endNodeID);

		// Loop over pavement junctions and select those touching the road link at the end of the planning horizon
		HashMap<String, List<Junction>> tacticalEndJunctions = new HashMap<String, List<Junction>>();
		tacticalEndJunctions.put("end", new ArrayList<Junction>());
		tacticalEndJunctions.put("outside", new ArrayList<Junction>());
		for (Junction j: pavementJunctions) {
			if (j.getv1rlID().contentEquals(rl.getPedRLID()) | j.getv2rlID().contentEquals(rl.getPedRLID())) {
				tacticalEndJunctions.get("end").add(j);
			}
		}
		
		return tacticalEndJunctions;
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
	
	public TacticalAlternative getTacticalPath() {
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

	public Junction getStartPavementJunction() {
		return startPavementJunction;
	}

	
	public Junction getDestPavementJunction() {
		return destPavementJunction;
	}	
	
}
