package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Stack;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.Transformer;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.LineString;

import repast.simphony.random.RandomHelper;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.OD;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.transformers.EdgeRoadLinkIDTransformer;
import repastInterSim.pathfinding.transformers.EdgeWeightTransformer;
import repastInterSim.environment.NetworkEdge;

public class PedPathFinder {
		
	private Ped ped;
    private Transformer<RepastEdge<Junction>,Integer> primaryCostHeuristic;
    private Transformer<RepastEdge<Junction>,Integer> secondaryCostHeuristic;
	
	private OD origin;
	private OD destination;
	
	private String fullStrategicPathString = ""; // Used to analyse crossings locations in relation to a pedestrian's route
	private List<RoadLink> strategicPath;
	private static int defaultLinksPerTacticalUpdate = 1;
	private boolean firstUpdateDone = false;
	private Junction startPavementJunction;
	private Junction destPavementJunction;
	private TacticalRoute tacticalPath = new TacticalRoute();
	private String fullTacticalPathString = "";
	private int nLinksPerTacticalUpdate;
	
	private Coordinate nextCrossingCoord;	

	public PedPathFinder(OD o, OD d, Geography<Junction> paveG, Network<Junction> paveNetwork, boolean minimiseCrossings) {
		init(o, d, paveG, paveNetwork, minimiseCrossings);
	}
	
	public PedPathFinder(Ped p, Geography<Junction> paveG, Network<Junction> paveNetwork, boolean minimiseCrossings) {
		this.ped = p;
		init(p.getOrigin(), p.getDestination(), paveG, paveNetwork, minimiseCrossings);
	}
	
	private void init(OD o, OD d, Geography<Junction> paveG, Network<Junction> paveNetwork, boolean minimiseCrossings) {
		
		this.origin = o;
		this.destination = d;
				
		planStrategicPath(this.origin.getGeom().getCoordinate(), this.destination.getGeom().getCoordinate(), paveG, paveNetwork);
		
		if (minimiseCrossings) {
			this.primaryCostHeuristic = new EdgeRoadLinkIDTransformer<Junction>();
			this.secondaryCostHeuristic = new EdgeWeightTransformer<Junction>();
		}
		else {
			this.primaryCostHeuristic = new EdgeWeightTransformer<Junction>();
			this.secondaryCostHeuristic = new EdgeRoadLinkIDTransformer<Junction>();
		}
		
		this.nLinksPerTacticalUpdate = PedPathFinder.defaultLinksPerTacticalUpdate;
	}
	
	public void step() {
		this.tacticalPath.step();
	}
	
	/**
	 * Initialise a new road link routing object that can be used to find a path along a topological road network.
	 * Use this to identify the shortest path through the network and assign this path to this classes' strategic path attribute.
	 * 
	 */
	public void planStrategicPath(Coordinate oC, Coordinate dC, Geography<Junction> paveG, Network<Junction> paveNetwork) {
		RoadNetworkRoute rnr = new RoadNetworkRoute(oC, dC);
		
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
		for (RoadLink rl: this.strategicPath) {
			this.fullStrategicPathString = this.fullStrategicPathString + ":" + rl.getFID();
		}
		
		
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
		// Only remove links from strategic path after first update
		if (firstUpdateDone==false) {
			firstUpdateDone = true;
		}
		else {
			for (int i = 0; i < this.nLinksPerTacticalUpdate; i++) {
				this.strategicPath.get(0).getPeds().remove(this.ped);
				this.strategicPath.remove(0);
			}
			this.nLinksPerTacticalUpdate = PedPathFinder.defaultLinksPerTacticalUpdate;
		}
		
		// Add ped to next road link it is walking along
		this.strategicPath.get(0).getPeds().add(this.ped);
		
		// If no tactical path has been set use the strategic path start junction, otherwise set the start junction as the end junction of previous tactical path
		Junction startJunction = null;
		if (this.tacticalPath.isBlank()) {
			startJunction = this.startPavementJunction;
		}
		else {
			startJunction = this.tacticalPath.getFinalJunction();
		}
		
		// Initialise Accumulator Route that agent will use to navigate along the planning horizon, and update the number of links in the tactical planning horizon
		// Calculate number of links in planning horizon
		int tacticalHorizonLinks = getNLinksWithinAngularDistance(this.strategicPath, this.ped.getpHorizon());
		Geography<Road> roadGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.ROAD_GEOGRAPHY);
		Geography<CrossingAlternative> caGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.CA_GEOGRAPHY);
		Network<Junction> pavementNetwork = SpaceBuilder.getNetwork(GlobalVars.CONTEXT_NAMES.PAVEMENT_NETWORK);
		this.tacticalPath.clear();
		this.tacticalPath = planTacticalPath(pavementNetwork, caGeography, roadGeography, tacticalHorizonLinks, this.ped, this.strategicPath, startJunction, this.destPavementJunction, this.primaryCostHeuristic, this.secondaryCostHeuristic);
    }
	
	/*
	 * Plan a tactical level path using the accumulator crossing choice path finding model.
	 */
	public TacticalRoute planTacticalPath(Network<Junction> pavementNetwork, Geography<CrossingAlternative> caG, Geography<Road> rG, int nTL, Ped p, List<RoadLink> sP, Junction currentJ, Junction destJ, Transformer<RepastEdge<Junction>,Integer> heuristic1, Transformer<RepastEdge<Junction>,Integer> heuristic2) {
		
		// NetworkPath object is used to find paths on the pavement network
		NetworkPathFinder<Junction> nP = new NetworkPathFinder<Junction>(pavementNetwork);
		
		boolean destInPlanningHorizon = false;
		if (nTL == sP.size()) {
			destInPlanningHorizon = true;
		}
		
		// Get road link ID of link at end of planning horizon and first strategic path road link outside of planning horizon
		List<Junction> endJunctions = null;
		if (destInPlanningHorizon==false) {
			RoadLink rlEndHorz = sP.get(nTL-1);
			RoadLink rlOutHorz = sP.get(nTL);
			HashMap<String, List<Junction>> horizonJunctions = tacticalHorizonJunctions(pavementNetwork, rlEndHorz, rlOutHorz);
			endJunctions = horizonJunctions.get("end");
		}
		else {
			endJunctions = new ArrayList<Junction>();
			endJunctions.add(destJ);
		}
		
		// Get filter to use for selecting pavement network nodes that lie in the planning horizon only
		List<Junction> horizonStrategicNodes = new ArrayList<Junction>();
		List<RepastEdge<Junction>> horizonTacticalEdges = new ArrayList<RepastEdge<Junction>>();
		for( int i=0; i<nTL; i++ ) {
			sP.get(i).getJunctions().stream().forEach(horizonStrategicNodes::add);
			
			List<Junction> orLinkJunctions = sP.get(i).getJunctions();
			
			List<Junction> roadLinkPaveJuncs = new ArrayList<Junction>();
			SpaceBuilder.orJuncToPaveJunc.get(orLinkJunctions.get(0)).stream().forEach(roadLinkPaveJuncs::add);
			SpaceBuilder.orJuncToPaveJunc.get(orLinkJunctions.get(1)).stream().forEach(roadLinkPaveJuncs::add);
			
			for (Junction paveJA: roadLinkPaveJuncs) {
				for (Junction paveJB: roadLinkPaveJuncs) {
					RepastEdge<Junction> paveEdge = pavementNetwork.getEdge(paveJA, paveJB);
					if (paveEdge != null) {
						horizonTacticalEdges.add(paveEdge);
					}
				}
			}
		}
		
		Predicate<RepastEdge<Junction>> tacticalEdgeFilter = e -> horizonTacticalEdges.stream().anyMatch( n -> n.equals(e));
						
		// Choose path to end of tactical horizon
		List<RepastEdge<Junction>> tacticalPath = chooseTacticalPathWithHorizon(nP, tacticalEdgeFilter, currentJ, endJunctions, heuristic1, heuristic2);
		
		// In cases where informal crossing is restricted it can be impossible to find a pavement network route within the planning horizon
		// In this case do not enforce planning horizon by not filtering any links 
	    if (tacticalPath==null) {
	    	tacticalPath = chooseTacticalPathNoHorizon(nP, currentJ, endJunctions, heuristic1, heuristic2);
			
	    	// Update strategic path here so that it matches tactical path. Only required where tactical path is planned without planning horizon filter
	    	updateStrategicPath(tacticalPath, nTL);
	    	sP = this.strategicPath;
	    }
		
		// Create tactical alternative from this path
		TacticalRoute tr = setupChosenTacticalAlternative(nP, sP, nTL, tacticalPath, currentJ, destJ, caG, rG, p);
		tr.updateTargetCoordiante();

		return tr;
	}
	
	/*
	 * Method that chooses the tactical path. Does by finding all simple paths to the junction(s) at the end of the
	 * planning horizon using a filtered version of the pavement graph. 
	 * 
	 * These paths are then ranked, first using heurisitc 1 then heuristic 2. The joint shortest paths by these measures are added to a list
	 * and the chosen path is selected at random from that list.
	 * 
	 * @param Network<Junction> pavementNetwork
	 * 		The network used for tactical routing
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
	public List<RepastEdge<Junction>> chooseTacticalPathWithHorizon(NetworkPathFinder<Junction> nP, Predicate<RepastEdge<Junction>> filter, Junction currentJ, Collection<Junction> targetJunctions, Transformer<RepastEdge<Junction>,Integer> heuristic1, Transformer<RepastEdge<Junction>,Integer> heuristic2) {

		List<Stack<RepastEdge<Junction>>> candidatePaths = nP.getAllShortestPaths(currentJ, targetJunctions, filter, heuristic1);
		
		// Due to removing crossing links from network can be unable to plan path using only planning horizon links.
		// In this case do not enforce planning horizon by not filtering any links
		if (candidatePaths.size()==0) {
			return null;
		}

		// If multiple paths have same length by this heuristic filter again using second heuristic
		if (candidatePaths.size()>1) {
			candidatePaths = nP.getShortestOfMultiplePaths(candidatePaths, heuristic2);
		}
			
		// Any paths in candidatePaths have equally low path length when measured using both heuristic 1 and heuristic 2.
		// To choose between these we choose at random
	    //int pathIndex = RandomHelper.nextIntFromTo(0, candidatePaths.size()-1);
	    List<RepastEdge<Junction>> chosenPath = candidatePaths.get(0);
	    
	    return chosenPath;
	}
	
	public List<RepastEdge<Junction>> chooseTacticalPathNoHorizon(NetworkPathFinder<Junction> nP, Junction currentJ, Collection<Junction> targetJunctions, Transformer<RepastEdge<Junction>,Integer> heuristic1, Transformer<RepastEdge<Junction>,Integer> heuristic2) {
		
		// Getting all shortest paths without a filter can take a long time when using the number of crossings transformer
		// To reduce computation identify an initial set of candidate paths based on link length, from which to filter based on heuristics. k=3 will return 3 shortest paths to each target node.
		Transformer<RepastEdge<Junction>, Integer> distanceTransformer = new EdgeWeightTransformer<Junction>();
		List<Stack<RepastEdge<Junction>>> candidatePaths = nP.getKShortestPaths(currentJ, targetJunctions, null, distanceTransformer, 3);

		// If multiple paths have same length by this heuristic filter again using second heuristic
		if (candidatePaths.size()>1) {
			candidatePaths = nP.getShortestOfMultiplePaths(candidatePaths, heuristic1);
		}
		
		if (candidatePaths.size()>1) {
			candidatePaths = nP.getShortestOfMultiplePaths(candidatePaths, heuristic1);
		}
			
		// Any paths in candidatePaths have equally low path length when measured using both heuristic 1 and heuristic 2.
		// To choose between these we choose at random
	    //int pathIndex = RandomHelper.nextIntFromTo(0, candidatePaths.size()-1);
	    List<RepastEdge<Junction>> chosenPath = candidatePaths.get(0);
	    
	    return chosenPath;
	}
	
	/*
	 * Sets the strategic path road links to match with the road links the tactical path traverses
	 */
	public void updateStrategicPath(List<RepastEdge<Junction>> chosenPath, int nTL) {
		
		List<RoadLink> sP = roadPathFromPavementPath(chosenPath);
		
		this.strategicPath.get(0).getPeds().remove(this.ped);
		// remove previous strategic links that lead up to tactical horizon
		for (int i=0; i<nTL; i++) {
			this.strategicPath.remove(0);
		}
		
		// Now prepend strategic path with new set of links that lead up to tactical horizon
		this.strategicPath = Stream.of(sP, this.strategicPath).flatMap(Collection::stream).collect(Collectors.toList());
		this.strategicPath.get(0).getPeds().add(this.ped);
	}
	
	private static List<RoadLink> roadPathFromPavementPath(List<RepastEdge<Junction>> pavementPath){
		Geography<Junction> orJunctionGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.OR_JUNCTION_GEOGRAPHY);
		Network<Junction> orNetwork = SpaceBuilder.getNetwork(GlobalVars.CONTEXT_NAMES.OR_ROAD_NETWORK);
		List<RoadLink> rp = new ArrayList<RoadLink>();
		
		// Loop over pavement network links, get the road network junctions they correspond to and from these get the road network links they travel along
		for (RepastEdge<Junction> paveEdge: pavementPath) {
			String orJ1 = paveEdge.getSource().getjuncNodeID();
			String orJ2 = paveEdge.getTarget().getjuncNodeID();
			
			Junction j1=null;
			Junction j2=null;
			for(Junction j: orJunctionGeography.getAllObjects()) {
				if (j.getFID().contentEquals(orJ1)) {
					j1=j;
				}
				if (j.getFID().contentEquals(orJ2)) {
					j2=j;
				}
			}
			
			RepastEdge<Junction> e = orNetwork.getEdge(j1, j2);
			if (e==null) {

				continue;
			}
			
			RoadLink rl = ( (NetworkEdge<Junction>) e).getRoadLink();
			rp.add(rl);
		}
		
		// Now remove duplicates so that road link route doesn't involve doubling back
		int i=0;
		while(i<rp.size()-1) {
			RoadLink rl1 = rp.get(i);
			RoadLink rl2 = rp.get(i+1);
			if (rl1.getFID().contentEquals(rl2.getFID())) {
				rp.remove(rl1);
				rp.remove(rl2);
				i=0; // Start searching from beginning again
			}
			i++;
		}
		
		return rp;
	}

	/*
	 * New method for setting up a tactical alternative. The method takes the chosen tactical path along with the strategic path and uses this to 
	 * set up the tactical alternative, which requires identifying at which points in the tactical path crossing locations need to be chosen and how to choose 
	 * crossing locations at those points
	 */
	public TacticalRoute setupChosenTacticalAlternative(NetworkPathFinder<Junction> nP, List<RoadLink> sP, int tacticalNLinks, List<RepastEdge<Junction>> tacticalPath, Junction currentJ, Junction destJ, Geography<CrossingAlternative> caG, Geography<Road> rG, Ped p) {
				
		// Need to split the chosen tactical path into three section sections
		List<RepastEdge<Junction>> initTacticalPath = new ArrayList<RepastEdge<Junction>>(); 
		List<RepastEdge<Junction>> firstLinkTacticalPath = new ArrayList<RepastEdge<Junction>>();
		List<RepastEdge<Junction>> remainderTacticalPath = new ArrayList<RepastEdge<Junction>>();
		
		// Need to get the junctions at the end of the first link in strategic path
		// It is possible that tactical path bypasses first n links of strategic path (occurs when planning horizon filter is removed) in which case find junctions of tactical path that first meet up with strategic path
		List<Junction> endFirstLinkJunctions = null;
		Integer indexEndFirstLinkPath = null;
		while ( (sP.size()>this.nLinksPerTacticalUpdate) & (indexEndFirstLinkPath==null) ) {
			HashMap<String, List<Junction>> tacticalHorzJuncs = tacticalHorizonJunctions(nP.getNet(),  sP.get(this.nLinksPerTacticalUpdate-1), sP.get(this.nLinksPerTacticalUpdate));
			endFirstLinkJunctions = tacticalHorzJuncs.get("end");
			indexEndFirstLinkPath = getIndexOfEdgeThatReachesTargetJunctions(tacticalPath, currentJ, endFirstLinkJunctions);
			
			// If end juncs don't match up with tacticalPath, try outside junctions
			if(indexEndFirstLinkPath==null) {
				endFirstLinkJunctions = tacticalHorzJuncs.get("outside");
				indexEndFirstLinkPath = getIndexOfEdgeThatReachesTargetJunctions(tacticalPath, currentJ, endFirstLinkJunctions);
			}
			
			// If still no match move onto next road link.
			if (indexEndFirstLinkPath==null) {
				this.nLinksPerTacticalUpdate++;
			}
		}
		
		// if index still null means that tactical path extends to the destination junction
		if (indexEndFirstLinkPath==null) {
			endFirstLinkJunctions = new ArrayList<Junction>();
			endFirstLinkJunctions.add(destJ);
			indexEndFirstLinkPath = getIndexOfEdgeThatReachesTargetJunctions(tacticalPath, currentJ, endFirstLinkJunctions);
			this.nLinksPerTacticalUpdate=sP.size()-1;
		}
		
		for (int i=0;i<tacticalPath.size();i++) {
			if (i<indexEndFirstLinkPath) {
				firstLinkTacticalPath.add(tacticalPath.get(i));
			}
			else {
				remainderTacticalPath.add(tacticalPath.get(i));
			}
		}
		
		// Initialise the tactical alternative - sets the path
		TacticalRoute tr = new TacticalRoute(p, sP, tacticalNLinks, initTacticalPath, firstLinkTacticalPath, remainderTacticalPath, currentJ, nP);
		
		return tr;
		
	}
	
	public static Integer getIndexOfEdgeThatReachesTargetJunctions(List<RepastEdge<Junction>> path, Junction startJunction, List<Junction> targetJunctions) {
		Junction prev = startJunction;
		boolean reachedEndJunction = false;
		Integer i = 0;
		while ( (reachedEndJunction == false) & (path.size()>0) ) {
			
			// It is possible that the target junction can't be found, in which case return null
			if (i>=path.size()) {
				i = null;
				break;
			}
			
			RepastEdge<Junction> e = path.get(i);
			
			Junction next = null;
			if (e.getSource().getGeom().equals(prev.getGeom())) {
				next = e.getTarget();
			}
			else {
				next = e.getSource();
			}
			
			final Junction prevFinal = prev;
			final Junction nextFinal = next;
			
			// Check whether previous node / next node is one of the end junctions
			boolean matchPrev = targetJunctions.stream().anyMatch(j -> j.getFID().contentEquals(prevFinal.getFID()));
			boolean matchNext = targetJunctions.stream().anyMatch(j -> j.getFID().contentEquals(nextFinal.getFID()));
			
			// Based on these matches decide whether to include the current edge in the path section that goes to the end of the first link
			if (matchPrev & !matchNext) {
				reachedEndJunction = true;
			}
			// If next junction matches target and the end of the path has been reached then include this link and end loop
			else if ( (i == path.size()-1) & matchNext) {
				i++;
				reachedEndJunction = true;
			}
			else {
				i++;
			}
			
			prev = next;
		}
		
		return i;
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
	
	public String getFullStrategicPathString(){
		return this.fullStrategicPathString;
	}
	
	public TacticalRoute getTacticalPath() {
		return this.tacticalPath;
	}
	
	public void addTacticalLinkToFullTacticalPathString(String rlID) {
		this.fullTacticalPathString += rlID + ":";
	}
	
	public String getFullTacticalPathString() {
		return this.fullTacticalPathString;
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
	
	public Transformer<RepastEdge<Junction>,Integer> getPrimaryCostHeuristic() {
		return this.primaryCostHeuristic;
	}
	
	public Transformer<RepastEdge<Junction>,Integer> getSecondaryCostHeuristic() {
		return this.secondaryCostHeuristic;
	}
	
	public void clear() {
		this.ped=null;
	    this.primaryCostHeuristic=null;
	    this.secondaryCostHeuristic=null;
		this.origin=null;
		this.destination=null;
		this.fullStrategicPathString = null;
		this.strategicPath=null;
		this.startPavementJunction=null;
		this.destPavementJunction=null;
		this.tacticalPath.clear();
		this.tacticalPath=null;
		this.nextCrossingCoord=null;	
	}
	
}
