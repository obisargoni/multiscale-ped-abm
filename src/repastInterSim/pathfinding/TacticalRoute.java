package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.UnmarkedCrossingAlternative;
import repastInterSim.main.SpaceBuilder;

public class TacticalRoute {
	
	private Ped ped;
	
	private List<RoadLink> strategicPath;
	
	private NetworkPathFinder<Junction> nP;
	private Junction currentJunction = null;
	private RepastEdge<Junction> currentEdge = null;
	private Junction endJunction = null;
	private Junction finalJunction = null;
	private LinkedList<RepastEdge<Junction>> routePath = null;
	private List<RepastEdge<Junction>> initPath; // Path that gets agent from the end of the previous road link they were on to the start of the next. Empty when last junction agent reached boarders previous and next road link.
	private List<RepastEdge<Junction>> pathToEnd; // Path that gets agent from start of tactical horizon to end of tactical horizon
	private List<RepastEdge<Junction>> pathRemainder; // Path that gets agent from first link outside tactical horizon to the end of their destination
	private boolean isBlank = false;
	private AccumulatorRoute accumulator = new AccumulatorRoute();
	
	public TacticalRoute() {
		// Blank constructor
		this.isBlank = true;
	}
	
	public TacticalRoute(Ped p, List<RoadLink> sP, int tNL, List<RepastEdge<Junction>> initTacticalPath, List<RepastEdge<Junction>> firstLinkTacticalPath, List<RepastEdge<Junction>> remainderTacticalPath, Junction startJunction, NetworkPathFinder<Junction> nP) {
		this.ped = p;
		
		this.strategicPath = sP;
		
		this.nP = nP;
		
		this.currentJunction = startJunction;
		this.endJunction = null; // Is this needed?
		this.initPath = initTacticalPath;
		this.pathToEnd = firstLinkTacticalPath;
		this.pathRemainder = remainderTacticalPath;
		
		List<RepastEdge<Junction>> rP = Stream.of(this.initPath, this.pathToEnd).flatMap(Collection::stream).collect(Collectors.toList());		
		this.routePath = rP.stream().collect(Collectors.toCollection(LinkedList::new)); 
		
	}
	
	public Junction getCurrentJunction() {
		return this.currentJunction;
	}
	
	/*
	 * Method for retrieving the coordinate the pedestrian agent should be walking towards.
	 * 
	 * Three possible cases. The tactical route of the pedestrian agent is set by the route junctions. Additional coordinate can be prepended to the route
	 * via the routeCoordinates list. This is used for adding the start and end coordinates of a crossing into the route. Also, when the pedestrian agent is as the end 
	 * of their journey, they will walk towards their destination coordinate rather than the final junction. Otherwise the agent walks towards the next junction in
	 * the route.
	 */
	public Coordinate getTargetCoordinate() {
		// route coordinates storde coordinates of crossing location
		if(this.getAccumulatorRoute().getCrossingCoordinates().size()>0) {
			return this.getAccumulatorRoute().getCrossingCoordinates().getLast();
		}
		// If a crossing is not required (no neeed to cross or ped has completed crossing) and ped is at end of journey (last strategic link and no remaining edges in route path) then walk towards dest coordinate
		else if ( (this.strategicPath.size() == 1) & (this.routePath.size()==0) & (this.accumulator.crossingRequired()==false) ) {
			return this.ped.getDestination().getGeom().getCoordinate();
		}
		else {
			return this.currentJunction.getGeom().getCoordinate();
		}
	}
	
	/*
	 * First tries to remove the top coordiante from the coordinates stack. If the coordinates stack is empty then the top junction is removed from the junction stack.
	 */
	public void updateTargetCoordiante() {		
		if (this.getAccumulatorRoute().getCrossingCoordinates().size()>0) {
			this.getAccumulatorRoute().removeCrossingCoordinate();
		}
		else {
			boolean dontUpdate = this.accumulator.crossingRequired() & (this.accumulator.caChosen() == false) & ( (this.strategicPath.size()==1) | this.accumulator.isDirectCrossing() ); 
			if ( dontUpdate==false ) {
				// do not update the current junction if crossing is required but crossing location is not chosen and:
				// - ped is at the end of their route OR
				// - ped is at a direct crossing
				// In these cases ped must continue to accumulate activation despite reaching their current junction
				updateCurrentJunction();
				
				// Update ped attributes that are recorded for analysis
				// If current edge is null, tactical path needs updating so do not update ped attributes yet.
				if (this.currentEdge != null) {
					this.ped.setChosenCrossingType("none"); // Looks like this happens at the wrong point. Crossing type set to None before crossing ends.
					NetworkEdge<Junction> ne = (NetworkEdge<Junction>) this.currentEdge; 
					this.ped.setCurrentPavementLinkID(ne.getRoadLink().getFID());
				}
			}
		}
	}
	
	/*
	 * Remove first entry from route junctions and then set current junction to new first entry
	 */
	public void updateCurrentJunction() {		
		// Update the current edge			
		this.currentEdge = this.routePath.poll();
		
		// initialise black accumulator initially. Ensures that TacticalRoute accumulators is specific to the 'currentEdge'
		this.accumulator = new AccumulatorRoute();
		
		// Identify the next junction
		Junction nextJunction = null;
		if (this.currentEdge == null) {
			this.finalJunction = this.currentJunction;
			nextJunction = null;
		}
		else {
			nextJunction = edgeAdjacentJunction(this.currentEdge, this.currentJunction);
			
			// Check if this edge requires crossing a road link. If it does, initialise an accumulator to choose crossing location
			// Crossing links only ever contains 
			RoadLink crossingLink = tacticalEdgeCrossingLink(this.currentEdge, SpaceBuilder.orRoadLinkGeography); 
			if (crossingLink != null) {
				
				// AccumulatorRoute requires a default edge and default junction the ped walks towards or stays at whilst choosing a crossing
				
				// For a diagonally crossing, the default jucntion is the junction on the same side of the road, at the otehr end of the road link
				// For a direct crossing it is the pedestrians current junction
				Junction defaultJunction = this.noCrossTargetJunction(this.currentJunction, nextJunction);
				
				// The default edge is the edge that connects the current junction to the default junction.
				// For direct crossing this will be null since current and default junctions are the same, in which case don't change current edge
				RepastEdge<Junction> defaultEdge = this.nP.getNet().getEdge(this.currentJunction, defaultJunction);
				RepastEdge<Junction> targetRouteEdge = this.currentEdge;
				boolean directCrossing; // Used to flag and AccumulatorRoute as choosing a crossing lcoation for a direct rather than diagonal crossing tactical link 
				if (defaultEdge!=null) {
					this.currentEdge = defaultEdge;
					directCrossing=false;						
				}
				else {
					directCrossing=true;
				}
				
				// Get road length - the length of the road that crossing alternatives are being considered for
				double roadLength = crossingLink.getGeom().getLength();
				
				// Identify crossing alternatives
				List<RoadLink> crossingLinks = new ArrayList<RoadLink>();
				crossingLinks.add(crossingLink);
				List<CrossingAlternative> cas = getCrossingAlternatives(SpaceBuilder.caGeography, crossingLinks, ped, SpaceBuilder.roadGeography);
				
				// Initialise accumulator crossing choice model
				this.accumulator = new AccumulatorRoute(this.ped, roadLength, defaultJunction, nextJunction, cas, targetRouteEdge, directCrossing);
				
				// Set target junction to be the default, no crossing, junction while agent chooses crossing location
				nextJunction = this.accumulator.getDefaultJunction();			
			}
		}
		
		this.currentJunction = nextJunction;
	}
	
	public void caChosenUpdateCurrentJunction() {
		// Set the current junction to be the target junction - this
		this.currentJunction = this.accumulator.getTargetJunction();
		this.currentEdge = this.accumulator.getTargetRouteEdge();
	}
	
	public Junction edgeAdjacentJunction(RepastEdge<Junction> e, Junction j) {
		Junction adj = null;
		if (e.getSource().equals(j)) {
			adj = e.getTarget();
		}
		else {
			adj = e.getSource();
		}
		return adj;
	}
	
	public List<RepastEdge<Junction>> getPathToEnd() {
		return pathToEnd;
	}
	
	public List<RepastEdge<Junction>> getRemainderPath() {
		return this.pathRemainder;
	}
	
	public List<RepastEdge<Junction>> getRoutePath(){
		return this.routePath;
	}
	
	public NetworkPathFinder<Junction> getNetworkPathFinder() {
		return this.nP;
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
	public static List<CrossingAlternative> getCrossingAlternatives(Geography<CrossingAlternative> caG, List<RoadLink> rls, Ped p, Geography<Road> rG){
		
		List<CrossingAlternative> cas = new ArrayList<CrossingAlternative>();
		
		// Agent identifies crossing locations on the road links passed in
		// Loop through these and get crossing alternatives that belong to these road links
		for (RoadLink rl: rls) {
			for (CrossingAlternative ca: caG.getAllObjects()) {
				if (ca.getRoadLinkID().contentEquals(rl.getFID())) {
					// Add to list
					cas.add(ca);
				}
			}
		}
		
		// Add unmarked crossing alternative to list
		// Don't set road for unmarked crossings
		UnmarkedCrossingAlternative caU = new UnmarkedCrossingAlternative();
		caU.setPed(p);
		cas.add(caU);
		
		return cas;		
	}
	
	/* 
	 * This method is used to identify the default junction a pedestrian walks towards whilst choosing a crossing location.
	 * 
	 * If the target junction is located at a different road node to the source junction, find the junction that belongs to the same road node
	 * as the target junction, is not the target junction itself, and does not require a road crossing to get to.
	 * 
	 * If the target junction is located at the same road node as the source junction, the no cross target junction is the source junction.
	 */
	public Junction noCrossTargetJunction(Junction sourceJunction, Junction targetJunction) {
		Junction noCrossJ = null;
		
		// get road node id
		String roadNodeID = targetJunction.getjuncNodeID();
		
		if (sourceJunction.getjuncNodeID().contentEquals(roadNodeID)) {
			noCrossJ = sourceJunction;
		}
		else {
			// Loop through junctions connected to the source junction
			for (Junction j:this.nP.getNet().getAdjacent(sourceJunction)) {
				
				if ( (j.getjuncNodeID().contentEquals(roadNodeID)) & !(j.getFID().contentEquals(targetJunction.getFID())) ) {
					noCrossJ = j;
					break;
				}
			}
		}
		
		return noCrossJ;
	}
	
	public Junction getEndJunction() {
		return this.endJunction;
	}
	
	public Junction getFinalJunction() {
		return this.finalJunction;
	}
	
	public AccumulatorRoute getAccumulatorRoute() {
		return this.accumulator;
	}
	
	public void step() {
		this.accumulator.step();
		
		// If a crossing has been chosen, update the tactical path to reflect this
		if (this.accumulator.getChosenCA() != null) {
			caChosenUpdateCurrentJunction();
			
			// Update the ped attributes that are collected for data analysis
			this.ped.setChosenCrossingType(this.accumulator.getChosenCA().getType());
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) this.currentEdge; 
			this.ped.setCurrentPavementLinkID(ne.getRoadLink().getFID());
		}
	}
	
	public boolean isBlank() {
		return this.isBlank;
	}
	
	static private List<RoadLink> tacticalPathPrimaryCrossingLinks(RepastEdge<Junction> edge, List<RoadLink> sP) {
		List<RepastEdge<Junction>> path = new ArrayList<RepastEdge<Junction>>();
		path.add(edge);
		return tacticalPathPrimaryCrossingLinks(path, sP);
	}
	
	static private List<RoadLink> tacticalPathPrimaryCrossingLinks(List<RepastEdge<Junction>> path, List<RoadLink> sP) {
		List<RoadLink> crossedLinks = new ArrayList<RoadLink>();
		for (RepastEdge<Junction> e: path) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) e;
			for (RoadLink rl : sP) {
				if (rl.getPedRLID().contentEquals(ne.getRoadLink().getPedRLID())) {
					crossedLinks.add(rl);
				}
			}
		}		
		return crossedLinks;
	}
	
	static private RoadLink tacticalEdgeCrossingLink(RepastEdge<Junction> edge, Geography<RoadLink> orRLG) {
		RoadLink crossedLink = null;
		NetworkEdge<Junction> ne = (NetworkEdge<Junction>) edge;
		for (RoadLink rl : orRLG.getAllObjects()) {
			if (rl.getPedRLID().contentEquals(ne.getRoadLink().getPedRLID())) {
				crossedLink = rl;
				break;
			}
		}
		return crossedLink;
	}

	public RepastEdge<Junction> getCurrentEdge() {
		return this.currentEdge;
	}
}