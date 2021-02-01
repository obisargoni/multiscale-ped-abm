package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.graph.RepastEdge;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;
import repastInterSim.environment.RoadLink;

public class TacticalAlternative {
	
	private NetworkPath<Junction> nP;
	private Junction currentJunction = null;
	private RepastEdge<Junction> currentEdge = null;
	private Junction endJunction = null;
	private Junction outsideJunction = null;
	private List<CrossingAlternative> crossingAlternatives;
	private LinkedList<Junction> routeJunctions = null;
	private LinkedList<Coordinate> routeCoordinates = new LinkedList<Coordinate>();
	private Coordinate destCoordinate = null;
	private List<RepastEdge<Junction>> initPath; // Path that gets agent from the end of the previous road link they were on to the start of the next. Empty when last junction agent reached boarders previous and next road link.
	private List<RepastEdge<Junction>> pathToEnd; // Path that gets agent from start of tactical horizon to end of tactical horizon
	private List<RepastEdge<Junction>> pathRemainder; // Path that gets agent from first link outside tactical horizon to the end of their destination
	private boolean recurringEndJunction = false;
	
	public TacticalAlternative() {
		// Blank constructor
	}
	
	public TacticalAlternative(List<RoadLink> tacticalStrategicPath, List<RepastEdge<Junction>> initTacticalPath, List<RepastEdge<Junction>> firstLinkTacticalPath, List<RepastEdge<Junction>> remainderTacticalPath, Junction startJunction) {
		this.currentJunction = startJunction;
		this.endJunction = ; // Is this needed?
		this.initPath = initTacticalPath;
		this.pathToEnd = firstLinkTacticalPath;
		this.pathRemainder = remainderTacticalPath;
		
	}
	
	public TacticalAlternative(NetworkPath<Junction> nP, Junction startJunction, Junction endJunction) {
		this.nP = nP;
		this.currentJunction = startJunction;
		this.endJunction = endJunction;
		this.pathToEnd = new ArrayList<RepastEdge<Junction>>();
		this.pathEndToOutside = new ArrayList<RepastEdge<Junction>>();
		this.pathRemainder = new ArrayList<RepastEdge<Junction>>();
	}
	
	private void setRouteJunctions() {
		// Get junctions from path
		List<RepastEdge<Junction>> path = this.getRoutePath();
		LinkedList<Junction> nodesToEnd = nP.nodePathFromEdges(path, this.currentJunction);
		this.routeJunctions = nodesToEnd;
	}
	
	public LinkedList<Junction> getRouteJunctions(){
		if (this.routeJunctions == null) {
			setRouteJunctions();
		}
		return this.routeJunctions;
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
		if(this.routeCoordinates.size()>0) {
			return routeCoordinates.getLast();
		}
		else if ((this.destCoordinate != null) & (this.getRouteJunctions().size()==1)) {
			return this.destCoordinate;
		}
		else {
			return this.currentJunction.getGeom().getCoordinate();
		}
	}
	
	/*
	 * First tries to remove the top coordiante from the coordinates stack. If the coordinates stack is empty then the top junction is removed from the junction stack.
	 */
	public void updateTargetCoordiante() {
		if (this.routeCoordinates.size()>0) {
			this.routeCoordinates.removeLast();
		}
		else {
			updateCurrentJunction();
		}
	}
	
	/*
	 * Remove first entry from route junctions and then set current junction to new first entry
	 */
	public void updateCurrentJunction() {
		if (this.recurringEndJunction & (this.getRouteJunctions().size()==1)) {
			// do not update the current junction
			
		}
		else {
			// Update the current junction
			updateCurrentEdge();
			
			Junction nextJunction = null;
			if (this.currentEdge.getSource().equals(this.currentJunction)) {
				nextJunction = this.currentEdge.getTarget();
			}
			else {
				nextJunction = this.currentEdge.getSource();
			}
			this.currentJunction = nextJunction;	
		}
	}
	
	/*
	 * Set the current edge to be either the next edge in the initial path or the next edge in the tactical path that goes up to the end of the 
	 * current strategic road link.
	 */
	public void updateCurrentEdge() {
		if (this.initPath.size() > 0) {
			this.currentEdge = this.initPath.get(0);
			this.initPath.remove(0);
		}
		else if (this.pathToEnd.size() > 0) {
			this.currentEdge = this.initPath.get(0);
			this.initPath.remove(0);
		}
	}
	
	public List<RepastEdge<Junction>> getPathToEnd() {
		return pathToEnd;
	}
	
	public void setPathToEnd(List<RepastEdge<Junction>> path) {
		this.pathToEnd = path;
	}
	
	public void setPathToEnd() {
		this.pathToEnd = nP.getShortestPath(this.currentJunction, this.endJunction);
	}
	
	public List<RepastEdge<Junction>> getPathEndToOutside() {
		return pathEndToOutside;
	}
	
	public void setPathEndToOutside(List<RepastEdge<Junction>> path) {
		this.pathEndToOutside = path;
	}
	
	public void setOutsideJunction(Junction j) {
		this.outsideJunction = j;
	}
	
	public Junction getOutsideJunction() {
		return this.outsideJunction;
	}
	
	public List<RepastEdge<Junction>> getAlternativeRemainderPath() {
		return this.pathRemainder;
	}
	
	public void setAlternativeRemainderPath(List<RepastEdge<Junction>> path) {
		this.pathRemainder = path;
	}
	
	/*
	 * Create the path of pavement edges along the tactical route on the fly by combining the path to the end of the tactical horizon and the
	 * path from the end to the start of the next tactical horizon.
	 */
	public List<RepastEdge<Junction>> getRoutePath() {
		return Stream.of(this.pathToEnd, this.pathEndToOutside).flatMap(Collection::stream).collect(Collectors.toList());
	}

	public void setCrossingAlternatives(List<CrossingAlternative> cas) {
		this.crossingAlternatives = cas;
	}
	
	public List<CrossingAlternative> getCrossingAlternatives() {
		return this.crossingAlternatives;
	}

	public void updatePathToEnd(Junction j) {
		// Get the earliest route junction that is connected to input junction j and belongs to the same road network node as j
		// Use this as the new starting junction.
		String roadNodeID = j.getjuncNodeID();
		boolean breakLoop = false;
		
		// Start at first route coord and continue searching along route until match found
		for (Junction rj:this.getRouteJunctions()) {
			
			// Get adjacent junctions and find the one that belongs to the same road node and matches the route junction
			for (Junction aj: this.nP.getNet().getAdjacent(j)) {
				if ( (aj.getjuncNodeID().contentEquals(roadNodeID)) & (aj.getFID().contentEquals(rj.getFID())) ) {
					this.currentJunction = aj;
					breakLoop = true;
					break;
				}
			}
			if (breakLoop) {
				break;
			}
		}
		
		// Get path from this junction to the end junction and set route junctions to null so that they will be re calculated
		this.setPathToEnd();
		this.routeJunctions = null;
	}

	void addCoordinate(Coordinate c) {
		this.routeCoordinates.push(c);
	}
	
	public Junction getEndJunction() {
		return this.endJunction;
	}

	public void setDestinationCoordinate(Coordinate dC) {
		this.destCoordinate = dC;
	}
	
	/*
	 * Used to control whether the end junction of this tactical alternative is recurring or not.
	 */
	public void setRecurringEndJunction(boolean b) {
		this.recurringEndJunction = b;
	}
}