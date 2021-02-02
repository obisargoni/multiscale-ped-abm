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
	private Junction finalJunction = null;
	private List<CrossingAlternative> crossingAlternatives;
	private LinkedList<RepastEdge<Junction>> routePath = null;
	private LinkedList<Coordinate> routeCoordinates = new LinkedList<Coordinate>();
	private Coordinate destCoordinate = null;
	private List<RepastEdge<Junction>> initPath; // Path that gets agent from the end of the previous road link they were on to the start of the next. Empty when last junction agent reached boarders previous and next road link.
	private List<RepastEdge<Junction>> pathToEnd; // Path that gets agent from start of tactical horizon to end of tactical horizon
	private List<RepastEdge<Junction>> pathRemainder; // Path that gets agent from first link outside tactical horizon to the end of their destination
	private boolean recurringEndJunction = false;
	private boolean isBlank = false;
	private AccumulatorRoute accumulator = new AccumulatorRoute();
	
	public TacticalAlternative() {
		// Blank constructor
		this.isBlank = true;
	}
	
	public TacticalAlternative(List<RoadLink> tacticalStrategicPath, List<RepastEdge<Junction>> initTacticalPath, List<RepastEdge<Junction>> firstLinkTacticalPath, List<RepastEdge<Junction>> remainderTacticalPath, Junction startJunction) {
		this.currentJunction = startJunction;
		this.endJunction = ; // Is this needed?
		this.initPath = initTacticalPath;
		this.pathToEnd = firstLinkTacticalPath;
		this.pathRemainder = remainderTacticalPath;
		this.routePath = (LinkedList<RepastEdge<Junction>>) Stream.of(this.initPath, this.pathToEnd).flatMap(Collection::stream).collect(Collectors.toList());
		
	}
	
	public TacticalAlternative(NetworkPath<Junction> nP, Junction startJunction, Junction endJunction) {
		this.nP = nP;
		this.currentJunction = startJunction;
		this.endJunction = endJunction;
		this.pathToEnd = new ArrayList<RepastEdge<Junction>>();
		this.pathEndToOutside = new ArrayList<RepastEdge<Junction>>();
		this.pathRemainder = new ArrayList<RepastEdge<Junction>>();
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
		else if ((this.destCoordinate != null) & (this.routePath.size()==1)) {
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
		if (this.recurringEndJunction & (this.routePath.size()==1)) {
			// do not update the current junction
			
		}
		
		else {
			// Update the current edge
			this.routePath.pollFirst();
			this.currentEdge = this.routePath.peekFirst();
			
			// Identify the next junction
			Junction nextJunction = null;
			if (this.currentEdge == null) {
				this.finalJunction = this.currentJunction;
				nextJunction = null;
			}
			else {
				if (this.currentEdge.getSource().equals(this.currentJunction)) {
					nextJunction = this.currentEdge.getTarget();
				}
				else {
					nextJunction = this.currentEdge.getSource();
				}
				
				// Check if this edge requires crossing a primary link. If it does, initialise an accumulator to choose crossing location
				int primaryCrossingParity = RoadNetworkRoute.calculateRouteParity(this.currentJunction.getGeom().getCoordinate(), nextJunction.getGeom().getCoordinate(), sP);
				if (primaryCrossingParity == 1) {
					
					Junction defaultJunction = this.noCrossTargetJunction(this.currentJunction, nextJunction);
					
					// Record the desired route path, and overwrite route path with the now default route path
					LinkedList<RepastEdge<Junction>> targetRoutePath = this.routePath;
					RepastEdge<Junction> defaultEdge = this.nP.getNet().getEdge(this.currentJunction, nextJunction);
					this.routePath = new LinkedList<RepastEdge<Junction>>();
					this.routePath.add(defaultEdge);
					
					// Identify crossing alternatives
					List<CrossingAlternative> cas = getCrossingAlternatives(caG, tSP, p, rG);
					
					// Initialise accumulator crossing choice model
					this.accumulator = new AccumulatorRoute(p, pHLength, tr, defaultDest, targetJunction, targetRoutePath);
					
					// Set target junction to be the default, no crossing, junction while agent chooses crossing location
					nextJunction = this.accumulator.getDefaultJunction();			
				}
			}
			
			this.currentJunction = nextJunction;
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
	/*
	 * Given the source and target junctions of an edge of the tactical network, return the junction that belongs to the same road node
	 * as the target junction, is not the target junction itself, and does not require a road crossing to get to. If no such junction is 
	 * available return null.
	 * 
	 * This method is used to identify the default junction a pedestrian walks towards whilst choosing a crossing location.
	 */
	public Junction noCrossTargetJunction(Junction sourceJunction, Junction targetJunction) {
		Junction noCrossJ = null;
		
		// get road node id
		String roadNodeID = targetJunction.getjuncNodeID();
		
		// Loop through junctions connected to the source junction
		for (Junction j:this.nP.getNet().getAdjacent(sourceJunction)) {
			
			if ( (j.getjuncNodeID().contentEquals(roadNodeID)) & !(j.getFID().contentEquals(targetJunction.getFID())) ) {
				noCrossJ = j;
				break;
			}
		}
		
		return noCrossJ;
	}

	void addCoordinate(Coordinate c) {
		this.routeCoordinates.push(c);
	}
	
	public Junction getEndJunction() {
		return this.endJunction;
	}
	
	public Junction getFinalJunction() {
		return this.finalJunction;
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
	
	public void step() {
		this.accumulator.step();
		
		// If a crossing has been chosen, update the tactical path to reflect this
		if (this.accumulator.getChosenCA() != null) {
			CrossingAlternative ca = this.accumulator.getChosenCA();
			
			// Add the coordinates of the start and end of the crossing to the route
			this.addCoordinate(ca.nearestCoord(this.ped.getLoc()));
			this.addCoordinate(ca.farthestCoord(this.ped.getLoc()));
			
			// Set the current junction to be the target junction - this
			this.currentJunction = this.accumulator.getTargetJunction();
			this.routePath = this.accumulator.getTargetRoutePath();			
			
			// Finally record which crossing type the pedestrian agent chose
			this.ped.setChosenCrossingType(ca.getType());
		}
	}
	
	public boolean isBlank() {
		return this.isBlank;
	}
}