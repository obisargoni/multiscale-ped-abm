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

public class TacticalAlternative {
	
	private NetworkPath<Junction> nP;
	private Junction currentJunction;
	private Junction endJunction;
	private List<CrossingAlternative> crossingAlternatives;
	private LinkedList<Junction> routeJunctions = null;
	private LinkedList<Coordinate> routeCoordinates = new LinkedList<Coordinate>();; 
	private List<RepastEdge<Junction>> pathToEnd; // Path that gets agent from start of tactical horizon to end of tactical horizon
	private List<RepastEdge<Junction>> pathEndToOutside; // Path that gets agent from end of tactical horizon to first link outside of tactical horizon
	private List<RepastEdge<Junction>> pathRemainder; // Path that gets agent from first link outside tactical horizon to the end of their destination	
	
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
	 * First tries to return a coordinate from the route coordinates stack, then if that is empty returns the coordinate of the next junction in the path.
	 */
	public Coordinate getTargetCoordinate() {
		if(this.routeCoordinates.size()>0) {
			return routeCoordinates.getFirst();
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
			this.routeCoordinates.removeFirst();
		}
		else {
			updateCurrentJunction();
		}
	}
	
	public void updateCurrentJunction() {
		this.currentJunction = this.getRouteJunctions().pollFirst();
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
		// Get the junction that is connected to junction j across the road, use this as the new starting junction
		for (Junction aj: this.nP.getNet().getAdjacent(j)) {
			NetworkEdge<Junction> e = (NetworkEdge<Junction>)this.nP.getNet().getEdge(j, aj);
			if (e.getRoadLink() != null) {
				this.currentJunction = aj;
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
}