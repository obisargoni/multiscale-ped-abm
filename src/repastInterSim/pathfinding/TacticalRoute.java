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

public class TacticalRoute {
	
	private NetworkPath<Junction> nP;
	private Junction currentJunction;
	private Junction endJunction;
	private List<CrossingAlternative> crossingAlternatives;
	private LinkedList<Junction> routeJunctions = null;
	private LinkedList<Coordinate> routeCoordinates = new LinkedList<Coordinate>();; 
	private List<RepastEdge<Junction>> pathToEnd; // Path that gets agent from start of tactical horizon to end of tactical horizon
	private List<RepastEdge<Junction>> pathEndToOutside; // Path that gets agent from end of tactical horizon to first link outside of tactical horizon
	private List<RepastEdge<Junction>> pathRemainder; // Path that gets agent from first link outside tactical horizon to the end of their destination
	private boolean routeCompleted;
	
	
	public TacticalRoute(NetworkPath<Junction> nP, Junction startJunction, Junction endJunction) {
		this.nP = nP;
		this.currentJunction = startJunction;
		this.endJunction = endJunction;
		this.pathToEnd = new ArrayList<RepastEdge<Junction>>();
		this.pathEndToOutside = new ArrayList<RepastEdge<Junction>>();
		this.pathRemainder = new ArrayList<RepastEdge<Junction>>();
		this.routeCompleted = false;
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
	
	public void setCurrentJunction(Junction j) {
		this.currentJunction = j;
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
	}
	}
	}
	}
	}

	}

	public void setRouteRemainderPath(List<RepastEdge<Junction>> path) {
		this.remainderPath = path;
	}
	
	public List<RepastEdge<Junction>> getRouteRemainderPath() {
		return this.remainderPath;
	}
	
	
}