package pathfinding;

import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.environment.Junction;

public class RoadNetworkRoute {
	
	protected Coordinate origin;
	protected Coordinate destination;
	
	public RoadNetworkRoute(Coordinate origin, Coordinate destination) {
		
		this.origin = origin;
		this.destination = destination;
		
	}
	
	public List<RepastEdge<Junction>> getShortestRoute(Network<Junction> net, Iterable<Junction> currentJunctions, Iterable<Junction> destJunctions,
			Junction[] routeEndpoints) throws Exception {

		List<RepastEdge<Junction>> shortestPath = null;
		return shortestPath;
	}

}
