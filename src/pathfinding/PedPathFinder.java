package pathfinding;

import repastInterSim.environment.OD;

public class PedPathFinder {
	
	private RoadNetworkRoute strategicPath;
	
	private OD origin;
	private OD destination;
	
	public PedPathFinder(OD o, OD d) {
		this.origin = o;
		this.destination = d;
		
		this.strategicPath = new RoadNetworkRoute(origin.getGeom().getCoordinate(), destination.getGeom().getCoordinate());
	}
	
	public RoadNetworkRoute getStrategicPath() {
		return this.strategicPath;
	}

}
