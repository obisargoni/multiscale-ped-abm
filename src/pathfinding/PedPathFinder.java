package pathfinding;

import java.util.HashMap;
import java.util.List;

import org.geotools.coverage.grid.GridCoverage2D;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.Geography;
import repastInterSim.agent.GridRoute;
import repastInterSim.agent.Ped;
import repastInterSim.environment.OD;
import repastInterSim.environment.RoadLink;
import repastInterSim.main.GlobalVars;

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
