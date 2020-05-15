package repastInterSim.pathfinding;

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
	
	private OD origin;
	private OD destination;
	private Geography<Object> geography;
	
	private List<RoadLink> strategicPath;
	private GridRoute tacticalPath = new GridRoute();
	
	private Coordinate nextTacticalPathCoord;

	private Coordinate nextCrossingCoord;
	private RoadLink currentRoadLink;

	public PedPathFinder(Geography<Object> g, OD o, OD d) {
		this.geography = g;
		this.origin = o;
		this.destination = d;
		
		planStrategicPath();
	}
	
	/**
	 * Initialise a new road link routing object that can be used to find a path along a topological road network.
	 * Use this to identify the shortest path through the network and assign this path to this classes' strategic path attribute.
	 * 
	 */
	public void planStrategicPath() {
		// Initialise road network route - needs to ne non-directed for pedestrians! fix this
		RoadNetworkRoute rnr = new RoadNetworkRoute(origin.getGeom().getCoordinate(), destination.getGeom().getCoordinate());
		
		// Find shortest path using road network route
		try {
			rnr.setRoadLinkRoute();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Get path of road links and set this as the strategic path
		this.strategicPath = rnr.getRoadsX();
		
		// The strategic path can be empty when origin and destination are very close
		if(this.strategicPath.isEmpty() == false) {
			this.currentRoadLink = this.strategicPath.get(0);
		}

	}
	
	public void updateTacticalPathCoordinate(HashMap<Integer, Double> gSPM, Coordinate tacticalOriginCoord) {	
		// If reached the end of one section of the route, or if route has just been created, need to produce next set of route coordinates.
		if(this.tacticalPath.getRouteX().size() == 0) {
			// Both use the roadLinkCoordX[0] to set, consider passing in as parameter?
			planTacticaPath(gSPM, tacticalOriginCoord);
		}
		
    	// If crossing coord is null, check the route for upcoming crossing coord
		// Importantly, do this before any coords are removed from the route
    	if(this.nextCrossingCoord == null) {
    		this.nextCrossingCoord = this.tacticalPath.getNextRouteCrossingCoord();
    	}
    	
		this.nextTacticalPathCoord = this.tacticalPath.getRouteX().get(0);
		this.tacticalPath.removeNextFromRoute();
    }
	
	public void planTacticaPath(HashMap<Integer, Double> gSPM, Coordinate tacticalOriginCoord) {
		
		Coordinate tacticalDestCoord;
		if (this.strategicPath.isEmpty() | this.strategicPath.size() == 1) {
			tacticalDestCoord = this.destination.getGeom().getCoordinate();
		}
		else {
			this.currentRoadLink = this.strategicPath.get(0);
			
			// Need to use the network edge associated to the road link to make sure the end of the road link is used as the destination
			tacticalDestCoord = this.currentRoadLink.getEdge().getTarget().getGeom().getCoordinate();
		}
		
		GridCoverage2D grid = this.geography.getCoverage(GlobalVars.CONTEXT_NAMES.BASE_COVERAGE);
		GridRoute tP = new GridRoute(grid, gSPM, tacticalOriginCoord, tacticalDestCoord, true);
		
    	// Get updated set of route coords to follow to next road link coordinate
		tP.setGroupedGridPath();
    	
    	// Adds coordinates to route from next section of grid path
		tP.setNextRouteSection();
		
		this.tacticalPath= tP;
		
		// Remove this road link from the strategic path, no longer the next road link
		if (this.strategicPath.isEmpty() == false) {
			this.strategicPath.remove(0);
		}
	}
	
	
	public List<RoadLink> getStrategicPath() {
		return this.strategicPath;
	}
	
	public GridRoute getTacticalPath() {
		return this.tacticalPath;
	}

	public Coordinate getNextTacticalPathCoord() {
		return nextTacticalPathCoord;
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
	
	public RoadLink getCurrentRoadLink() {
		return this.currentRoadLink;
	}

}
