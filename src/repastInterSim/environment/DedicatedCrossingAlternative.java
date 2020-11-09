package repastInterSim.environment;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;

import repast.simphony.space.gis.Geography;
import repastInterSim.agent.Ped;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
import repastInterSim.pathfinding.RoadNetworkRoute;

public class DedicatedCrossingAlternative implements CrossingAlternative {
	
	// Coordinates at which the ca meets pavement
	private Geometry geom = null;
	private Coordinate c1 = null;
	private Coordinate c2 = null;
	
	// Destination coordinate. This is the destination this crossing alternative leads to
	private Coordinate destination = null;
	
	// Default type is unmarked
	private String type = "unmarked";
	
	// id of the road link this crossing is located on
	private String roadLinkID;
	private Road road;
	
	// pedestrian agent perceiving this crossing
	private Ped ped = null;
	
	// Road geography containing the pavement and carriageway polygons
	private Geography<Road> roadGeography = null;
	private List<RoadLink> strategicPathsection;

	public DedicatedCrossingAlternative(){
				
	}
	
	/*
	 * Calculate the distance to the input coordinate. Distance return is the distance to the nearest
	 * crossing alternative coordiante
	 * 
	 * @param Coordinate loc
	 * 		The coordinate to calculate the distance from
	 */
	public Double distanceTo(Coordinate loc) {
		double d1 = getC1().distance(loc);
		double d2 = getC2().distance(loc);
		return Math.min(d1, d2);
	}
	
	public Coordinate nearestCoord(Coordinate loc) {
		Coordinate[] coords = {getC1(), getC2()};
		
		Coordinate cNear = Arrays.stream(coords).min(Comparator.comparing(c->c.distance(loc))).get();
		
		return cNear;
	}

	public Coordinate farthestCoord(Coordinate loc) {
		Coordinate[] coords = {getC1(), getC2()};
		
		Coordinate cFar = Arrays.stream(coords).max(Comparator.comparing(c->c.distance(loc))).get();
		
		return cFar;
	}	
	
	/*
	 * Vehicles are assumed to yield at dedicated crossings resulting in zero vehicle flow
	 */
	public Integer getvFlow() {
		return 0;
	}
	
	public void setRoadLinkID(String rlID) {
		this.roadLinkID = rlID;
	}
	
	public String getRoadLinkID() {
		if (this.roadLinkID == null) {
			return this.getRoad().getRoadLinkID();
		}
		else {
			return this.roadLinkID;
		}
	}
	
	public Road getRoad() {
		if (this.road == null) {
			return this.ped.getCurrentRoad();
		}
		else {
			return this.road;
		}
	}

	public void setRoad() {
		this.road = this.ped.getCurrentRoad();
	}

	/*
	 * Get the nearest coordinate on the pavement on the opposite side of the road to the input coordinate.
	 * Used to identify the end point of unmarked crossing alternatives.
	 */
	public Coordinate nearestOppositePedestrianCoord(Coordinate c, String roadLinkID, Geography<Road> rG, List<RoadLink> sps) {
		
		// Get pedestrian road objects on this road link
		List<Road> caPedRoads = null;
		try {
			caPedRoads = RoadNetworkRoute.getRoadLinkPedestrianRoads(roadLinkID, rG);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Loop through ped roads and find nearest coordinate on each
		Double minDist = Double.MAX_VALUE;
		Coordinate nearestOpCoord = null;
		for (Road rd:caPedRoads) {
			// Get nearest coord
			Coordinate nearC = GISFunctions.xestGeomCoordinate(c, rd.getGeom(), false);
			
			// Check parity to make sure coordinate is on opposite side of the road
			int p = RoadNetworkRoute.calculateRouteParity(c, nearC, sps);
			if (p==0) {
				continue;
			}
			
			// Check if nearer
			double d = c.distance(nearC);
			if (d < minDist) {
				minDist = d;
				nearestOpCoord = nearC;
			}
		}
		
		return nearestOpCoord;
	}

	public Coordinate getC1() {
		if (c1==null) {
			return this.ped.getLoc();
		}
		else {
			return c1;
		}
	}

	public void setC1(Coordinate c1) {		
		this.c1 = c1;
	}

	public Coordinate getC2() {
		if(c2==null) {
			return nearestOppositePedestrianCoord(this.ped.getLoc(), this.getRoadLinkID(), this.roadGeography, this.strategicPathsection);
		}
		else {
			return c2;
		}
	}

	public void setC2(Coordinate c2) {
		this.c2 = c2;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public Coordinate getDestination() {
		return destination;
	}

	public void setDestination(Coordinate destination) {
		this.destination = destination;
	}

	public Ped getPed() {
		return ped;
	}

	public void setPed(Ped ped) {
		this.ped = ped;
	}

	public Geography<Road> getRoadGeography() {
		return roadGeography;
	}

	public void setRoadGeography(Geography<Road> roadGeography) {
		this.roadGeography = roadGeography;
	}

	@Override
	public Geometry getGeom() {
		return this.geom;
	}

	@Override
	public void setGeom(Geometry g) {
		this.geom = g;
		
		int ncoords = g.getCoordinates().length;
		this.c1 = g.getCoordinates()[0];
		this.c2 = g.getCoordinates()[ncoords-1];
	}

	public void setStrategicPathSection(List<RoadLink> sps) {
		this.strategicPathsection = sps;		
	}

}
