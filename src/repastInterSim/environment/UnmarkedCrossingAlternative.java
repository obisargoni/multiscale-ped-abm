package repastInterSim.environment;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.space.gis.Geography;
import repastInterSim.agent.Ped;
import repastInterSim.pathfinding.RoadNetworkRoute;

public class UnmarkedCrossingAlternative extends CrossingAlternative {
	
	private Ped ped;
	
	private String type = "unmarked";
	private Coordinate destination;
		
	// Road geography containing the pavement and carriageway polygons
	private Geography<Road> roadGeography = null;
	private List<RoadLink> strategicPathsection;


	public UnmarkedCrossingAlternative() {
		// TODO Auto-generated constructor stub
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

	@Override
	/*
	 * Get the number of vehicles on the road link. 
	 * Ideally will calculate exactly the number of cars that would pass through the crossing in a given time period
	 */
	public Integer getvFlow() {
		int vehicleNumber = this.getRoad().getRoadLinksVehicleCount();
		return vehicleNumber;
	}
	
	public String getRoadLinkID() {
		return this.getRoad().getRoadLinkID();
	}
	
	public Road getRoad() {
		return this.ped.getCurrentRoad();
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
		return this.ped.getLoc();
	}

	public Coordinate getC2() {
		return nearestOppositePedestrianCoord(this.ped.getLoc(), this.getRoadLinkID(), this.roadGeography, this.strategicPathsection);
	}
	
	public String getType() {
		return type;
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

	public void setStrategicPathSection(List<RoadLink> rls) {
		this.strategicPathsection = rls;
	}

	@Override
	/*
	 * Unmarked crossing doesn't have a geometry
	 */
	public Geometry getGeom() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	/*
	 * Unmarked crossing doesn't have a geometry
	 */
	public void setGeom(Geometry c) {
		// TODO Auto-generated method stub
		
	}

}
