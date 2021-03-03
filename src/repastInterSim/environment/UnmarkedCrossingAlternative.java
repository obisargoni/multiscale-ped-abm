package repastInterSim.environment;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;

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
		
		// Get rays perpendicular to agent's bearing to find crossing coordinate
		LineString ray1 = GISFunctions.linestringRay(c, this.ped.getBearing() - (Math.PI/2), 50.0);
		LineString ray2 = GISFunctions.linestringRay(c, this.ped.getBearing() + (Math.PI/2), 50.0);
		
		LineString[] rays = {ray1, ray2};
		
		// Loop through ped roads and find nearest coordinate on each
		Double minDist = Double.MAX_VALUE;
		Coordinate nearestOpCoord = null;
		int iRay = 0;
		boolean oppRoadSideUnknown = true;
		while (oppRoadSideUnknown) {		
			for (Road rd:caPedRoads) {

				// First check whether both rays have been tried
				// If not all rays tried use ray intersection method
				// Otherwise use nearest coord method
				Coordinate[] intersectingCoords = null;
				if (iRay<rays.length) {
					LineString ray = rays[iRay];
					// Get intersection between ray and this polygon
					intersectingCoords = rd.getGeom().intersection(ray).getCoordinates();
				}
				else {
					Coordinate nearC = GISFunctions.xestGeomCoordinate(c, rd.getGeom(), false);
					intersectingCoords = new Coordinate[1];
					intersectingCoords[0] = nearC;
				}

				// Now loop through these intersecting coords to find one that is
				// - on the opposite side of the road
				// - nearest to ped
				for (int j=0; j<intersectingCoords.length; j++) {
					
					Coordinate intC = intersectingCoords[j];
					
					// Check parity to make sure coordinate is on opposite side of the road
					int p = RoadNetworkRoute.calculateRouteParity(c, intC, sps);
					if (p==0) {
						continue;
					}
					else {
						oppRoadSideUnknown = false;
						
						// Check if nearer
						double d = c.distance(intC);
						if (d < minDist) {
							minDist = d;
							nearestOpCoord = intC;
						}
					}
				}
			}
			iRay++;
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
