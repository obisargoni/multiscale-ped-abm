package repastInterSim.environment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;

import repast.simphony.space.gis.Geography;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.main.SpaceBuilder;
import repastInterSim.pathfinding.RoadNetworkRoute;

public class UnmarkedCrossingAlternative extends CrossingAlternative {
	
	private Ped ped;
	
	private String type = "unmarked";

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

	/*
	 * Get the number of vehicles on the road link. 
	 * Ideally will calculate exactly the number of cars that would pass through the crossing in a given time period
	 */
	@Override
	public Integer getvFlow() {
		double crossingTime = this.getC1().distance(this.getC2()) / this.ped.getSpeed();
		
		// Loop through vehicles on the road links this crossing covers, count the number that will pass crossing point in crossing time
		int vehicleCount = 0;
		for (int i=0; i<this.getRoad().getRoadLinks().size(); i++){
			RoadLink rl = this.getRoad().getRoadLinks().get(i);
			for(int j = 0; j<rl.getQueue().count(); j++){
				int vi = rl.getQueue().readPos() + j;
				if (vi>=rl.getQueue().capacity()) {
					vi = vi-rl.getQueue().capacity();
				}
				Vehicle v = rl.getQueue().elements[vi];
				
				// Check if crossing is in front of vehicle, if not continue to next vehicle
				if (!GISFunctions.coordInFront(v.getLoc(), v.getBearing(), getC1())) {
					continue;
				}
				
				// Get expected future location of vehicle
				double travelDist = v.getSpeed() * crossingTime;
				Coordinate futureLoc = new Coordinate(v.getLoc().x + Math.sin(v.getBearing())*travelDist, v.getLoc().y + Math.cos(v.getBearing())*travelDist);
				
				// Check if crossing location is in front of future location or not
				boolean crossingAhead = GISFunctions.coordInFront(futureLoc, v.getBearing(), getC1());
							
				if (!crossingAhead){
					vehicleCount++;
				}
			}
		}
		return vehicleCount;
	}
	
	public String getRoadLinkID() {
		return this.ped.getPathFinder().getStrategicPath().get(0).getFID();
	}
	
	public Road getRoad() {
		return this.ped.getCurrentRoad();
	}
	
	/*
	 * Get the nearest coordinate on the pavement on the opposite side of the road to the input coordinate.
	 * Used to identify the end point of unmarked crossing alternatives.
	 * 
	 * @param Coordinate c
	 * 		The location of the pedestrian agent
	 * @param String roadLinkID
	 * 		The ID of the road link the pedestrian agent is currently walking beside
	 * @param Geography<Road> rG
	 * 		The Road geography. Contains pavement and carriageway polygons objects
	 * @param Geography<PedObstruction> poG
	 * 		Geography containing the ped obstructions.
	 * @param List<RoadLink> fsp
	 * 		The strategic path road links. Contains all strategic path links, even those the pedestrian agent has passed.
	 * 
	 * @returns
	 * 		Coordinate
	 */
	public Coordinate nearestOppositePedestrianCoord(Coordinate c, String roadLinkID, Geography<Road> rG, Geography<PedObstruction> poG, List<RoadLink> fsp) {
		
		// Get pedestrian roads via attributes of road link and road objects
		
		// Identify geometries that demark the edge of the road. First try pedestrian pavement polygons. Failing that try obstruction geometries.
		List<Geometry> roadEdgeGeoms = new ArrayList<Geometry>();
		List<Road> caPedRoads = this.getRoad().getORRoadLink().getRoads().stream().filter(r -> r.getPriority().contentEquals("pedestrian")).collect(Collectors.toList());
		if (caPedRoads.size()==0) {
			Geometry nearby = GISFunctions.pointGeometryFromCoordinate(c).buffer(30);
			roadEdgeGeoms = SpatialIndexManager.searchGeoms(poG, nearby);
		}
		else {
			for (Road rd: caPedRoads) {
				roadEdgeGeoms.add(rd.getGeom());
			}
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
			for (Geometry g: roadEdgeGeoms) {

				// First check whether both rays have been tried
				// If not all rays tried use ray intersection method
				// Otherwise use nearest coord method
				Coordinate[] intersectingCoords = null;
				if (iRay<rays.length) {
					LineString ray = rays[iRay];
					// Get intersection between ray and this polygon
					intersectingCoords = g.intersection(ray).getCoordinates();
				}
				else {
					Coordinate nearC = GISFunctions.xestGeomCoordinate(c, g, false);
					intersectingCoords = new Coordinate[1];
					intersectingCoords[0] = nearC;
				}

				// Now loop through these intersecting coords to find one that is
				// - on the opposite side of the road
				// - nearest to ped
				for (int j=0; j<intersectingCoords.length; j++) {
					
					Coordinate intC = intersectingCoords[j];
					
					// Check parity to make sure coordinate is on opposite side of the road
					int p = RoadNetworkRoute.calculateRouteParity(c, intC, fsp);
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
			
			// Break loop if agent can't find opposite coord
			if (iRay>10) {
				break;
			}
		}
		
		return nearestOpCoord;
	}
	
	/*
	 * Find the coordiante on the opposite side of the road to the pedestrians current position.
	 * 
	 * Opposite side of the road defined as in the direction perpendicular to the bearing of the road link the pedestrian is walking beside
	 * and at the far edge of the carriadgeway from the pedestrian.
	 * 
	 * @param Coordinate c
	 * 		The location of the pedestrian agent
	 * @param RoadLink roadLink
	 * 		The road link the pedestrian agent is currently walking beside
	 * @param Geography<Road> rG
	 * 		The Road geography. Contains pavement and carriageway polygons objects
	 * @param Geography<PedObstruction> poG
	 * 		Geography containing the ped obstructions.
	 * @param List<RoadLink> fsp
	 * 		The strategic path road links. Contains all strategic path links, even those the pedestrian agent has passed.
	 * 
	 * @returns
	 * 		Coordinate
	 */
	public Coordinate oppositeSideOfRoadCoord(Coordinate c, RoadLink roadLink, Geography<Road> rG, Geography<PedObstruction> poG, List<RoadLink> fsp) {
		
		// Opposite side of the road is in direction perpendicular to road link. Find the bearing to the opposite side of the road
		Coordinate[] rlCoords = roadLink.getGeom().getCoordinates(); 
		double rlBearing = GISFunctions.bearingBetweenCoordinates(rlCoords[0], rlCoords[rlCoords.length-1]);
		double perp1 = rlBearing - Math.PI / 2;
		double perp2 = rlBearing + Math.PI / 2;
		
		// Find which of these bearings points to opp side of road to ped
		Coordinate rlCent = roadLink.getGeom().getCentroid().getCoordinate();
		double rlToPedBearing = GISFunctions.bearingBetweenCoordinates(rlCent, c);
		
		double range1 = Vector.acuteRangeBetweenAngles(rlToPedBearing, perp1);
		double range2 = Vector.acuteRangeBetweenAngles(rlToPedBearing, perp2);
		
		// Bearing to opposite side will be more than 90 deg from bearing to ped
		double oppRoadAngle;
		if (range1>Math.PI/2) {
			oppRoadAngle = perp1;
		} 
		else {
			oppRoadAngle = perp2;
		}
		
		// Identify geometries that demark the edge of the road. Try carriadgeway polygons. Failing that try obstruction geometries.
		List<Geometry> roadEdgeGeoms = new ArrayList<Geometry>();
		List<Road> caPedRoads = this.getRoad().getORRoadLink().getRoads().stream().filter(r -> r.getPriority().contentEquals("vehicle")).collect(Collectors.toList());
		if (caPedRoads.size()==0) {
			Geometry nearby = GISFunctions.pointGeometryFromCoordinate(c).buffer(30);
			roadEdgeGeoms = SpatialIndexManager.searchGeoms(poG, nearby);
		}
		else {
			for (Road rd: caPedRoads) {
				roadEdgeGeoms.add(rd.getGeom());
			}
		}
		
		// Get ray that points in direction of opposite side of the road from ped
		LineString ray = GISFunctions.linestringRay(c, oppRoadAngle, 30.0);
				
		// Loop through ped roads and find nearest coordinate on each
		Double maxDist = Double.MIN_VALUE;
		Coordinate nearestOpCoord = null;
		while (nearestOpCoord==null) {		
			for (Geometry g: roadEdgeGeoms) {

				// Intersect ray with geometry that demarks where the other side of the road is
				Coordinate[] intersectingCoords = null;
				intersectingCoords = g.intersection(ray).getCoordinates();

				// Now loop through these intersecting coords to find one that is
				// - on the opposite side of the road
				// - nearest to ped
				for (int j=0; j<intersectingCoords.length; j++) {
					
					Coordinate intC = intersectingCoords[j];
					
					// Check if nearer
					double d = c.distance(intC);
					if (d > maxDist) {
						maxDist = d;
						nearestOpCoord = intC;
					}
				}
			}
			
			// If unable to find an opp road coord, increase length of ray.
			if ( (nearestOpCoord==null) & (ray.getLength()<50.0) ) {
				ray = GISFunctions.linestringRay(c, oppRoadAngle, 50.1);
			}
			else {
				break;
			}
		}
		
		return nearestOpCoord;
	}
	
	public Coordinate getC1() {
		return this.ped.getLoc();
	}

	public Coordinate getC2() {
		List<RoadLink> fullStrategicPath = this.ped.getPathFinder().getFullStrategicPath();
		RoadLink currentRL = this.ped.getPathFinder().getStrategicPath().get(0);
		return oppositeSideOfRoadCoord(this.ped.getLoc(), currentRL, SpaceBuilder.roadGeography, SpaceBuilder.pedObstructGeography, fullStrategicPath);
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
