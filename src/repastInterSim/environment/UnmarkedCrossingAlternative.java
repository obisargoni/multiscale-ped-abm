package repastInterSim.environment;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.space.gis.Geography;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class UnmarkedCrossingAlternative extends CrossingAlternative {
	
	private Ped ped;
	
	private String type = "unmarked";
	
	private CrossingCoordsCache ccc;

	public UnmarkedCrossingAlternative() {
		this.ccc = new CrossingCoordsCache();
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
	public double getvFlow() {
		double crossingTime = this.getC1().distance(this.getC2()) / this.ped.getSpeed();
		
		// Loop through vehicles on the road links this crossing covers, count the number that will pass crossing point in crossing time
		double vehicleCount = 0;
		List<RoadLink> itnLinks = this.getCurrentVehicleRoadLinks(); 
		for (int i=0; i<itnLinks.size(); i++){
			RoadLink rl = itnLinks.get(i);
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
	
	/*
	 * Find the nearest coordinate on the opposite side of the road to the pedestrians current position.
	 * 
	 * Opposite side of the road defined as in the direction perpendicular to the bearing of the road link the pedestrian is walking beside
	 * and at the far edge of the carriadgeway from the pedestrian.
	 * 
	 * @param Coordinate c
	 * 		The location of the pedestrian agent
	 * @param Geography<PedObstruction> poG
	 * 		Geography containing the ped obstructions.
	 * 
	 * @returns
	 * 		Coordinate
	 */
	public Coordinate oppositeSideOfRoadCoord(Coordinate c, Geography<PedObstruction> poG) {
		RoadLink orRoadLink = this.getORRoadLink();
		Coordinate rlCent = orRoadLink.getGeom().getCentroid().getCoordinate();
		double oppRoadAngle = GISFunctions.oppositeSideOfRoadAngle(c, orRoadLink, rlCent);
		Coordinate nearestOpCoord = GISFunctions.xSideOfRoadCoord(c, poG, "opposite", orRoadLink, rlCent, oppRoadAngle);
		return nearestOpCoord;
	}
	
	/*
	 * Find the nearest carriageway coordinate on the same side of the road to the pedestrians current position.
	 * 
	 * Same side of the road defined as being in the direction perpendicular to the bearing of the road link the pedestrian is walking beside
	 * and at the near edge of the carriadgeway from the pedestrian.
	 * 
	 * @param Coordinate c
	 * 		The location of the pedestrian agent
	 * @param Geography<PedObstruction> poG
	 * 		Geography containing the ped obstructions.
	 * 
	 * @returns
	 * 		Coordinate
	 */
	public Coordinate sameSideOfRoadCoord(Coordinate c, Geography<PedObstruction> poG) {
		RoadLink orRoadLink = this.getORRoadLink();
		Coordinate rlCent = orRoadLink.getGeom().getCentroid().getCoordinate();
		double oppRoadAngle = GISFunctions.oppositeSideOfRoadAngle(c, orRoadLink, rlCent);
		Coordinate nearestSameCoord = GISFunctions.xSideOfRoadCoord(c, poG, "same", orRoadLink, rlCent, oppRoadAngle);
		return nearestSameCoord;
	}
	
	public Coordinate getC1() {
		Geography<PedObstruction> pedObstructGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.PED_OBSTRUCTION_GEOGRAPHY);
		RoadLink orRoadLink = this.getORRoadLink();
		boolean caChosen = this.ped.getPathFinder().getTacticalPath().getAccumulatorRoute().caChosen();
		
		return ccc.getC1(this.ped.getLoc(), pedObstructGeography, orRoadLink, caChosen);
	}

	public Coordinate getC2() {
		Geography<PedObstruction> pedObstructGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.PED_OBSTRUCTION_GEOGRAPHY);
		RoadLink orRoadLink = this.getORRoadLink();
		boolean caChosen = this.ped.getPathFinder().getTacticalPath().getAccumulatorRoute().caChosen();
		
		return ccc.getC2(this.ped.getLoc(), pedObstructGeography, orRoadLink, caChosen);
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
	
	@Override
	public double getCrossingBearing() {
		return GISFunctions.bearingBetweenCoordinates(this.getC1(), this.getC2());
	}
	
	@Override
	public double getCrossingDistance() {
		return this.getC1().distance(this.getC2());
	}

}

/**
 * Caches the start and end coordinates of the crossing. Cache is recreated when new input coordinate is given (ie the pedestrian's location changes)
 * 
 * @author Obi Thompson Sargoni
 */
class CrossingCoordsCache {
	
	// Used to control when to update the cache
	private boolean caChosen;
	private Coordinate cacheCoord;
	private RoadLink cacheRL;

	private Coordinate[] theCache = new Coordinate[2]; // Records the crossing coordinates
	private Double oppRoadBearing;
	private Coordinate rlCentroid;

	CrossingCoordsCache() {
		
	}

	public void clear() {
		this.theCache[0]=null;
		this.theCache[1]=null;
		this.oppRoadBearing=null;
		this.rlCentroid=null;
	}

	
	private void populateCache(Coordinate c, Geography<PedObstruction> poG, RoadLink orRoadLink, boolean caChosen) {
		
		this.caChosen = caChosen;
		this.cacheCoord=c;
		this.cacheRL=orRoadLink;
		
		this.rlCentroid = orRoadLink.getGeom().getCentroid().getCoordinate();
		this.oppRoadBearing = GISFunctions.oppositeSideOfRoadAngle(c, orRoadLink, this.rlCentroid);
		
		Coordinate nearestOppCoord = GISFunctions.xSideOfRoadCoord(c, poG, "opposite", orRoadLink, this.rlCentroid, this.oppRoadBearing);
		Coordinate nearestSameCoord = GISFunctions.xSideOfRoadCoord(c, poG, "same", orRoadLink, this.rlCentroid, this.oppRoadBearing);
		
		this.theCache[0] = nearestSameCoord;
		this.theCache[1] = nearestOppCoord;
	}

	/**
	 * 
	 * @param c
	 * @return
	 * @throws Exception
	 */
	public Coordinate getC1(Coordinate c, Geography<PedObstruction> poG, RoadLink orRoadLink, boolean caChosen) {
		
		boolean populate = decidePopulate(c, orRoadLink, caChosen);
		
		if (populate) {
			populateCache(c, poG, orRoadLink, caChosen);
		}
		
		return this.theCache[0];
	}
	
	/**
	 * 
	 * @param c
	 * @return
	 * @throws Exception
	 */
	public Coordinate getC2(Coordinate c, Geography<PedObstruction> poG, RoadLink orRoadLink, boolean caChosen) {
		
		boolean populate = decidePopulate(c, orRoadLink, caChosen);
		
		if (populate) {
			populateCache(c, poG, orRoadLink, caChosen);
		}
		
		return this.theCache[1];
	}
	
	private boolean decidePopulate(Coordinate c, RoadLink orRoadLink, boolean caChosen) {
		boolean populate;
		if ( this.cacheCoord==null ) {
			populate=true;
		}
		else if (caChosen==true) {
			populate = false;
		}
		else if(c.equals2D(this.cacheCoord)) {
			populate = false;
		}
		else {
			populate = true;
		}
		return populate;
	}

}
