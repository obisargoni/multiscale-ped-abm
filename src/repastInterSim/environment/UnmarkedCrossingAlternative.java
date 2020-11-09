package repastInterSim.environment;

import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.space.gis.Geography;
import repastInterSim.agent.Ped;

public class UnmarkedCrossingAlternative implements CrossingAlternative {
	
	private Ped ped;
	
	// Road geography containing the pavement and carriageway polygons
	private Geography<Road> roadGeography = null;
	private List<RoadLink> strategicPathsection;


	public UnmarkedCrossingAlternative() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public Geometry getGeom() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setGeom(Geometry c) {
		// TODO Auto-generated method stub

	}

	@Override
	public Coordinate getC1() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Coordinate getC2() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Double distanceTo(Coordinate loc) {
		// TODO Auto-generated method stub
		return null;
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

	@Override
	public Road getRoad() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setRoad() {
		// TODO Auto-generated method stub

	}

	@Override
	public String getType() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setType(String type) {
		// TODO Auto-generated method stub

	}

	@Override
	public Coordinate getDestination() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setDestination(Coordinate d) {
		// TODO Auto-generated method stub

	}

	@Override
	public Coordinate nearestCoord(Coordinate loc) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Coordinate farthestCoord(Coordinate loc) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void setPed(Ped p) {
		this.ped = p;
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

}
