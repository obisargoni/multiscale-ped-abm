package repastInterSim.environment;

import com.vividsolutions.jts.geom.Coordinate;

public interface CrossingAlternative extends FixedGeography {
	
	public Coordinate getC1();
	public Coordinate getC2();
	public Double distanceTo(Coordinate loc);
	public Integer getvFlow();
	public Road getRoad();
	public void setRoad();
	public String getType();
	public void setType(String type);
	public Coordinate getDestination();
	public void setDestination(Coordinate d);
	public Coordinate nearestCoord(Coordinate loc);
	public Coordinate farthestCoord(Coordinate loc);

}
