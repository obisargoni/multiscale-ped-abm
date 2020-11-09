package repastInterSim.environment;

import com.vividsolutions.jts.geom.Coordinate;

public interface CrossingAlternative {
	
	public Double distanceTo(Coordinate loc);
	public Coordinate nearestCoord(Coordinate loc);
	public Coordinate farthestCoord(Coordinate loc);
	public Integer getvFlow();
	public String getRoadLinkID();
	public Road getRoad();
	public Coordinate getC1();
	public Coordinate getC2();
	public String getType();
	public Coordinate getDestination();
	public void setDestination(Coordinate d);


}
