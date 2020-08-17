package repastInterSim.environment;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.Geography;

public class CrossingAlternative {
	
	// Coordinates at which the ca meets pavement
	private Coordinate c1 = null;
	private Coordinate c2 = null;
	
	// Destination coordinate. This is the destination this crossing alternative leads to
	private Coordinate destination = null;
	
	// Default type is unmarked
	private String type = "unmarked";
	
	// id of the road link this crossing is located on
	private RoadLink roadLink;
	

	public CrossingAlternative(){
				
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
	
	public Integer getvFlow(Geography<Road> rG) {
		// Get the number of vehicles on the road link
		Road r = null;
    	for (Road ri: rG.getAllObjects()) {
    		if (ri.getRoadLinkID().contentEquals(getRoadLinkID())){
    			r = ri;
    		}
    	}
    	int vehicleNumber = r.getRoadLinksVehicleCount();
		return vehicleNumber;
	}
	
	public String getRoadLinkID() {
		return this.roadLink.getFID();
	}
	
	public RoadLink getRoadLink() {
		return this.roadLink;
	}

	public void setRoadLink(RoadLink rL) {
		this.roadLink = rL;
	}


	public Coordinate getC1() {
		return c1;
	}

	public void setC1(Coordinate c1) {
		this.c1 = c1;
	}

	public Coordinate getC2() {
		return c2;
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

}
