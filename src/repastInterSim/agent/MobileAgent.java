package repastInterSim.agent;

import java.util.HashMap;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.Geography;
import repastInterSim.environment.Destination;
import repastInterSim.environment.Route;

public interface MobileAgent {
	
	Destination getDestination();
	void setLoc();
	Coordinate getLoc();
	Route getRoute();
	HashMap<Integer, Integer> getGridPrioritySummandMap();
	Geography<Object> getGeography();

}
