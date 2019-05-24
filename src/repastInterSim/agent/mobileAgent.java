package repastInterSim.agent;

import com.vividsolutions.jts.geom.Coordinate;

import repastInterSim.environment.Destination;
import repastInterSim.environment.Route;

public interface mobileAgent {
	
	Destination getDestination();
	void setLoc();
	Coordinate getLoc();
	Route getRoute();

}
