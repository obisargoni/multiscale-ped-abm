package repastInterSim.agent;

import com.vividsolutions.jts.geom.Coordinate;

import repastInterSim.environment.Destination;

public interface mobileAgent {
	
	Destination getDestination();
	void setLoc();
	Coordinate getLoc();

}
