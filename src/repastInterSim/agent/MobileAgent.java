package repastInterSim.agent;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.Geography;
import repastInterSim.environment.Destination;

public interface MobileAgent {
	
	Destination getDestination();
	void setLoc();
	Coordinate getLoc();
	Geography<Object> getGeography();
	void tidyForRemoval();

}
