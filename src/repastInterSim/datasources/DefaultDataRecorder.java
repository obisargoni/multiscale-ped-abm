package repastInterSim.datasources;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

import repastInterSim.environment.GISFunctions;

public class DefaultDataRecorder {
	
	protected Point p;
	
	public DefaultDataRecorder() {
		this.p = GISFunctions.pointGeometryFromCoordinate(new Coordinate(0,0));
	}
	
	public Geometry getGeom() {
		return this.p;
	}

}
