package repastInterSim.environment;

import com.vividsolutions.jts.geom.Coordinate;

public class CrossingAlternative {
	
	// Coordinates at which the ca meets pavement
	private Coordinate c1 = null;
	private Coordinate c2 = null;
	
	// Default type is unmarked
	private String type = "unmarked";

	public CrossingAlternative(){
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

}
