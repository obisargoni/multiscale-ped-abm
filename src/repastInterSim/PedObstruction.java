package repastInterSim;

import com.vividsolutions.jts.geom.Coordinate;

public class PedObstruction implements FixedGeography {
	
	private Coordinate coord;
	
	public PedObstruction() {
		
	}
	
	public Coordinate getCoords() {
		return this.coord;
	}
	
	public void setCoords(Coordinate c) {
		this.coord = c;
	}

}
