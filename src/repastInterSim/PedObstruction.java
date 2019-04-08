package repastInterSim;

import com.vividsolutions.jts.geom.Geometry;

public class PedObstruction implements FixedGeography {
	
	private Geometry geom;
    
	public PedObstruction() {
		
	}
	
	public Geometry getGeom() {
		return this.geom;
	}
	
	public void setGeom(Geometry g) {
		this.geom = g;
	}
}
