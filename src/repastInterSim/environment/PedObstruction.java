package repastInterSim.environment;

import com.vividsolutions.jts.geom.Geometry;

public class PedObstruction implements FixedGeography {
	
	private Geometry geom;
	private String priority = ""; // Priority information comes from GIS data
    
	public PedObstruction() {
		
	}
	
	public Geometry getGeom() {
		return this.geom;
	}
	
	public void setGeom(Geometry g) {
		this.geom = g;
	}
	
	public String getPriority() {
		return this.priority;
	}
	
	public void setPriority(String p) {
		this.priority = p;
	}
}
