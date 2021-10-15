package repastInterSim.environment;

import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.context.Context;
import repast.simphony.space.gis.Geography;

public class OD implements FixedGeography{
	
	private Geometry geom;
	private Long id;
	private String fid;
		
	public OD() {

	}

	@Override
	public Geometry getGeom() {
		return this.geom;
	}

	@Override
	public void setGeom(Geometry g) {
		this.geom = g;
	}

	public Long getId() {
		return id;
	}

	public void setId(Long id) {
		this.id = id;
	}
	
	public String getFID() {
		return this.fid;
	}
	
	public void setFID(String s) {
		this.fid = s;
	}
	
	public void clear() {
		this.geom=null;
		this.fid=null;
		this.id=null;
		
	}
}
