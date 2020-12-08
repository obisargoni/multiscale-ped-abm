package repastInterSim.environment;

import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.context.Context;
import repast.simphony.space.gis.Geography;

public class OD implements FixedGeography{
	
	private Context<Object> context;
	private Geography<Object> geography;
	private Geography<OD> destinationGeography;
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
	
	public void setRootGeography(Geography<Object> G) {
		this.geography = G;
	}
	
	public void setDestinationGeography(Geography<OD> G) {
		this.destinationGeography = G;
	}
	
	public void setRootContext(Context<Object> C) {
		this.context = C;
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
}
