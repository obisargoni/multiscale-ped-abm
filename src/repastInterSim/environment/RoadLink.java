package repastInterSim.environment;

import java.util.ArrayList;
import java.util.List;

import com.vividsolutions.jts.geom.Geometry;

public class RoadLink implements FixedGeography {
	
	private Geometry geom;
	private List<Junction> junctions; // The Roads connected to this Junction, used in GIS road network
	private NetworkEdge<Junction> edge;

	
	public RoadLink() {
		this.junctions = new ArrayList<Junction>();
		this.edge = null;
	}

	@Override
	public Geometry getGeom() {
		return this.geom;
	}

	@Override
	public void setGeom(Geometry g) {
		this.geom = g;
	}
	
	public void addJunction(Junction j) {
		this.junctions.add(j);
	}
	
	public void setEdge(NetworkEdge<Junction> e) {
		this.edge = e;
	}
	
	public NetworkEdge<Junction> getEdge(){
		return this.edge;
	}
	
	
	

}
