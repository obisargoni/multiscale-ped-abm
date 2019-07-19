package repastInterSim.environment;

import java.util.ArrayList;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.PrecisionModel;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;

public class RoadLink implements FixedGeography {
	
	private Geometry geom;
	private List<Junction> junctions; // The Roads connected to this Junction, used in GIS road network
	private NetworkEdge<Junction> edge;
	private String direction = null;
	
	/**
	 * The null road represents Road objects that do not actually exist, preventing NullPointerExceptions. This is
	 * necessary for routes that include transport networks as these wont necessarily have a Road object associated with
	 * them (e.g. train lines).
	 */
	public static RoadLink nullRoad;
	static {
		RoadLink.nullRoad = new RoadLink();
		Coordinate[] c = {new Coordinate(), new Coordinate()};
		CoordinateArraySequence cs = new CoordinateArraySequence(c);
		RoadLink.nullRoad.setGeom(new LineString(cs, new GeometryFactory()));
	}

	
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
	
	public List<Junction> getJunctions(){
		return this.junctions;
	}
	
	public void setEdge(NetworkEdge<Junction> edge2) {
		this.edge = edge2;
	}
	
	public NetworkEdge<Junction> getEdge(){
		return this.edge;
	}
	
	
	public String getDirection() {
		return this.direction;
	}
	
	public void setDirection(String dir) {
		
		String d = (String)dir;
		
		// Only allow "-" or "+" as orientations
		if (d.equals("-")) {
			this.direction = "-";
		}
		
		else if (d.equals("+")) {
			this.direction = "-";
		}
		
	}

}
