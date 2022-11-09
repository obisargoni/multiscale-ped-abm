package repastInterSim.environment;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.apache.velocity.runtime.RuntimeInstance;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;

import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.environment.RunState;
import repast.simphony.engine.schedule.ScheduledMethod;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.main.GlobalVars;
import repastInterSim.util.RingBufferFillCount;

public class RoadLink implements FixedGeography, Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private Geometry geom;
	private String priority = "";
	private List<Junction> junctions; // The Roads connected to this Junction, used in GIS road network
	private NetworkEdge<Junction> edge;
	private String fid = null;
	private String pedRLID = null; // ID of the strategic road network link this object corresponds to or cuts across in the case of pavement links
	private String pedRoadID = null; // The ID of the pedestrian pavement polygon this link corresponds to in the case of pavement links
	private String direction = null;
	private String MNodeFID = null;
	private String PNodeFID = null;
	private String infCross = "true"; // Controls whether or not informal crossing is permitted on this road link
	private RingBufferFillCount<Vehicle> queue;
	private List<Ped> peds = new ArrayList<Ped>();
	private List<Road> roads = new ArrayList<Road>();
	private double averageVehicleCount=0;
	
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
	
	private void initQueue(double roadLength) {
		int capacity = Math.max(1, (int) (roadLength / GlobalVars.vehicleLength)); 
		Vehicle[] queueArr = new Vehicle[capacity];
		this.queue = new RingBufferFillCount<Vehicle>(queueArr);
	}
	
	@ScheduledMethod(start = 1, interval = 1, priority = 3)
	public void updateAverageVehicleCount() {
		int vehCount = this.queue.count();
		this.averageVehicleCount+= (vehCount - this.averageVehicleCount) / RunEnvironment.getInstance().getCurrentSchedule().getTickCount();
	}

	@Override
	public Geometry getGeom() {
		return this.geom;
	}

	@Override
	public void setGeom(Geometry g) {
		this.geom = g;
		initQueue(this.geom.getLength());
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
			this.direction = "+";
		}
		
	}
	
	public String getMNodeFID() {
		return this.MNodeFID;
	}
	
	public void setMNodeFID(String mnfid) {
		this.MNodeFID = mnfid;
	}
	
	public String getPNodeFID() {
		return this.PNodeFID;
	}
	
	public void setPNodeFID(String pnfid) {
		this.PNodeFID = pnfid;
	}
	
	public String getFID() {
		return this.fid;
	}
	
	public void setFID(String fid) {
		this.fid = fid;
	}
	
	public int getVehicleCount() {
		return this.queue.count();
	}
	
	public Integer addVehicleToQueue(Vehicle v) {
		return this.queue.put(v);
	}
	
	public Vehicle removeVehicleFromQueue() {
		return this.queue.take();
	}
	
	public Vehicle removeVehicleFromQueue(int i) {
		return this.queue.take(i);
	}
	
	public String getPriority() {
		return this.priority;
	}
	
	public void setPriority(String pri) {
		this.priority = pri;
	}

	public String getPedRLID() {
		return pedRLID;
	}

	public void setPedRLID(String pedRoadLinkID) {
		this.pedRLID = pedRoadLinkID;
	}
	

	public String getPedRoadID() {
		return pedRoadID;
	}

	public void setPedRoadID(String pedRoadID) {
		this.pedRoadID = pedRoadID;
	}
	
	public RingBufferFillCount<Vehicle> getQueue() {
		return this.queue;
	}
	
	public List<Ped> getPeds() {
		return this.peds;
	}
	
	public List<Road> getRoads(){
		return this.roads;
	}
	
	public double getAverageVehicleCount() {
		return this.averageVehicleCount;
	}
	
	public String infCross() {
		return this.infCross;
	}
	
	public void setInfCross(String ic) {
		this.infCross = ic;
	}
	
	public void clear() {
		this.geom = null;
		this.priority = null;
		this.junctions = null;
		
		this.edge.clear();
		this.edge = null;
		
		this.fid = null;
		this.pedRLID = null;
		this.pedRoadID = null;
		this.direction = null;
		this.MNodeFID = null;
		this.PNodeFID = null;
		this.queue = null;
		this.peds = null;
		this.roads = null;
	}

}
