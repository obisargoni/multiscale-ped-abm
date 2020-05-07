/*
 * This class draws heavily on Nick Malleson's Road.java class
 * from the RepastCity project
 */

package repastInterSim.environment;

import java.util.ArrayList;
import java.util.List;

import com.vividsolutions.jts.geom.Geometry;

import repastInterSim.main.SpaceBuilder;

/*
 * Not sure what attributes and methods are needed for the road object.
 */
public class Road implements FixedGeography {
	
	private String roadLinkID; // The id of the road link this road object should be linked with
	private Geometry geom;	
	private String priority = ""; // Priority information comes from GIS data
	
	// Allows road agents to be joined with a particular RoadLink agent
	private List<RoadLink> roadLinks = new ArrayList<RoadLink>();
	
	
	/*
	 * Instance method
	 */
	public Road() {
	}
	
	public Road(String pri) {
		this.priority = pri;
	}
	
	public int getRoadLinksVehicleCount() {
		int count = 0;
		for(RoadLink rl:this.roadLinks) {
			count += rl.getVehicleCount();
		}
		return count;
	}
	
	/**
	 * Estimate number of leader vehicles by increasing count by one for each road link with
	 * a non-zero vehicle count. If a road link has a non-zero vehicle count then it has a leader vehicle
	 * @return
	 * 		Count of leader vehicles
	 */
	public int getNumberLeadVehicles() {
		int leaderCount = 0;
		for(RoadLink rl: this.roadLinks) {
			if(rl.getVehicleCount() > 0) {
				leaderCount++;
			}
		}
		return leaderCount;
	}
	
	public List<RoadLink> getRoadLinks() {
		return this.roadLinks;
	}
	
	public void setPriority(String pri) {
		this.priority = pri;
	}
	
	public String getPriority() {
		return this.priority;
	}
	
	public Geometry getGeom() {
		return this.geom;
	}
	
	public void setGeom(Geometry g) {
		this.geom = g;
	}
	
	public void setRoadLinkID(String rlFID) {
		this.roadLinkID = rlFID;
		
		// Once road link FID is set, can set the Road Link object
		setRoadLinks();
	}
	
	public String getRoadLinkID() {
		return this.roadLinkID;
	}
	
	/**
	 * To model two way roads, two RoadLink objects are created with the same id. Therefore, when getting the RoadLink objects
	 * that this Road object is associated, allow for multiple RoadLink objects to be associated.
	 */
	public void setRoadLinks() {
		for(RoadLink rl: SpaceBuilder.roadLinkContext.getObjects(RoadLink.class)) {
			if (rl.getFID().contentEquals(this.roadLinkID)) {
				this.roadLinks.add(rl);
			}
		}
	}

}
