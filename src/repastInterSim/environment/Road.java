/*
 * This class draws heavily on Nick Malleson's Road.java class
 * from the RepastCity project
 */

package repastInterSim.environment;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import com.vividsolutions.jts.geom.Geometry;

import repastInterSim.main.SpaceBuilder;

/*
 * Not sure what attributes and methods are needed for the road object.
 */
public class Road implements FixedGeography, Serializable {
	
	private static final long serialVersionUID = 1L;
	private String roadLinkID; // The id of the road link this road object should be linked with
	private Geometry geom;	
	private String priority = ""; // Priority information comes from GIS data
	
	// Allows road agents to be joined with a particular RoadLink agent - ITN links and OR links
	private List<RoadLink> roadLinks = null;
	private RoadLink orRoadLink = null;
	
	
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
		for(RoadLink rl:getRoadLinks()) {
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
		for(RoadLink rl: getRoadLinks()) {
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
	}
	
	public String getRoadLinkID() {
		return this.roadLinkID;
	}

	public void setRoadLinks(List<RoadLink> rLs) {
		this.roadLinks = rLs;
	}
	
	public void setORRoadLink(RoadLink orRL) {
		this.orRoadLink = orRL;
	}
	
	public RoadLink getORRoadLink() {
		return this.orRoadLink;
	}

}
