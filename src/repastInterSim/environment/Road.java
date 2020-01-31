/*
 * This class draws heavily on Nick Malleson's Road.java class
 * from the RepastCity project
 */

package repastInterSim.environment;

import com.vividsolutions.jts.geom.Geometry;

import repastInterSim.main.SpaceBuilder;

/*
 * Not sure what attributes and methods are needed for the road object.
 */
public class Road implements FixedGeography {
	
	private String roadLinkFI; // The id of the road link this road object should be linked with
	private Geometry geom;	
	private String priority = ""; // Priority information comes from GIS data
	
	// Allows road agents to be joined with a particular RoadLink agent
	private RoadLink roadLink = null;
	
	
	/*
	 * Instance method
	 */
	public Road() {
	}
	
	public Road(String pri) {
		this.priority = pri;
	}
	
	public RoadLink getRoadLink() {
		return this.roadLink;
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
	
	public void setRoadLinkFI(String rlFID) {
		this.roadLinkFI = rlFID;
		
		// Once road link FID is set, can set the Road Link object
		setRoadLink();
	}
	
	public String getRoadLinkFI() {
		return this.roadLinkFI;
	}
	
	public void setRoadLink() {
		for(RoadLink rl: SpaceBuilder.roadLinkContext.getObjects(RoadLink.class)) {
			if (rl.getFID().contentEquals(this.roadLinkFI)) {
				this.roadLink = rl;
			}
		}
	}

}
