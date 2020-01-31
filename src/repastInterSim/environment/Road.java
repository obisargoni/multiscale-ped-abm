/*
 * This class draws heavily on Nick Malleson's Road.java class
 * from the RepastCity project
 */

package repastInterSim.environment;

import com.vividsolutions.jts.geom.Geometry;

/*
 * Not sure what attributes and methods are needed for the road object.
 */
public class Road implements FixedGeography {
	
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
		setPriority(pri);
	}
	
	public Road(String pri, RoadLink rl) {
		setPriority(pri);
		setRoadLink(rl);
	}
	
	public void setRoadLink(RoadLink rl) {
		this.roadLink = rl;
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

}
