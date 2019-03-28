/*
 * This class draws heavily on Nick Malleson's Road.java class
 * from the RepastCity project
 */

package repastInterSim;

import com.vividsolutions.jts.geom.Coordinate;

/*
 * Not sure what attributes and methods are needed for the road object.
 */
public class Road implements FixedGeography {
	
	private Coordinate coord;
	
	private String priority = "";
	
	
	/*
	 * Instance method
	 */
	public Road() {
		
	}
	
	public Road(String pri) {
		setPriority(pri);
	}
	
	
	public void setPriority(String pri) {
		this.priority = pri;
	}
	
	public String getPriority() {
		return this.priority;
	}
	
	public Coordinate getCoords() {
		return this.coord;
	}
	
	public void setCoords(Coordinate c) {
		this.coord = c;
	}

}
