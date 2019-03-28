package repastInterSim;

public class PedestrianRoad extends Road {


	private Integer capacity;
	
	/* Instance method
	 * 
	 */
	public PedestrianRoad() {
		super("pedestrian");
	}
	
	
	// Might want to include some other attributes to do with junctions and the road network
	
	/*
	 * Getters and setters for capacity, speedLimit and lanes attributes.
	 */
	public void setCapacity(Integer cap) {
		this.capacity = cap;
	}
	
	public Integer getCapacity() {
		return this.capacity;
	}

}
