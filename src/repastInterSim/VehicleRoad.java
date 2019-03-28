package repastInterSim;

public class VehicleRoad extends Road {


	private Integer capacity;
	private Integer speedLimit; // speed limit in kms-1
	private Integer nLanes;	
	
	/*
	 * Instance method
	 */
	public VehicleRoad() {
		super("vehicle");
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
	
	public void setSpeedLimit(Integer spdLim) {
		this.speedLimit = spdLim;
	}
	
	public Integer getSpeedLimit() {
		return this.speedLimit;
	}
	
	public void setNumLanes(Integer nL) {
		this.nLanes = nL;
	}
	
	public Integer getNumLanes() {
		return this.nLanes;
	}

}
