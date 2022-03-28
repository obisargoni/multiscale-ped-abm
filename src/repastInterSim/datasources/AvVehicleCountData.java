package repastInterSim.datasources;

public class AvVehicleCountData extends DefaultDataRecorder {
	
	private String fid;
	private double avVehicleCount;

	public AvVehicleCountData(String fid, double aVC) {
		super();
		this.fid = fid;
		this.avVehicleCount = aVC;
	}
	
	public String getFID() {
		return this.fid;
	}
	
	public double getAvVehCount() {
		return this.avVehicleCount;
	}

}
