package repastInterSim.datasources;

public class VehicleRouteData extends DefaultDataRecorder {
	
	private int id;
	private String startJunctionID;
	private String destJunctionID;
	private String fullStrategicPathString;
	private String journeyDistance;
	private String journeyDuration;
	

	public VehicleRouteData(int vID, String sjID, String djID, String fsp, double jD, int jDr) {
		super();
		this.id = vID;
		this.startJunctionID = sjID;
		this.destJunctionID = djID;
		this.fullStrategicPathString = fsp;
		this.journeyDistance = String.valueOf(jD);
		this.journeyDuration = String.valueOf(jDr);
	}
	
	public int getID() {
		return this.id;
	}
	
	public String getStartJunctionID() {
		return this.startJunctionID;
	}

	public String getDestJunctionID() {
		return destJunctionID;
	}


	public String getFullStrategicPathString() {
		return fullStrategicPathString;
	}

	public String getJourneyDistance() {
		return journeyDistance;
	}

	public String getJourneyDuration() {
		return journeyDuration;
	}

}
