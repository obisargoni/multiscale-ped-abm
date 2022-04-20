package repastInterSim.datasources;

public class PedRouteData extends DefaultDataRecorder {
	
	private int id;
	private String startPavementJunctionID;
	private String destPavementJunctionID;
	private String fullStrategicPathString;
	private String fullTacticalRouteString;
	private String journeyDistance;
	private String journeyDuration;
	

	public PedRouteData(int pID, String spjID, String dpjID, String fsp, String ftp, double jD, int jDr) {
		super();
		this.id = pID;
		this.startPavementJunctionID = spjID;
		this.destPavementJunctionID = dpjID;
		this.fullStrategicPathString = fsp;
		this.fullTacticalRouteString = ftp;
		this.journeyDistance = String.valueOf(jD);
		this.journeyDuration = String.valueOf(jDr);
	}
	
	public int getID() {
		return this.id;
	}
	
	public String getStartPavementJunctionID() {
		return this.startPavementJunctionID;
	}

	public String getDestPavementJunctionID() {
		return destPavementJunctionID;
	}


	public String getFullStrategicPathString() {
		return fullStrategicPathString;
	}

	public String getFullTacticalRouteString() {
		return fullTacticalRouteString;
	}

	public String getJourneyDistance() {
		return journeyDistance;
	}

	public String getJourneyDuration() {
		return journeyDuration;
	}

}
