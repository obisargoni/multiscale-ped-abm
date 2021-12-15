package repastInterSim.datasources;

public class PedRouteData {
	
	private int id;
	private String startPavementJunctionID;
	private String destPavementJunctionID;
	private String fullStrategicPathString;
	private String fullTacticalRouteString;
	

	public PedRouteData(int pID, String spjID, String dpjID, String fsp, String ftp) {
		// TODO Auto-generated constructor stub
		this.id = pID;
		this.startPavementJunctionID = spjID;
		this.destPavementJunctionID = dpjID;
		this.fullStrategicPathString = fsp;
		this.fullTacticalRouteString = ftp;
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

}
