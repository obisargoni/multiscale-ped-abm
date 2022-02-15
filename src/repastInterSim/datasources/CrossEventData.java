package repastInterSim.datasources;

import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;

public class CrossEventData extends DefaultDataRecorder {
	
	private int pedID;
	private String tacticalEdgeID;
	private String crossingType;
	private String crossingCoordinatesString;
	private Double ttc;

	public CrossEventData(int pID, String tEID, String cT, List<Coordinate> cCs) {
		super();
		this.pedID=pID;
		this.tacticalEdgeID=tEID;
		this.crossingType=cT;
		
    	String ccString = "";
    	for (Coordinate c: cCs) {
    		ccString += c.toString();
    		ccString += ",";
    	}
		this.crossingCoordinatesString=ccString;
		
		this.ttc=null;
	}

	public int getPedID() {
		return pedID;
	}

	public String getTacticalEdgeID() {
		return tacticalEdgeID;
	}

	public String getCrossingType() {
		return crossingType;
	}

	public String getCrossingCoordinatesString() {
		return crossingCoordinatesString;
	}

	public String getTTCString() {
    	if (this.ttc==null) {
    		return "";
    	}
    	else {
    		return this.ttc.toString();
    	}
	}

	public void setTTC(Double ttc) {
		this.ttc = ttc;
	}

}
