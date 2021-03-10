package repastInterSim.environment;

import java.util.Arrays;
import java.util.Comparator;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

public class CrossingAlternative extends Signal implements FixedGeography  {
	
	// Coordinates at which the ca meets pavement
	private long id;
	private Geometry geom = null;
	private Coordinate c1 = null;
	private Coordinate c2 = null;
	
	private String type;
	
	// id of the road link this crossing is located on
	private String roadLinkID;

	public CrossingAlternative(){
		super();
	}
	
	/*
	 * Calculate the distance to the input coordinate. Distance return is the distance to the nearest
	 * crossing alternative coordiante
	 * 
	 * @param Coordinate loc
	 * 		The coordinate to calculate the distance from
	 */
	public Double distanceTo(Coordinate loc) {
		double d1 = getC1().distance(loc);
		double d2 = getC2().distance(loc);
		return Math.min(d1, d2);
	}
	
	public Coordinate nearestCoord(Coordinate loc) {
		Coordinate[] coords = {getC1(), getC2()};
		
		Coordinate cNear = Arrays.stream(coords).min(Comparator.comparing(c->c.distance(loc))).get();
		
		return cNear;
	}

	public Coordinate farthestCoord(Coordinate loc) {
		Coordinate[] coords = {getC1(), getC2()};
		
		Coordinate cFar = Arrays.stream(coords).max(Comparator.comparing(c->c.distance(loc))).get();
		
		return cFar;
	}	
	
	/*
	 * Vehicles are assumed to yield at dedicated crossings resulting in zero vehicle flow
	 */
	public Integer getvFlow() {
		return 0;
	}
	
	public void setRoadLinkID(String rlID) {
		this.roadLinkID = rlID;
	}
	
	public String getRoadLinkID() {
		return this.roadLinkID;
	}
	
	/*
	 * This method uses the shape file data to set the list of road link ID this crossing alternative's signal controls
	 */
	public void setITNLinkIDs(String ids) {
		String[] itnRLIDs = ids.split(",");
		this.itnLinkIDs = itnRLIDs;
	}
	
	public void setSigPhases(String phases) {
		String[] allPhases = phases.split(",");
		
		// Initialise complete phase array. Assume that each phase is of same length as first. 
		char[][] finalAllPhasesArray = new char[allPhases.length][allPhases[0].length()];
		
		for (int i=0; i<allPhases.length; i++) {
			char[] phaseArray = new char[allPhases[i].length()];
			for(int j=0; j<allPhases[i].length(); j++) {
				phaseArray[j] = allPhases[i].charAt(j);
			}
			finalAllPhasesArray[i] = phaseArray;
		}
		
		this.phases = finalAllPhasesArray;
	}
	
	public void setPhaseDurs(String dursString) {
		String[] dursArray = dursString.split(",");
		int[] finalDursArray = new int[dursArray.length];
		for(int i=0; i<dursArray.length; i++) {
			finalDursArray[i] = Integer.parseInt(dursArray[i]);
		}
		this.phaseDurations = finalDursArray;
	}

	public Coordinate getC1() {
		return c1;
	}

	public void setC1(Coordinate c1) {		
		this.c1 = c1;
	}

	public Coordinate getC2() {
		return c2;
	}

	public void setC2(Coordinate c2) {
		this.c2 = c2;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	@Override
	public Geometry getGeom() {
		return this.geom;
	}

	@Override
	public void setGeom(Geometry g) {
		this.geom = g;
		
		int ncoords = g.getCoordinates().length;
		this.c1 = g.getCoordinates()[0];
		this.c2 = g.getCoordinates()[ncoords-1];
	}
	
	public void setID(long id) {
		this.id = id;
	}
	
	public long getID() {
		return this.id;
	}
}
