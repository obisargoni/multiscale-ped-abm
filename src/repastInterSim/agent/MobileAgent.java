package repastInterSim.agent;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.Geography;
import repastInterSim.environment.OD;

public class MobileAgent {
    private static int uniqueID = 0;

    protected int id;
    protected Geography<Object> geography; // Space the agent exists in
    protected Coordinate maLoc; // The coordinate of the centroid of the agent.
    protected OD origin; // The origin agent this agent starts at
    protected OD destination; // The destination agent that this agents is heading towards.
    
    MobileAgent(Geography<Object> geography, OD o, OD d){
    	this.id = MobileAgent.uniqueID++;
    	this.geography = geography;
    	this.origin = o;
    	this.destination = d;
    }
	
	public OD getOrigin() {
		return this.origin;
	}
	
	public OD getDestination() {
		return this.destination;
	}
	
	public void setLoc(Coordinate c) {
		this.maLoc = c;
	}
	
	public Coordinate getLoc() {
		return this.maLoc;
	}
	
    public String getLocString() {
    	return this.maLoc.toString();
    }
    
    public String getLocXString() {
    	return String.valueOf(this.maLoc.x);
    }
    
    public String getLocYString() {
    	return String.valueOf(this.maLoc.y);
    }
    
    public String getOriginXString() {
    	return String.valueOf(this.origin.getGeom().getCentroid().getCoordinate().x);
    }
    
    public String getOriginYString() {
    	return String.valueOf(this.origin.getGeom().getCentroid().getCoordinate().y);
    }
    
    public String getDestinationXString() {
    	return String.valueOf(this.destination.getGeom().getCentroid().getCoordinate().x);
    }
    
    public String getDestinationYString() {
    	return String.valueOf(this.destination.getGeom().getCentroid().getCoordinate().y);
    }
    
	public void setLoc() {
		// TODO Auto-generated method stub
		
	}
	
	public Geography<Object> getGeography(){
		return this.geography;
	}
	
	public int getID() {
		return this.id;
	}
	
	public void tidyForRemoval() {
		;
	}
}
