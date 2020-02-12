package repastInterSim.agent;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.Geography;
import repastInterSim.environment.Destination;

public class MobileAgent {
    private static int uniqueID = 0;

    protected int id;
    protected Geography<Object> geography; // Space the agent exists in
    protected Coordinate maLoc; // The coordinate of the centroid of the agent.
    protected Destination destination; // The destination agent that this agents is heading towards.
    
    MobileAgent(Geography<Object> geography, Destination d){
    	this.id = MobileAgent.uniqueID++;
    	this.geography = geography;
    	this.destination = d;
    }
	
	public Destination getDestination() {
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
