package repastInterSim.agent;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.Geography;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.OD;
import repastInterSim.environment.Road;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.SpaceBuilder;

public class MobileAgent {
    private static int uniqueID = 0;

    protected int id;
    protected Coordinate maLoc; // The coordinate of the centroid of the agent.
    protected OD origin; // The origin agent this agent starts at
    protected OD destination; // The destination agent that this agents is heading towards.
    protected Coordinate defaultDestination;
    protected double bearing;
    
    MobileAgent(OD o, OD d){
    	this.id = MobileAgent.uniqueID++;
    	this.origin = o;
    	this.destination = d;
    	
    	// Initialise default destination as the actual destination
    	this.defaultDestination = d.getGeom().getCoordinate();
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
	
	public Road getCurrentRoad() {
		Road r = null;
		try {
			r = GISFunctions.getCoordinateRoad(this.getLoc(), SpaceBuilder.roadGeography);
		} catch (RoutingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return r;
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
		return SpaceBuilder.geography;
	}
	
	public int getID() {
		return this.id;
	}
	
	public void tidyForRemoval() {
		;
	}
	
	public Coordinate getDefaultDestination() {
		return this.defaultDestination;
	}
	
	public void setDefaultDestination(Coordinate c) {
		this.defaultDestination = c;
	}
	
	public double getBearing() {
		return this.bearing;
	}
	
	public void setBearing(double b) {
		this.bearing = b;
	}
}
