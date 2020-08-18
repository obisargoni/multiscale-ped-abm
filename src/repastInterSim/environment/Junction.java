/*
©Copyright 2012 Nick Malleson
This file is part of RepastCity.

RepastCity is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RepastCity is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RepastCity.  If not, see <http://www.gnu.org/licenses/>.
*/

package repastInterSim.environment;

import java.util.ArrayList;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

public class Junction implements FixedGeography{
	
	public static int UniqueID = 0;
	private int id ;
	private Geometry geom;
	private List<RoadLink> roadLinks; // The Roads connected to this Junction, used in GIS road network
	private String fid = null;
	
	// Attributes specific to pedestrian junctions
	private String p1pID;
	private String p2pID;
	private String p1rlID;
	private String p2rlID;
	
	private String v1pID;
	private String v2pID;
	private String v1rlID;
	private String v2rlID;
	
	
	public Junction() {
		this.id = UniqueID++;
		this.roadLinks = new ArrayList<RoadLink>();
	}
	
	public Junction(String fid) {
		this.id = UniqueID++;
		this.roadLinks = new ArrayList<RoadLink>();
		this.fid = fid;
	}
	
	
	/**
	 * Get the junction's unique ID (these are assigned in incremental order as
	 * junctions are created.
	 * 
	 * Is this made redundant by the inclusion of the junctions fid from the gis data?
	 */
	public int getId() {
		return id;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		Coordinate coord = this.geom.getCoordinate();
		return "Junction "+this.id+" ("+coord.x+","+coord.y+")";
	}
	
	public List<RoadLink> getRoads() {
		return this.roadLinks;
	}
	
	public void addRoadLink(RoadLink rl) {
		this.roadLinks.add(rl);
	}
	
	/**
	 * Tests if Junctions are equal by comparing the coorinates.
	 * @param j The junction to be compared with this one
	 * @return True if their coordinates are equal, false otherwise
	 */
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof Junction)) {
			return false;
		}
		Junction j = (Junction) obj;
		return this.getGeom().equals(j.getGeom());
	}

	/**
	 * Get the coordinate of this junction
	 */
	public Geometry getGeom() {
		return geom;
	}
	
	@Override
	public void setGeom(Geometry g) {
		this.geom = g;
		
	}
	
	public String getFID() {
		return this.fid;
	}

	public String getP1pID() {
		return p1pID;
	}

	public void setP1pID(String p1pID) {
		this.p1pID = p1pID;
	}

	public String getP2pID() {
		return p2pID;
	}

	public void setP2pID(String p2pID) {
		this.p2pID = p2pID;
	}

	public String getP1rlID() {
		return p1rlID;
	}

	public void setP1rlID(String p1rlID) {
		this.p1rlID = p1rlID;
	}

	public String getP2rlID() {
		return p2rlID;
	}

	public void setP2rlID(String p2rlID) {
		this.p2rlID = p2rlID;
	}

	public String getV1pID() {
		return v1pID;
	}

	public void setV1pID(String v1pID) {
		this.v1pID = v1pID;
	}

	public String getV2pID() {
		return v2pID;
	}

	public void setV2pID(String v2pID) {
		this.v2pID = v2pID;
	}

	public String getV1rlID() {
		return v1rlID;
	}

	public void setV1rlID(String v1rlID) {
		this.v1rlID = v1rlID;
	}

	public String getV2rlID() {
		return v2rlID;
	}

	public void setV2rlID(String v2rlID) {
		this.v2rlID = v2rlID;
	}
	
}
