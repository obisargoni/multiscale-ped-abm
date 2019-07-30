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
	
}
