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

package repastInterSim.main;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.vividsolutions.jts.geom.Geometry;

/**
 * 
 * @author nick
 *
 */
public abstract class GlobalVars {
	
	private static Logger LOGGER = Logger.getLogger(GlobalVars.class.getName());
	
	/* These are strings that match entries in the repastcity.properties file.*/
	public static final String GISDataDirectory = "GISDataDirectory";
	public static final String BuildingShapefile = "BuildingShapefile";
	public static final String RoadShapefile = "RoadShapefile";
	public static final String BuildingsRoadsCoordsCache = "BuildingsRoadsCoordsCache";
	public static final String BuildingsRoadsCache = "BuildingsRoadsCache";
	
	public static double tStep = 1;
	public static double pedVsd = 0.1; // standard dev of ped speeds
	public static double pedVavg = 0.8; // average pedestrian speed
	public static double pedMassAv = 60; // 60kg average mass
	public static double pedMasssd = 10; // 60kg average mass
	public static double interactionForceConstant = 100;
	
	// Locations of the GIS data
	public static String GISDataDir = ".//data//";
	public static String VehicleRoadShapefile = "topographicAreaVehicle_EPSG27700_Clipped_Single.shp";
	public static String PedestrianRoadShapefile = "topographicAreaPedestrain_EPSG27700_Clipped_Single.shp";
	public static String RoadLinkShapefile = "mastermap-itn RoadLink Clipped.shp";
	public static String PedestrianObstructionShapefile = "topographicLineObstructing_EPSG27700_Clipped_RoadPathIntersection_Single.shp";
	public static String StartingZonesFile = "StartZones.shp";
	public static String DestinationsFile = "destCoords.shp";
	
	public static final class GEOGRAPHY_PARAMS {
		
		/**
		 * Different search distances used in functions that need to find objects that are
		 * close to them. A bigger buffer means that more objects will be analysed (less
		 * efficient) but if the buffer is too small then no objects might be found. 
		 * The units represent a lat/long distance so I'm not entirely sure what they are,
		 * but the <code>Route.distanceToMeters()</code> method can be used to roughly 
		 * convert between these units and meters.
		 * @see Geometry
		 * @see Route
		 */
		public enum BUFFER_DISTANCE {
			/** The smallest distance, rarely used. Approximately 0.001m*/
			SMALL(0.00000001, "0.001"),
			/** Most commonly used distance, OK for looking for nearby houses or roads.
			 * Approximatey 110m */
			MEDIUM(0.001,"110"),
			/** Largest buffer, approximately 550m. I use this when doing things that
			 * don't need to be done often, like populating caches.*/
			LARGE(0.005,"550");
			/**
			 * @param dist The distance to be passed to the search function (in lat/long?)
			 * @param distInMeters An approximate equivalent distance in meters.
			 */
			BUFFER_DISTANCE(double dist, String distInMeters) {
				this.dist = dist;
				this.distInMeters = distInMeters;
			}
			public double dist;
			public String distInMeters;
		}

		public static final double TRAVEL_PER_TURN = 1; // TODO Make a proper value for this
	}
	
	/** Names of contexts and projections. These names must match those in the
	 * parameters.xml file so that they can be displayed properly in the GUI. */
	public static final class CONTEXT_NAMES {
		
		public static final String MAIN_CONTEXT = "repastInterSim";
		public static final String MAIN_GEOGRAPHY = "Geography";
		
		public static final String BUILDING_CONTEXT = "BuildingContext";
		public static final String BUILDING_GEOGRAPHY = "BuildingGeography";
		
		public static final String ROAD_LINK_CONTEXT = "RoadLinkContext";
		public static final String ROAD_LINK_GEOGRAPHY = "RoadLinkGeography";
		
		public static final String JUNCTION_CONTEXT = "JunctionContext";
		public static final String JUNCTION_GEOGRAPHY = "JunctionGeography";
		
		public static final String ROAD_NETWORK = "RoadNetwork";
		
		public static final String AGENT_CONTEXT = "AgentContext";
		public static final String AGENT_GEOGRAPHY = "AgentGeography";
	
	}
	
	// Parameters used by transport networks
	public static final class TRANSPORT_PARAMS {
		
		// Used to synchronise routing code blocks
		public static Object currentBurglarLock = new Object();

	}

}
