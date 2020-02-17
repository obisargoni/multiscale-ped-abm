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

import java.util.HashMap;
import java.util.logging.Logger;

import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.parameter.Parameters;

/**
 * 
 * @author nick
 *
 */
public abstract class GlobalVars {
	
	private static Logger LOGGER = Logger.getLogger(GlobalVars.class.getName());
	
	public static Parameters params = RunEnvironment.getInstance ().getParameters();
	
	// Use to manage transformations between the CRS used in the geography and the CRS used for spatial calculations
	static String geographyCRSString = "EPSG:27700";
	
	/* These are strings that match entries in the repastcity.properties file.*/
	public static final String GISDataDirectory = "GISDataDirectory";
	public static final String BuildingShapefile = "BuildingShapefile";
	public static final String RoadShapefile = "RoadShapefile";
	public static final String BuildingsRoadsCoordsCache = "BuildingsRoadsCoordsCache";
	public static final String BuildingsRoadsCache = "BuildingsRoadsCache";
	
	public static double spaceScale = 1;
	public static double stepToTimeRatio = 1;
	public static double[] north = {0,1}; // Defines north, against which bearings are taken
	
	public static double tStep = 1;
	public static double pedVsd = 0.1; // standard dev of ped speeds
	public static double pedVavg = 0.8; // average pedestrian speed
	public static double pedMassAv = 60; // 60kg average mass
	public static double pedMasssd = 10; // 60kg average mass
	public static double interactionForceConstant = 100;
	public static int lookAheadTimeSteps = 3; // The number of timesteps to use when calculating an agents lookahead coordinate, used to identifying upcoming crossings.
	
	public static double defaultVehicleAcceleration = 0.1;
	public static double initialVehicleSpeed = 0.5;
	public static Integer maxVehicleSpeed = 3;
	
	public static final class GEOGRAPHY_PARAMS {
		
		/**
		 * Different search distances used in functions that need to find objects that are
		 * close to them. A bigger buffer means that more objects will be analysed (less
		 * efficient) but if the buffer is too small then no objects might be found. 
		 * The units represent distance in meters to match the CRS given in the GlobalVars.geographyCRSString.
		 * @see Geometry
		 * @see Route
		 */
		public enum BUFFER_DISTANCE {
			/** The smallest distance, rarely used. Approximately 0.001m*/
			SMALL(0.001, "0.001"),
			/** Intermediate distance between small and medium. Used for checking pedestrian proximity to
			 * crossing locations. */
			SMALLPLUS(2, "2"),
			/** Most commonly used distance, OK for looking for nearby houses or roads.
			 * Approximatey 110m */
			MEDIUM(110,"110"),
			/** Largest buffer, approximately 550m. I use this when doing things that
			 * don't need to be done often, like populating caches.*/
			LARGE(550,"550");
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
	}
	
	/** Names of contexts and projections. These names must match those in the
	 * parameters.xml file so that they can be displayed properly in the GUI. */
	public static final class CONTEXT_NAMES {
		
		public static final String MAIN_CONTEXT = "repastInterSim";
		public static final String MAIN_GEOGRAPHY = "Geography";
		
		public static final String ROAD_CONTEXT = "roadContext";
		public static final String ROAD_GEOGRAPHY = "roadGeography";
		
		public static final String PED_OBSTRUCTION_CONTEXT = "pedObstructContext";
		public static final String PED_OBSTRUCTION_GEOGRAPHY = "pedObstructGeography";
		
		public static final String ROAD_LINK_CONTEXT = "RoadLinkContext";
		public static final String ROAD_LINK_GEOGRAPHY = "RoadLinkGeography";
		
		public static final String JUNCTION_CONTEXT = "JunctionContext";
		public static final String JUNCTION_GEOGRAPHY = "JunctionGeography";
		
		public static final String VEHICLE_DESTINATION_CONTEXT = "DestinationContext";
		public static final String VEHICLE_DESTINATION_GEOGRAPHY = "DestinationGeography";
		public static final String PEDESTRIAN_DESTINATION_CONTEXT = "DestinationContext";
		public static final String PEDESTRIAN_DESTINATION_GEOGRAPHY = "DestinationGeography";
		
		public static final String ROAD_NETWORK = "RoadNetwork";
		
		public static final String AGENT_CONTEXT = "AgentContext";
		public static final String AGENT_GEOGRAPHY = "AgentGeography";
		
		public static final String BASE_COVERAGE = "baseGrid";
	}
	
	// Parameters used by transport networks
	public static final class TRANSPORT_PARAMS {
		
		// Used to synchronise routing code blocks
		public static Object currentBurglarLock = new Object();
		
		public static String routeDefaultDescription = "route";
		public static String routeCrossingDescription = "crossing";
		public static String routeRoadLinkChangeDescription = "link";
		
		// Proportion by which to expand bounding box of start-end coordinates for partial routing
		public static double partialBoundingBoxIncrease = 0.15;
	}
	
	public static final class GRID_PARAMS {
		
		private static HashMap<String, Integer> priorityValueMap = new HashMap<String, Integer> ();
		
		public static Integer defaultGridValue = 0;
		
		public static HashMap<String, Integer> getPriorityValueMap() {
			if (priorityValueMap.isEmpty()) {
				priorityValueMap.put("pedestrian", 1);
				priorityValueMap.put("vehicle", 2);
				priorityValueMap.put("pedestrian_obstruction", defaultGridValue);
			}
			return GRID_PARAMS.priorityValueMap;
		}
	}
	
	public static final class MOBILE_AGENT_PARAMS {
		
		public static double cautiousPriorityCostRatio = 10; // A pedestrian initialised with this perceives the cost of moving through vehicle priority space as 10x more costly than moving through ped priority space
		public static double boldPriorityCostRatio = 1; // Perceives cost of moving through vehicle space as the same as moving through pedestrian space
		
		// Parameter used to convert from proportion of road space occupied by vehicles to perception of grid cell cost
		public static double gridCellCostParam = params.getDouble("cellCostUpdate");
		
		public static double destinationArrivalDistance = 1.0; // Distance from destination at which agent is considered to ahve arrived
		
		// Constants used to estimate space taken up by vehicles 
		public static double vehicleWidth = 1.8;
		public static double vehicleLength = 4.0;
		public static double vehicleMPH = 20;
		public static double vehicleSpeed = vehicleMPH * 1609.34 * (1.0/(60*60)); // Convert mph to ms-1
		public static double vehicleReactionTime = 0.85; // Taken from trl report https://trl.co.uk/sites/default/files/PPR313_new.pdf
		public static double laneWidth = 3.65;
	}
}
