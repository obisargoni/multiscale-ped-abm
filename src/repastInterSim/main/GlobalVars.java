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
import java.util.Properties;

import com.vividsolutions.jts.geom.Geometry;
/**
 * 
 * @author nick
 *
 */
public abstract class GlobalVars {
	
	// Directories for data exports
	public static String outputDir = ".\\output\\";
	public static String exportDir = outputDir + "\\export\\";
	
	//private static Logger LOGGER = Logger.getLogger(GlobalVars.class.getName());
	
	// Store properties in Global Vars. Used to access data file names
	public static Properties properties;
	
	
	// Use to manage transformations between the CRS used in the geography and the CRS used for spatial calculations
	public static String geographyCRSString = "EPSG:27700";
	
	/* These are strings that match entries in the repastcity.properties file.*/
	public static final String GISDataDirectory = "GISDataDir";
	public static final String TestDataDir = "TestDataDir";
	public static final String VehicleRoadShapefile = "VehicleRoadShapefile";
	public static final String PedestrianRoadShapefile = "PedestrianRoadShapefile";
	public static final String RoadLinkShapefile = "RoadLinkShapefile";
	public static final String ORRoadLinkShapefile = "ORRoadLinkShapefile";
	public static final String PedestrianObstructionShapefile = "PedestrianObstructionShapefile";
	public static final String VehicleDestinationsFile = "VehicleDestinationsFile";
	public static final String PedestrianDestinationsFile = "PedestrianDestinationsFile";
	public static final String vehicleODFlowsFile = "vehicleODFlowsFile";
	public static final String pedestrianODFlowsFile = "pedestrianODFlowsFile";
	public static final String PavementJunctionShapeFile = "PavementJunctionsShapefile";
	
	public static final String ODRoadLinkCoordsCache = "ODRoadLinkCoordsCache";
	public static final String ODORRoadLinkCoordsCache = "ODORRoadLinkCoordsCache";
	public static final String BuildingsRoadsCache = "BuildingsRoadsCache";
	public static final String RoadLinkRoadsCache = "RoadLinkRoadsCache";
	public static final String ODPavementJunctionCache = "ODPavementJunctionCache";
			
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
	public static double deafultTacticalPlanningHorizon = 20.0; // Degrees
	
	public static double defaultVehicleAcceleration = 0.8;
	public static double defaultVehicleDecceleration = 4.5; // These two values taken from SUMO car following model https://sumo.dlr.de/pdf/KraussDiss.pdf
	public static double initialVehicleSpeed = 0.5;
	public static Integer maxVehicleSpeed = 3;
	public static double vehicleLength = 3.5; // 3.5m
	public static double vehicleWidth = 1.7; // 1.7m
	public static double obstacleYieldDistance = 2; // Keep 2m distance from obstacles such as peds and traffic signals
	
	
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
		
		public static final String OR_ROAD_LINK_CONTEXT = "ORRoadLinkContext";
		public static final String OR_ROAD_LINK_GEOGRAPHY = "ORRoadLinkGeography";
		
		public static final String JUNCTION_CONTEXT = "JunctionContext";
		public static final String JUNCTION_GEOGRAPHY = "JunctionGeography";
		
		public static final String OR_JUNCTION_CONTEXT = "ORJunctionContext";
		public static final String OR_JUNCTION_GEOGRAPHY = "ORJunctionGeography";
		
		public static final String PAVEMENT_JUNCTION_CONTEXT = "PavementJunctionContext";
		public static final String PAVEMENT_JUNCTION_GEOGRAPHY = "PavementJunctionGeography";
		
		public static final String PAVEMENT_LINK_CONTEXT = "PavementLinkContext";
		public static final String PAVEMENT_LINK_GEOGRAPHY = "PavementLinkGeography";

		
		public static final String VEHICLE_DESTINATION_CONTEXT = "DestinationContext";
		public static final String VEHICLE_DESTINATION_GEOGRAPHY = "DestinationGeography";
		public static final String PEDESTRIAN_DESTINATION_CONTEXT = "DestinationContext";
		public static final String PEDESTRIAN_DESTINATION_GEOGRAPHY = "DestinationGeography";
		
		public static final String CA_CONTEXT = "CAContext";
		public static final String CA_GEOGRAPHY = "CAGeography";
		
		public static final String ROAD_NETWORK = "RoadNetwork";
		public static final String OR_ROAD_NETWORK = "ORRoadNetwork";
		public static final String GRID_NETWORK = "GridNetwork";
		
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
		public static double partialBoundingBoxIncrease = 0.5;
	}
	
	public static final class GRID_PARAMS {
		
		private static HashMap<String, Integer> priorityValueMap = new HashMap<String, Integer> ();
		
		public static Integer defaultGridValue = 0;
		
		public static HashMap<String, Integer> getPriorityValueMap() {
			if (priorityValueMap.isEmpty()) {
				priorityValueMap.put("pedestrian", 1);
				priorityValueMap.put("pedestrian_crossing", 2);
				priorityValueMap.put("vehicle", 3);
				priorityValueMap.put("road_link", 4);
				priorityValueMap.put("pedestrian_obstruction", defaultGridValue);
			}
			return GRID_PARAMS.priorityValueMap;
		}
	}
	
	public static class MOBILE_AGENT_PARAMS {		
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
