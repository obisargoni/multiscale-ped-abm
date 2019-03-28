package repastInterSim;

/*
 * A placeholder until the real user panel is set up.
 */
public class UserPanel {
	static double tStep = 1;
	static double pedVsd = 0.1; // standard dev of ped speeds
	static double pedVavg = 0.8; // average pedestrian speed
	static double pedMassAv = 60; // 60kg average mass
	static double pedMasssd = 10; // 60kg average mass
	static double interactionForceConstant = 1;
	
	public static String MAIN_CONTEXT = "repastInterSim";
	public static String ROAD_CONTEXT = "RoadContext";
	public static String DESTINATION_CONTEXT = "DestinationContext";
	public static String PEDESTRIAN_CONTEXT = "PedestrianContext";
	
	public static String ROAD_GEOGRAPHY = "RoadGeography";
	public static String DESTINATION_GEOGRAPHY = "DestinationGeography";
	public static String PEDESTRIAN_GEOGRAPHY = "PedestrianGeography";
	
	// Locations of the GIS data
	public static String GISDataDir = ".//data//";
	public static String VehicleRoadShapefile = "mastermap-topo_2903032_0 TopographicArea Roads Paths//topographicAreaVehicle_EPSG4277.shp";
	public static String PedestrianRoadShapefile = "mastermap-topo_2903032_0 TopographicArea Roads Paths//topographicAreaPedestrian_EPSG4277.shp";
	public static String DestinationsFile = "destCoordsEPSG4277.shp";
	
}
