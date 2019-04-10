package repastInterSim.main;

/*
 * A placeholder until the real user panel is set up.
 */
public class UserPanel {
	public static double tStep = 1;
	public static double pedVsd = 0.1; // standard dev of ped speeds
	public static double pedVavg = 0.8; // average pedestrian speed
	public static double pedMassAv = 60; // 60kg average mass
	public static double pedMasssd = 10; // 60kg average mass
	public static double interactionForceConstant = 100;
	
	// Locations of the GIS data
	public static String GISDataDir = ".//data//";
	public static String VehicleRoadShapefile = "mastermap-topo_2903032_0 TopographicArea Roads Paths//topographicAreaVehicle_EPSG4277_Clipped_Single.shp";
	public static String PedestrianRoadShapefile = "mastermap-topo_2903032_0 TopographicArea Roads Paths//topographicAreaPedestrain_EPSG4277_Clipped_Single.shp";
	public static String PedestrianObstructionShapefile = "mastermap-topo_2903032_0 TopographicArea Roads Paths//topographicLineObstructing_EPSG4277_Clipped_RoadPathIntersection_Single.shp";
	public static String StartingZonesFile = "StartZones.shp";
	public static String DestinationsFile = "destCoordsEPSG4277.shp";
	
	public static String MAIN_CONTEXT = "repastInterSim";
	public static String MAIN_GEOGRAPHY = "Geography";
	
}
