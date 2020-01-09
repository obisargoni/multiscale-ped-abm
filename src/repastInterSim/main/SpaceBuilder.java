package repastInterSim.main;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.geometry.Envelope2D;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.operation.TransformException;

import com.opencsv.CSVReader;
import com.opencsv.exceptions.CsvException;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ISchedule;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.gis.util.GeometryUtil;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.space.gis.RepastCoverageFactory;
import repast.simphony.space.gis.WritableGridCoverage2D;
import repast.simphony.space.graph.Network;
import repast.simphony.util.collections.IndexedIterable;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.environment.Destination;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.VehicleDestinationContext;
import repastInterSim.environment.contexts.RoadContext;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.PedObstructionContext;
import repastInterSim.environment.contexts.PedestrianDestinationContext;
import repastInterSim.environment.contexts.RoadLinkContext;

public class SpaceBuilder extends DefaultContext<Object> implements ContextBuilder<Object> {
	
	private static Properties properties;
	private static Boolean isDirected = true; // Indicates whether the road network is directed ot not. 
	
	private static Integer pDI = 0; // Pedestrian destination index. Used to select which destination to assign to pedestrians
	private static Integer vDI = 0; // Vehicle destination index. Used to select which destination to assign to pedestrians
	
	private static Context<Object> context;
	private static Geography<Object> geography; 
	
	public static Context<Road> roadContext;
	public static Geography<Road> roadGeography;
	
	public static Context<PedObstruction> pedObstructContext;
	public static Geography<PedObstruction> pedObstructGeography;
	
	public static Context<Destination> vehicleDestinationContext;
	public static Geography<Destination> vehicleDestinationGeography;
	
	public static Context<Destination> pedestrianDestinationContext;
	public static Geography<Destination> pedestrianDestinationGeography;
	
	public static Context<RoadLink> roadLinkContext;
	public static Geography<RoadLink> roadLinkGeography;
	
	public static Context<Junction> junctionContext;
	public static Geography<Junction> junctionGeography;
	public static Network<Junction> roadNetwork;
	
	private static ArrayList<Geography> fixedGeographies = new ArrayList<Geography>();
	
	public static GeometryFactory fac;
	
	public static Map<String, String> values;
	
	/*
	 * A logger for this class. Note that there is a static block that is used to configure all logging for the model
	 * (at the bottom of this file).
	 */
	private static Logger LOGGER = Logger.getLogger(SpaceBuilder.class.getName());
	
	    /* (non-Javadoc)
	 * @see repast.simphony.dataLoader.ContextBuilder#build(repast.simphony.context.Context)
	 * 
	 */
	@Override
	public Context<Object> build(Context<Object> c) {
	    
		context = c;
		context.setId(GlobalVars.CONTEXT_NAMES.MAIN_CONTEXT);
		
		fac = new GeometryFactory();
		
		// Read in the model properties
		try {
			readProperties();
		} catch (IOException ex) {
			throw new RuntimeException("Could not read model properties,  reason: " + ex.toString(), ex);
		}
	   
		// Initiate geographic spaces
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, context, geoParams);
		geography.setCRS(GlobalVars.geographyCRSString);
		context.add(geography);
		
		// Road Geography stores polygons representing road and pavement surfaces
		roadContext = new RoadContext();
		roadGeography = createTypedGeography(Road.class, roadContext, GlobalVars.CONTEXT_NAMES.ROAD_GEOGRAPHY);
		context.addSubContext(roadContext);
		fixedGeographies.add(roadGeography);
		
		// Ped Obstruction context stores GIS linestrings representing barriers to pedestrian movement
		pedObstructContext = new PedObstructionContext();
		pedObstructGeography = createTypedGeography(PedObstruction.class, pedObstructContext, GlobalVars.CONTEXT_NAMES.PED_OBSTRUCTION_GEOGRAPHY);
		context.addSubContext(pedObstructContext);
		fixedGeographies.add(pedObstructGeography);

		// Road link geography is used to create the road network projection
		roadLinkContext = new RoadLinkContext();
		roadLinkGeography = createTypedGeography(RoadLink.class, roadLinkContext, GlobalVars.CONTEXT_NAMES.ROAD_LINK_GEOGRAPHY);
		fixedGeographies.add(roadLinkGeography);

		// Junction geography also used to create the road network
		junctionContext = new JunctionContext();
		junctionGeography = createTypedGeography(Junction.class, junctionContext, GlobalVars.CONTEXT_NAMES.JUNCTION_GEOGRAPHY);
		context.addSubContext(junctionContext);
		fixedGeographies.add(junctionGeography);
		
		// Destinations geography used for creating cache of destinations and their nearest road coordinates		
		vehicleDestinationContext = new VehicleDestinationContext();
		vehicleDestinationGeography = createTypedGeography(Destination.class, vehicleDestinationContext, GlobalVars.CONTEXT_NAMES.VEHICLE_DESTINATION_GEOGRAPHY);
		context.addSubContext(vehicleDestinationContext);
		fixedGeographies.add(vehicleDestinationGeography);

		pedestrianDestinationContext = new PedestrianDestinationContext();
		pedestrianDestinationGeography = createTypedGeography(Destination.class, pedestrianDestinationContext, GlobalVars.CONTEXT_NAMES.PEDESTRIAN_DESTINATION_GEOGRAPHY);
		context.addSubContext(pedestrianDestinationContext);
		fixedGeographies.add(pedestrianDestinationGeography);
		
		
	    // Load agents from shapefiles
		String GISDataDir = getProperty("GISDataDir");
		try {
			
			// Build the fixed environment
			
			// 1. Load vehicle and pedestrian destinations
			String vehicleDestinationsFile = GISDataDir + getProperty("VehicleDestinationsFile");
			GISFunctions.readShapefile(Destination.class, vehicleDestinationsFile, vehicleDestinationGeography, vehicleDestinationContext);
			
			String pedestrianDestinationsFile = GISDataDir + getProperty("PedestrianDestinationsFile");
			GISFunctions.readShapefile(Destination.class, pedestrianDestinationsFile, pedestrianDestinationGeography, pedestrianDestinationContext);
			
			// 2. Load roads
			String vehicleRoadFile = GISDataDir + getProperty("VehicleRoadShapefile");
			String pedestrianRoadFile = GISDataDir + getProperty("PedestrianRoadShapefile");
			GISFunctions.readShapefile(Road.class, vehicleRoadFile, roadGeography, roadContext);
			GISFunctions.readShapefile(Road.class, pedestrianRoadFile, roadGeography, roadContext);
			
			// 3. Load pedestrian obstruction boundaries
			String pedObstructionFile = GISDataDir + getProperty("PedestrianObstructionShapefile");
			GISFunctions.readShapefile(PedObstruction.class, pedObstructionFile, pedObstructGeography, pedObstructContext);
			
			// Build the road network
			
			// 1. Load the road links
			String roadLinkFile = GISDataDir + getProperty("RoadLinkShapefile");
			GISFunctions.readShapefile(RoadLink.class, roadLinkFile, roadLinkGeography, roadLinkContext);
			SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);

			
			// 2. roadNetwork
			NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, isDirected);
			builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
			roadNetwork = builder.buildNetwork();
			GISFunctions.buildGISRoadNetwork(roadLinkGeography, junctionContext,junctionGeography, roadNetwork);
			
		} catch (MalformedURLException | FileNotFoundException | MismatchedDimensionException e1 ) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		// Set the internal context and geography attributes of the destination agents
		IndexedIterable<Destination> vehicleDestinations = vehicleDestinationContext.getObjects(Destination.class);
		for (Destination d : vehicleDestinations) {
			d.setRootContext(context);
			d.setRootGeography(geography);
			d.setDestinationGeography(vehicleDestinationGeography);
		}
		
		IndexedIterable<Destination> pedestrianDestinations = pedestrianDestinationContext.getObjects(Destination.class);
		for (Destination d : pedestrianDestinations) {
			d.setRootContext(context);
			d.setRootGeography(geography);
			d.setDestinationGeography(pedestrianDestinationGeography);
		}
		
		
		// Now that shapefiles have been read and loaded as agents get the envelope that covers 
		// all fixed geography agents and use to create grid coverage 
		ReferencedEnvelope fixedGeographyEnvelope = GISFunctions.getMultipleGeographiesEnvelope(fixedGeographies);
		
		// Creates GIS grid with 1mx1m cells and adds to the geography projection
		int width = (int) Math.ceil(fixedGeographyEnvelope.getWidth())*2;
		int height = (int) Math.ceil(fixedGeographyEnvelope.getHeight())*2;
		// 2-category coverage (pedestrian priority areas and vehicle priority areas)
		 Category[] categories	= new Category[] {	
	        new Category("No data", Color.BLACK, 0),
	        new Category("Level 1", Color.GREEN, 1),
	        new Category("Level 3", Color.RED, 10)
	    };

		WritableGridCoverage2D pedGrid = RepastCoverageFactory.createWritableByteIndexedCoverage("pedGrid", width, height, fixedGeographyEnvelope, categories, null, 0);
		geography.addCoverage("pedGrid", pedGrid);
		
		// Initialise map from attribute value to numeric grid cell value
		Map<String,Integer> priorityGridValueMap = new HashMap<String, Integer> ();
		priorityGridValueMap.put("pedestrian", 1);
		priorityGridValueMap.put("vehicle", 10);
		// Loop over coverage grid cells to check values and number of cells
		GISFunctions.setGridCoverageValuesFromGeography(pedGrid, Road.class, roadGeography, "priority", priorityGridValueMap);
    	
    	// Get the number of pedestrian agents to add to the space from the parameters
    	Parameters params = RunEnvironment.getInstance().getParameters();
    	int nP = (int)params.getInteger("nPeds");		
		
		// Read in OD matrix data for vehicles from CSV
		List<String[]> vehicleFlows = readCSV(GISDataDir + getProperty("vehicleODFlowsFile"));
		
		// Read in OD matrix data for pedestrians from CSV
		List<String[]> pedestrianFlows = readCSV(GISDataDir + getProperty("pedestrianODFlowsFile"));

		
		// Schedule the creation of vehicle agents - tried doing this with annotations but it didnt work
		ISchedule schedule = RunEnvironment.getInstance().getCurrentSchedule();
	    ScheduleParameters vehicleScheduleParams = ScheduleParameters.createRepeating(1,50);
	    schedule.schedule(vehicleScheduleParams, this, "addVehicleAgents", vehicleFlows);
	    
		// Schedule the creation of pedestrian agents
	    ScheduleParameters pedestrianScheduleParams = ScheduleParameters.createRepeating(1,100);
	    schedule.schedule(pedestrianScheduleParams, this, "addPedestrianAgents", pedestrianFlows);
		
		return context;
		
	}
	
	/*
	 * Scheduled method that adds vehicle agents to the simulation. Each method call vehicle agents
	 * are initialised with origins and destinations taken from an OD matrix. The OD matrix values 
	 * are used to control the frequency of vehicles initialised with each OD pair, therefore controlling
	 * the flow of vehicles 
	 * 
	 */
	public void addVehicleAgents(List<String[]> odData) {
		
		// Get number of possible origins/destinations
		int nOD = vehicleDestinationGeography.size();
		
		// Use to generate random number and use to determine which OD pairs to use when creating a vehicle agent - check this is valid
		Random rn = new Random();
		
		// Iterate through all OD pairs and initialise vehicle moving between these two if prob is above threshold
		for (int iO = 0; iO<nOD; iO++) {
			int iD = vDI % nOD;
			
			// Get the OD matrix entry
			Float flow = Float.parseFloat(odData.get(iO)[iD]);
			float threshold = rn.nextFloat();

			
			// Create vehicle instance probabilistically according to flow rates
			if (flow > threshold) {
				Coordinate o = vehicleDestinationContext.getObjects(Destination.class).get(iO).getGeom().getCentroid().getCoordinate();
				Destination d = vehicleDestinationContext.getObjects(Destination.class).get(iD);
				addVehicle(o, d);
			}
		}
		
		// Increment the vehicle destination index so that a different destination is considered next time.
		vDI++;
	}
	
    /*
     * Scheduled method that adds pedestrian agents to the simulation. Origin and destination of the pedestrian
     * is set probabilistically based on OD matrix flow values.
     * 
     * @param odData The OD flow data used to create pedestrian agents with origins and destinations that match the flow data
     */
	public void addPedestrianAgents(List<String[]> odData) {
		
		// Get number of possible origins/destinations
		int nOD = pedestrianDestinationGeography.size();
		
		// Use to generate random number and use to determine which OD pairs to use when creating a vehicle agent - check this is valid
		Random rn = new Random();
		
		// Iterate through all OD pairs and initialise vehicle moving between these two if prob is above threshold
		for (int iO = 0; iO<nOD; iO++) {
			
			int iD = pDI % nOD;
			
			// Get the OD matrix entry
			Float flow = Float.parseFloat(odData.get(iO)[iD]);
			float threshold = rn.nextFloat();
			
			// Create vehicle instance probabilistically according to flow rates
			if (flow > threshold) {
				Coordinate o = pedestrianDestinationContext.getObjects(Destination.class).get(iO).getGeom().getCentroid().getCoordinate();
				Destination d = pedestrianDestinationContext.getObjects(Destination.class).get(iD);
				addPed(o, d);
			}
		}
		
		// Increment the vehicle destination index so that a different destination is considered next time.
		pDI++; 
	}
	
	/*
	 * Create a destination agent and add the agent to the context. Generate a random coordinate that lies within a 
	 * boundary in the geography and moves the agent to that coordinate.
	 * 
	 * @param context
	 * 			The context to add the destination agent to
	 * @param geography
	 * 			The geography to move the destination agent to
	 * @param gF
	 * 			The geometry factory used to generate the geometry to move the destination agent to
	 * @param bndry
	 * 			The geometry to select a random coordinate from within
	 * @returns d
	 * 			The destination agent that was added to the context and geography
	 */
	public Destination addRandomDestination(Geometry bndry) {
		
		Destination d = new Destination();
		context.add(d);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = GeometryUtil.generateRandomPointsInPolygon(bndry, 1).get(0);
		
		Geometry destGeom = fac.createPoint(destCoord);
		geography.move(d, destGeom);
				
		return d;
	}
	
	/*
	 * Create a destination agent and add the agent to the context. Creates a coordinate from input x and y values
	 * and moves the agent to that coordinate.
	 * 
	 * @param context
	 * 			The context to add the destination agent to
	 * @param geography
	 * 			The geography to move the destination agent to
	 * @param gF
	 * 			The geometry factory used to generate the geometry to move the destination agent to
	 * @param paramX
	 * 			The x coordinate for the destination
	 * @param paramY
	 * 			The y coordinate for the destination
	 * @returns d
	 * 			The destination agent that was added to the context and geography
	 */
	public Destination addUserDestination( String paramX, String paramY) {
		
		Destination d = new Destination();
		context.add(d);
		
		// Get the x&y coords for the destination set by the user
		Parameters  params = RunEnvironment.getInstance().getParameters();
		double xCoord = (double)params.getInteger(paramX);
		double yCoord = (double)params.getInteger(paramY);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = new Coordinate(xCoord, yCoord);
		
		Geometry destGeom = fac.createPoint(destCoord);
		
		geography.move(d, destGeom);
		
		return d;
		
	}
	
	/*
	 * Adds pedestrian agents to the context and geography. Creates a new pedestrian agent and adds the pedestrian to the context.
	 * Creates a circle geometry in the geography and moves the pedestrian agent to this geometry. Assigns the pedestrian agent
	 * its destination and initial direction.
	 * 
	 * @param context
	 * 			The context to add the pedestrian agent to
	 * @param geography
	 * 			The geography to move the pedestrian agent to
	 * @param 
	 *			The geometry factory used to generate the geometry to move the pedestrian agent to
	 * @param coord
	 * 			The coordinate to move the centroid of the pedestrian to in the geography
	 */
    public Ped addPed(Coordinate coord, Destination d)  {
        
        // Instantiate a new pedestrian agent and add the agent to the context

    	Ped newPed = new Ped(geography, pedestrianDestinationGeography, fac, d);
        context.add(newPed);
        
        // Create a new point geometry.
		Point pt = GISFunctions.pointGeometryFromCoordinate(coord);
		
		// Transform the coordinate so that the circle can be created using a radius in metres
		Geometry circle = pt.buffer(newPed.getRad());
		
		// Move the pedestrian to this geometry
		moveAgentToGeometry(geography, circle, newPed);
		
		// Set the location attribute of the pedestrian agent to be its current location. Simplifies future calculations
		newPed.setLoc();
		
		// Once pedestrian location has been set, can set the coordinates to travel along
		try {
			newPed.getRoute().setPedestrianRoute();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Now set the initial bearing of the pedestrian to be the direction of the first coordinate on the route
        Coordinate routeCoord = newPed.getRoute().getRouteXCoordinate(0);
		double ang = newPed.setBearingToDestinationCoord(routeCoord);
		newPed.setPedestrianBearing(ang);
        
        return newPed;
    }
    
    /*
     * Initialise a vehicle agent and add to to the context and projection
     */
    private Vehicle addVehicle(Coordinate o, Destination d) {
		Vehicle V = new Vehicle(geography, vehicleDestinationGeography, GlobalVars.maxVehicleSpeed, GlobalVars.defaultVehicleAcceleration, GlobalVars.initialVehicleSpeed, d);
		context.add(V);
		Point pt = fac.createPoint(o);
		Geometry vehicleCircle = pt.buffer(2);
		moveAgentToGeometry(geography, vehicleCircle, V);
		V.setLoc();
		
		return V;
    }
	
	/*
	 * Return the geometry associated to an agent.
	 * 
	 * @param geography
	 * 			The geography the agent belongs to
	 * @param agent
	 * 			The agent to get the associated geography of
	 */
	public static <T> Geometry getAgentGeometry(Geography<T> geography, Object agent) {
		Geometry geom = geography.getGeometry(agent);
		return geom;
	}
	
	/*
	 * Move an agent to the input geometry in the input geography.
	 * 
	 * @param geography
	 * 			The geography to add the agent to
	 * @param geom
	 * 			The geometry to move the agent to
	 * @param agent
	 * 			The agent to move to the geometry
	 */
	public static <T> void moveAgentToGeometry(Geography<T> geography, Geometry geom, T agent) {
		geography.move(agent, geom);
	}	
	
	/*
	 * Move an agent to the input geometry in the input geography.
	 * 
	 * @param geography
	 * 			The geography to add the agent to
	 * @param agent
	 * 			The agent to move to the geometry
	 * @param distance
	 * 			The distance to move the agent by
	 * @param angle
	 * 			The direction to move the agent in
	 */
	public static <T> void moveAgentByVector(Geography<T> geography, T agent, double distance, double angle) {
		geography.moveByVector(agent, distance, angle);
	}
	
	/*
	 * Calculate the size of an iterable
	 * 
	 * @param i
	 * 			The iterable to calculate the size of
	 * 
	 * @returns 
	 * 			The size of the iterable
	 */
	public static int sizeOfIterable(Iterable<?> i) {
		int size = 0;
		Iterator<?> it = i.iterator();
		while (it.hasNext()) {
			size++;
			it.next();
		}
		return size;
	}
	
	/**
	 * Get the value of a property in the properties file. If the input is empty or null or if there is no property with
	 * a matching name, throw a RuntimeException.
	 * 
	 * @param property
	 *            The property to look for.
	 * @return A value for the property with the given name.
	 */
	public static String getProperty(String property) {
		if (property == null || property.equals("")) {
			throw new RuntimeException("getProperty() error, input parameter (" + property + ") is "
					+ (property == null ? "null" : "empty"));
		} else {
			String val = SpaceBuilder.properties.getProperty(property);
			if (val == null || val.equals("")) { // No value exists in the
													// properties file
				throw new RuntimeException("checkProperty() error, the required property (" + property + ") is "
						+ (property == null ? "null" : "empty"));
			}
			return val;
		}
	}

	/**
	 * Read the properties file and add properties. Will check if any properties have been included on the command line
	 * as well as in the properties file, in these cases the entries in the properties file are ignored in preference
	 * for those specified on the command line.
	 * 
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	private void readProperties() throws FileNotFoundException, IOException {

		File propFile = new File("./respastInterSim.properties");
		if (!propFile.exists()) {
			throw new FileNotFoundException("Could not find properties file in the default location: "
					+ propFile.getAbsolutePath());
		}

		LOGGER.log(Level.FINE, "Initialising properties from file " + propFile.toString());

		SpaceBuilder.properties = new Properties();

		FileInputStream in = new FileInputStream(propFile.getAbsolutePath());
		SpaceBuilder.properties.load(in);
		in.close();

		// See if any properties are being overridden by command-line arguments
		for (Enumeration<?> e = properties.propertyNames(); e.hasMoreElements();) {
			String k = (String) e.nextElement();
			String newVal = System.getProperty(k);
			if (newVal != null) {
				// The system property has the same name as the one from the
				// properties file, replace the one in the properties file.
				LOGGER.log(Level.INFO, "Found a system property '" + k + "->" + newVal
						+ "' which matches a NeissModel property '" + k + "->" + properties.getProperty(k)
						+ "', replacing the non-system one.");
				properties.setProperty(k, newVal);
			}
		} // for
		return;
	} // readProperties
	
	/*
	 * Read CSV data into List of String arrays. Catch exceptions.
	 * 
	 * @param String filePath. Path of the CSV file to read
	 * 
	 * @returns List<String[]> The csv data as a list of string arrays
	 */
	private List<String[]> readCSV(String filePath) {
		// Read in OD matrix data for pedestrians from CSV
		List<String[]> csvData = null;
		try {
		     CSVReader reader = new CSVReader(new FileReader(filePath));
		     try {
				csvData = reader.readAll();
			} catch (CsvException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		     reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return csvData;
	}
	
	/**
	 * Initialise a geography projection. Set the geography name and coordinate reference system.
	 * 
	 * @param <T> 
	 * 			The type of object to be contained in the context and geography
	 * @param cl
	 * 			The class of object to be contained in the context and geography
	 * @param context
	 * 			Context to use when creating geography projection. Must be type Context<T>
	 * @param geographyName
	 * 			The name to use for the geography
	 * @return
	 * 			Geography project of type T
	 */
	private <T> Geography<T> createTypedGeography(Class<T> cl, Context<T> context, String geographyName) {
		GeographyParameters<T> GeoParams = new GeographyParameters<T>();
		Geography<T> g = GeographyFactoryFinder.createGeographyFactory(null).createGeography(geographyName, context, GeoParams);
		g.setCRS(GlobalVars.geographyCRSString);
		
		return g;
	}
}
