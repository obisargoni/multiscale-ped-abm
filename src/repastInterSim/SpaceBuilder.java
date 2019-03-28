package repastInterSim;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.factory.Hints;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Point;

import contexts.DestinationContext;
import contexts.PedestrianContext;
import contexts.RoadContext;
import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.gis.util.GeometryUtil;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.space.gis.ShapefileLoader;
import repast.simphony.space.gis.SimpleAdder;
import repast.simphony.util.collections.IndexedIterable;

public class SpaceBuilder extends DefaultContext<Object> implements ContextBuilder<Object> {
	
	public static double spaceScale = 1;
	public static double[] north = {0,1}; // Defines north, against which bearings are taken
	
	// Use to manage transformations between the CRS used in the geography and the CRS used for spatial calculations
	public static String geographyCRSString = "EPSG:4277";
	public static String calculationCRSString = "EPSG:27700";
	public static MathTransform transformToGeog;
	public static MathTransform transformToCalc;
	
	// Initialise contexts and geographies	
	private static Context<Object> mainContext;
	
	public static Context<Road> roadContext;
	public static Geography<Road> roadGeography;
	private static Context<Ped> pedContext;
	private static Geography<Ped> pedGeography;
	private static Context<Destination> destContext;
	private static Geography<Destination> destGeography;
	
	
	
	 /* (non-Javadoc)
	 * @see repast.simphony.dataLoader.ContextBuilder#build(repast.simphony.context.Context)
	 * 
	 */
	@Override
	public Context<Object> build(Context<Object> context) {
		
		mainContext = context;
		mainContext.setId(UserPanel.MAIN_CONTEXT);
	   
		// Initiate geography parameters
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>(new SimpleAdder<Object>());
		GeographyParameters<Road> geoParamsRoad = new GeographyParameters<Road>(new SimpleAdder<Road>());
		GeographyParameters<Ped> geoParamsPed = new GeographyParameters<Ped>(new SimpleAdder<Ped>());
		GeographyParameters<Destination> geoParamsDest = new GeographyParameters<Destination>(new SimpleAdder<Destination>());
		
		// Use GB Coordinate projection, also define a transform between degree and metre projections
		geoParams.setCrs(geographyCRSString);
		geoParamsRoad.setCrs(geographyCRSString);
		geoParamsPed.setCrs(geographyCRSString);
		geoParamsDest.setCrs(geographyCRSString);
		
		
		Geography<Object> geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("Geography", mainContext, geoParams);
		mainContext.add(geography);
		

		// Not sure what this line does and whether it is required for coordinate transformation
		Hints.putSystemDefault(Hints.FORCE_LONGITUDE_FIRST_AXIS_ORDER, Boolean.TRUE);
		CoordinateReferenceSystem geographyCRS = null;
		CoordinateReferenceSystem calculationCRS = null;
		try {
			geographyCRS = CRS.decode(geographyCRSString);
			calculationCRS = CRS.decode(calculationCRSString);
			transformToGeog = CRS.findMathTransform(calculationCRS, geographyCRS);
			transformToCalc = CRS.findMathTransform(geographyCRS, calculationCRS);

		} catch (FactoryException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
		GeometryFactory fac = new GeometryFactory();
	    
	    // Code below taken from the 'Geography' repast example
		// Create an area in which to create agents.  This border is loaded from a shapefile.
		String boundaryFilename = ".//data//JunctClipEPSG4277.shp";
		List<SimpleFeature> features = loadFeaturesFromShapefile(boundaryFilename);
		Geometry boundary = (MultiPolygon)features.iterator().next().getDefaultGeometry();
		

		
		try {
			// Add the pedestrian and vehicle roadways to the simulation - need to figure out how to assign priority to road when being read in from file
			roadContext = new RoadContext();
			roadGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(
					UserPanel.ROAD_GEOGRAPHY, roadContext, geoParamsRoad);
			String vehicleRoadFile = UserPanel.GISDataDir + UserPanel.VehicleRoadShapefile;
			String pedestrianRoadFile = UserPanel.GISDataDir + UserPanel.PedestrianRoadShapefile;
			readShapefile(Road.class, vehicleRoadFile, roadGeography, roadContext);
			readShapefile(Road.class, pedestrianRoadFile, roadGeography, roadContext);
			// Not sure what this line does
			//SpatialIndexManager.createIndex(roadGeography, Road.class);
			mainContext.addSubContext(roadContext);
			
			
			// Add destinations to the model
			destContext = new DestinationContext();
			destGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(
					UserPanel.DESTINATION_GEOGRAPHY, destContext, geoParamsDest);
			String destinationsFile = UserPanel.GISDataDir + UserPanel.DestinationsFile;
			readShapefile(Destination.class, destinationsFile, destGeography, destContext);
			// Not sure what this line does
			//SpatialIndexManager.createIndex(destGeography, Destination.class);
			mainContext.addSubContext(destContext);
			
			// Create the pedestrian context but do not yet add pedestrians to the model
			pedContext = new PedestrianContext();
			pedGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(
					UserPanel.PEDESTRIAN_GEOGRAPHY, pedContext, geoParamsPed);
			mainContext.addSubContext(pedContext);
			
			
		} catch (MalformedURLException | FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		           
	    // Get the destinations so that they can be assigned to pedestrian agents
		IndexedIterable<Destination> destinations = destContext.getObjects(Destination.class); 
		
		// Comment out anything to do with adding pedestrians to the simulation
		
	    
    	// Get the number of pedetrian agent to add to the space from the parameters
    	Parameters params = RunEnvironment.getInstance().getParameters();
    	int nP = (int)params.getInteger("nPeds");
    	nP= 0;
    	
		// Generate random points in the area to create agents.
		List<Coordinate> agentCoords = GeometryUtil.generateRandomPointsInPolygon(boundary, nP);
		
		// Create the agents from the collection of random coords.
		Random randCoord = new Random();
		int i = 0;
		int destinationIndex = 0;
		for (Coordinate coord : agentCoords) {
			
			// Generate a random initial direction for the pedestrian
    		double randBearing = randCoord.nextFloat() * 2 * Math.PI;
    		
    		// Crude way to assign different destinations
    		if (i > nP / 2) {
    			destinationIndex = 1; 
    		}
			
    		try {
				Ped newPed = addPed(context, geography, fac, randBearing, coord, destinations.get(destinationIndex), Color.BLUE);
			} catch (MismatchedDimensionException | TransformException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
    		i+=1;
		}
		
		return mainContext;
	}
	
	public Destination addRandomDestination(Context<Object> context, Geography<Object> geography, GeometryFactory gF, Geometry bndry, double destExtent, Color c, MathTransform ttM, MathTransform ttD) {
		
		Destination d = new Destination(geography, c);
		context.add(d);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = GeometryUtil.generateRandomPointsInPolygon(bndry, 1).get(0);
		
		// Get the coordinate and buffer it by the extent of the destination to create a 
		// circle that defines the destination
		Geometry destGeom = gF.createPoint(destCoord);//.buffer(destExtent);
		geography.move(d, destGeom);
				
		return d;
	}
	
	public Destination addUserDestination(Context<Object> context, Geography<Object> geography,GeometryFactory gF, String paramX, String paramY, int destExtent, Color c, MathTransform ttM, MathTransform ttD) {
		
		Destination d = new Destination(geography, c);
		context.add(d);
		
		// Get the x&y coords for the destination set by the user
		Parameters  params = RunEnvironment.getInstance().getParameters();
		double xCoord = (double)params.getInteger(paramX);
		double yCoord = (double)params.getInteger(paramY);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = new Coordinate(xCoord, yCoord);
		
		// Get the coordinate and buffer it by the extent of the destination to create a 
		// circle that defines the destination
		Geometry destGeom = gF.createPoint(destCoord).buffer(destExtent);
		
		geography.move(d, destGeom);
		
		return d;
		
	}
	
    public Ped addPed(Context context, Geography geography, GeometryFactory gF, double pedAngle, Coordinate coord, Destination d, Color c) throws MismatchedDimensionException, TransformException {
        
        // Instantiate a new pedestrian agent and add the agent to the context
        Ped newPed = new Ped(geography, gF, pedAngle, d, c);
        context.add(newPed);
        
        // Create a new point geometry. Move the pedestrian to this point. In doing so this 
        // pedestrian agent becomes associated with this geometry.
		Point pt = gF.createPoint(coord);
		
		// Transform the coordinate so that the circle can be created using a radius in metres
		Point ptCalc = (Point)JTS.transform(pt, transformToCalc);
		Geometry circle = ptCalc.buffer(newPed.getRad());
		moveAgentToCalculationGeometry(geography, circle, newPed);
		//geography.move(newPed, circle);
		
		// Set the angle to the destination and point the pedestrian in the direction of that direction.
		double a0 = newPed.seta0FromDestinationCoord();
		newPed.setaP(a0);
        	
        return newPed;
    }
	
	/*
	 * Taken from the Geography RS example. Returns a list of SimpleFeature
	 * objects from the shapefile path passed to the function.
	 */
	private List<SimpleFeature> loadFeaturesFromShapefile(String filename){
		URL url = null;
		try {
			url = new File(filename).toURL();
		} catch (MalformedURLException e1) {
			e1.printStackTrace();
		}

		List<SimpleFeature> features = new ArrayList<SimpleFeature>();
		
		// Try to load the shapefile
		SimpleFeatureIterator fiter = null;
		ShapefileDataStore store = null;
		store = new ShapefileDataStore(url);

		try {
			fiter = store.getFeatureSource().getFeatures().features();

			while(fiter.hasNext()){
				features.add(fiter.next());
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		finally{
			fiter.close();
			store.dispose();
		}
		
		return features;
	}
	
	/**
	 * Loads features from the specified shapefile.  The appropriate type of agents
	 * will be created depending on the geometry type in the shapefile (point, 
	 * line, polygon).
	 * 
	 * @param filename the name of the shapefile from which to load agents
	 * @param context the context
	 * @param geography the geography
	 */
	private List<Destination> loadFeatures (String filename, Context context, Geography geography){

		List<SimpleFeature> features = loadFeaturesFromShapefile(filename);
		List<Destination> destinations = new ArrayList<Destination>();
		
		// For each feature in the file
		for (SimpleFeature feature : features){
			Geometry geom = (Geometry)feature.getDefaultGeometry();
			Object agent = null;

			if (!geom.isValid()){
				System.out.println("Invalid geometry: " + feature.getID());
			}
			

			// For Points, create Destination agents
			if (geom instanceof Point){
				geom = (Point)feature.getDefaultGeometry();		
				
				agent = new Destination(geography, Color.RED);
				
				destinations.add((Destination)agent);
								
			}

			if (agent != null){
				context.add(agent);
				geography.move(agent, geom);
			}
			else{
				System.out.println("Error creating agent for  " + geom);
			}
		}
		
		return destinations;
	}
	
	static Geometry getGeometryForCalculation(Geography G, Object agent) throws MismatchedDimensionException, TransformException {
		Geometry geom = G.getGeometry(agent);
		
		return JTS.transform(geom, transformToCalc);
	}
	
	static void moveAgentToCalculationGeometry(Geography G, Geometry geomCalc, Object agent) throws MismatchedDimensionException, TransformException {
		G.move(agent, JTS.transform(geomCalc, transformToGeog));
	}
	
	/**
	 * Taken from Nick Malleson's RepastCity.
	 * 
	 * Nice generic function :-) that reads in objects from shapefiles.
	 * <p>
	 * The objects (agents) created must extend FixedGeography to guarantee that they will have a setCoords() method.
	 * This is necessary because, for simplicity, geographical objects which don't move store their coordinates
	 * alongside the projection which stores them as well. So the coordinates must be set manually by this function once
	 * the shapefile has been read and the objects have been given coordinates in their projection.
	 * 
	 * @param <T>
	 *            The type of object to be read (e.g. PecsHouse). Must exted
	 * @param cl
	 *            The class of the building being read (e.g. PecsHouse.class).
	 * @param shapefileLocation
	 *            The location of the shapefile containing the objects.
	 * @param geog
	 *            A geography to add the objects to.
	 * @param context
	 *            A context to add the objects to.
	 * @throws MalformedURLException
	 *             If the location of the shapefile cannot be converted into a URL
	 * @throws FileNotFoundException
	 *             if the shapefile does not exist.
	 * @see FixedGeography
	 */
	public static <T extends FixedGeography> void readShapefile(Class<T> cl, String shapefileLocation,
			Geography<T> geog, Context<T> context) throws MalformedURLException, FileNotFoundException {
		File shapefile = null;
		
		// ShapefileLoader is useful when the attributes of the features of the shapefile are to be assigned to agents
		// Might need to take care to ensure the agents has an attribute/method that matches the attributes of the shapefile features.
		ShapefileLoader<T> loader = null;
		shapefile = new File(shapefileLocation);
		if (!shapefile.exists()) {
			throw new FileNotFoundException("Could not find the given shapefile: " + shapefile.getAbsolutePath());
		}
		loader = new ShapefileLoader<T>(cl, shapefile.toURI().toURL(), geog, context);
		while (loader.hasNext()) {
			loader.next();
		}
		for (T obj : context.getObjects(cl)) {
			obj.setCoords(geog.getGeometry(obj).getCentroid().getCoordinate());
		}
}

}
