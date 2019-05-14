package repastInterSim.main;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

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
import com.vividsolutions.jts.geom.Polygon;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.gis.util.GeometryUtil;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.space.graph.Network;
import repast.simphony.util.collections.IndexedIterable;
import repastInterSim.agent.Ped;
import repastInterSim.environment.Destination;
import repastInterSim.environment.FixedGeography;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.ShapefileLoader;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.RoadLinkContext;

public class SpaceBuilder extends DefaultContext<Object> implements ContextBuilder<Object> {
	
	public static double spaceScale = 1;
	public static double[] north = {0,1}; // Defines north, against which bearings are taken
	
	public static Context<RoadLink> roadLinkContext;
	public static Geography<RoadLink> roadLinkGeography;
	
	public static Context<Junction> junctionContext;
	public static Network<Junction> roadNetwork;
	
	    /* (non-Javadoc)
	 * @see repast.simphony.dataLoader.ContextBuilder#build(repast.simphony.context.Context)
	 * 
	 */
	@Override
	public Context<Object> build(Context<Object> context) {
	    
		context.setId(GlobalVars.CONTEXT_NAMES.MAIN_CONTEXT);
	   
		// Initiate geographic spaces
		GeographyParameters<Object> geoParams = new GeographyParameters<Object>();
		GeometryFactory fac = new GeometryFactory();
		Geography<Object> geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY, context, geoParams);
		geography.setCRS(GlobalVars.geographyCRSString);
		context.add(geography);
		
		
		// Road link geography is used to create the road network projection
		GeographyParameters<RoadLink> roadLinkGeoParams = new GeographyParameters<RoadLink>();
		roadLinkContext = new RoadLinkContext();
		roadLinkGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.ROAD_LINK_GEOGRAPHY, roadLinkContext, roadLinkGeoParams);
		roadLinkGeography.setCRS(GlobalVars.geographyCRSString);

		// Junction geography also used to create the road network
		GeographyParameters<Junction> junctionGeoParams = new GeographyParameters<Junction>();
		junctionContext = new JunctionContext();
		Geography<Junction> junctionGeography = GeographyFactoryFinder.createGeographyFactory(null).createGeography(GlobalVars.CONTEXT_NAMES.JUNCTION_GEOGRAPHY, junctionContext, junctionGeoParams);
		junctionGeography.setCRS(GlobalVars.geographyCRSString);
		context.addSubContext(junctionContext);
		
		
	    // Load agents from shapefiles
		try {
			
			// Build the fixed environment
			
			// 1. Load destinations
			String destinationsFile = GlobalVars.GISDataDir + GlobalVars.DestinationsFile;
			readShapefile(Destination.class, destinationsFile, geography, context);
			
			// 2. Load roads
			String vehicleRoadFile = GlobalVars.GISDataDir + GlobalVars.VehicleRoadShapefile;
			String pedestrianRoadFile = GlobalVars.GISDataDir + GlobalVars.PedestrianRoadShapefile;
			readShapefile(Road.class, vehicleRoadFile, geography, context);
			readShapefile(Road.class, pedestrianRoadFile, geography, context);
			
			// 3. Load pedestrian obstruction boundaries
			String pedObstructionFile = GlobalVars.GISDataDir + GlobalVars.PedestrianObstructionShapefile;
			readShapefile(PedObstruction.class, pedObstructionFile, geography, context);
			
			
			// Build the road network
			
			// 1. Load the road links
			String roadLinkFile = GlobalVars.GISDataDir + GlobalVars.RoadLinkShapefile;
			readShapefileWithType(RoadLink.class, roadLinkFile, roadLinkGeography, roadLinkContext);
			SpatialIndexManager.createIndex(roadLinkGeography, RoadLink.class);

			
			// 2. roadNetwork
			NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>(GlobalVars.CONTEXT_NAMES.ROAD_NETWORK,junctionContext, false);
			builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
			roadNetwork = builder.buildNetwork();
			GISFunctions.buildGISRoadNetwork(roadLinkGeography, junctionContext,junctionGeography, roadNetwork);
			
		} catch (MalformedURLException | FileNotFoundException | MismatchedDimensionException e1 ) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		// Set the internal context and geography attributes of the destination agents
		IndexedIterable<Object> destinations = context.getObjects(Destination.class);
		for (Object d : destinations) {
			((Destination)d).setContext(context);
			((Destination)d).setGeography(geography);
		}
		

		
    	// Get the number of pedestrian agent to add to the space from the parameters
    	Parameters params = RunEnvironment.getInstance().getParameters();
    	int nP = (int)params.getInteger("nPeds");
    	
    	String startingZonesFile = GlobalVars.GISDataDir + GlobalVars.StartingZonesFile;
		List<Coordinate> agentCoords = getRandomCoordinatesWithinStartingZones(startingZonesFile,  fac,  nP);
		
		
		// Create the pedestrian agents
		int i = 0;
		int destinationIndex = 0;
		for (Coordinate coord : agentCoords) {
    		
    		// Crude way to assign different destinations
    		if (i > nP / 2) {
    			destinationIndex = 1; // i
    		}
			
    		Ped newPed = addPed(context, geography, fac, coord, (Destination)destinations.get(destinationIndex));
    		i+=1;
		}
		
		return context;
		
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
	public Destination addRandomDestination(Context<Object> context, Geography<Object> geography, GeometryFactory gF, Geometry bndry) {
		
		Destination d = new Destination();
		context.add(d);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = GeometryUtil.generateRandomPointsInPolygon(bndry, 1).get(0);
		
		Geometry destGeom = gF.createPoint(destCoord);
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
	public Destination addUserDestination(Context<Object> context, Geography<Object> geography, GeometryFactory gF, String paramX, String paramY) {
		
		Destination d = new Destination();
		context.add(d);
		
		// Get the x&y coords for the destination set by the user
		Parameters  params = RunEnvironment.getInstance().getParameters();
		double xCoord = (double)params.getInteger(paramX);
		double yCoord = (double)params.getInteger(paramY);
		
		// Initialize random coordinates for the destination
		Coordinate destCoord = new Coordinate(xCoord, yCoord);
		
		Geometry destGeom = gF.createPoint(destCoord);
		
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
    public <T> Ped addPed(Context<Object> context, Geography<Object> geography, GeometryFactory gF, Coordinate coord, Destination d)  {
        
        // Instantiate a new pedestrian agent and add the agent to the context
        Ped newPed = new Ped(geography, gF, d);
        context.add(newPed);
        
        // Create a new point geometry.
		Point pt = gF.createPoint(coord);
		
		// Transform the coordinate so that the circle can be created using a radius in metres
		Geometry circle = pt.buffer(newPed.getRad());
		
		// Move the pedestrian to this geometry
		moveAgentToGeometry(geography, circle, newPed);
		
		// Set the location attribute of the pedestrian agent to be its current location. Simplifies future calculations
		newPed.setLoc();
		
		// Set the angle to the destination and point the pedestrian in the direction of that direction.
		Coordinate dLoc = getAgentGeometry(geography, d).getCoordinate();
		double a0 = newPed.setDirectionFromDestinationCoord(dLoc);
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
				
				agent = new Destination();
				
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
	
	public static <T> Geometry getAgentGeometry(Geography<T> G, Object agent) {
		Geometry geom = G.getGeometry(agent);
		return geom;
	}
	
	public static <T> void moveAgentToGeometry(Geography<T> G, Geometry geomCalc, T agent) {
		G.move(agent, geomCalc);
	}
	
	/**
	 * This function was taken from Nick Malleson. I have edited it so that generic type contexts and geographies can 
	 * be passed to the function. His comments below.
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
	 * @throws TransformException 
	 * @throws MismatchedDimensionException 
	 * @see FixedGeography
	 */
	public static <T extends FixedGeography> void readShapefile(Class<T> cl, String shapefileLocation,
		Geography<Object> geog, Context<Object> context) throws MalformedURLException, FileNotFoundException  {
		File shapefile = null;
		ShapefileLoader<T> loader = null;
		shapefile = new File(shapefileLocation);
		if (!shapefile.exists()) {
			throw new FileNotFoundException("Could not find the given shapefile: " + shapefile.getAbsolutePath());
		}
		loader = new ShapefileLoader<T>(cl, shapefile.toURI().toURL(), geog, context);
		while (loader.hasNext()) {
			loader.next();
		}
		for (Object obj : context.getObjects(cl)) {
			// Warning of unchecked type cast below should be ok since only objects of this type were selected from the context
			((T)obj).setGeom(getAgentGeometry(geog, obj));
		}
	}
	
	/**
	 * This function was taken from Nick Malleson.
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
	 * @throws TransformException 
	 * @throws MismatchedDimensionException 
	 * @see FixedGeography
	 */
	public static <T extends FixedGeography> void readShapefileWithType(Class<T> cl, String shapefileLocation,
		Geography<T> geog, Context<T> context) throws MalformedURLException, FileNotFoundException {
		File shapefile = null;
		ShapefileLoader<T> loader = null;
		shapefile = new File(shapefileLocation);
		if (!shapefile.exists()) {
			throw new FileNotFoundException("Could not find the given shapefile: " + shapefile.getAbsolutePath());
		}
		loader = new ShapefileLoader<T>(cl, shapefile.toURI().toURL(), geog, context);
		while (loader.hasNext()) {
			loader.next();
		}
		for (Object obj : context.getObjects(cl)) {
			// Warning of unchecked type cast below should be ok since only objects of this type were selected from the context
			((T)obj).setGeom(getAgentGeometry(geog, obj));
		}
	}
	
	public List<Coordinate> getRandomCoordinatesWithinRoads(Context<Object> c, Geography<Object> g, GeometryFactory fac, Integer nPoints){
		
		IndexedIterable<Object> agents = c.getObjects(Road.class);
		Polygon[] roadPolygons = new Polygon[agents.size()];

		
		// Iterate over the agents and get their polygon geometry
		int i = 0;
		for (Object a : agents) {
			Polygon p = (Polygon)g.getGeometry(a);
			roadPolygons[i] = p;
			i++;
		}
		
		// Create single MultiPolygon that includes all road polygons and generate random coordinates that are within this multipolygon
	    MultiPolygon combined = new MultiPolygon(roadPolygons, fac);
		List<Coordinate> randCoords = GeometryUtil.generateRandomPointsInPolygon(combined, nPoints);
		
		return randCoords;
	}
	
	public List<Coordinate> getRandomCoordinatesWithinStartingZones(String startingZonesFile, GeometryFactory fac, Integer nPoints){
		
		// Load the starting zones
		List<SimpleFeature> startingZones = loadFeaturesFromShapefile(startingZonesFile);
		Polygon[] startingPolygons = new Polygon[startingZones.size()];
		
		// Iterate over the agents and get their polygon geometry
		int i = 0;
		for (SimpleFeature sf: startingZones) {
			MultiPolygon mp = (MultiPolygon)sf.getDefaultGeometry();
			startingPolygons[i] = (Polygon)mp.getGeometryN(0);
			i++;
		}
		
		// Create single MultiPolygon that includes all road polygons and generate random coordinates that are within this multipolygon
	    MultiPolygon combined = new MultiPolygon(startingPolygons, fac);
		List<Coordinate> randCoords = GeometryUtil.generateRandomPointsInPolygon(combined, nPoints);
		
		return randCoords;
	}
	
	public static int sizeOfIterable(Iterable<?> i) {
		int size = 0;
		Iterator<?> it = i.iterator();
		while (it.hasNext()) {
			size++;
			it.next();
		}
		return size;
	}

}
