package repastInterSim.environment;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateSequence;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;

import repast.simphony.context.Context;
import repast.simphony.gis.util.GeometryUtil;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.Network;
import repast.simphony.util.collections.IndexedIterable;
import repastInterSim.main.SpaceBuilder;


public class GISFunctions {
	
	/**
	 * Create the road network. Runs through the roads in the <code>roadGeography</code> and, for each one, will create
	 * <code>Junction</code> objects at their end points and an edge linking them. The <code>Junction</code> objects are
	 * added to the given <code>Geography</code> (so that we know where they are spatially) and they are also
	 * added, along with the edge between them, to the <code>junctionNetwork</code> so that topographical relationships
	 * can be established. (The <code>junctionNetwork</code> is part of the <code>Context</code>
	 * 
	 * @param roadGeography
	 * @param Context
	 * @param Geography
	 * @param roadNetwork
	 * @throws TransformException 
	 * @throws MismatchedDimensionException 
	 */
	public static void buildGISRoadNetwork(Geography<RoadLink> roadLinkGeography, Context<Junction> junctionContext,
			Geography<Junction> junctionGeography, Network<Junction> roadNetwork) {

		// Create a GeometryFactory so we can create points/lines from the junctions and roads
		// (this is so they can be displayed on the same display to check if the network has been created successfully)
		GeometryFactory geomFac = new GeometryFactory();

		// Create a cache of all Junctions and coordinates so we know if a junction has already been created at a
		// particular coordinate
		Map<Coordinate, Junction> coordMap = new HashMap<Coordinate, Junction>();
		
		// Iterate through all roads
		// Iterate through all roads
		Iterable<RoadLink> roadIt = roadLinkGeography.getAllObjects();
		for (RoadLink roadLink : roadIt) {
			
			// Create a LineString from the road so we can extract coordinates
			Geometry roadGeom = roadLink.getGeom();
			Coordinate c1 = roadGeom.getCoordinates()[0]; // First coord
			Coordinate c2 = roadGeom.getCoordinates()[roadGeom.getNumPoints() - 1]; // Last coord

			// Create Junctions from these coordinates and add them to the Geography (if they haven't been
			// created already)
			Junction junc1, junc2;
			if (coordMap.containsKey(c1)) {
				// A Junction with those coordinates (c1) has been created, get it so we can add an edge to it
				junc1 = coordMap.get(c1);
			} else { // Junction does not exit
				junc1 = new Junction();
				Point p1 = geomFac.createPoint(c1);
				junc1.setGeom(p1);
				junctionContext.add(junc1);
				coordMap.put(c1, junc1);
				SpaceBuilder.moveAgentToGeometry(junctionGeography, p1, junc1);
			}
			if (coordMap.containsKey(c2)) {
				junc2 = coordMap.get(c2);
			} else { // Junction does not exit
				junc2 = new Junction();
				Point p2 = geomFac.createPoint(c2);
				junc2.setGeom(p2);
				junctionContext.add(junc2);
				coordMap.put(c2, junc2);
				SpaceBuilder.moveAgentToGeometry(junctionGeography, p2, junc2);
			}
			// Tell the road object who it's junctions are
			roadLink.addJunction(junc1);
			roadLink.addJunction(junc2);
			// Tell the junctions about this roadLink
			junc1.addRoadLink(roadLink);
			junc2.addRoadLink(roadLink);

			// Create an edge between the two junctions, assigning a weight equal to it's length
			String direction = roadLink.getDirection();
			NetworkEdge<Junction> edge = null;
			if (direction.equals("-")) {
				edge = new NetworkEdge<Junction>(junc1, junc2, true, roadGeom.getLength(), null);
			}
			else if (direction.equals("+")) {
				edge = new NetworkEdge<Junction>(junc2, junc1, true, roadGeom.getLength(), null);
			}
				

			// Tell the roadLink and the Edge about each other
			roadLink.setEdge(edge);
			edge.setRoad(roadLink);

			// Add the edge to the network
			if (!roadNetwork.containsEdge(edge)) {
				roadNetwork.addEdge(edge);
			} else {
				//LOGGER.severe("CityContext: buildRoadNetwork: for some reason this edge that has just been created "
					//	+ "already exists in the RoadNetwork.");
			}

		} // for roadLink:
	}
	
	/*
	 * Get a list of coordinates that are randomly distributed within the geometries associated with Road agents.
	 * 
	 * @param context
	 * 			The context the Road agents belong to
	 * @param geography
	 * 			The geography containing the road geometries
	 * @param fac
	 * 			The geometry factory to use when generating coordinates
	 * @param nPoints
	 * 			The number of coordinates to generate and return
	 * @returns
	 * 			A list of coordinates 
	 */
	public static List<Coordinate> getRandomCoordinatesWithinRoads(Context<Object> context, Geography<Object> geography, GeometryFactory fac, Integer nPoints){
		
		IndexedIterable<Object> agents = context.getObjects(Road.class);
		Polygon[] roadPolygons = new Polygon[agents.size()];

		
		// Iterate over the agents and get their polygon geometry
		int i = 0;
		for (Object a : agents) {
			Polygon p = (Polygon)geography.getGeometry(a);
			roadPolygons[i] = p;
			i++;
		}
		
		// Create single MultiPolygon that includes all road polygons and generate random coordinates that are within this multipolygon
	    MultiPolygon combined = new MultiPolygon(roadPolygons, fac);
		List<Coordinate> randCoords = GeometryUtil.generateRandomPointsInPolygon(combined, nPoints);
		
		return randCoords;
	}
	
	/*
	 * Taken from the Geography RS example. Returns a list of SimpleFeature
	 * objects from the shapefile path passed to the function.
	 * 
	 * @param filename
	 * 			The path to the shapefile to load features from
	 */
	private static List<SimpleFeature> loadFeaturesFromShapefile(String filename){
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
	private static List<Destination> loadFeatures (String filename, Context context, Geography geography){

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
	
	/*
	 * Get a list of coordinates that are randomly distributed within the geometries in a shapefile.
	 * 
	 * @param shapeFilePath
	 * 			The path of the shape file containing the geometries to generate random coordinates within
	 * @param fac
	 * 			The geometry factory to use when generating coordinates
	 * @param nPoints
	 * 			The number of coordinates to generate and return
	 * @returns
	 * 			A list of coordinates 
	 */
	public static List<Coordinate> getRandomCoordinatesWithinShapeFileGeometries(String shapeFilePath, GeometryFactory fac, Integer nPoints){
		
		// Load the starting zones
		List<SimpleFeature> startingZones = loadFeaturesFromShapefile(shapeFilePath);
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
	 * @param geography
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
		Geography<Object> geography, Context<Object> context) throws MalformedURLException, FileNotFoundException  {
		File shapefile = null;
		ShapefileLoader<T> loader = null;
		shapefile = new File(shapefileLocation);
		if (!shapefile.exists()) {
			throw new FileNotFoundException("Could not find the given shapefile: " + shapefile.getAbsolutePath());
		}
		loader = new ShapefileLoader<T>(cl, shapefile.toURI().toURL(), geography, context);
		while (loader.hasNext()) {
			loader.next();
		}
		for (Object obj : context.getObjects(cl)) {
			// Warning of unchecked type cast below should be ok since only objects of this type were selected from the context
			((T)obj).setGeom(SpaceBuilder.getAgentGeometry(geography, obj));
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
	 * @param geography
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
		Geography<T> geography, Context<T> context) throws MalformedURLException, FileNotFoundException {
		File shapefile = null;
		ShapefileLoader<T> loader = null;
		shapefile = new File(shapefileLocation);
		if (!shapefile.exists()) {
			throw new FileNotFoundException("Could not find the given shapefile: " + shapefile.getAbsolutePath());
		}
		loader = new ShapefileLoader<T>(cl, shapefile.toURI().toURL(), geography, context);
		while (loader.hasNext()) {
			loader.next();
		}
		for (Object obj : context.getObjects(cl)) {
			// Warning of unchecked type cast below should be ok since only objects of this type were selected from the context
			((T)obj).setGeom(SpaceBuilder.getAgentGeometry(geography, obj));
		}
	}
	
	
	/*
	 * Function to calculate the angle in radians between a vector that points north and the
	 * vector from one coordinate to another.
	 * 
	 * @param Coordinate c1
	 * @param Coordinate c2
	 * 
	 * @return double angle in radians between north and the vector from c1 to c2
	 */
	public static double bearingBetweenCoordinates(Coordinate c1, Coordinate c2) {
		
		// Get the vector between the coordinates and convert to a unit vector
		double [] u = new double[2];
		u[0] = c2.x - c1.x;
		u[1] = c2.y - c1.y;
		
		u = Vector.unitV(u);
		
		double a = Vector.angleBetweenNorthAndUnitVector(u);
		
		return a;
	}
	
	/*
	 * Creates a point geometry from a coordinate.
	 * 
	 * @param coord
	 * 			The coordinate to move the agent to
	 */
	public static Point pointGeometryFromCoordinate(Coordinate coord) {
	    // Now create new geometry at the location of the new centroid
	    Coordinate[] cArray = {coord};
	    GeometryFactory fac = new GeometryFactory();
	    CoordinateSequence cs = fac.getCoordinateSequenceFactory().create(cArray);
	    Point pt = new Point(cs, fac);
	    
	    return pt;
	}

}
