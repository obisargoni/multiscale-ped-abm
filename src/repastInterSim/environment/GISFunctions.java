package repastInterSim.environment;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.geometry.Envelope2D;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
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
import repast.simphony.space.gis.WritableGridCoverage2D;
import repast.simphony.space.graph.Network;
import repast.simphony.util.collections.IndexedIterable;
import repastInterSim.main.SpaceBuilder;


public class GISFunctions {
	
	private static volatile GridEnvelopeGeometryCache gridEnvelopeGeometryCache = null;
	
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
			Coordinate cFirst = roadGeom.getCoordinates()[0]; // First coord
			Coordinate cLast = roadGeom.getCoordinates()[roadGeom.getNumPoints() - 1]; // Last coord

			// Create Junctions from these coordinates and add them to the Geography (if they haven't been
			// created already)
			Junction juncFirst, juncLast;
			if (coordMap.containsKey(cFirst)) {
				// A Junction with those coordinates (cFirst) has been created, get it so we can add an edge to it
				juncFirst = coordMap.get(cFirst);
			} else { // Junction does not exit
				juncFirst = new Junction(roadLink.getMNodeFID()); // The minus node fid corresponds to the junction at the first coord of the road link
				Point pF = geomFac.createPoint(cFirst);
				juncFirst.setGeom(pF);
				junctionContext.add(juncFirst);
				coordMap.put(cFirst, juncFirst);
				SpaceBuilder.moveAgentToGeometry(junctionGeography, pF, juncFirst);
			}
			if (coordMap.containsKey(cLast)) {
				juncLast = coordMap.get(cLast);
			} else { // Junction does not exit
				juncLast = new Junction(roadLink.getPNodeFID()); // The plus node fid corresponds to the junction at the last coord of the road link
				Point pL = geomFac.createPoint(cLast);
				juncLast.setGeom(pL);
				junctionContext.add(juncLast);
				coordMap.put(cLast, juncLast);
				SpaceBuilder.moveAgentToGeometry(junctionGeography, pL, juncLast);
			}
			// Tell the road object who it's junctions are
			roadLink.addJunction(juncFirst);
			roadLink.addJunction(juncLast);
			// Tell the junctions about this roadLink
			juncFirst.addRoadLink(roadLink);
			juncLast.addRoadLink(roadLink);

			// Create an edge between the two junctions, assigning a weight equal to it's length
			String direction = roadLink.getDirection();
			NetworkEdge<Junction> edge = null;
			if (direction.equals("+")) {
				edge = new NetworkEdge<Junction>(juncFirst, juncLast, true, roadGeom.getLength(), null);
			}
			else if (direction.equals("-")) {
				edge = new NetworkEdge<Junction>(juncLast, juncFirst, true, roadGeom.getLength(), null);
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
		Geography<T> geography, Context<T> context) throws MalformedURLException, FileNotFoundException  {
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
		for (T obj : context.getObjects(cl)) {
			// Warning of unchecked type cast below should be ok since only objects of this type were selected from the context
			(obj).setGeom(SpaceBuilder.getAgentGeometry(geography, obj));
		}
	}
	
	/**
	 * Iterates over the geometries in a Geography Projection and expands an envelope such that 
	 * the envelope matches the extent of the geometries in the Geography.
	 * 
	 * Restricted to use with FixedGeography type Geographies since these contain object that do not 
	 * change geometry.
	 * 
	 * @param <T> 
	 * 			The class of objects contained in the geography
	 * @param geography 
	 * 			The repast simphony geography projection containing the geometries to return an envelope over
	 * @return Envelope matching the extent of the geometries in the geography.
	 */
	public static <T extends FixedGeography> ReferencedEnvelope getGeographyEnvelope(Geography<T> geography) {
		ReferencedEnvelope rEnv = new ReferencedEnvelope();

		for (T obj : geography.getAllObjects()) {
			Geometry g = obj.getGeom();
			rEnv.expandToInclude(g.getEnvelopeInternal());
		}
		
		return rEnv;
	}
	
	/**
	 * Iterates over the geometries in a Geography Projection and expands an envelope such that 
	 * the envelope matches the extent of the geometries in the Geography.
	 * 
	 * Restricted to use with FixedGeography type Geographies since these contain object that do not 
	 * change geometry.
	 * 
	 * @param <T> 
	 * 			The class of objects contained in the geography
	 * @param geography 
	 * 			The repast simphony geography projection containing the geometries to return an envelope over
	 * @param Envelope
	 * 			Envelope to expand to include the geometries in the geography
	 * @return Envelope matching the extent of the geometries in the geography
	 */
	public static <T extends FixedGeography> ReferencedEnvelope expandEnvelopeToGeography(Geography<T> geography, ReferencedEnvelope rEnv) {
		
		for (T obj : geography.getAllObjects()) {
			Geometry g = obj.getGeom();
			rEnv.expandToInclude(g.getEnvelopeInternal());
		}
		
		return rEnv;
	}
	
	/**
	 * Iterate over multiple geographies to produce an envelope that matches the extent of the geometries in all geographies.
	 * 
	 * Geographies must have matching Coordinate Reference Systems.
	 * 
	 * @param <T>
	 * 			The class of objects contained in the geography
	 * @param geographies
	 * 			IndexedIterable of repast simphony geography projections
	 * @return Envelope
	 * 			Envelope matching the extent of the geometries in the geography
	 */
	public static <T extends FixedGeography> ReferencedEnvelope getMultipleGeographiesEnvelope(ArrayList<Geography> geographies) {
		ReferencedEnvelope rEnv = new ReferencedEnvelope();
		
		CoordinateReferenceSystem firstCRS = geographies.iterator().next().getCRS();
		
		for(Geography<T> g : geographies) {
			// Throw exception if geography coordinate reference systems don't all match
			assert firstCRS.equals(g.getCRS());
			
			rEnv = expandEnvelopeToGeography(g, rEnv);
		}
		return rEnv;
		
	}
	
	/**
	 * Sets the values of grid coverage cells according to an attribute of objects contained within a geography. Grid cells that
	 * intersect with a geography object have their value set according to the objects attribute, the name of which is given as an input to 
	 * the method.
	 * 
	 * @param <T> The type of the agent in the geography
	 * @param grid The grid coverage whose cell values are to be set
	 * @param cl The class of agent in the geography 
	 * @param geography The geography projection containing the agents used to set grid cell values
	 * @param attributeName The name of the agent attribute to use to set grid cell values
	 * @param agentAttributeValueMap A map between attribute value and grid cell value
	 */
	public static <T extends FixedGeography> void setGridCoverageValuesFromGeography(WritableGridCoverage2D grid, List<GridEnvelope2D> gridEnvelopeList, Class<T> cl, Geography<T> geography, String attributeName, Map<String,Integer> agentAttributeValueMap) {
		
	    // Get class attributeMethodMap read attribute methods
	    Map<String, Method> readAttributeMethodMap = ShapefileLoader.getAttributeMethodMap(cl, "r");
	    
	    Method readAttributeMethod = readAttributeMethodMap.get(attributeName);

		Iterable<T> Obs = geography.getAllObjects();
		
		for(GridEnvelope2D gridEnv: gridEnvelopeList) {
			
			GridCoordinates2D gridPos = new GridCoordinates2D(gridEnv.x,gridEnv.y);
			Polygon worldPoly = null;
			try {
				worldPoly = getGridEnveloplePolygon(grid, gridEnv);
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			
			for(T Ob: Obs) {
				Object attributeValue = null;
				String strAttrVal = null;
				try {
					attributeValue = readAttributeMethod.invoke(Ob);
					strAttrVal = attributeValue.toString();
				} catch (IllegalAccessException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalArgumentException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (InvocationTargetException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				// Choose GIS method based on attribute value
				if (strAttrVal.contentEquals("pedestrian")) {
					if(worldPoly.within((Ob.getGeom()))) {
						Integer gridValue = agentAttributeValueMap.get(strAttrVal);
						grid.setValue(gridPos, gridValue);
					}
				}
				else {
					if(worldPoly.intersects((Ob.getGeom()))) {
						Integer gridValue = agentAttributeValueMap.get(strAttrVal);
						grid.setValue(gridPos, gridValue);
					}	
				}
			}
		}
	}
	
	private static Polygon getGridEnveloplePolygon(WritableGridCoverage2D grid, GridEnvelope2D gridEnv) throws Exception {		
		if (gridEnvelopeGeometryCache == null) {
			gridEnvelopeGeometryCache = GridEnvelopeGeometryCache.getInstance(grid);
		} // if not cached
		return GridEnvelopeGeometryCache.get(grid, gridEnv);
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


/* ********************************************************************************** */
/**
 * Caches the lookup from a grid envelope to the corresponding geometry the envelope represents in the
 * world crs
 * 
 * <p>
 * This class can be serialised so that if the GIS data doesn't change it doesn't have to be re-calculated each time.
 * 
 * @author Nick Malleson
 */
class GridEnvelopeGeometryCache implements Serializable {

	private static Logger LOGGER = Logger.getLogger(GridEnvelopeGeometryCache.class.getName());

	private static final long serialVersionUID = 1L;
	private static Hashtable<GridEnvelope2D, Polygon> theCache; // The actual cache

	// The location that the serialised object might be found.
	private File serialisedLoc;
	// The time that this cache was created, can be used to check data hasn't
	// changed since
	private long createdTime;

	private GridEnvelopeGeometryCache(WritableGridCoverage2D grid)
			throws Exception {


		//this.serialisedLoc = serialisedLoc;
		GridEnvelopeGeometryCache.theCache = new Hashtable<GridEnvelope2D, Polygon>();

		LOGGER.log(Level.FINE, "GridEnvelopeGeometryCache() creating new cache from grid (and modification date):\n\t"
				//+ this.serialisedLoc.getAbsolutePath()
				);
		
		// Don't populate the cache for now, instead add to on the fly
		//populateCache(grid);
		this.createdTime = new Date().getTime();
		//serialise(); don't worry about serialising for now
	}

	public void clear() {
		this.theCache.clear();
	}

	
	private void populateCache(WritableGridCoverage2D grid)
			throws Exception {
		double time = System.nanoTime();
		theCache = new Hashtable<GridEnvelope2D, Polygon>();
		
		// Cast to int. 
		int width = grid.getRenderedImage().getWidth();
		int height = grid.getRenderedImage().getHeight();
		
		// Loop over coverage grid cells
		for(int i=0;i<width;i++) {
			for (int j=0;j<height;j++) {
				GridEnvelope2D gridEnv = new GridEnvelope2D(i, j,1, 1);
				Polygon gridPoly = getWorldPolygonFromGridEnvelope(grid, gridEnv);
				theCache.put(gridEnv, gridPoly);
			}
		}
		LOGGER.log(Level.FINER, "Finished caching grid envelope polygons (" + (0.000001 * (System.nanoTime() - time)) + "ms)");
	} // if nearestRoadCoordCache = null;

	/**
	 * 
	 * @param c
	 * @return
	 * @throws Exception
	 */
	public static Polygon get(WritableGridCoverage2D grid, GridEnvelope2D gE) throws Exception {
		if (gE == null) {
			throw new Exception("Route.GridEnvelopeGeometryCache.get() error: the given grid envelope is null.");
		}
		double time = System.nanoTime();
		Polygon gP = theCache.get(gE);
		if (gP != null) {
			LOGGER.log(Level.FINER, "GridEnvelopeGeometryCache.get() (using cache) - ("
					+ (0.000001 * (System.nanoTime() - time)) + "ms)");
			return gP;
		}
		// If get here then the grid envelope is not in the cache, get polygon on the fly
		gP = getWorldPolygonFromGridEnvelope(grid, gE);
		if (gP != null) {
			LOGGER.log(Level.FINER, "GridEnvelopeGeometryCache.get() (not using cache) - ("
					+ (0.000001 * (System.nanoTime() - time)) + "ms)");
			theCache.put(gE, gP);
			return gP;
		}
		/* IF HERE THEN ERROR, PRINT DEBUGGING INFO */
		StringBuilder debugIntro = new StringBuilder(); // Add in some extra infor for debugging
		debugIntro.append("Route.GridEnvelopeGeometryCache.get() error: couldn't find a coordinate to return.\\n");
		debugIntro.append(gE.x + " " + gE.y);
		throw new Exception(debugIntro.toString());
	}
	
	/**
	 * Takes an envelope expressed in grid coordinates and transforms it to the equivalent envelope 
	 * expressed in the coordinate reference system the grid maps to - the 'real world' space.
	 * @param grid
	 * 			The grid coverage the grid envelope belongs to
	 * @param gridEnvelope
	 * 			The grid envelope to transform
	 * @return
	 *			The polygon (a square) representing the same envelop in a geographical space
	 */
	public static Polygon getWorldPolygonFromGridEnvelope(WritableGridCoverage2D grid, GridEnvelope2D gridEnvelope) {
		
		Envelope2D worldEnv = null;
		
		// Transform grid envelope into envelope with coordinates in gis reference system
		try {
			worldEnv = grid.getGridGeometry().gridToWorld(gridEnvelope);
		} catch (TransformException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
					
		Coordinate[] coords = {
				new Coordinate(worldEnv.getMinX(), worldEnv.getMinY()),
				new Coordinate(worldEnv.getMinX(), worldEnv.getMaxY()),
				new Coordinate(worldEnv.getMaxX(), worldEnv.getMaxY()),
				new Coordinate(worldEnv.getMaxX(), worldEnv.getMinY()),
				new Coordinate(worldEnv.getMinX(), worldEnv.getMinY())
		};
		
		Polygon wEPoly = new GeometryFactory().createPolygon(coords);
		
		return wEPoly;
	}
	
	
	// Not sure how to make use of serialising
	private void serialise() throws IOException {
		double time = System.nanoTime();
		FileOutputStream fos = null;
		ObjectOutputStream out = null;
		try {
			if (!this.serialisedLoc.exists())
				this.serialisedLoc.createNewFile();
			fos = new FileOutputStream(this.serialisedLoc);
			out = new ObjectOutputStream(fos);
			out.writeObject(this);
			out.close();
		} catch (IOException ex) {
			if (serialisedLoc.exists()) {
				// delete to stop problems loading incomplete file next time
				serialisedLoc.delete();
			}
			throw ex;
		}
		LOGGER.log(Level.FINE, "... serialised GridEnvelopeGeometryCache to " + this.serialisedLoc.getAbsolutePath()
				+ " in (" + 0.000001 * (System.nanoTime() - time) + "ms)");
	}

	/**
	 * Used to create a new BuildingsOnRoadCache object. This function is used instead of the constructor directly so
	 * that the class can check if there is a serialised version on disk already. If not then a new one is created and
	 * returned.
	 * 
	 * @param grid
	 * @return
	 * @throws Exception
	 */
	public synchronized static GridEnvelopeGeometryCache getInstance(WritableGridCoverage2D grid) throws Exception {
		double time = System.nanoTime();
		// See if there is a cache object on disk. - not sure how this works atm so commenting out
		/*
		if (serialisedLoc.exists()) {
			FileInputStream fis = null;
			ObjectInputStream in = null;
			GridEnvelopeGeometryCache ncc = null;
			try {

				fis = new FileInputStream(serialisedLoc);
				in = new ObjectInputStream(fis);
				ncc = (GridEnvelopeGeometryCache) in.readObject();
				in.close();

				// Check that the cache is representing the correct data and the
				// modification dates are ok
				return ncc;
			} catch (IOException ex) {
				if (serialisedLoc.exists())
					serialisedLoc.delete(); // delete to stop problems loading incomplete file next tinme
				throw ex;
			} catch (ClassNotFoundException ex) {
				if (serialisedLoc.exists())
					serialisedLoc.delete();
				throw ex;
			}
		}
		*/
		// No serialised object, or got an error when opening it, just create a new one
		return new GridEnvelopeGeometryCache(grid);
	}

}
