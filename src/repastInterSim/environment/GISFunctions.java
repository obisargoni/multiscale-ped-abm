package repastInterSim.environment;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridCoverage2D;
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
import repastInterSim.exceptions.RoutingException;
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
	public static <T extends FixedGeography> void setGridCoverageValuesFromGeography(WritableGridCoverage2D grid, Class<T> cl, Geography<T> geography, String attributeName, Map<String,Integer> agentAttributeValueMap) {
		
	    // Get class attributeMethodMap read attribute methods
	    Map<String, Method> readAttributeMethodMap = ShapefileLoader.getAttributeMethodMap(cl, "r");
	    
	    Method readAttributeMethod = readAttributeMethodMap.get(attributeName);

		Iterable<T> Obs = geography.getAllObjects();
		
		// Cast to int. 
		int width = grid.getRenderedImage().getWidth();
		int height = grid.getRenderedImage().getHeight();
		
		// Loop over coverage grid cells
		for(int i=0;i<width;i++) {
			for (int j=0;j<height;j++) {
						
				GridCoordinates2D gridPos = new GridCoordinates2D(i,j);
				GridEnvelope2D gridEnv = new GridEnvelope2D(gridPos.x, gridPos.y,1, 1);
				Polygon worldPoly = getWorldPolygonFromGridEnvelope(grid, gridEnv);
				
				List<T> intersectingObs = SpatialIndexManager.findIntersectingObjects(geography, worldPoly);
				if (intersectingObs.size() == 0) {
					continue;
				}
				
				// For each of the objects that intersect the geometry get the
				// string attribute values. Set grid cell value based on these
				List<String> obsAttributeValues = new ArrayList<String>();
				for(T Ob: intersectingObs) {
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
					
					if(!obsAttributeValues.contains(strAttrVal)) {
						obsAttributeValues.add(strAttrVal);
					}
				}
				
				// Choose GIS method based on attribute value
				Integer gridValue = null;
				
				// First check for ped obstructions, set value based on this
				if (obsAttributeValues.contains("pedestrian_obstruction")) {
					gridValue = agentAttributeValueMap.get("pedestrian_obstruction");
				}
				// Next check if grid cell intersects with a vehicle priority geom
				else if (obsAttributeValues.contains("vehicle")){
					gridValue = agentAttributeValueMap.get("vehicle");
				}
				// Next check if grid cell intersects with a pedestrian crossing geom
				else if (obsAttributeValues.contains("pedestrian_crossing")){
					gridValue = agentAttributeValueMap.get("pedestrian_crossing");
				}
				// Finally, can set grid value to be pedestrian if no other types intersect
				else {
					gridValue = agentAttributeValueMap.get("pedestrian");
				}
				grid.setValue(gridPos, gridValue);
			} // j
		}//i
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
		GridCoverage2D grd = (GridCoverage2D)grid;
		
		Polygon wEPoly = getWorldPolygonFromGridEnvelope(grd, gridEnvelope);
		
		return wEPoly;
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
	public static Polygon getWorldPolygonFromGridEnvelope(GridCoverage2D grid, GridEnvelope2D gridEnvelope) {
		
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
	
    /**
     *  Gets list of roads the input polygon intersects with.
     *  
     *  List should have one, where road is one-way, or two, where road is two way, RoadLink objects in.
     * @return
     * 		Road Link the agent is on
     * @throws RoutingException 
     */
	public static List<Road> getGridPolygonRoads(Polygon p) throws RoutingException {
		Road r = null;
		
		List<Road> containingRoads = SpatialIndexManager.findIntersectingObjects(SpaceBuilder.roadGeography, p, "contains");
    	
    	return containingRoads;
	}
	
	/**
	 * Gets list of unique road link IDs the input polygon intersects.
	 * 
	 * Typically set will contain only one entry, but polygon can overlap Roads that have different road link IDs
	 * @param p
	 * 		Polygon to get intersecting road links for
	 * @return
	 * 		List of unique intersecting road link IDs
	 * @throws RoutingException 
	 */
	public static List<String> getGridPolygonRoadLinkIDs(Polygon p) throws RoutingException{
		List<String> roadLinkIDs = new ArrayList<String>();
		
		List<Road> intersectingRoads = getGridPolygonRoads(p);
		
		// Gather the ID of road links the current and previous grid cells overlap with
		for(Road r:intersectingRoads) {
			if (!roadLinkIDs.contains(r.getRoadLinkFI())) {
				roadLinkIDs.add(r.getRoadLinkFI());
			}
		}
		return roadLinkIDs;
	}
    /**
     *  Gets Road the input coordinate intersects with.
     *  
     *  List should have one, where road is one-way, or two, where road is two way, RoadLink objects in.
     * @return
     * 		Road Link the agent is on
     * @throws RoutingException 
     */
	public static Road getCoordinateRoad(Coordinate c) throws RoutingException {
		Road r = null;
		
		List<Road> intersectingRoads = SpatialIndexManager.findIntersectingObjects(SpaceBuilder.roadGeography, c);
    	
    	if(intersectingRoads.size() == 0) {
    		// Method returns default value, null, if there are no intersecting roads
    	}
    	else if (intersectingRoads.size() == 1) {
        	r = intersectingRoads.get(0);
    	}
    	else {
    		throw new RoutingException("Input coordinate intersects with multiple road objects. Unexpected");
    	}
    	
    	return r;
	}
	
	/**
	 * Get the gis coordinate that corresponds to the location of the input Grid Coordinate
	 * in the coordinate reference system used by the grid coverage
	 * 
	 * @param grid
	 * 			The grid in which the grid cell sits
	 * @param cell
	 * 			The grid coordinate to get the gis coordinate of
	 * @return
	 * 			Coordinate. The gis coordinate
	 */
	public static Coordinate gridCellToCoordinate(GridCoverage2D grid, GridCoordinates2D cell) {
		double[] cellCoord = null;
		try {
			cellCoord = grid.getGridGeometry().gridToWorld(cell).getCoordinate();
		} catch (TransformException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Coordinate c = new Coordinate(cellCoord[0], cellCoord[1]);
		
		return c;
	}

}
