package repastSocialForce;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.opengis.feature.simple.SimpleFeature;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;

import repast.simphony.context.Context;
import repast.simphony.context.DefaultContext;
import repast.simphony.context.space.gis.GeographyFactoryFinder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduleParameters;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.gis.util.GeometryUtil;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.gis.GeographyParameters;
import repast.simphony.util.ContextUtils;

public class SpaceBuilder extends DefaultContext<Object> implements ContextBuilder<Object> {
	
	static double spaceScale = 1;
	static double[] north = {0,1}; // Defines north, against which bearings are taken
	
	    /* (non-Javadoc)
	 * @see repast.simphony.dataLoader.ContextBuilder#build(repast.simphony.context.Context)
	 * 

	 */
	@Override
	public Context<Object> build(Context<Object> context) {
	    
		context.setId("repastSocialForce");
	   
		// Initiate geographic spaces
		GeographyParameters geoParams = new GeographyParameters();
		
		// Use GB Coordinate projection
		geoParams.setCrs("EPSG:27700");
		Geography<Object> geography = GeographyFactoryFinder.createGeographyFactory(null).createGeography("Geography", context, geoParams);
		
		String checkCRS = geoParams.getCrs();
		context.add(geography);
	    
		GeometryFactory fac = new GeometryFactory();
	    
	    // Code below taken from the 'Geography' repast example
		// Create an area in which to create agents.  This border is loaded from a shapefile.
		String boundaryFilename = ".//data//JunctClip.shp";
		List<SimpleFeature> features = loadFeaturesFromShapefile(boundaryFilename);
		Geometry boundary = (MultiPolygon)features.iterator().next().getDefaultGeometry();
		
	    
	    // A separate class is used to handle the creation of pedestrians
	    //Destination d1 = addRandomDestination(context, geography, fac, boundary, 0.01, Color.RED);
	    //Destination d2 = addRandomDestination(context, geography, fac, boundary, 0.01, Color.RED);
	    //Destination d1 = addUserDestination(context, geography, fac, "destX1", "destY1", 5);
	    //Destination d2 = addUserDestination(context, geography, fac, "destX2", "destY2", 5);
		List<Destination> destinations = loadFeatures (".//data//destCoords.shp", context, geography);
	    
	    // Calculate the distance between destinations to better understand units
	    Geometry g1 = (Point)geography.getGeometry(destinations.get(0));
	    Geometry g2 = (Point)geography.getGeometry(destinations.get(1));
	    
	    double dist = g2.distance(g1);
	    
	    
    	// Get the number of pedetrian agent to add to the space from the parameters
    	Parameters params = RunEnvironment.getInstance().getParameters();
    	int nP = (int)params.getInteger("nPeds");

		// Generate random points in the area to create agents.
		List<Coordinate> agentCoords = GeometryUtil.generateRandomPointsInPolygon(boundary, nP);
		
		// Create the agents from the collection of random coords.
		Random randCoord = new Random();
		for (Coordinate coord : agentCoords) {
			
			// Generate a random initial direction for the pedestrian
    		double randBearing = randCoord.nextFloat() * FastMath.PI * 2;
    		double[] dir = {FastMath.sin(randBearing), FastMath.cos(randBearing)};
			
    		Ped newPed = addPed(context, geography, fac, dir,coord, destinations.get(0), Color.BLUE);
		}
		
		return context;
	}
	
	public Destination addRandomDestination(Context<Object> context, Geography<Object> geography, GeometryFactory gF, Geometry bndry, double destExtent, Color c) {
		
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
	
	public Destination addUserDestination(Context<Object> context, Geography<Object> geography,GeometryFactory gF, String paramX, String paramY, int destExtent, Color c) {
		
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
	
    public Ped addPed(Context context, Geography geography, GeometryFactory gF, double[] direction, Coordinate coord, Destination d, Color c) {
        
        // Instantiate a new pedestrian agent and add the agent to the context
        Ped newPed = new Ped(geography, direction, d, c);
        context.add(newPed);
        
        // Create a new point geometry. Move the pedestrian to this point. In doing so this 
        // pedestrian agent becomes associated with this geometry.
		Point geom = gF.createPoint(coord);
		geography.move(newPed, geom);
        	
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
	
}
