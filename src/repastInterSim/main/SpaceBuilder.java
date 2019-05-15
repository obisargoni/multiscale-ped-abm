package repastInterSim.main;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.util.Iterator;
import java.util.List;

import org.opengis.geometry.MismatchedDimensionException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

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
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.environment.contexts.JunctionContext;
import repastInterSim.environment.contexts.RoadLinkContext;

public class SpaceBuilder extends DefaultContext<Object> implements ContextBuilder<Object> {
	
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
			GISFunctions.readShapefile(Destination.class, destinationsFile, geography, context);
			
			// 2. Load roads
			String vehicleRoadFile = GlobalVars.GISDataDir + GlobalVars.VehicleRoadShapefile;
			String pedestrianRoadFile = GlobalVars.GISDataDir + GlobalVars.PedestrianRoadShapefile;
			GISFunctions.readShapefile(Road.class, vehicleRoadFile, geography, context);
			GISFunctions.readShapefile(Road.class, pedestrianRoadFile, geography, context);
			
			// 3. Load pedestrian obstruction boundaries
			String pedObstructionFile = GlobalVars.GISDataDir + GlobalVars.PedestrianObstructionShapefile;
			GISFunctions.readShapefile(PedObstruction.class, pedObstructionFile, geography, context);
			
			
			// Build the road network
			
			// 1. Load the road links
			String roadLinkFile = GlobalVars.GISDataDir + GlobalVars.RoadLinkShapefile;
			GISFunctions.readShapefileWithType(RoadLink.class, roadLinkFile, roadLinkGeography, roadLinkContext);
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
		List<Coordinate> agentCoords = GISFunctions.getRandomCoordinatesWithinShapeFileGeometries(startingZonesFile,  fac,  nP);
		
		
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

}
