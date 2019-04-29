package repastInterSim.environment;

import java.util.HashMap;
import java.util.Map;

import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.context.Context;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.Network;
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
			Geography<Junction> junctionGeography, Network<Junction> roadNetwork) throws MismatchedDimensionException, TransformException {

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
				SpaceBuilder.moveAgentToCalculationGeometry(junctionGeography, p1, junc1);
			}
			if (coordMap.containsKey(c2)) {
				junc2 = coordMap.get(c2);
			} else { // Junction does not exit
				junc2 = new Junction();
				Point p2 = geomFac.createPoint(c2);
				junc2.setGeom(p2);
				junctionContext.add(junc2);
				coordMap.put(c2, junc2);
				SpaceBuilder.moveAgentToCalculationGeometry(junctionGeography, p2, junc2);
			}
			// Tell the road object who it's junctions are
			roadLink.addJunction(junc1);
			roadLink.addJunction(junc2);
			// Tell the junctions about this roadLink
			junc1.addRoadLink(roadLink);
			junc2.addRoadLink(roadLink);

			// Create an edge between the two junctions, assigning a weight equal to it's length
			NetworkEdge<Junction> edge = new NetworkEdge<Junction>(junc1, junc2, false, roadGeom.getLength(), null);

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

}
