package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.stream.Collectors;

import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.InvalidGridGeometryException;
import org.geotools.geometry.DirectPosition2D;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.LineString;

import repast.simphony.random.RandomHelper;
import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repast.simphony.space.graph.ShortestPath;
import repast.simphony.space.graph.UndirectedJungNetwork;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class GridRoute {
	
	protected GridCoverage2D grid;

	protected Coordinate origin;
	protected Coordinate destination;

	/*
	 * The route consists of a list of coordinates which describe how to get to the destination. Each coordinate might
	 * have an attached 'speed' which acts as a multiplier and is used to indicate whether or not the agent is
	 * travelling along a transport route (i.e. if a coordinate has an attached speed of '2' the agent will be able to
	 * get to the next coordinate twice as fast as they would do if they were walking). The current position incicate
	 * where in the lists of coords the agent is up to. Other attribute information about the route can be included as
	 * separate arrays with indices that match those of the 'route' array below.
	 */
	private int currentPosition;
	protected List<Coordinate> routeX;
	protected List<Double> routeSpeedsX;
	/*
	 * This maps route coordinates to their containing Road, used so that when travelling we know which road/community
	 * the agent is on. private
	 */
	protected List<RoadLink> roadsX;

	// Record which function has added each coord, useful for debugging
	protected List<String> routeDescriptionX;
	
	//private static Logger LOGGER = Logger.getLogger(GridRoute.class.getName());
	
	// Used to get grid cell summand value when running flood fill algorithm for routing. Single agent can produce routes from different costs, reflecting agent's changing perceptions of costs.
	private HashMap<Integer, Double> gridSummandPriorityMap;
	
	// The grid coordiantes of the agents' route
	private List<GridCoordinates2D> gridPath;
	
	// The filtered list of grid coordinates, with unrequired coordiantes removed
	private List<GridCoordinates2D> prunedGridPath;
	
	// The grid coordinates of crosisng points - used for visualisation
	private List<GridCoordinates2D> gridPathCrossings;	
	
	// Record cells/coordinates at which agent enters new road link
	private Map<Coordinate, GridCoordinates2D> routeCoordMap;
	private List<Coordinate> primaryRouteX;
	
	// Used to recrode the route of grid cells or coordinates, grouped by the road link they belong to
	private Map<GridCoordinates2D, List<GridCoordinates2D>> groupedGridPath;
	
	// Record the value of grid cells following flood fill (used when routing via a grid)
	private double[][] floodFillValues = null;
	
	// Sets whether to run flood fill on full grid or just a partial section of it
	private boolean partialFF = false;
	
	/**
	 * Create a blank grid route object
	 *
	 */
	public GridRoute() {
		this.initialiseRouteVectors();
	}

	
	/**
	 * Create a new grid route object
	 * 
	 * @param grid
	 * 		The grid coverage to use for path finding
	 * @param mA
	 * 		The mobile agent this route belongs to
	 * @param gSPM
	 * 		The map from integers used to indicate the road user priority of grid cells to the agents perceived cost of moving through those grid cells. 
	 * Used for routing on a grid.
	 * @param origin
	 * 		The origin coordinate of the route
	 * @param destination
	 * 		The destination coordinate of the route
	 * @param
	 * 		A boolean value indicating whether to consider only a partial area of the grid when producing the route.
	 */
	public GridRoute(GridCoverage2D grd, HashMap<Integer, Double> gSPM, Coordinate origin, Coordinate destination, boolean partial) {

		this.grid = grd;
		this.origin = origin;
		this.destination = destination;
		// TODO Auto-generated constructor stub
		this.gridSummandPriorityMap = gSPM;
		this.partialFF = partial;
		
		this.initialiseRouteVectors();

	}
	
	/**
	 * This method produces a path of grid cells from the agents current position to the destination set when creating a GridRoute instance. 
	 * The grid path is grouped into sections according to the road link the grid cells belong to. 
	 */
	public void setGroupedGridPath() {
		
		this.initialiseRouteVectors();
		
		// sequence of grid cell coordinates leading from agent' current position to end destination
		this.gridPath = getDijkstraGridPath(this.grid, this.origin, this.destination);
				
		// First grid cell coordinate is agent's first coordinate along that road link, so use to index first group of coordinates
		GridCoordinates2D primaryRouteCell = null;
		Coordinate primaryRouteCoord = null;
		String prevCellRoadLinkFID = null;
		String cellRoadLinkFID = null;
		
		// Identify the primary route coordinates as those coordinate along the path where the road link changes
		// Use each primary grid cell to index a new group of grid cell coordinates (the path along a particular road link)
		for (int i = 0; i < gridPath.size(); i++) {			
			GridCoordinates2D gridCell = gridPath.get(i);
			Coordinate cellCoord = GISFunctions.gridCellToCoordinate(grid, gridCell);
			
			Road cellRoad = null;
			try {
				cellRoad = GISFunctions.getCoordinateRoad(cellCoord);
			} catch (RoutingException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			try {
				cellRoadLinkFID = cellRoad.getRoadLinkID();
			} catch (NullPointerException e) {
				// The coordinate that a grid cell transforms to does not intersect a road
				// Grid cell can be included in path because this the coord not intersecting with a road shouldnt affect the pruning of grid cells to coords
				// But don't consider for primary route coordinate
				/*
		    	String msg = "GridRoute.setGroupedGridPath(): Grid cell coordinate does not intersect Road \n\r" + 
						"GridCell: " + gridCell.toString()+ "\n\r" +
						"cellCoord: " + cellCoord.toString() + "\n\r";
				//LOGGER.log(Level.SEVERE,msg);
				System.out.print(msg);
				*/
			}
			
			// Check for change in road link id
			// IDs will be the same unless cellRoadLinkID has been updated. Ensures that only coordinates that unambiguously belong to a road linj
			// are added as primary route coordinates
			if (!cellRoadLinkFID.equals(prevCellRoadLinkFID)) {
				primaryRouteCell = gridCell;
				primaryRouteCoord = cellCoord;
				this.routeCoordMap.put(primaryRouteCoord, primaryRouteCell);
				this.primaryRouteX.add(primaryRouteCoord);
				this.groupedGridPath.put(primaryRouteCell, new ArrayList<GridCoordinates2D>());
				prevCellRoadLinkFID = cellRoadLinkFID;
			}
			this.groupedGridPath.get(primaryRouteCell).add(gridCell);
		}
		
		checkGridRoute();
		
		// Add the destination to the list of primary route coordinates (those related to changing road link)
		this.primaryRouteX.add(this.destination);
	}
	
	private void addCoordinatesToRouteFromGridPath(List<GridCoordinates2D> gP) {
		
		Set<Integer> routeIndices = new HashSet<Integer>();
		Map<Integer, String> descriptionMap = new HashMap<Integer, String>();
		
		int[] cellValue = new int[1];
		
		// Search for primary crossings
		// Identify grid cells in the route that are located on primary road links
		for (int i = 1; i < gP.size(); i++) {
			GridCoordinates2D gridCell = gP.get(i);
			cellValue = grid.evaluate(gridCell, cellValue);
			Integer val = cellValue[0];
			
			// If grid cell is a primary road link, a primary crossing has occured in the route
			// Search forward and back in the route to find the crossing entry and exit points
			if (val.compareTo(GlobalVars.GRID_PARAMS.getPriorityValueMap().get("road_link")) == 0) {
				
				// Search forward
				for (int j = i; j<gP.size(); j++) {
					GridCoordinates2D crossingCell = gP.get(j);
					cellValue = grid.evaluate(crossingCell, cellValue);
					val = cellValue[0];
					if (val.compareTo(GlobalVars.GRID_PARAMS.getPriorityValueMap().get("pedestrian")) == 0) {
						routeIndices.add(j);
						descriptionMap.put(j, GlobalVars.TRANSPORT_PARAMS.routeDefaultDescription);
						break;
					}
				}
				
				// Search backwards
				for (int j = i; j>1; j--) {
					GridCoordinates2D crossingCell = gP.get(j);
					cellValue = grid.evaluate(crossingCell, cellValue);
					val = cellValue[0];
					if (val.compareTo(GlobalVars.GRID_PARAMS.getPriorityValueMap().get("pedestrian")) == 0) {
						routeIndices.add(j);
						descriptionMap.put(j, GlobalVars.TRANSPORT_PARAMS.routeDefaultDescription);
						break;
					}
				}
			}
		}
		
		
		// Now prune path coordinates that are redundant.
		// These are defined as those which lay between coordinates which are not separated by a ped obstruction
		// or change in road priority
		int startCellIndex = 0;
		GridCoordinates2D startCell = gP.get(startCellIndex);
		Coordinate startCoord = GISFunctions.gridCellToCoordinate(grid, startCell);
		
		// Given a fixed starting path coordinate, loop through path coordinates, create line between pairs
		// if line intersects with an obstacle, set the route coord to be the previous path coord for which there 
		// was no intersection
		for (int i = startCellIndex + 1; i < gP.size(); i++) {
			int gridPathIndexToIncludeInRoute;
			GridCoordinates2D gridCell = gP.get(i);
			Coordinate gridCoord = GISFunctions.gridCellToCoordinate(grid, gridCell);
			
			if (checkForObstaclesBetweenRouteCoordinates(startCoord, gridCoord)) {
				// Handle the case where lines between neighbouring cells intersect with obstrucctions in order to avoid infinite loop
				if (i-1 == startCellIndex) {
					gridPathIndexToIncludeInRoute = i;
				}
				else {
					gridPathIndexToIncludeInRoute = i-1;
				}
				
				// Only add this grid cell to the route if not already included
				if (!routeIndices.contains(gridPathIndexToIncludeInRoute)) {
					routeIndices.add(gridPathIndexToIncludeInRoute);
					descriptionMap.put(gridPathIndexToIncludeInRoute, GlobalVars.TRANSPORT_PARAMS.routeDefaultDescription);
				}

				startCellIndex = gridPathIndexToIncludeInRoute;
				
				// Update start cell
				startCell = gP.get(startCellIndex);
				startCoord = GISFunctions.gridCellToCoordinate(grid, startCell);
				
				// Loop continues two steps ahead from starting cell but this doesn't change route outcome
			}
		}
		
		// Order final set of route indicies
		List<Integer> routeIndicesSorted = routeIndices.stream().sorted().collect(Collectors.toList());
		for (int i:routeIndicesSorted) {
			GridCoordinates2D routeCell = gP.get(i);
			prunedGridPath.add(routeCell);
			Coordinate routeCoord = GISFunctions.gridCellToCoordinate(grid, routeCell);
			addToRoute(routeCoord, RoadLink.nullRoad, 1, descriptionMap.get(i));
		}
	}
	
	private boolean checkForObstaclesBetweenRouteCoordinates(Coordinate startCoord, Coordinate endCoord) {
		boolean isObstructingObjects = false;
		
		Coordinate[] lineCoords = {startCoord, endCoord};
		LineString pathLine = GISFunctions.lineStringGeometryFromCoordinates(lineCoords);

		
		// Check if line passes through a ped obstruction
		Geography<PedObstruction> pedObstructGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.PED_OBSTRUCTION_GEOGRAPHY);
		isObstructingObjects = GISFunctions.doesIntersectGeographyObjects(pathLine, pedObstructGeography);
		
		// check if line passes through two different priority road objects
		Geography<Road> roadGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.ROAD_GEOGRAPHY);
		List<Road> intersectingRoads = SpatialIndexManager.findIntersectingObjects(roadGeography, pathLine);
		if (intersectingRoads.size()>0) {
			String priority = intersectingRoads.get(0).getPriority();
			intersectingRoads.remove(0);
			for (Road intersectingR: intersectingRoads) {
				if (!intersectingR.getPriority().contentEquals(priority)) {
					isObstructingObjects = true;
				}
			}
		}
		
		return isObstructingObjects;
	}
	/**
	 * Find a path through a grid coverage layer by using the flood fill algorithm to calculate cell 'costs'
	 * and acting greedily to identify a path.
	 * 
	 * @param gridCoverageName
	 * 			The name of the coverage layer to use
	 * @return
	 * 			List<GridCoordinates2D> The grid coordinates path
	 */
	private List<GridCoordinates2D> getFloodFillGridPath(GridCoverage2D grid, Coordinate o, Coordinate d){
		
		List<GridCoordinates2D> gP = new ArrayList<GridCoordinates2D>();

		DirectPosition2D dpStart = new DirectPosition2D(o.x, o.y);
		DirectPosition2D dpEnd = new DirectPosition2D(d.x, d.y);
		GridCoordinates2D start = null;
		GridCoordinates2D end = null;
		try {
			start = grid.getGridGeometry().worldToGrid(dpStart);
			end = grid.getGridGeometry().worldToGrid(dpEnd);
		} catch (InvalidGridGeometryException | TransformException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Map<String, Integer> gridBounds = gridBounds(grid, start, end, this.partialFF);
		int mini = gridBounds.get("mini");
		int minj = gridBounds.get("minj");
		int maxi = gridBounds.get("maxi");
		int maxj = gridBounds.get("maxj");
		
		double[][] cellValues = gridCoverageFloodFill(grid, end, mini, minj, maxi, maxj);
		boolean atEnd = false;
		
		GridCoordinates2D next = start;
		gP.add(next);
		while(!atEnd) {
			next = greedyMooreNeighbour(next, cellValues, gP, mini, minj, maxi, maxj);
			gP.add(next);
			if (next.equals(end)) {
				atEnd = true;
			}
		}
		
		return gP;
	}
	
	/**
	 * Runs the flood fill algorithm on the grid coverage with the name given as an input parameter. This algorithm assigns 
	 * each cell a value based on its distance from the mobile agent's end destination. THe process stops once the cell 
	 * containing the mobile agents current position is reached.
	 * 
	 * The grid coverage cell values are used in the calculation of the flood fill values, they represent the 'cost' of travelling through that cell.
	 * @param gridCoverageName
	 * 			The name of the grid coverage to use when calculating cell flood fill values
	 * @return
	 * 			2D double array of cell values
	 */
	private double[][] gridCoverageFloodFill(GridCoverage2D grid, GridCoordinates2D end, int mini, int minj, int maxi, int maxj) {

		int width = grid.getRenderedImage().getTileWidth();
		int height = grid.getRenderedImage().getTileHeight();
		
			floodFillValues = new double[height][width]; // Initialised with zeros
		int [][] n = new int[height][width]; // Use to log number of times a cell is visited. All cells should get visited once.
		List<GridCoordinates2D> q = new ArrayList<GridCoordinates2D>();
		
		// Staring at the destination, flood fill cell values using distance measure between cells
		GridCoordinates2D thisCell;
		double thisCellValue;
		double nextCellValue;
		int[] cellValue = new int[1];
		
		int i = end.x;
		int j = end.y;
		n[j][i] = 1; // Make sure the end cell value doesn't get updated
		q.add(end);
		while(q.size() > 0) {
			// Get random cell from queue - random is important
			int nextCellIndex = RandomHelper.nextIntFromTo(0, q.size()-1);
			thisCell = q.get(nextCellIndex);
			q.remove(nextCellIndex);
			
			thisCellValue = floodFillValues[thisCell.y][thisCell.x];
			for (GridCoordinates2D nextCell: xNeighbours(thisCell, "moore", mini, minj, maxi, maxj)) {
				
				cellValue = grid.evaluate(nextCell, cellValue);
				i = nextCell.x;
				j = nextCell.y;
				
				// If cell with default value, assign value the max int value and exclude from further computation
				if (cellValue[0] == GlobalVars.GRID_PARAMS.defaultGridValue) {
					floodFillValues[j][i] = Integer.MAX_VALUE;
					n[j][i] += 1;
					continue;
				}
				// Ensure the next cell doesn't already have a value
				if (n[j][i] == 0) {
					// Get the cost of moving through this cell for the mobile agent by mapping from cell value using agents priority map
					double summand = this.gridSummandPriorityMap.get(cellValue[0]);
					nextCellValue = thisCellValue + summand;
					floodFillValues[j][i] = nextCellValue;
					n[j][i] += 1;
					q.add(nextCell);
				}
			}
		}
		return floodFillValues;
	}
	
	private GridCoordinates2D greedyManhattanNeighbour(GridCoordinates2D cell, double[][] cellValues, List<GridCoordinates2D> path, int mini, int minj, int maxi, int maxj) {
		GridCoordinates2D greedyNeighbour = greedyXNeighbour(cell, cellValues, path, "manhattan", mini, minj, maxi, maxj);
		return greedyNeighbour;
	}
	
	private GridCoordinates2D greedyMooreNeighbour(GridCoordinates2D cell, double[][] cellValues, List<GridCoordinates2D> path, int mini, int minj, int maxi, int maxj) {
		GridCoordinates2D greedyNeighbour = greedyXNeighbour(cell, cellValues, path, "moore", mini, minj, maxi, maxj);
		return greedyNeighbour;
	}
	
	private GridCoordinates2D greedyXNeighbour(GridCoordinates2D cell, double[][] cellValues, List<GridCoordinates2D> path, String neighbourType, int mini, int minj, int maxi, int maxj) {
		
		List<GridCoordinates2D> neighbours = xNeighbours(cell, neighbourType, mini, minj, maxi, maxj);
		
		// Initialise greedy options
		List<Double> minVal = new ArrayList<Double>();
		List<GridCoordinates2D> greedyNeighbours = new ArrayList<GridCoordinates2D>();
		
		minVal.add((double) Integer.MAX_VALUE);

		
		for(GridCoordinates2D neighbour:neighbours) {
			// Don't consider cells already in the path
			if (path.contains(neighbour)) {
				continue;
			}
			double val = cellValues[neighbour.y][neighbour.x];
			
			// If cell value equal to current minimum include in greedy option
			if (Math.abs(val - minVal.get(0)) < 0.0000000001) {
				minVal.add(val);
				greedyNeighbours.add(neighbour);
			}
			
			// Else  clear the current min values and replace with new min
			else if (val < minVal.get(0)) {
				// Replace old values with new ones
				minVal.clear();
				greedyNeighbours.clear();
				
				minVal.add(val);
				greedyNeighbours.add(neighbour);
			}
			else {
				continue;
			}
		}
		
	    GridCoordinates2D greedyNeighbour = null;
	    try {
	    	int neighbourIndex = RandomHelper.nextIntFromTo(0, greedyNeighbours.size()-1);
		    greedyNeighbour =  greedyNeighbours.get(neighbourIndex);
	    } catch (IllegalArgumentException e) {
	    	String msg = "GridRoute.greedyXNeighbour(): Problem getting greedy neighbour \n\r" + 
	    			"neighbour type: " + neighbourType + "\n\r" +
					"n greedy neighbours: " + String.valueOf(greedyNeighbours.size()) + "\n\r" +
					"n neighbourghs: " + String.valueOf(neighbours.size()) + "\n\r" +
					"origin coord: " + this.origin.toString() + "\n\r" +
					"destination coord: " + this.destination.toString() + "\n\r";
			//LOGGER.log(Level.SEVERE,msg);
			System.out.print(msg);
			throw e;
	    }
	    return greedyNeighbour;
	}
	
	
	/**
	 * Find a path through a grid coverage layer by using the A star algorithm.
	 * 
	 * @param grid
	 * 		The grid coverage object
	 * @param o
	 * 		The coordinate of the origin
	 * @param d
	 * 		The coordinate of the destination
	 * @return
	 * 			List<GridCoordinates2D> The grid coordinates path
	 */
	private List<GridCoordinates2D> getAStarGridPath(GridCoverage2D grid, Coordinate o, Coordinate d){
		
		List<GridCoordinates2D> gP = new ArrayList<GridCoordinates2D>();

		DirectPosition2D dpStart = new DirectPosition2D(o.x, o.y);
		DirectPosition2D dpEnd = new DirectPosition2D(d.x, d.y);
		GridCoordinates2D start = null;
		GridCoordinates2D end = null;
		try {
			start = grid.getGridGeometry().worldToGrid(dpStart);
			end = grid.getGridGeometry().worldToGrid(dpEnd);
		} catch (InvalidGridGeometryException | TransformException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Map<String, Integer> gridBounds = gridBounds(grid, start, end, this.partialFF);
		int mini = gridBounds.get("mini");
		int minj = gridBounds.get("minj");
		int maxi = gridBounds.get("maxi");
		int maxj = gridBounds.get("maxj");
		
		/*
		 * Do A Star Algo
		 * 
		 * Taken, in part, from Justin Wetherell's Java implementation of the A Star algorithm:
		 * https://github.com/phishman3579/java-algorithms-implementation/blob/master/src/com/jwetherell/algorithms/graph/AStar.java
		 * 
		 * 
		 */
		
		// Initialise records of grid cell distance measures used to choose path
		Set<GridCoordinates2D> closedSet = new HashSet<GridCoordinates2D>();
		List<GridCoordinates2D> openSet = new ArrayList<GridCoordinates2D>();
		Map<GridCoordinates2D, GridCoordinates2D> cameFrom = new HashMap<GridCoordinates2D, GridCoordinates2D>();
		Map<GridCoordinates2D, Double> gScore = new HashMap<GridCoordinates2D, Double>();
		Map<GridCoordinates2D, Double> fScore = new HashMap<GridCoordinates2D, Double>();
		
		GridCoordinates2D thisCell;
		int[] cellValue = new int[1];
		
	    Comparator<GridCoordinates2D> fComparator = new Comparator<GridCoordinates2D>() {
	        /**
	         * {@inheritDoc}
	         */
	        @Override
	        public int compare(GridCoordinates2D c1, GridCoordinates2D c2) {
	            if (fScore.get(c1) < fScore.get(c2))
	                return -1;
	            if (fScore.get(c2) < fScore.get(c1))
	                return 1;
	            return 0;
	        }
	    };
		
		// Initialise f score values
		for (int j = minj; j<maxj;j++) {
			for (int i = mini; i<maxi;i++) {
				fScore.put(new GridCoordinates2D(i, j), (double) Integer.MAX_VALUE);
			}
		}
		
		// Add stating coordinate to the openSet
		openSet.add(start);
        gScore.put(start, 0.0);
		
		// Run algorithm until destination is reached
		while(!openSet.isEmpty()) {
			thisCell = openSet.get(0);
			
			// If this is the destination, job's done
			if (thisCell.equals(end)) {
				// Unpack and return the path
				gP = reconstructPath(cameFrom, thisCell);
				return gP;
			}
			
			openSet.remove(0);
			closedSet.add(thisCell);
			
			// Get neighbours and compute fScores for each of these
			for (GridCoordinates2D nextCell: xNeighbours(thisCell, "moore", mini, minj, maxi, maxj)) {
				if(closedSet.contains(nextCell)) {
					continue;
				}
				
				cellValue = grid.evaluate(nextCell, cellValue);
				double summand;
				if (cellValue[0] == GlobalVars.GRID_PARAMS.defaultGridValue) {
					summand = Integer.MAX_VALUE;
				}
				else {
					summand = this.gridSummandPriorityMap.get(cellValue[0]);
				}
				double tentativeGScore = gScore.get(thisCell) + summand;
				if(!openSet.contains(nextCell)) {
					openSet.add(nextCell); // New candidate cell identified
				}
				else if (tentativeGScore >= gScore.get(nextCell)) {
					continue; // Don't bother if neighbour not going to have good f score
				}
				
				// thisCell is best path so far, record it
				cameFrom.put(nextCell, thisCell);
				gScore.put(nextCell, tentativeGScore);
				double tentativeFScore = gScore.get(nextCell) + heuristicCostEstimate(nextCell, end);
				fScore.put(nextCell, tentativeFScore);
				Collections.sort(openSet, fComparator);
			}	
		}
		return null;
	}
	
	/**
	 * Find a path through a grid coverage layer by using the Dijkstra's algorithm.
	 * 
	 * To do this, need to first create a network object from the grid , which enables the use of
	 * existing implementation of Dijkstra's.
	 * 
	 * @param grid
	 * 		The grid coverage object
	 * @param o
	 * 		The coordinate of the origin
	 * @param d
	 * 		The coordinate of the destination
	 * @return
	 * 			List<GridCoordinates2D> The grid coordinates path
	 */
	private List<GridCoordinates2D> getDijkstraGridPath(GridCoverage2D grid, Coordinate o, Coordinate d){
		
		// Initialise the grid path to be produced
		List<GridCoordinates2D> gP = new ArrayList<GridCoordinates2D>();
		
		DirectPosition2D dpStart = new DirectPosition2D(o.x, o.y);
		DirectPosition2D dpEnd = new DirectPosition2D(d.x, d.y);
		GridCoordinates2D start = null;
		GridCoordinates2D end = null;
		try {
			start = grid.getGridGeometry().worldToGrid(dpStart);
			end = grid.getGridGeometry().worldToGrid(dpEnd);
		} catch (InvalidGridGeometryException | TransformException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Network<GridCoordinates2D> gridNetwork = new UndirectedJungNetwork<GridCoordinates2D>(GlobalVars.CONTEXT_NAMES.GRID_NETWORK, new NetworkEdgeCreator<GridCoordinates2D>());
		
		// Initialise records of grid cell distance measures used to choose path
		Set<GridCoordinates2D> closedSet = new HashSet<GridCoordinates2D>();
		List<GridCoordinates2D> openSet = new ArrayList<GridCoordinates2D>();
		
		GridCoordinates2D thisCell;
		int[] cellValue = new int[1];
		
		Map<String, Integer> gridBounds = gridBounds(grid, start, end, this.partialFF);
		int mini = gridBounds.get("mini");
		int minj = gridBounds.get("minj");
		int maxi = gridBounds.get("maxi");
		int maxj = gridBounds.get("maxj");
		
		// Create the network
		openSet.add(start);
		while(!openSet.isEmpty()) {
			thisCell = openSet.get(0);			
			openSet.remove(0);
			closedSet.add(thisCell);
			
			// Get neighbours and compute fScores for each of these
			for (GridCoordinates2D nextCell: xNeighbours(thisCell, "moore", mini, minj, maxi, maxj)) {
				if(closedSet.contains(nextCell)) {
					continue;
				}
				
				cellValue = grid.evaluate(nextCell, cellValue);
				double perceivedValue;
				if (cellValue[0] == GlobalVars.GRID_PARAMS.defaultGridValue) {
					perceivedValue = Integer.MAX_VALUE;
				}
				else {
					perceivedValue = this.gridSummandPriorityMap.get(cellValue[0]);
				}
		
				// Create an edge between the two junctions, assigning a weight equal to it's length
				RepastEdge<GridCoordinates2D> edge = new RepastEdge<GridCoordinates2D>(thisCell, nextCell, false, perceivedValue);
		
				// Add the edge to the network
				if (!gridNetwork.containsEdge(edge)) {
					gridNetwork.addEdge(edge);
				} else {
					//LOGGER.severe("CityContext: buildRoadNetwork: for some reason this edge that has just been created "
						//	+ "already exists in the RoadNetwork.");
				}
				
				if(!openSet.contains(nextCell)) {
					openSet.add(nextCell); // New candidate cell identified
				}
			}
		}
		
		// Get path using Dijkstras
		ShortestPath<GridCoordinates2D> p = new ShortestPath<GridCoordinates2D>(gridNetwork);
		List<RepastEdge<GridCoordinates2D>> shortestPath = p.getPath(start, end);
		
		// Iterate over the path to get the grid cells that make up the path
		for (RepastEdge<GridCoordinates2D> edge : shortestPath) {
			gP.add(edge.getSource());
		}

		return gP;
	}
	
	private double heuristicCostEstimate(GridCoordinates2D cell, GridCoordinates2D destinationCell) {
		// Calculate as the crow flies distance
		double dist = Math.sqrt(Math.pow(cell.x - destinationCell.x, 2) +	Math.pow(cell.y - destinationCell.y,2));
		return dist;
	}
	
    private List<GridCoordinates2D> reconstructPath(Map<GridCoordinates2D,GridCoordinates2D> cameFrom, GridCoordinates2D thisCell) {
        List<GridCoordinates2D> path = new ArrayList<GridCoordinates2D>();

        while (thisCell != null) {
            GridCoordinates2D previousCell = thisCell;
            thisCell = cameFrom.get(thisCell);
            if (thisCell != null) {
            	path.add(thisCell);
            }
        }
        Collections.reverse(path);
        return path;
    }
	
	/**
	 * Given a grid coordinate return a list of the Manhattan neighbours of this coordinate (N, E, S, W)
	 * @param cell
	 * 			The grid coordinate to get the neighbours of
	 * @return
	 * 			List of GridCoordinates2D objects
	 */
	private List<GridCoordinates2D> manhattanNeighbourghs(GridCoordinates2D cell){
		
		List<GridCoordinates2D> mN = new ArrayList<GridCoordinates2D>();
		
		int[] range = {-1,1};
		int i = cell.x;
		int j = cell.y;
		for (int dx: range) {
			mN.add(new GridCoordinates2D(i + dx, j));
		}
		
		for (int dy: range) {
			mN.add(new GridCoordinates2D(i, j + dy));
		}
		
		return mN;
	}
	
	/**
	 * Given a grid coordinate return a list of the Manhattan neighbours of this coordinate (N, E, S, W)
	 * 
	 * Exclude coordinates that are lower than the minimum i and j parameters or greater than or equal to the
	 * maximum i and j parameters.
	 * 
	 * @param cell
	 * 			The grid coordinate to get the neighbours of
	 * @param mini
	 * 			Minimum i value
	 * @param minj
	 * 			Minimum j value
	 * @param maxi
	 * 			Maximum i value
	 * @param maxj
	 * 			Maximum j value
	 * @return
	 * 			List of GridCoordinates2D objects
	 */
	private List<GridCoordinates2D> manhattanNeighbourghs(GridCoordinates2D cell, int mini, int minj, int maxi, int maxj){
		
		List<GridCoordinates2D> mN = new ArrayList<GridCoordinates2D>();
		
		int[] range = {-1,1};
		int i = cell.x;
		int j = cell.y;
		for (int dx: range) {
			if ((i + dx >= mini) & (i + dx < maxi) & (j >= minj) & (j < maxj)) {
				mN.add(new GridCoordinates2D(i + dx, j));
			}
		}
		
		for (int dy: range) {
			if ((i >= mini) & (i < maxi) & (j + dy >= minj) & (j + dy < maxj)) {
				mN.add(new GridCoordinates2D(i, j + dy));
			}
		}
		
		return mN;
	}
	
	/**
	 * Given a grid coordinate return a list of the Moore neighbours of this coordinate.
	 * 
	 * Exclude coordinates that are lower than the minimum i and j parameters or greater than or equal to the
	 * maximum i and j parameters.
	 * 
	 * @param cell
	 * 			The grid coordinate to get the neighbours of
	 * @param mini
	 * 			Minimum i value
	 * @param minj
	 * 			Minimum j value
	 * @param maxi
	 * 			Maximum i value
	 * @param maxj
	 * 			Maximum j value
	 * @return
	 * 			List of GridCoordinates2D objects
	 */
	private List<GridCoordinates2D> mooreNeighbourghs(GridCoordinates2D cell, int mini, int minj, int maxi, int maxj){
		
		List<GridCoordinates2D> mN = new ArrayList<GridCoordinates2D>();
		
		int[] range = {-1,1,0};
		int i = cell.x;
		int j = cell.y;
		for (int dx: range) {
			for (int dy: range) {
				if ((dx==0) & (dy==0)) {
					continue;
				}
				if ((i + dx >= mini) & (i + dx < maxi) & (j + dy >= minj) & (j + dy < maxj)) {
					mN.add(new GridCoordinates2D(i + dx, j + dy));
				}
			}
		}
		
		return mN;
	}
	
	private List<GridCoordinates2D> xNeighbours(GridCoordinates2D cell, String neighbourType, int mini, int minj, int maxi, int maxj){
		List<GridCoordinates2D> neighbours;
		if (neighbourType.contentEquals("manhattan")) {
			neighbours = manhattanNeighbourghs(cell, mini, minj, maxi, maxj);
		}
		else {
			neighbours = mooreNeighbourghs(cell, mini, minj, maxi, maxj);
		}
		return neighbours;
	}
	
	private Map<String, Integer> gridBounds(GridCoverage2D grid, GridCoordinates2D start, GridCoordinates2D end, Boolean isPartial){
		Map<String, Integer> bounds = new HashMap<String, Integer>();
		
		// Set the bounds of the flood fill search based on whether this is a partial route search or not
		int width = grid.getRenderedImage().getTileWidth();
		int height = grid.getRenderedImage().getTileHeight();
		int mini = 0;
		int minj = 0;
		int maxi = width;
		int maxj = height;
		if (this.partialFF == true) {
			// Bounds set to bounding box of start-destination +- partialBoundingBoxIncrease % in x and y direction
			int dx = Math.abs(start.x-end.x);
			int dy = Math.abs(start.y-end.y);
			
			int xIncrease =  Math.max(5, (int) Math.floor(dx*GlobalVars.TRANSPORT_PARAMS.partialBoundingBoxIncrease));
			int yIncrease =  Math.max(5, (int) Math.floor(dy*GlobalVars.TRANSPORT_PARAMS.partialBoundingBoxIncrease));
			
			mini = Math.max(Math.min(start.x, end.x) - xIncrease,0);
			minj = Math.max(Math.min(start.y, end.y) - yIncrease,0);
			maxi = Math.min(Math.max(start.x, end.x) + xIncrease, width);
			maxj = Math.min(Math.max(start.y, end.y) + yIncrease, height);
		}
		bounds.put("mini", mini);
		bounds.put("minj", minj);
		bounds.put("maxi", maxi);
		bounds.put("maxj", maxj);
		return bounds;
	}
	
	private void checkGridRoute() {
		// Each primary route coordinate should be in the route map
		assert this.routeCoordMap.keySet().containsAll(this.primaryRouteX) : 
				"Mismatch between list of primary route coordinates and the coordiante:grid cell map";
		
		// Each grid cell the primary route coordinate looks up to should be used as an index in the grouped grid path
		assert this.routeCoordMap.values().containsAll(this.groupedGridPath.keySet()) : 
				"Mismatch between the grid cells used to index sections of the path and those in the coordiante:grid cell map";
	}	
	
	public void setNextRouteSection() {
		Coordinate nextRoadLinkCoord = this.primaryRouteX.get(0);
		List<GridCoordinates2D> nextPathSection = this.groupedGridPath.get(this.routeCoordMap.get(nextRoadLinkCoord));
		addCoordinatesToRouteFromGridPath(nextPathSection);
		
		// Finish by adding next road link route coord
		addToRoute(this.primaryRouteX.get(1), RoadLink.nullRoad, 1, GlobalVars.TRANSPORT_PARAMS.routeRoadLinkChangeDescription);
		
		// Remove the current road link coord since this section of the path has been added to the route
		this.primaryRouteX.remove(0);
	}
	
	private void initialiseRouteVectors() {
		
		this.routeX = new Vector<Coordinate>();
		this.roadsX = new Vector<RoadLink>();
		this.routeDescriptionX = new Vector<String>();
		this.routeSpeedsX = new Vector<Double>();
		this.gridPathCrossings = new Vector<GridCoordinates2D>();
		this.prunedGridPath = new Vector<GridCoordinates2D>();
		this.primaryRouteX = new Vector<Coordinate>();
		this.routeCoordMap = new HashMap<Coordinate, GridCoordinates2D>();
		this.groupedGridPath = new HashMap<GridCoordinates2D, List<GridCoordinates2D>>();
		
	}
	
	/**
	 * Convenience function that can be used to add details to the route. This should be used rather than updating
	 * individual lists because it makes sure that all lists stay in sync
	 * 
	 * @param coord
	 *            The coordinate to add to the route
	 * @param nullRoad
	 *            The road that the coordinate is part of
	 * @param speed
	 *            The speed that the road can be travelled along
	 * @param description
	 *            A description of why the coordinate has been added
	 */
	protected void addToRoute(Coordinate coord, RoadLink nullRoad, double speed, String description) {
		this.routeX.add(coord);
		this.roadsX.add(nullRoad);
		this.routeSpeedsX.add(speed);
		this.routeDescriptionX.add(description);
	}
	
	public void removeNextFromRoute() {
		removeFromRoute(0);
	}
	
	protected void removeFromRoute(int index) {
		this.routeX.remove(index);
		this.roadsX.remove(index);
		this.routeSpeedsX.remove(index);
		this.routeDescriptionX.remove(index);
	}
	
	public Coordinate getNextRouteCrossingCoord() {
		Coordinate crossingC = null;
		for (int i = 0; i< this.routeDescriptionX.size(); i++) {
			if (this.routeDescriptionX.get(i).contentEquals(GlobalVars.TRANSPORT_PARAMS.routeCrossingDescription)) {
				crossingC = this.routeX.get(i);
				break;
			}
		}
		// Could be null if there is not a crossing in the upcoming section of the route.
		return crossingC;
	}
	
	public List<GridCoordinates2D> getNextGroupedGridPathSection(){
		
		Coordinate thisRoadLinkCoord = this.primaryRouteX.get(0);
		GridCoordinates2D thisRoadLinkCell = this.routeCoordMap.get(thisRoadLinkCoord);
		List<GridCoordinates2D> pathSection = this.groupedGridPath.get(thisRoadLinkCell);
		return pathSection;
	}
	
	public double[][] getFloodFillGridValues() {
		return this.floodFillValues;
	}
	
	public List<GridCoordinates2D> getGridPath(){
		return this.gridPath;
	}
	
	public List<GridCoordinates2D> getPrunedGridPath(){
		return this.prunedGridPath;
	}
	
	public List<GridCoordinates2D> getGridPathCrossings(){
		return this.gridPathCrossings;
	}
	
	public List<Coordinate> getPrimaryRouteX(){
		return this.primaryRouteX;
	}

	public Map<GridCoordinates2D, List<GridCoordinates2D>> getGroupedGridPath() {
		return groupedGridPath;
	}

	public Map<Coordinate, GridCoordinates2D> getRouteCoordMap() {
		return routeCoordMap;
	}

	public List<Coordinate> getRouteX() {
		return routeX;
	}
}
