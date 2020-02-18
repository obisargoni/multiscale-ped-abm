package repastInterSim.agent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.coverage.grid.InvalidGridGeometryException;
import org.geotools.geometry.DirectPosition2D;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Polygon;

import repast.simphony.space.gis.Geography;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class GridRoute extends Route {
	
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
	 * Create a new route object
	 * 
	 * @param geography
	 * 		The geography projection that the mobile agent this route belongs to is in
	 * @param mA
	 * 		The mobile agent this route belongs to
	 * @param gSPM
	 * 		The map from integers used to indicate the road user priority of grid cells to the agents perceived cost of moving through those grid cells. 
	 * Used for routing on a grid.
	 * @param destination
	 * 		The destination coordinate of the route
	 */
	public GridRoute(Geography<Object> geography, MobileAgent mA,  HashMap<Integer, Double> gSPM, Coordinate destination) {
		super(geography, mA, destination);
		// TODO Auto-generated constructor stub
		this.gridSummandPriorityMap = gSPM;
	}
	
	/**
	 * Create a new route object
	 * 
	 * @param geography
	 * 		The geography projection that the mobile agent this route belongs to is in
	 * @param mA
	 * 		The mobile agent this route belongs to
	 * @param gSPM
	 * 		The map from integers used to indicate the road user priority of grid cells to the agents perceived cost of moving through those grid cells. 
	 * Used for routing on a grid.
	 * @param destination
	 * 		The destination coordinate of the route
	 * @param
	 * 		A boolean value indicating whether to consider only a partial area of the grid when producing the route.
	 */
	public GridRoute(Geography<Object> geography, MobileAgent mA, HashMap<Integer, Double> gSPM, Coordinate destination, boolean partial) {
		super(geography, mA, destination);
		// TODO Auto-generated constructor stub
		this.gridSummandPriorityMap = gSPM;
		this.partialFF = partial;
	}
	
	/**
	 * Create a new route object
	 * 
	 * @param geography
	 * 		The geography projection that the mobile agent this route belongs to is in
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
	public GridRoute(Geography<Object> geography, MobileAgent mA, HashMap<Integer, Double> gSPM, Coordinate origin, Coordinate destination, boolean partial) {
		super(geography, mA, origin, destination);
		// TODO Auto-generated constructor stub
		this.gridSummandPriorityMap = gSPM;
		this.partialFF = partial;
	}
	
	/**
	 * This method produces a path of grid cells from the agents current position to the destination set when creating a GridRoute instance. 
	 * The grid path is grouped into sections according to the road link the grid cells belong to. 
	 */
	public void setGroupedGridPath() {
		this.routeX = new Vector<Coordinate>();
		this.roadsX = new Vector<RoadLink>();
		this.routeDescriptionX = new Vector<String>();
		this.routeSpeedsX = new Vector<Double>();
		this.gridPathCrossings = new Vector<GridCoordinates2D>();
		this.prunedGridPath = new Vector<GridCoordinates2D>();
		this.primaryRouteX = new Vector<Coordinate>();
		this.routeCoordMap = new HashMap<Coordinate, GridCoordinates2D>();
		this.groupedGridPath = new HashMap<GridCoordinates2D, List<GridCoordinates2D>>();
		
		GridCoverage2D grid = geography.getCoverage(GlobalVars.CONTEXT_NAMES.BASE_COVERAGE);
		
		if (this.origin == null) {
			this.origin = this.mA.getLoc();
		}
		
		// sequence of grid cell coordinates leading from agent' current position to end destination
		this.gridPath = getGridCoveragePath(grid, this.origin, this.destination);
				
		// First grid cell coordinate is agent's first coordinate along that road link, so use to index first group of coordinates
		GridCoordinates2D primaryRouteCell = null;
		Coordinate primaryRouteCoord = null;
		String prevCellRoadLinkFID = null;
		String cellRoadLinkFID = null;
		
		// Identify the primary route coordinates as those coordinate along the path where the road link changes
		// Use each primary grid cell to index a new group of grid cell coordinates (the path along a particular road link)
		for (int i = 0; i < gridPath.size(); i++) {			
			GridCoordinates2D gridCell = gridPath.get(i);
			GridEnvelope2D gridEnv = new GridEnvelope2D(gridCell.x, gridCell.y, 1, 1);
			Polygon cellPoly = GISFunctions.getWorldPolygonFromGridEnvelope(grid, gridEnv);
			
			List<String> cellRoadLinkFIDs = null;
			try {
				cellRoadLinkFIDs = GISFunctions.getGridPolygonRoadLinkIDs(cellPoly);
			} catch (RoutingException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			// Only update the current grid cell road link ID if there is only one road link ID associated with this cell
			// Avoids ambiguity about which road link the primary route coordinates link to
			if (cellRoadLinkFIDs.size()==1) {
				cellRoadLinkFID = cellRoadLinkFIDs.get(0);
			}
			
			// Check for change in road link id
			// IDs will be the same unless cellRoadLinkID has been updated. Ensures that only coordinates that unambiguously belong to a road linj
			// are added as primary route coordinates
			if (!cellRoadLinkFID.equals(prevCellRoadLinkFID)) {
				primaryRouteCell = gridCell;
				primaryRouteCoord = gridCellToCoordinate(grid, primaryRouteCell);
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
	
	private void addCoordinatesToRouteFromGridPath(List<GridCoordinates2D> gridPath) {
		
		Set<Integer> routeIndices = new HashSet<Integer>();
		Map<Integer, String> descriptionMap = new HashMap<Integer, String>();
				
		GridCoverage2D grid = geography.getCoverage(GlobalVars.CONTEXT_NAMES.BASE_COVERAGE);
		
		double[] prevCellValue = new double[1];
		double[] cellValue = new double[1];
		GridCoordinates2D prevCell = gridPath.get(0);
		
		// Get indices of grid cells that are at location where road priority changes (crossing points)
		for (int i = 1; i < gridPath.size(); i++) {
			GridCoordinates2D gridCell = gridPath.get(i);
			
			// Get grid cell value of this and previous coord. If values differ this means they are located in
			// road space with different priority and therefore the previous grid cell should be included in the route
			prevCellValue = grid.evaluate(prevCell, prevCellValue);
			cellValue = grid.evaluate(gridCell, cellValue);
			Double prevVal = prevCellValue[0];
			Double val = cellValue[0];
			
			// If grid cell value increases, priority has decreased for this agent. Indicates crossing point where yielding is possible
			if (val.compareTo(prevVal) > 0) {
				routeIndices.add(i);
				descriptionMap.put(i, GlobalVars.TRANSPORT_PARAMS.routeCrossingDescription);
			}
			// If grid cell value decreases, this indicates this agents' priority is greater. Also crossing point but not one where yielding required
			else if (val.compareTo(prevVal) < 0) {
				routeIndices.add(i);
				descriptionMap.put(i, GlobalVars.TRANSPORT_PARAMS.routeDefaultDescription);

			}
			prevCell = gridCell;
		}
		
		
		// Now prune path coordinates that are redundant.
		// These are defined as those which lay between coordinates which are not separated by a ped obstruction
		// or change in road priority
		int startCellIndex = 0;
		GridCoordinates2D startCell = gridPath.get(startCellIndex);
		Coordinate startCoord = gridCellToCoordinate(grid, startCell);
		
		// Given a fixed starting path coordinate, loop through path coordinates, create line between pairs
		// if line intersects with an obstacle, set the route coord to be the previous path coord for which there 
		// was no intersection
		for (int i = startCellIndex + 1; i < gridPath.size(); i++) {
			int gridPathIndexToIncludeInRoute;
			GridCoordinates2D gridCell = gridPath.get(i);
			Coordinate gridCoord = gridCellToCoordinate(grid, gridCell);
			
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
				startCell = gridPath.get(startCellIndex);
				startCoord = gridCellToCoordinate(grid, startCell);
				
				// Loop continues two steps ahead from starting cell but this doesn't change route outcome
			}
		}
		
		// Order final set of route indicies
		List<Integer> routeIndicesSorted = routeIndices.stream().sorted().collect(Collectors.toList());
		for (int i:routeIndicesSorted) {
			GridCoordinates2D routeCell = gridPath.get(i);
			prunedGridPath.add(routeCell);
			Coordinate routeCoord = gridCellToCoordinate(grid, routeCell);
			addToRoute(routeCoord, RoadLink.nullRoad, 1, descriptionMap.get(i));
		}
	}
	
	private boolean checkForObstaclesBetweenRouteCoordinates(Coordinate startCoord, Coordinate endCoord) {
		boolean isObstructingObjects = false;
		
		Coordinate[] lineCoords = {startCoord, endCoord};
		LineString pathLine = new GeometryFactory().createLineString(lineCoords);
		
		// Check if line passes through a ped obstruction
		// If it does add the previous index to the pruned path list
		List<PedObstruction> intersectingObs = SpatialIndexManager.findIntersectingObjects(SpaceBuilder.pedObstructGeography, pathLine);
		if (intersectingObs.size() > 0){
			isObstructingObjects = true;
		}
		
		List<Road> intersectingRoads = SpatialIndexManager.findIntersectingObjects(SpaceBuilder.roadGeography, pathLine);
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
	private List<GridCoordinates2D> getGridCoveragePath(GridCoverage2D grid, Coordinate o, Coordinate d){
		
		List<GridCoordinates2D> gridPath = new ArrayList<GridCoordinates2D>();

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
		
		// Set the bounds of the flood fill search based on whether this is a partial route search or not
		int width = grid.getRenderedImage().getTileWidth();
		int height = grid.getRenderedImage().getTileHeight();
		int mini = 0;
		int minj = 0;
		int maxi = width;
		int maxj = height;
		if (this.partialFF == true) {
			// Bounds set to bounding box of start-destination +- 30% in x and y direction
			int dx = Math.abs(start.x-end.x);
			int dy = Math.abs(start.y-end.y);
			
			mini = Math.max(Math.min(start.x, end.x) - (int) Math.floor(dx*GlobalVars.TRANSPORT_PARAMS.partialBoundingBoxIncrease),0);
			minj = Math.max(Math.min(start.y, end.y) - (int) Math.floor(dy*GlobalVars.TRANSPORT_PARAMS.partialBoundingBoxIncrease),0);
			maxi = Math.min(Math.max(start.x, end.x) + (int) Math.floor(dx*GlobalVars.TRANSPORT_PARAMS.partialBoundingBoxIncrease), width);
			maxj = Math.min(Math.max(start.y, end.y) + (int) Math.floor(dy*GlobalVars.TRANSPORT_PARAMS.partialBoundingBoxIncrease), height);
		}
		
		double[][] cellValues = gridCoverageFloodFill(grid, end, mini, minj, maxi, maxj);
		boolean atEnd = false;
		
		GridCoordinates2D next = start;
		while(!atEnd) {
			if (next.equals(end)) {
				atEnd = true;
			}
			gridPath.add(next);
			next = greedyManhattanNeighbour(next, cellValues, gridPath, mini, minj, maxi, maxj);
		}
		
		return gridPath;
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
			thisCell = q.get(0);
			q.remove(0);
			
			thisCellValue = floodFillValues[thisCell.y][thisCell.x];
			for (GridCoordinates2D nextCell: manhattanNeighbourghs(thisCell, mini, minj, maxi, maxj)) {
				
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
	
	private GridCoordinates2D greedyManhattanNeighbour(GridCoordinates2D cell, double[][] cellValues, List<GridCoordinates2D> path, int mini, int minj, int maxi, int maxj) {
		List<GridCoordinates2D> manhattanNeighbours = manhattanNeighbourghs(cell, mini, minj, maxi, maxj);
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
				if ((i + dx >= mini) & (i + dx < maxi) & (j >= minj) & (j < maxj)) {
					mN.add(new GridCoordinates2D(i + dx, j + dy));
				}
			}
		}
		
		return mN;
	}
		
		// Initialise greedy options
		List<Double> minVal = new ArrayList<Double>();
		List<GridCoordinates2D> greedyNeighbours = new ArrayList<GridCoordinates2D>();
		
		minVal.add((double) Integer.MAX_VALUE);

		
		for(GridCoordinates2D neighbour:manhattanNeighbours) {
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
		
	    Random rand = new Random();
	    GridCoordinates2D greedyNeighbour = null;
	    try {
		    greedyNeighbour =  greedyNeighbours.get(rand.nextInt(greedyNeighbours.size()));
	    } catch (IllegalArgumentException e) {
	    	String msg = "GridRoute.greedyManhattanNeighbour(): Problem getting greedy neighbour \n\r" + 
					"n greedy neighbours: " + String.valueOf(greedyNeighbours.size()) + "\n\r" +
					"n neighbourghs: " + String.valueOf(manhattanNeighbours.size()) + "\n\r" +
					"origin coord: " + this.origin.toString() + "\n\r" +
					"destination coord: " + this.destination.toString() + "\n\r";
			//LOGGER.log(Level.SEVERE,msg);
			System.out.print(msg);
			throw e;
	    }
	    return greedyNeighbour;
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
	private Coordinate gridCellToCoordinate(GridCoverage2D grid, GridCoordinates2D cell) {
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
}
