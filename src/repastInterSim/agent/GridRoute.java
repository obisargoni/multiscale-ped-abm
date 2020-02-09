package repastInterSim.agent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.Vector;
import java.util.stream.Collectors;

import org.geotools.coverage.grid.GridCoordinates2D;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.InvalidGridGeometryException;
import org.geotools.geometry.DirectPosition2D;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;

import repast.simphony.space.gis.Geography;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.PedObstruction;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.SpatialIndexManager;
import repastInterSim.exceptions.RoutingException;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.IO;
import repastInterSim.main.SpaceBuilder;

public class GridRoute extends Route {
	
	// Used to get grid cell summand value when running flood fill algorithm for routing. Single agent can produce routes from different costs, reflecting agent's changing perceptions of costs.
	private HashMap<Integer, Double> gridSummandPriorityMap;
	
	// The grid coordiantes of the agents' route
	private List<GridCoordinates2D> gridPath;
	
	// The filtered list of grid coordinates, with unrequired coordiantes removed
	private List<GridCoordinates2D> prunedGridPath;
	
	// The grid coordinates of crosisng points - used for visualisation
	private List<GridCoordinates2D> gridPathCrossings;	
	
	// Record the coordinates of route points that correspond to crossing locations
	private List<Coordinate> routeCrossingsX;
	
	// Record coordinates at which agent enters new road link
	private List<Coordinate> routeRoadLinkX;
	
	private Map<Coordinate, List<Coordinate>> routeSectionsX;
	
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
	 * Set the route of a mobile agent using the agent's corresponding grid coverage
	 * layer and the flood fill algorithm. This produces a path of coordinates that each correspond
	 * to a grid cell in the coverage. The full path coordinates are pruned to retain
	 * only those that indicate where to turn and points of transition between road priority.
	 * @throws RoutingException 
	 * 
	 */
	public void setPedestrianGridRoute() {
				
		// Initialise class attributes
		this.routeX = new Vector<Coordinate>();
		this.roadsX = new Vector<RoadLink>();
		this.routeDescriptionX = new Vector<String>();
		this.routeSpeedsX = new Vector<Double>();
		this.routeSectionsX = new HashMap<Coordinate, List<Coordinate>>();
		this.routeCrossingsX = new Vector<Coordinate>();
		this.gridPathCrossings = new Vector<GridCoordinates2D>();
		this.prunedGridPath = new Vector<GridCoordinates2D>();
		
		GridCoverage2D grid = geography.getCoverage(GlobalVars.CONTEXT_NAMES.BASE_COVERAGE);

		gridPath = getGridCoveragePath(grid);
		IO.gridCoordiantesIterableToCSV(gridPath, ".\\data\\export\\export_grid_coverage_path.csv");
		Set<Integer> routeIndices = new HashSet<Integer>();
		Set<Integer> crossingIndices = new HashSet<Integer>();
		Set<Integer> roadLinkChangeIndices = new HashSet<Integer>();
		
		double[] prevCellValue = new double[1];
		double[] cellValue = new double[1];
		
		String prevCellRoadLinkFID = null;
		String cellRoadLinkFID = null;
		
		// The first cell in the path corresponds to the agents starting position,
		// therefore don't need to include this coordinate in the route
		GridCoordinates2D prevCell = gridPath.get(0);
		Coordinate prevCellCoord = gridCellToCoordinate(grid, prevCell);
		
		// Include first coord in the roadLinkChangeIndices so agents can update its route at the first road Link (does this mean first coord is on route?)
		roadLinkChangeIndices.add(0);
		
		// Get indices of grid cells that are at location where road priority changes (crossing points)
		for (int i = 1; i < gridPath.size(); i++) {
			GridCoordinates2D gridCell = gridPath.get(i);
			Coordinate cellCoord = gridCellToCoordinate(grid, gridCell);

			// Get grid cell value of this and previous coord. If values differ this means they are located in
			// road space with different priority and therefore the previous grid cell should be included in the route
			prevCellValue = grid.evaluate(prevCell, prevCellValue);
			cellValue = grid.evaluate(gridCell, cellValue);
			Double prevVal = prevCellValue[0];
			Double val = cellValue[0];
			
			try {
				prevCellRoadLinkFID = GISFunctions.getCoordinateRoad(prevCellCoord).getRoadLinkFI();
				cellRoadLinkFID = GISFunctions.getCoordinateRoad(cellCoord).getRoadLinkFI();
			} catch (RoutingException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (NullPointerException e) {
				e.printStackTrace();
			}
			
			// If grid cell value increases, priority has decreased for this agent. Indicates crossing point where yielding is possible
			if (val.compareTo(prevVal) > 0) {
				routeIndices.add(i);
				crossingIndices.add(i);
			}
			// If grid cell value decreases, this indicates this agents' priority is greater. Also crossing point but not one where yielding required
			else if (val.compareTo(prevVal) < 0) {
				routeIndices.add(i);
			}
			
			// Check for change in road link id, if road link id changes add to route
			if (!cellRoadLinkFID.contentEquals(prevCellRoadLinkFID)) {
				roadLinkChangeIndices.add(i);
				routeIndices.add(i);
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
				routeIndices.add(gridPathIndexToIncludeInRoute);
				startCellIndex = gridPathIndexToIncludeInRoute;
				
				// Update start cell
				startCell = gridPath.get(startCellIndex);
				startCoord = gridCellToCoordinate(grid, startCell);
				
				// Loop continues two steps ahead from starting cell but this doesn't change route outcome
			}
		}
		
		// Order final set of route indicies
		List<Integer> routeIndicesSorted = routeIndices.stream().sorted().collect(Collectors.toList());
		List<Integer> routeRoadLinkChangeIndicesSorted = roadLinkChangeIndices.stream().sorted().collect(Collectors.toList());
		Coordinate roadLinkChangeCoord = null;
		for (int i:routeIndicesSorted) {
			GridCoordinates2D routeCell = gridPath.get(i);
			prunedGridPath.add(routeCell);
			Coordinate routeCoord = gridCellToCoordinate(grid, routeCell);
			addToRoute(routeCoord, RoadLink.nullRoad, 1, "grid coverage path");
		}
		
		// Similarly order crossing point indices and get these coordinates
		List<Integer> crossingIndicesSorted = crossingIndices.stream().sorted().collect(Collectors.toList());
		for (int i:crossingIndicesSorted) {
			GridCoordinates2D crossingCell = gridPath.get(i);
			Coordinate crossingCoord = gridCellToCoordinate(grid, crossingCell);
			routeCrossingsX.add(crossingCoord);
			gridPathCrossings.add(crossingCell);
		}
		
		// Finally add the destination as a route coordinate, roadLinkChangeCoord at this point should be the last route change coord
		addToRoute(this.destination, RoadLink.nullRoad, 1, "grid coverage path");
	}
	
	public boolean checkForObstaclesBetweenRouteCoordinates(Coordinate startCoord, Coordinate endCoord) {
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
	public List<GridCoordinates2D> getGridCoveragePath(GridCoverage2D grid){
		
		List<GridCoordinates2D> gridPath = new ArrayList<GridCoordinates2D>();

		DirectPosition2D dpStart = new DirectPosition2D(this.mA.getLoc().x, this.mA.getLoc().y);
		DirectPosition2D dpEnd = new DirectPosition2D(this.destination.x, this.destination.y);
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
		int mini = 0;
		int minj = 0;
		int maxi = grid.getRenderedImage().getTileWidth();
		int maxj = grid.getRenderedImage().getTileHeight();
		if (this.partialFF == true) {
			mini = Math.min(start.x, end.x);
			minj = Math.min(start.y, end.y);
			maxi = Math.max(start.x, end.x);
			maxj = Math.max(start.y, end.y);
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
	public double[][] gridCoverageFloodFill(GridCoverage2D grid, GridCoordinates2D end, int mini, int minj, int maxi, int maxj) {

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
	
	public GridCoordinates2D greedyManhattanNeighbour(GridCoordinates2D cell, double[][] cellValues, List<GridCoordinates2D> path, int mini, int minj, int maxi, int maxj) {
		int width = cellValues[0].length;
		int height = cellValues.length;
		List<GridCoordinates2D> manhattanNeighbours = manhattanNeighbourghs(cell, mini, minj, maxi, maxj);
		
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
	    GridCoordinates2D greedyNeighbour =  greedyNeighbours.get(rand.nextInt(greedyNeighbours.size()));
	    return greedyNeighbour;
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
	private void addToRoute(Coordinate coord, RoadLink nullRoad, double speed, String description) {
		this.routeX.add(coord);
		this.roadsX.add(nullRoad);
		this.routeSpeedsX.add(speed);
		this.routeDescriptionX.add(description);
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
	
	public Coordinate getNextRouteCoord() {
		Coordinate nextRoadLinkCoord = this.routeRoadLinkX.get(0);
		List<Coordinate> nextCoordSection = this.routeSectionsX.get(nextRoadLinkCoord);
		Coordinate nextCoord = nextCoordSection.get(0);
		return nextCoord;
	}
	
	public void removeNextRouteCoord() {
		Coordinate nextRoadLinkCoord = this.routeRoadLinkX.get(0);
		List<Coordinate> nextCoordSection = this.routeSectionsX.get(nextRoadLinkCoord);
		nextCoordSection.remove(0);
		if (nextCoordSection.size() == 0) {
			this.routeRoadLinkX.remove(0);
		}
	}
	
	public List<Coordinate> getRouteCrossings(){
		return this.routeCrossingsX;
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
}
