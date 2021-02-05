package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;
import repastInterSim.environment.Road;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.UnmarkedCrossingAlternative;

public class TacticalAlternative {
	
	private Ped ped;
	private Geography<CrossingAlternative> caG;
	private Geography<Road> rG;
	
	private List<RoadLink> strategicPath;
	
	private NetworkPath<Junction> nP;
	private Junction currentJunction = null;
	private RepastEdge<Junction> currentEdge = null;
	private Junction endJunction = null;
	private Junction finalJunction = null;
	private LinkedList<RepastEdge<Junction>> routePath = null;
	private LinkedList<Coordinate> routeCoordinates = new LinkedList<Coordinate>();
	private Coordinate destCoordinate = null;
	private List<RepastEdge<Junction>> initPath; // Path that gets agent from the end of the previous road link they were on to the start of the next. Empty when last junction agent reached boarders previous and next road link.
	private List<RepastEdge<Junction>> pathToEnd; // Path that gets agent from start of tactical horizon to end of tactical horizon
	private List<RepastEdge<Junction>> pathRemainder; // Path that gets agent from first link outside tactical horizon to the end of their destination
	private boolean recurringEndJunction = false;
	private boolean isBlank = false;
	private AccumulatorRoute accumulator = new AccumulatorRoute();
	
	public TacticalAlternative() {
		// Blank constructor
		this.isBlank = true;
	}
	
	public TacticalAlternative(Ped p, List<RoadLink> sP, int tNL, List<RepastEdge<Junction>> initTacticalPath, List<RepastEdge<Junction>> firstLinkTacticalPath, List<RepastEdge<Junction>> remainderTacticalPath, Junction startJunction, Geography<CrossingAlternative> caG, Geography<Road> rG) {
		this.ped = p;
		this.caG = caG;
		this.rG = rG;
		
		this.strategicPath = sP;
		
		this.currentJunction = startJunction;
		this.endJunction = null; // Is this needed?
		this.initPath = initTacticalPath;
		this.pathToEnd = firstLinkTacticalPath;
		this.pathRemainder = remainderTacticalPath;
		
		this.routePath = (LinkedList<RepastEdge<Junction>>) Stream.of(this.initPath, this.pathToEnd).flatMap(Collection::stream).collect(Collectors.toList());
		
	}
	
	public Junction getCurrentJunction() {
		return this.currentJunction;
	}
	
	/*
	 * Method for retrieving the coordinate the pedestrian agent should be walking towards.
	 * 
	 * Three possible cases. The tactical route of the pedestrian agent is set by the route junctions. Additional coordinate can be prepended to the route
	 * via the routeCoordinates list. This is used for adding the start and end coordinates of a crossing into the route. Also, when the pedestrian agent is as the end 
	 * of their journey, they will walk towards their destination coordinate rather than the final junction. Otherwise the agent walks towards the next junction in
	 * the route.
	 */
	public Coordinate getTargetCoordinate() {
		if(this.routeCoordinates.size()>0) {
			return routeCoordinates.getLast();
		}
		else if ((this.destCoordinate != null) & (this.routePath.size()==1)) {
			return this.destCoordinate;
		}
		else {
			return this.currentJunction.getGeom().getCoordinate();
		}
	}
	
	/*
	 * First tries to remove the top coordiante from the coordinates stack. If the coordinates stack is empty then the top junction is removed from the junction stack.
	 */
	public void updateTargetCoordiante() {
		if (this.routeCoordinates.size()>0) {
			this.routeCoordinates.removeLast();
		}
		else {
			updateCurrentJunction();
		}
	}
	
	/*
	 * Remove first entry from route junctions and then set current junction to new first entry
	 */
	public void updateCurrentJunction() {
		if (this.recurringEndJunction & (this.routePath.size()==1)) {
			// do not update the current junction
			
		}
		
		else {
			// Update the current edge
			this.routePath.pollFirst();
			this.currentEdge = this.routePath.peekFirst();
			
			// If ped is now at the end of their journey, set the destination coordinate attribute and set recurring end junction to true
			if ((this.strategicPath.size() == 1) & (this.routePath.size()==1)) {
				this.recurringEndJunction = true;
				this.destCoordinate = this.ped.getDestination().getGeom().getCoordinate();
			}
			
			// Identify the next junction
			Junction nextJunction = null;
			if (this.currentEdge == null) {
				this.finalJunction = this.currentJunction;
				nextJunction = null;
			}
			else {
				if (this.currentEdge.getSource().equals(this.currentJunction)) {
					nextJunction = this.currentEdge.getTarget();
				}
				else {
					nextJunction = this.currentEdge.getSource();
				}
				
				// Check if this edge requires crossing a primary link. If it does, initialise an accumulator to choose crossing location				
				List<RoadLink> crossingLinks = tacticalPathPrimaryCrossingLinks(this.currentEdge, this.strategicPath); 
				if (crossingLinks.size()>0) {
					
					Junction defaultJunction = this.noCrossTargetJunction(this.currentJunction, nextJunction);
					
					// Record the desired route path, and overwrite route path with the now default route path
					LinkedList<RepastEdge<Junction>> targetRoutePath = this.routePath;
					RepastEdge<Junction> defaultEdge = this.nP.getNet().getEdge(this.currentJunction, nextJunction);
					this.routePath = new LinkedList<RepastEdge<Junction>>();
					this.routePath.add(defaultEdge);
					
					// Get road length - the length of the road that crossing alternatives are being considered for
					double roadLength = 0;
					for (RoadLink rl: crossingLinks) {
						roadLength += rl.getGeom().getLength();
					}
					
					// Identify crossing alternatives
					List<CrossingAlternative> cas = getCrossingAlternatives(caG, crossingLinks, ped, rG);
					
					// Initialise accumulator crossing choice model
					this.accumulator = new AccumulatorRoute(this.ped, roadLength, defaultJunction, nextJunction, cas, targetRoutePath);
					
					// Set target junction to be the default, no crossing, junction while agent chooses crossing location
					nextJunction = this.accumulator.getDefaultJunction();			
				}
			}
			
			this.currentJunction = nextJunction;
		}
	}
	
	public List<RepastEdge<Junction>> getPathToEnd() {
		return pathToEnd;
	}
	
	public List<RepastEdge<Junction>> getAlternativeRemainderPath() {
		return this.pathRemainder;
	}
	
	/*
	 * Identify the crossing alternatives that lie on the given road links. Prepare these crossing alternatives
	 * for use in the accumulator choice model.
	 * 
	 * @param Geography<CrossingAlternative> caG
	 * 		The geography containing all crossing alternative objects
	 * @param List<RoadLink> rls
	 * 		The road links to identify crossing alternatives along, typically the tactical planning horizon section of the strategic path.
	 * @param Ped p
	 * 		The pedestrian agent perceiving the crossing alternatives
	 * @param Geography<Road> rG
	 * 		The Geography containing Road objects
	 */
	public static List<CrossingAlternative> getCrossingAlternatives(Geography<CrossingAlternative> caG, List<RoadLink> rls, Ped p, Geography<Road> rG){
		
		List<CrossingAlternative> cas = new ArrayList<CrossingAlternative>();
		
		// Agent identifies crossing locations on the road links passed in
		// Loop through these and get crossing alternatives that belong to these road links
		for (RoadLink rl: rls) {
			for (CrossingAlternative ca: caG.getAllObjects()) {
				if (ca.getRoadLinkID().contentEquals(rl.getFID())) {
					// Add to list
					cas.add(ca);
				}
			}
		}
		
		// Add unmarked crossing alternative to list
		// Don't set road for unmarked crossings
		UnmarkedCrossingAlternative caU = new UnmarkedCrossingAlternative();
		caU.setPed(p);
		caU.setRoadGeography(rG);
		caU.setStrategicPathSection(rls);
		cas.add(caU);
		
		return cas;		
	}
	
	/*
	 * Given the source and target junctions of an edge of the tactical network, return the junction that belongs to the same road node
	 * as the target junction, is not the target junction itself, and does not require a road crossing to get to. If no such junction is 
	 * available return null.
	 * 
	 * This method is used to identify the default junction a pedestrian walks towards whilst choosing a crossing location.
	 */
	public Junction noCrossTargetJunction(Junction sourceJunction, Junction targetJunction) {
		Junction noCrossJ = null;
		
		// get road node id
		String roadNodeID = targetJunction.getjuncNodeID();
		
		// Loop through junctions connected to the source junction
		for (Junction j:this.nP.getNet().getAdjacent(sourceJunction)) {
			
			if ( (j.getjuncNodeID().contentEquals(roadNodeID)) & !(j.getFID().contentEquals(targetJunction.getFID())) ) {
				noCrossJ = j;
				break;
			}
		}
		
		return noCrossJ;
	}

	void addCoordinate(Coordinate c) {
		this.routeCoordinates.push(c);
	}
	
	public Junction getEndJunction() {
		return this.endJunction;
	}
	
	public Junction getFinalJunction() {
		return this.finalJunction;
	}

	public void setDestinationCoordinate(Coordinate dC) {
		this.destCoordinate = dC;
	}
	
	public void step() {
		this.accumulator.step();
		
		// If a crossing has been chosen, update the tactical path to reflect this
		if (this.accumulator.getChosenCA() != null) {
			CrossingAlternative ca = this.accumulator.getChosenCA();
			
			// Add the coordinates of the start and end of the crossing to the route
			this.addCoordinate(ca.nearestCoord(this.ped.getLoc()));
			this.addCoordinate(ca.farthestCoord(this.ped.getLoc()));
			
			// Set the current junction to be the target junction - this
			this.currentJunction = this.accumulator.getTargetJunction();
			this.routePath = this.accumulator.getTargetRoutePath();			
			
			// Finally record which crossing type the pedestrian agent chose
			this.ped.setChosenCrossingType(ca.getType());
		}
	}
	
	public boolean isBlank() {
		return this.isBlank;
	}
	
	private List<RoadLink> tacticalPathPrimaryCrossingLinks(RepastEdge<Junction> edge, List<RoadLink> sP) {
		List<RepastEdge<Junction>> path = new ArrayList<RepastEdge<Junction>>();
		path.add(edge);
		return tacticalPathPrimaryCrossingLinks(path, sP);
	}
	
	private List<RoadLink> tacticalPathPrimaryCrossingLinks(List<RepastEdge<Junction>> path, List<RoadLink> sP) {
		List<RoadLink> crossedLinks = new ArrayList<RoadLink>();
		for (RepastEdge<Junction> e: path) {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) e;
			for (RoadLink rl : sP) {
				if (rl.getPedRLID().contentEquals(ne.getRoadLink().getPedRLID())) {
					crossedLinks.add(rl);
				}
			}
		}		
		return crossedLinks;
	}
}