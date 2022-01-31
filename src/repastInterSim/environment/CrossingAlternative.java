package repastInterSim.environment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

import repast.simphony.space.gis.Geography;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class CrossingAlternative extends Signal implements FixedGeography  {
	
	// Coordinates at which the ca meets pavement
	private long id;
	private Geometry geom = null;
	private Coordinate c1 = null;
	private Coordinate c2 = null;
	private Double crossingBearing = null;
	
	private String type;
	
	// id of the road link this crossing is located on
	protected String roadLinkID;
	protected RoadLink orRoadLink;

	public CrossingAlternative(){
		super();
	}
	
	/*
	 * Calculate the distance to the input coordinate. Distance return is the distance to the nearest
	 * crossing alternative coordiante
	 * 
	 * @param Coordinate loc
	 * 		The coordinate to calculate the distance from
	 */
	public Double distanceTo(Coordinate loc) {
		double d1 = getC1().distance(loc);
		double d2 = getC2().distance(loc);
		return Math.min(d1, d2);
	}
	
	public Coordinate nearestCoord(Coordinate loc) {
		Coordinate[] coords = {getC1(), getC2()};
		
		Coordinate cNear = Arrays.stream(coords).min(Comparator.comparing(c->c.distance(loc))).get();
		
		return cNear;
	}

	public Coordinate farthestCoord(Coordinate loc) {
		Coordinate[] coords = {getC1(), getC2()};
		
		Coordinate cFar = Arrays.stream(coords).max(Comparator.comparing(c->c.distance(loc))).get();
		
		return cFar;
	}	
	
	/*
	 * Vehicles are assumed to yield at dedicated crossings resulting in zero vehicle flow
	 */
	public double getvFlow() {
		return 0.0;
	}
	
	/*
	 * Returns 1 if a conflict would occur between the input pedestrian and vehicles if the pedestrian were 
	 * to cross at this crossing alternative.
	 * 
	 * By definition, pedestrian expects no conflicts at marked crossing alternatives.
	 * 
	 * @param Ped p
	 */
	public int wouldConflictOccur(Ped p) {
		return 0;
	}
	
	/*
	 * Return the number of vehicles passing through the crossing alternative divided by the total vehicle capacity of the road links that pass through the crossing alternative.
	 */
	public double getNormalisedVFlow() {
		double vFlow = this.getvFlow();
		double totalCapacity = 0;
		for (RoadLink rl: this.getCurrentVehicleRoadLinks()) {
			totalCapacity += rl.getQueue().capacity();
		}
		if (vFlow==0) {
			return 0;
		}
		else {
			double normVFlow = vFlow / totalCapacity;
			return normVFlow;	
		}
	}
	
	public HashMap<Vehicle, Double> vehicleTTCs(double[] pLoc, double[] pV) {
		List<RoadLink> itnLinks = this.getCurrentVehicleRoadLinks();
		return vehicleTTCs(pLoc, pV, itnLinks);
	}
	
	/*
	 * Get the time to collision between the vehicles travelling on this road link and and a pedestrian agent.
	 * 
	 * @param double[] pLoc
	 * 		the position of pedestrian agent to calculate ttc to
	 * @param double[] pV
	 * 		The velocity of the pedestrian agent 
	 * @param List<RoadLink> itnLinks
	 * 		The ITN Road Links to look for vehicles on.
	 */
	public HashMap<Vehicle, Double> vehicleTTCs(double[] pLoc, double[] pV, List<RoadLink> itnLinks) {
		
		HashMap<Vehicle, Double> vehicleTTCs = new HashMap<Vehicle, Double>();

		// Loop through all vehicles on road links assocated with this crossing and calculate time to collision to pedestrian, either using peds current location and velocaity,
		// or assuming ped walks from start to end coordinate on crossing.
		for (int i=0; i<itnLinks.size(); i++){
			RoadLink rl = itnLinks.get(i);
			for(int j = 0; j<rl.getQueue().count(); j++){
				int vi = rl.getQueue().readPos() + j;
				if (vi>=rl.getQueue().capacity()) {
					vi = vi-rl.getQueue().capacity();
				}
				Vehicle v = rl.getQueue().elements[vi];
				
				Double ttc = v.TTC(pLoc, pV);
				vehicleTTCs.put(v, ttc);
			}
		}
		
		return vehicleTTCs;
	}
	
	public List<RoadLink> getCurrentVehicleRoadLinks() {
		// Get or link this crossing is on and the neighbournig links
		Network<Junction> orRoadNetwork = SpaceBuilder.getNetwork(GlobalVars.CONTEXT_NAMES.OR_ROAD_NETWORK);
		List<RoadLink> rls = new ArrayList<RoadLink>();
		
		rls.add(this.getORRoadLink());
				
		for (Junction j: this.getORRoadLink().getJunctions()) {
			for (RepastEdge<Junction> e: orRoadNetwork.getEdges(j)) {
				NetworkEdge<Junction> ne = (NetworkEdge<Junction>) e;
				if (rls.contains(ne.getRoadLink())==false) {
					rls.add(ne.getRoadLink());
				}
			}
		}

		List<RoadLink> itnLinks = new ArrayList<RoadLink>();
		for (RoadLink orLink: rls) {
			for (RoadLink itnLink: SpaceBuilder.orToITN.get(orLink)){
				itnLinks.add(itnLink);
			}
		}
		return itnLinks;
	}
	
	public void setRoadLinkID(String rlID) {
		this.roadLinkID = rlID;
	}
	
	public String getRoadLinkID() {
		return this.roadLinkID;
	}
	
	public RoadLink getORRoadLink() {
		if (this.orRoadLink == null) {
			Geography<RoadLink> orRoadLinkGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.OR_ROAD_LINK_GEOGRAPHY);
			for(RoadLink rl: orRoadLinkGeography.getAllObjects()) {
				if (rl.getFID().contentEquals(this.roadLinkID)) {
					this.orRoadLink = rl;
					break;
				}
			}
		}
		return this.orRoadLink;
	}
	
	/*
	 * This method uses the shape file data to set the list of road link ID this crossing alternative's signal controls
	 */
	public void setITNLinkIDs(String ids) {
		String[] itnRLIDs = ids.split(",");
		this.itnLinkIDs = itnRLIDs;
	}
	
	public void setSigPhases(String phases) {
		if(!phases.contentEquals("")) {
			String[] allPhases = phases.split(",");
			
			// Initialise complete phase array. Assume that each phase is of same length as first. 
			char[][] finalAllPhasesArray = new char[allPhases.length][allPhases[0].length()];
			
			for (int i=0; i<allPhases.length; i++) {
				char[] phaseArray = new char[allPhases[i].length()];
				for(int j=0; j<allPhases[i].length(); j++) {
					phaseArray[j] = allPhases[i].charAt(j);
				}
				finalAllPhasesArray[i] = phaseArray;
			}
			
			this.phases = finalAllPhasesArray;
		}
	}
	
	public void setPhaseDurs(String dursString) {
		if(!dursString.contentEquals("")) {
			String[] dursArray = dursString.split(",");
			int[] finalDursArray = new int[dursArray.length];
			for(int i=0; i<dursArray.length; i++) {
				finalDursArray[i] = Integer.parseInt(dursArray[i]);
			}
			this.phaseDurations = finalDursArray;
		}
	}

	public Coordinate getC1() {
		return c1;
	}

	public void setC1(Coordinate c1) {		
		this.c1 = c1;
	}

	public Coordinate getC2() {
		return c2;
	}

	public void setC2(Coordinate c2) {
		this.c2 = c2;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	@Override
	public Geometry getGeom() {
		return this.geom;
	}

	@Override
	public void setGeom(Geometry g) {
		this.geom = g;
		
		int ncoords = g.getCoordinates().length;
		this.c1 = g.getCoordinates()[0];
		this.c2 = g.getCoordinates()[ncoords-1];
		
		this.signalLoc = GISFunctions.midwayBetweenTwoCoordinates(this.c1, this.c2);
	}
	
	public void setID(long id) {
		this.id = id;
	}
	
	public long getID() {
		return this.id;
	}
	
	/*
	 * Used in cleanup to break links between agents
	 */
	public void clear() {
		this.geom = null;
		this.orRoadLink=null;
		this.c1 = null;
		this.c2 = null;
		this.type=null;
		this.roadLinkID=null;
		this.orRoadLink=null;
	}
	
	/*
	 * Calculates and retreives the bearing between C1 and C2
	 */
	public double getCrossingBearing() {
		if (this.crossingBearing==null) {
			this.crossingBearing = GISFunctions.bearingBetweenCoordinates(this.getC1(), this.getC2());
		}
		return this.crossingBearing;
	}
}
