package repastInterSim.agent;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunState;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.Geography;
import repastInterSim.environment.OD;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.RoadLink;
import repastInterSim.environment.Vector;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class Vehicle extends MobileAgent {

	private int maxSpeed; // The distance from a vehicle ahead at which the agent adjusts speed to follow
	
	private double speed;
	private double acc;
	private double dmax;	    
	private RoadLink currentRoadLink; // Used for identifying when the vehicle moves from one road link to another
	private Integer queuePos;
	private Route route;


	public Vehicle(int mS, double a, double s, OD o, OD d) {
		super(o, d);
		this.maxSpeed = mS;
		this.tau = 1;
		this.acc = a;
		this.speed = s;
		this.dmax = 20/GlobalVars.spaceScale; // Assuming vehicle drivers adjust their driving according to what's happening 20m in front.
		
		this.destination = d;
		Coordinate dCoord = this.destination.getGeom().getCentroid().getCoordinate();
		
		Context context = RunState.getInstance().getMasterContext();
		Geography<Object> geography = (Geography<Object>) context.getProjection(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY);
		this.route = new Route(geography, this, dCoord);
	}

	/*
	 * Default behaviour, shuffle = true, randomised the scheduling of collections
	 * of agents. In the case of a space, the process may not be random since the
	 * vehicle in the front might have priority.
	 */
	@ScheduledMethod(start = 1, interval = 1, shuffle = false)
	public void step() throws Exception {
		
    	// Check that a route has been generated
    	if (this.route.getRouteX() == null) {
    		this.route.setRoute();
    	}
    	
    	if (this.queuePos==null) {
    		// Set queue pos. If no capacity on link will return null and vehicle will not be able to drive
    		this.setCurrentRoadLinkAndQueuePos(this.route.getRoadsX().get(0));
		}

		// Drive only if vehicle has been added to road link
    	if(this.queuePos!=null) {
    		drive();
    	}
	}
	
	/*
	 * Drive the vehicle agent. 
	 * 
	 * Identify obstacles and set vehicle accelaration with respect to these obstacles.
	 * 
	 * Get the displacement of the vehicle this time step. Update the vehicle's position and speed.  
	 * 
	 */
	public void drive() {
		Context context = RunState.getInstance().getMasterContext();
		Geography<Object> geography = (Geography<Object>) context.getProjection(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY);
		
		// Check for nearby cars, pedestrians and signals
		Vehicle vehicleInFront = getVehicleInFront();
		List<Ped> crossingPeds = this.getCrossingPedestrians();
		List<CrossingAlternative> cas = this.getRoadLinkCrossingAlterantives(this.route.getRoadsX().get(0).getFID());
		
		// Set speed based on vehicle ahead, crossing pedestrians and traffic signal.
		double perturbation = 0;
		double vDes = getDesiredSpeed(vehicleInFront, crossingPeds, cas); 
		this.speed = Math.max(0, vDes-perturbation);
		double distanceToTravel = this.speed * GlobalVars.stepToTimeRatio;
		Coordinate vehicleLoc = this.maLoc;
		
		// get the next coordinate along the route
		double distanceTraveled = 0;
		boolean isFinal = false;
		
		while (Double.compare(distanceToTravel, distanceTraveled) > 0) {
			// Get next coordinate along the route
	        Coordinate routeCoord = this.route.routeX.get(0);
	        RoadLink nextRoadLink = this.route.getRoadsX().get(0);
	        
	        // Is this the final destination?
	        Coordinate destC = this.destination.getGeom().getCentroid().getCoordinate();
	        isFinal = (routeCoord.equals2D(destC));
	        
	        // Calculate the distance to this coordinate
			double distToCoord = vehicleLoc.distance(routeCoord);
			
			// Calculate the angle
			double newBearing = GISFunctions.bearingBetweenCoordinates(maLoc, routeCoord);
			if (Double.isNaN(newBearing)==false) {
				this.bearing = newBearing;
			}
			
			// If vehicle travel distance is too small to get to the next route coordinate move towards next coordinate
			if (Double.compare(distToCoord, distanceToTravel) > 0) {
				// Move agent in the direction of the route coordinate the amount it is able to travel
				vehicleLoc = new Coordinate(vehicleLoc.x + distanceToTravel*Math.sin(this.bearing), vehicleLoc.y + distanceToTravel*Math.cos(this.bearing));
				distanceTraveled += distanceToTravel;
				distanceToTravel = 0;
			}
			// The vehicle is able to travel up to or beyond its next route coordinate
			else {
				vehicleLoc = routeCoord;
				
				// If vehicle has moved distance such that it can enter the next road link, check if the next link has capacity and progress as appropriate
				boolean routeCoordIsOnNextLink = !nextRoadLink.getFID().contentEquals(currentRoadLink.getFID()); 
				boolean canProgressToNextLink=false;
				Integer newQPos = null;
		        if (routeCoordIsOnNextLink) {
		        	// Check if vehicle can move onto the next link. Can't if there is no capacity
		        	// If successfully added will return true
		        	newQPos = nextRoadLink.addVehicleToQueue(this);
		        	if (newQPos!=null) {
		        		canProgressToNextLink = true;
		        	}
		        }
				

		        if (canProgressToNextLink) {
		        	boolean posOK = currentRoadLink.getQueue().readPos() == this.queuePos; // Check that the vehicle that will be removed from the queue is this vehicle
		        	assert posOK;
		        	currentRoadLink.removeVehicleFromQueue();
		        	this.queuePos = newQPos;
		        }
		        else if (routeCoordIsOnNextLink & !canProgressToNextLink) {
		        	isFinal = true; // If can't progress to next link must stop here, for now. Can resume next tick.
		        }
		        
				// If vehicle can't go any further set distanceToTravel to zero
				if (isFinal) {
					// NOTE: this means the distanceAlongRoute isn't the actual distance moved by the vehicle since it was moved up to its final coordinate only and not beyond
					distanceTraveled = distanceToTravel;
					distanceToTravel = 0;
				}
				else {
					distanceTraveled += distToCoord;
					distanceToTravel-=distToCoord;
				}
				
				if( (routeCoordIsOnNextLink==false) | (canProgressToNextLink==true) ) {
					this.route.routeX.remove(routeCoord);
					currentRoadLink = nextRoadLink;
					this.route.getRoadsX().remove(0); // Every route coordinate has its corresponding road link added to roadsX. Removing a link doesn't necessarily mean the vehicle has progressed to the next link.
				}
			}
			
		}
		
		Polygon vehicleGeom = vehicleRectangePolygon(vehicleLoc, this.bearing);
		GISFunctions.moveAgentToGeometry(geography, vehicleGeom, this);
		
		setLoc();
		
	}
	
	/*
	 * Get the desired speed of the of the vehicle with respect to any vehicle, pedestrian or signal obstacles in the immediate vicinity. 
	 */
	public double getDesiredSpeed(Vehicle vif, List<Ped> cPeds, List<CrossingAlternative> cas) {
		
		// Get car following acceleration
		double cfs = carFollowingSafeSpeed(vif);
		
		// Get pedestrian yielding acceleration
		double pys = pedYieldingSafeSpeed(cPeds);
		
		// Get traffic signal acceleration
		double tss = crossingAlternativeSafeSpeed(cas);
		
		// Choose the lowest of these as the vehicle's acceleration		
		double safeSpeed = Math.min(cfs, Math.min(pys, tss));
		
		// Now choose desired speed as the minimum speed between the speed limit, safe following speed and speed if vehicle increased speed by default acceleration.
		double vDes = Math.min(this.maxSpeed, Math.min(safeSpeed, this.speed + GlobalVars.defaultVehicleAcceleration*GlobalVars.stepToTimeRatio));
		
		return vDes;
	}

	/*
	 * Calculate the safe speed of the vehicle agent in response to the traffic signal given by the nearest crossing alternative in front
	 * of the vehicle agent.
	 * 
	 * Doesn't account for leaving space for other cars.
	 */
	public double crossingAlternativeSafeSpeed(List<CrossingAlternative> cas) {
				
		// Get nearest ca in front of vehicle. Crossing alternative has to be closer than dmax to be perceived
		double nearestD = this.dmax;
		CrossingAlternative nearestCAInFront = null;
		for(int i=0; i<cas.size();i++) {
			Coordinate signalLoc = cas.get(i).getSignalLoc();
			if (GISFunctions.coordInFront(this.maLoc, this.bearing, signalLoc)) {
				double d = this.maLoc.distance(signalLoc);
				if (d<nearestD) {
					nearestCAInFront = cas.get(i);
					nearestD = d;
				}
			}
		}
		
		// If no crossing alternative signal in front then return max double value, this prevents vehicle from choosing speed due to crossing alternative
		if (nearestCAInFront==null) {
			return  Double.MAX_VALUE;
		}
		
		// Check for a traffic signal
		char signalState = nearestCAInFront.getState(this.route.getRoadsX().get(0).getFID());
		
		// If signal is green, also return max value so that vehicle ignores signal in its speed choice 
		if (signalState == 'g') {
			return Double.MAX_VALUE;
		} 
		// If signal state is red vehicle must yield to it
		else if (signalState == 'r') {
			double vSafe = safeFollowingSpeedObstacle(nearestCAInFront.getGeom().getCoordinate());
			return vSafe;
		}
		else {
			return Double.MAX_VALUE;
		}
	}
	
	/*
	 * Get safe speed required to avoid collision with crossing pedestrians, assuming pedestrians are not going to yield.
	 * 
	 * @param List<Ped> cPeds
	 * 		A list of pedestrian agents that are crossing the road the vehicle is on.
	 */
	public double pedYieldingSafeSpeed(List<Ped> cPeds) {
		
		// Initialise safe speed as max value, since vehicle chooses minimum between speed limit and safe speed this ensures default value is not chosen.
		double vSafe = Double.MAX_VALUE;
		
		// Check for pedestrians that are in front of the vehicle and within perception distance
		Ped nearestPed = null;
		double pedDist = Double.MAX_VALUE;
		for (int i=0; i<cPeds.size(); i++) {
			
			// If pedestrian is not in front of vehicle then continue
			if (GISFunctions.coordInFront(maLoc, this.bearing, cPeds.get(i).getLoc())==false) {
				continue;
			}
			
			// Find nearest ped within vehicle's perception distance
			double pDist = this.maLoc.distance(cPeds.get(i).getLoc());
			if ( (pDist<pedDist) & (pDist<this.dmax) ) {
				nearestPed = cPeds.get(i);
				pedDist = pDist;
			}
		}
		
		// Finally if crossing ped within perception distance identified get acceleration required to avoid collision with this ped
		if (nearestPed != null) {

			
			// Speed set to zero since vehicle must come to complete stop to avoid collision with peds
			vSafe = safeFollowingSpeedObstacle(nearestPed.getLoc());
		}
		
		return vSafe;
	}

	/*
	 * Calculates safe following speed to the vehicle in front
	 * 
	 * @param vehicleInFront The vehicle immediately in front of the ego vehicle
	 * 
	 * @return The updated acceleration
	 */
	public double carFollowingSafeSpeed(Vehicle vehicleInFront) {
		// initialise safe following speed as Max values, this means that this speed won;t get chosen since it exceeds speed limit
		double vSafe = Double.MAX_VALUE;
	
		// Only do this if there is a vehicle in front to follow
		if (vehicleInFront != null) {
			
			// Get the location and speed of the vehicle in front and use to calculate safe following speed
			vSafe = safeFollowingSpeedVehicle(vehicleInFront.getSpeed(), vehicleInFront.getLoc());
		}

		return vSafe;
	}
	
    public Vehicle getVehicleInFront()  {
    	
    	// Use road link queue to check for any vehicles in front on the current link
    	Vehicle vInFront = this.currentRoadLink.getQueue().getElementAhead(this.queuePos);
    	
    	if (vInFront == null) {
    		// Get vehicle at the back of the road link ahead, returns null if there are no road links ahead.
    		vInFront = getVehicleAtEndOfNextRoadLink();
    	}
    	
    	// If still no vehicle in front return null, otherwise check distance to vehicle in front
    	if (vInFront==null) {
    		return null;
    	}
    	else if (maLoc.distance(vInFront.getLoc()) < this.dmax) {
    		return vInFront;
    	}
    	else {
    		return null;
    	}
    }
	
	/*
	 * Get vehicle at the end of the next road link  by looping over roads in route until the 
	 * next road link is reached and check this queue.
	 */
	public Vehicle getVehicleAtEndOfNextRoadLink() {
		int i = 0;
		String nextRoadID = this.route.getRoadsX().get(i).getFID();
		while ( (this.currentRoadLink.getFID() == nextRoadID) & (i<this.route.getRoadsX().size())) {
			nextRoadID = this.route.getRoadsX().get(i).getFID();
			i++;
		}
		
		// Get vehicle at the back of the road link ahead if there is a link ahead
		if (nextRoadID.contentEquals(this.currentRoadLink.getFID())==false) {
			return this.route.getRoadsX().get(i).getQueue().getEndElement();
		}
		else {
			return null;
		}
	}
	
    
    /*
     * Gets the pedestrians on the current road link that are crossing the road.
     * 
     * Change this to get OR link from ITN link via SpaceBuilder lookups.
     */
    public List<Ped> getCrossingPedestrians() {
    	// Get peds on current road by getting the OR road link associated to the vehicle's current ITN road link
    	RoadLink currentITNLink = this.route.getRoadsX().get(0);
    	List<Ped> crossingPedsOnRoad = new ArrayList<Ped>();
    	
    	// Get crossing peds via OR Link this ITN link is associated to
    	RoadLink orLink = SpaceBuilder.itnToOR.get(currentITNLink);
    	crossingPedsOnRoad = orLink.getPeds().stream().filter(p -> p.isCrossing()).collect(Collectors.toList());
    	return crossingPedsOnRoad;
    }


	/*
	 * Get the crossings that are located on this road link. These objects are also the traffic lights for vehicles.
	 * 
	 * @param String itnRoadLinkID
	 * 		The id of the road link the vehicle is travelling along.
	 *  
	 * @return List<CrossingAlternative> cas.
	 * 		The crossing alternatives that control the traffic flow of this road link.
	 */
	public List<CrossingAlternative> getRoadLinkCrossingAlterantives(String itnRoadLinkID) {
		List<CrossingAlternative> cas = new ArrayList<CrossingAlternative>();
		
		// Agent identifies crossing locations on the road links passed in
		// Loop through these and get crossing alternatives that belong to these road links
		Geography<CrossingAlternative> caGeography = SpaceBuilder.getGeography(GlobalVars.CONTEXT_NAMES.CA_GEOGRAPHY);
		for (CrossingAlternative ca: caGeography.getAllObjects()) {
			if (ca.getITNRoadLinkIDs()==null) continue;
			
			for(int i=0; i< ca.getITNRoadLinkIDs().length; i++) {
				if(ca.getITNRoadLinkIDs()[i].contentEquals(itnRoadLinkID)) {
					cas.add(ca);
				}
			}
		}
		return cas;
	}
	
	/**
	 * Method to be run when removing the agent from the context. 
	 * 
	 * In this case make sure to reduce the count of vehicles on the current road link
	 */
	@Override
	public void tidyForRemoval() {
		this.currentRoadLink.removeVehicleFromQueue();
	}
	
	/*
	 * Calculate the safe following speed for a stationary obstacle such as a pedestrian or traffic signal.
	 * 
	 * In this case the desired distance is fixed to a value given in GlobalVars. The speed of the obstacle is set to zero.
	 * 
	 * @param Coordinate obsLoc
	 * 		The location of the obstacle the vehicle is yielding to
	 * 
	 * @return double
	 * 		The safe following speed
	 */
	public Double safeFollowingSpeedObstacle(Coordinate obsLoc) {
		
		// Calculate distance to obstacle as by resolving the direct distance to the obstacle in the direction of the vehicle
		double obsBearing = GISFunctions.bearingBetweenCoordinates(maLoc, obsLoc);
		double angleBetweenVehicleObstacle = this.bearing - obsBearing;
		
		double d = this.maLoc.distance(obsLoc) * Math.cos(angleBetweenVehicleObstacle);
		
		double safeSpeed = safeFollowingSpeed(GlobalVars.obstacleYieldDistance, d, 0);
		return safeSpeed;	
	}
	
	/*
	 * Calculate the safe following speed from the input vehicle using the simple car following model in the SUMO documentation. 
	 * 
	 * https://sumo.dlr.de/pdf/KraussDiss.pdf
	 * 
	 * @param double leaderSpeed
	 * 		The speed of the agent the vehicle is following
	 * @param Coordinate leaderLoc
	 * 		The location of the agent the vehicle is following
	 * 
	 * @return double
	 * 		The safe following speed
	 */
	public Double safeFollowingSpeedVehicle(double leaderSpeed, Coordinate leaderLoc) {
		// Get desired distance from vehicle in front - driver's reaction time * leader speed
		double dDesired = this.tau * leaderSpeed;
		double d = this.maLoc.distance(leaderLoc) - GlobalVars.vehicleLength;
		double safeSpeed = safeFollowingSpeed(dDesired, d, leaderSpeed);
		return safeSpeed;
	}
	
	/*
	 * A more general method for calculating the safe following speed that takes desired distance to be as an input.   
	 * https://sumo.dlr.de/pdf/KraussDiss.pdf
	 * 
	 * @param double dDesired
	 * 		The desired distance to maintain from agent vehicle is following
	 * @param double leaderSpeed
	 * 		The speed of the agent the vehicle is following
	 * @param Coordinate leaderLoc
	 * 		The location of the agent the vehicle is following
	 * 
	 * @return double
	 * 		The safe following speed
	 */
	public Double safeFollowingSpeed(double dDesired, double d, double leaderSpeed) {
		
		Double vSafe = null;			
		
		// get characteristic time scale used in model
		double tauB = ((this.speed + leaderSpeed) / 2.0) / GlobalVars.defaultVehicleDecceleration;
		
		vSafe = leaderSpeed + (d - dDesired) / (tauB + this.tau);

		return vSafe;
	}
	
	public void setSpeed(double s) {
		this.speed = s;
	}
	
	public double getSpeed() {
		return this.speed;
	}
    
    /*
     * Set the location attribute of the agent to be the coordinate of its 
     * centroid, in the coordinate reference frame used by the agent for GIS calculations. 
     */
	@Override
    public void setLoc()  {
		Context context = RunState.getInstance().getMasterContext();
		Geography<Object> geography = (Geography<Object>) context.getProjection(GlobalVars.CONTEXT_NAMES.MAIN_GEOGRAPHY);
		
    	// Get centroid coordinate of this agent
    	Coordinate vL = GISFunctions.getAgentGeometry(geography, this).getCentroid().getCoordinate();
    	DecimalFormat newFormat = new DecimalFormat("#.#######");
    	vL.x = Double.valueOf(newFormat.format(vL.x));
    	vL.y = Double.valueOf(newFormat.format(vL.y));
    	this.maLoc = vL;
    }
    
    /*
     * Get the destination of this vehicle
     * 
     *  @returns
     *  	The Destination object of this vehicle
     */
	@Override
    public OD getDestination() {
    	return this.destination;
    }
	
    /*
     * Getter for the route
     * 
     * @returns Route of the vehicle
     * 
     */
	public Route getRoute() {
    	return this.route;
    }
    
    public void setCurrentRoadLinkAndQueuePos(RoadLink rl) {
    	this.currentRoadLink = rl;
		this.queuePos = currentRoadLink.addVehicleToQueue(this);
    }
    
    public double getDMax() {
    	return this.dmax;
    }
    
    public double getAcc() {
    	return this.acc;
    }
    
    public String getCurrentRoadLinkID() {
    	return this.currentRoadLink.getFID();
    }
    
    public static Coordinate[] vehicleRectangleCoordiantes(Coordinate loc, double bearing) {
    	// Bearing is clockwise angle from north so use -ve bearing since rotation matrix is for counter clockwise turn
    	double[][] rotationMatrix = {{Math.cos(-bearing), -Math.sin(-bearing)},{Math.sin(-bearing), Math.cos(-bearing)}};
    	
    	// Get points of north facing vehicle rectangle, relative to centre
    	double[] p1 = {-GlobalVars.vehicleWidth/2, GlobalVars.vehicleLength/2};
    	double[] p2 = {GlobalVars.vehicleWidth/2, GlobalVars.vehicleLength/2};
    	double[] p3 = {GlobalVars.vehicleWidth/2, -GlobalVars.vehicleLength/2};
    	double[] p4 = {-GlobalVars.vehicleWidth/2, -GlobalVars.vehicleLength/2};
    	
    	double[][] corners = {p1,p2,p3,p4};
    	double[][] rotatedCorners = new double[4][2];
    	
    	// Rotate corners
    	for(int p=0; p<corners.length; p++) {
    		rotatedCorners[p][0] = corners[p][0];
    		rotatedCorners[p][1] = corners[p][1];
    		for (int i=0; i<rotationMatrix.length; i++) {
    			double[] mRow = rotationMatrix[i];
    			double newc = corners[p][0]*mRow[0] + corners[p][1]*mRow[1];
    			rotatedCorners[p][i] = newc;
    		}
    	}
    	
    	// Get coordinates of corners of vehicle rectangle
    	Coordinate[] recCoords = new Coordinate[rotatedCorners.length+1];
    	for (int i=0; i<corners.length; i++) {
    		Coordinate c = new Coordinate(loc.x+rotatedCorners[i][0], loc.y+rotatedCorners[i][1]);
    		recCoords[i] = c;
    	}
    	recCoords[rotatedCorners.length] = recCoords[0];
    	
    	return recCoords;
    }
    
    public static Polygon vehicleRectangePolygon(Coordinate loc, double bearing) {
    	
    	// Get coordinates of corners of vehicle rectangle
    	Coordinate[] recCoords = vehicleRectangleCoordiantes(loc, bearing);
    	
    	Polygon vehicleRec = null;
    	try {
        	vehicleRec = GISFunctions.polygonGeometryFromCoordinates(recCoords);
    	} catch (IllegalArgumentException e) {
    		e.printStackTrace();
    	}
    	
    	return vehicleRec;
    }
    
    /*
     * Calculate the time to collision to an object at location pLoc travelling at velocity pV.
     */
    public Double TTC(double[] pLoc, double[] pV) {
    	
    	// Get coordinates of the edges of the vehicle geometry
    	Coordinate[] recCorners = vehicleRectangleCoordiantes(this.maLoc, this.bearing);
    	double[] v = {this.speed*Math.sin(this.bearing), this.speed*Math.cos(this.bearing)}; 
    	
    	// Get TTC value for ped and each edge of vehicle, find lowest TTC
    	Double minTTC = null;
    	for (int i=0; i<recCorners.length-1;i++) {
    		double[] e1 = {recCorners[i].x, recCorners[i].y};
    		double[] e2 = {recCorners[i+1].x, recCorners[i+1].y};
    		
    		Double ttc = Vector.edgeTTC(pLoc, pV, e1, e2, v);
    		
    		if (ttc==null) {
    			continue;
    		}
    		else if (ttc<0) {
    			// Only care about conflicts that occur in the future
    			continue;
    		}
    		else if (minTTC==null) {
    			minTTC = ttc;
    		}
    		else if (minTTC>ttc) {
    			minTTC = ttc;
    		}
    	}
    	
    	return minTTC;
    }

}
