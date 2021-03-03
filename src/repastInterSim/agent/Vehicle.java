package repastInterSim.agent;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Point;

import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.space.gis.Geography;
import repastInterSim.environment.OD;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.RoadLink;
import repastInterSim.main.GlobalVars;
import repastInterSim.main.SpaceBuilder;

public class Vehicle extends MobileAgent {

	private int maxSpeed, followDist; // The distance from a vehicle ahead at which the agent adjusts speed to follow
	private double speed;
	private double acc;
	private double dmax;	    
	private RoadLink currentRoadLink; // Used for identifying when the vehicle moves from one road link to another
	private int queuePos;
	private Route route;


	public Vehicle(int mS, double a, double s, OD o, OD d) {
		super(o, d);
		this.maxSpeed = mS;
		this.acc = a;
		this.speed = s;
		this.dmax = 20/GlobalVars.spaceScale; // Assuming vehicle drivers adjust their driving according to what's happening 20m in front.
		
		this.destination = d;
		Coordinate dCoord = this.destination.getGeom().getCentroid().getCoordinate(); 
		this.route = new Route(SpaceBuilder.geography, this, dCoord);
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
    		this.setCurrentRoadLinkAndQueuePos(this.route.getRoadsX().get(0));
		}
    	
		// Check for nearby cars
		Vehicle vehicleInFront = getVehicleInFront();

		// Drive
		drive(vehicleInFront);
		// moveForward();
	}
	
	/*
	 * Drive the vehicle agent. Set the vehicle agent's speed and update the
	 * x-coordinate of the vehicle using its current speed and acceleration and
	 * preventing overtaking. Move the vehicle agent to its new location.
	 * 
	 * Prevention of overtaking is not currently working.
	 * 
	 * @param vehicleInFront Vehicle. The vehicle in front of the agent vehicle
	 * 
	 */
	public void drive(Vehicle vehicleInFront) {

		// Check for a traffic signal
		boolean sigState = true;
		double disp = 0; // Initialise the amount to move the vehicle by

		if (sigState == true) {
			// Continue driving by following car ahead
			// Update acceleration. travel for time step at this acceleration, leading to an updated speed
			setAccFollowing(vehicleInFront);
			disp = this.speed * GlobalVars.stepToTimeRatio + 0.5 * this.acc * Math.pow(GlobalVars.stepToTimeRatio, 2);
			updateSpeed();

			// setAcc(vehicleInFront);
		} 
		/*
		else if (sigState == false) {
			// Set speed based on distance from signal
			// In this case signal will be within a certain distance of the vehicle
			Signal sig = getSignal();
			setAccSignal(sig, vehicleInFront);
			disp = this.speed * GlobalVars.stepToTimeRatio + 0.5 * this.acc * Math.pow(GlobalVars.stepToTimeRatio, 2);
			setSpeed();
		}
		*/
		
		// get the next coordinate along the route
		double distanceAlongRoute = 0;
		
		while (disp > distanceAlongRoute) {
			// Get next coordinate along the route
	        Coordinate routeCoord = this.route.routeX.get(0);
	        RoadLink nextRoadLink = this.route.getRoadsX().get(0);
	        
	        // Is this the final destination?
	        Coordinate destC = this.destination.getGeom().getCentroid().getCoordinate();
	        boolean isFinal = (routeCoord.equals2D(destC));
	        
	        // Calculate the distance to this coordinate
			double distToCoord = maLoc.distance(routeCoord);
			
			// Calculate the angle
			double angleToCoord = GISFunctions.bearingBetweenCoordinates(maLoc, routeCoord);
			
			// If vehicle travel distance is too small to get to the next route coordinate move towards next coordinate
			if (distToCoord > disp) {
				// Move agent in the direction of the route coordinate the amount it is able to travel
				Coordinate newCoord = new Coordinate(maLoc.x + disp*Math.sin(angleToCoord), maLoc.y + disp*Math.cos(angleToCoord));
				Point p = GISFunctions.pointGeometryFromCoordinate(newCoord);
				Geometry g = p.buffer(1); // For now represent cars by 1m radius circles. Later will need to change to rectangles
				GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, g, this);
				distanceAlongRoute += disp;
			}
			// The vehicle is able to travel up to or beyond its next route coordinate
			else {
				// Move to the coordinate and repeat with next coordinate along
				Point p = GISFunctions.pointGeometryFromCoordinate(routeCoord);
				Geometry g = p.buffer(1);
				GISFunctions.moveAgentToGeometry(SpaceBuilder.geography, g, this);
				
				// If vehicle has been moved onto a different road link update the road link queues
		        if (!nextRoadLink.getFID().contentEquals(currentRoadLink.getFID())) {
		        	this.queuePos = nextRoadLink.getQueue().writePos();
		        	assert nextRoadLink.addVehicleToQueue(this); // If successfully added will return true
		        	assert currentRoadLink.getQueue().readPos() == this.queuePos; // Check that the vehicle that will be removed from the queue is this vehicle
		        	currentRoadLink.removeVehicleFromQueue();
		        }
				
				// If this is the final coordinate in the vehicle's route set distance travelled to be the vehicle displacement
				// since the vehicle has now reached the destination and can't go any further
				if (isFinal) {
					// NOTE: this means the distanceAlongRoute isn't the actual distance moved by the vehicle since it was moved up to its final coordinate only and not beyond
					distanceAlongRoute = disp;
				}
				else {
					distanceAlongRoute += distToCoord;	
				}
				
				this.route.routeX.remove(routeCoord);
				currentRoadLink = nextRoadLink;
				this.route.getRoadsX().remove(0); // Every route coordinate has its corresponding road link added to roadsX. Removing a link doesn't necessarily mean the vehicle has progressed to the next link.
			}
			
		}
		
		setLoc();
		
	}
	
	/*
	 * Simply sets speed to be the speed of the vehicle in front accounting for the
	 * acceleration or deceleration required to get to that speed in the next
	 * timestep. This assumes a fixed acceleration.
	 * 
	 * @param vehicleInFront Vehicle. The vehicle agent in front of this vehicle
	 * agent
	 * 
	 * @return Double. The speed set for this vehicle
	 */
	public double setSpeedFollowing(Vehicle vehicleInFront) {

		/*
		 * Set speed so that after one time step it will be the same as the car in front
		 * The vehicle in front is not null only if it is within the following distance
		 * of this agent vehicle.
		 */
		if (vehicleInFront != null) {
			// Get speed of vehicle in front
			double vifSpeed = vehicleInFront.getSpeed();
			this.speed = vifSpeed - (this.acc * GlobalVars.stepToTimeRatio);

		}
		// If there is no vehicle in front just speed up
		else {
			this.speed = this.speed + (this.acc * GlobalVars.stepToTimeRatio);
		}

		enforceSpeedLimit();
		return this.speed;
	}

	public double enforceSpeedLimit() {
		// Enforce speed limits
		if (this.speed > this.maxSpeed) {
			this.speed = this.maxSpeed;
		}
		// Min speed is zero (no-reversing)
		if (this.speed < 0) {
			this.speed = 0;
		}
		return this.speed;
	}

	/*
	 * Set the speed and acceleration of the vehicle agent such that it will come to
	 * a complete stop at the signal.
	 * 
	 * Doesn't account for leaving space for other cars.
	 */
	/*
	public double setAccSignal(Signal s, Vehicle vehicleInFront) {
		double d; // initialise the distance the vehicle must stop in
		double sigX = this.geography.getLocation(s).getX();
		double vX = this.geography.getLocation(this).getX();
		
		// How to ensure there is some sort of buffer
		if (vehicleInFront == null) {
			// Follow the signal
			setAccFollowing(s);
		} else {
			double vifX = this.geography.getLocation(vehicleInFront).getX();

			// Depending on whether the vehicle in front or the signal is closer, set the
			// stopping distance
			if (vifX < sigX) {
				//d = vifX - vX - this.buffer;
				// Follow the vehicle in front
				setAccFollowing(vehicleInFront);
			} else {
				//d = sigX - vX - this.buffer;
				setAccFollowing(s);
	/*
	 * Get acceleration required to avoid collision with crossing pedestrians, assuming pedestrians are not going to yield.
	 * 
	 * @param List<Ped> cPeds
	 * 		A list of pedestrian agents that are crossing the road the vehicle is on.
	 */
	public double pedYieldingAcceleration(List<Ped> cPeds) {
		
		// Initialise acceleration as max value, since vehicle chooses minimum acceleration this ensure default acceleration is not chosen.
		double pedYieldAcc = Double.MAX_VALUE;
		
		// Check for pedestrians that are in front of the vehicle and within perception distance
		Ped nearestPed = null;
		double pedDist = Double.MAX_VALUE;
		for (int i=0; i<cPeds.size(); i++) {
			
			// Get bearing to pedestrian in order to see if it is in front or behind vehicle
			double bearingToPed = GISFunctions.bearingBetweenCoordinates(maLoc, cPeds.get(i).getLoc());
			
			// If difference between bearings is > 90 degs then ped is behind vehicle, move onto next ped
			if (Math.abs(this.bearing - bearingToPed) > Math.PI/2) {
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
			int alpha, m, l;
			alpha = 1;
			m = 0;
			l = 0; // Parameters for the car following model. Needs refactor.
			
			double objectiveVelocity = 0;
			this.acc = (((alpha * Math.pow(this.speed,m)) / Math.pow(maLoc.distance(nearestPed.getLoc()),l)) * (objectiveVelocity - this.speed));
		}
		
		return pedYieldAcc;
	}

	/*
	 * Updates the vehicle's acceleration using the General Motors car following
	 * model described here: {@link
	 * https://nptel.ac.in/courses/105101008/downloads/cete_14.pdf}. This model
	 * might not be suitable for this simple exercise.
	 * 
	 * @param vehicleInFront The vehicle immediately in front of the ego vehicle
	 * 
	 * @return The updated acceleration
	 */
	public double carFollowingAcceleration(Vehicle vehicleInFront) {
		// Update acceleration based on the position and velocity of the vehicle in
		// front.
		
		double objV = 0;
	
		// Only do this if there is a vehicle in front to follow
		if (vehicleInFront != null) {
			int alpha, m, l;
			alpha = 1;
			m = 0;
			l = 0; // Parameters for the car following model. Needs refactor.
			
			Coordinate vifPt = GISFunctions.getAgentGeometry(SpaceBuilder.geography, vehicleInFront).getCentroid().getCoordinate();

			// Acceleration is negative since in order to have caught up to car in front
			// will have been travelling faster
			this.acc = (((alpha * Math.pow(this.speed,m)) / Math.pow(maLoc.distance(vifPt),l)) * (objV - this.speed));
		} else {
			this.acc = GlobalVars.defaultVehicleAcceleration; // Default acceleration
		}

		return this.acc;
	}
	
    public Vehicle getVehicleInFront()  {
    	
    	// Use road link queue to check for any vehicles in front on the current link
    	Vehicle vInFront = this.currentRoadLink.getQueue().getElementAhead(this.queuePos);
    	
    	if (vInFront == null) {
    		// Get vehicle at the back of the road link ahead
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
		while (this.currentRoadLink.getFID() == nextRoadID) {
			i++;
			nextRoadID = this.route.getRoadsX().get(i).getFID();
		}
		
		// Get vehicle at the back of the road link ahead
		return this.route.getRoadsX().get(i).getQueue().getEndElement();
	}
	
    
    /*
     * Gets the pedestrians on the current road link that are crossing the road.
     */
    public List<Ped> getCrossingPedestrians() {
    	// Get peds on current road by getting the OR road link associated to the vehicle's current ITN road link
    	RoadLink currentITNLink = this.route.getRoadsX().get(0);
    	List<Ped> crossingPedsOnRoad = currentITNLink.getRoads().get(0).getORRoadLink().getPeds().stream().filter(p -> p.isCrossing()).collect(Collectors.toList());
    	return crossingPedsOnRoad;
    }


	/*
	 * Get the signal agent in the space continuous space
	 * 
	 * @return Signal. The signal agent
	 */
	/*
	public Signal getSignal() {
		Signal sig = null;
		for (Object o : this.geography.getObjects()) {
			if (o instanceof Signal) {
				sig = (Signal) o;
			}
		}

		return sig;
	}
	*/

	/*
	 * Identify the signal agent in the space and how far away it is. If signal
	 * agent is within a threshold distance get the state of the signal
	 * 
	 * @return Boolean. True if there is no need to adjust driving. False if signal
	 * means stop.
	 */
	/*
	public boolean checkSignal() {
		double threshDist = 5; // Distance at which drivers will alter behaviour depending on signal
		Signal sig = getSignal();

		// First check if vehicle has past the signal, in which case there is no signal
		// to check
		double sigX = this.geography.getLocation(sig).getX();
		double vX = this.geography.getLocation(this).getX();
		if (vX > sigX) {
			return true;
		} else if ((sigX - vX) < threshDist) {
			return sig.getState();
		} else {
			return true;
		}

	}
	*/
	
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
	 * Updates the vehicle's speed using the General Motors car following model
	 * described here: {@link
	 * https://nptel.ac.in/courses/105101008/downloads/cete_14.pdf} In future this
	 * will be revised to ensure a good academic car following model is used
	 * 
	 * @param vehcileInFront The vehicle immediately in front of the ego vehicle
	 * 
	 * @return The new speed
	 */
	public double updateSpeed() {
		// Update velocity
		this.speed = this.speed + this.acc * GlobalVars.stepToTimeRatio;
		
		enforceSpeedLimit();

		return this.speed;

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
    	// Get centroid coordinate of this agent
    	Coordinate vL = GISFunctions.getAgentGeometry(SpaceBuilder.geography, this).getCentroid().getCoordinate();
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
    
    @Override
    public Geography<Object> getGeography() {
    	return SpaceBuilder.geography;
    }
    
    public void setCurrentRoadLinkAndQueuePos(RoadLink rl) {
    	this.currentRoadLink = rl;
		this.queuePos = currentRoadLink.getQueue().writePos();
		currentRoadLink.addVehicleToQueue(this);
    }

}
