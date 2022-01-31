/*
 * The AccumulatorRoute class is used to model a pedestrian agent's choice of crossing location, given their origin and destination
 */

package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.environment.RunState;
import repast.simphony.parameter.Parameters;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.agent.Ped;
import repastInterSim.agent.Vehicle;
import repastInterSim.datasources.CrossEventData;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdge;

public class AccumulatorRoute {
	
	private Ped ped;
	
	private double roadLength;
	private List<CrossingAlternative> cas = new ArrayList<CrossingAlternative>();
	private double[] caActivations;
	private int nAccumulations = 0; // Used to trigger choice based on time limit 
	
	private RepastEdge<Junction> targetRouteEdge = null;
	private Junction targetJunction;
	private Junction defaultJunction;
	private Coordinate targetCoordinate;
	
	private boolean endOfRoute;
	private boolean directCrossing;
	private boolean caChosen;
	private boolean crossingRequired;
	private boolean isCrossing = false;
	private CrossingAlternative sampledCA = null;
	private double sampledCAve;
	private double sampledCAwt;
	private CrossingAlternative chosenCA = null;
	private LinkedList<Coordinate> crossingCoordinates = new LinkedList<Coordinate>();
	
	private boolean isBlank = false;
	
	private CrossEventData ced = null;
	
	public AccumulatorRoute() {
		// Blank constructor allows ped agent to be initialised with an AccumulatorRoute object which returns a null initial coordinate
		this.isBlank = true;
		this.caChosen = false;
		this.crossingRequired = false;
	}
	
	public AccumulatorRoute(Ped p, double rL, Junction dJ, Junction tJ, Coordinate tC, List<CrossingAlternative> cas, RepastEdge<Junction> tRE, boolean dC) {
		this.ped = p;
		this.roadLength = rL;
		
		this.cas = cas;
		this.caActivations = new double[this.cas.size()];
		for (int i=0;i<this.caActivations.length; i++) {
			this.caActivations[i]=0;
		}
		
		this.targetRouteEdge = tRE;
		this.targetJunction = tJ;
		this.targetCoordinate=tC;
		this.defaultJunction = dJ;
		this.caChosen = false;
		
		// Expect that cas will always contain elements but adding in this check just in case
		if(cas.size()>0) {
			this.crossingRequired = true;
		}
		else {
			this.crossingRequired = false;
		}
		
		this.directCrossing = dC;
		
		if (this.ped.getDestination().getGeom().getCoordinate().equals2D(targetCoordinate)) {
			this.endOfRoute=true;
		}
		else {
			this.endOfRoute=false;
		}
	}
	/*
	 * Calculate the probability of sampling each crossing alternative using the softmax function
	 * with the argument given by the distance to the crossing alternative from the pedestrian multiplied by the
	 * pedestrian's lambda factor.
	 */
	public List<Double> caSamplingProbabilities(){
		List<Double> probs = new ArrayList<Double>();
		
		double lambda = this.ped.getLambda();
		
		// Loop through crossing alternatives and find distance from ped to each one
		List<Double> distances = new ArrayList<Double>();
		for (CrossingAlternative ca: this.cas) {
			Double dTo = ca.distanceTo(ped.getLoc());
			
			Double dSalience = (this.roadLength - Math.abs(dTo)) / this.roadLength;
			
			// Multiply distance by ped attribute lambda to control effect of distance on probability
			distances.add(dSalience*lambda);
		}
		
		// Calculate sampling probability using softmax function
	    double total = distances.stream().map(Math::exp).mapToDouble(f->f).sum();
	    probs = distances.stream().map(f-> (Math.exp(f) / total)).collect(Collectors.toList());
		
		return probs;
	}
	
	/*
	 * Sample a single crossing alternative at random according to the probabilities given by the crossing alterantives'
	 * salience, defined by their distance to the pedestrian agent.
	 */
	public int sampleCA() {
		
		int nCAs = this.cas.size();
		List<Double> probs = caSamplingProbabilities();
		
		// sample random number between 0 and 1 and use to sample a CA
		Double r = RandomHelper.getDistribution("caSampleDistribution").nextDouble();
		
		Integer sampledi = null;
		double sum_prob = 0;
		for (int i=0; i<nCAs; i++) {
			sum_prob += probs.get(i);
			if (sum_prob > r) {
				sampledi = i;
				break;
			}
		}
		
		// If sampledCA is still null means last CA should have been chosen but wasn't due to double comparisons
		if (sampledi == null) {
			sampledi = nCAs-1;
		}
		return sampledi;
	}
	
	/*
	 * Binary indicator of whether or not crossing at this crossing alternative now will result in a vehicle conflict if
	 * if pedestrian agent and vehicle agents do not alter course. Returns 1 if no conflicts would occur, 0 otherwise.
	 * 
	 * @param CrossingAlternative ca
	 */
	public int isCAConfictsFree(CrossingAlternative ca) {
		return 1 - ca.wouldConflictOccur(this.ped);
	}
	
	/*
	 * Get the crossing exposure indicator value for the input crossing alternative
	/*
	 * Returns 1 if there is zero vehicle flow through the crossing line during the time it would take
	 * the pedestrian to cross the road. 0 otherwise, indicating the presence of vehicles.
	 * 
	 * @param CrossingAlternative ca
	 */
	public double caVehicleExposureIndicator(CrossingAlternative ca) {
		
		double vFlow = ca.getvFlow();
		if (vFlow==0) {
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	
	/*
	 * Get the roadside walk time indicator for the input crossing alternative.
	 * 
	 * @param CrossingAlternative ca
	 */
	public double caRoadsideDetourIndicator(CrossingAlternative ca) {
		
		// If the input crossing alternative is unmarked then return 1, as by definition there is no detour using this alternative
		if (ca.getType().contentEquals("unmarked")) {
			return 1;
		}
		
		// Get walk time to crossing alternative and from crossing alternative to destination
		Double dToCA = ca.distanceTo(this.ped.getLoc());
		Double dFromCAToDest = ca.distanceTo(this.targetCoordinate);
				
		// Get the unmarked crossing alternative and calculate the distance to dest using it
		CrossingAlternative umCA = null;
		for (CrossingAlternative ca_: this.cas) {
			if (ca_.getType().contentEquals("unmarked")) {
				umCA = ca_;
			}
		}
		
		Double umDToCA = umCA.distanceTo(this.ped.getLoc());
		Double umDFromCAToDest = umCA.distanceTo(this.targetCoordinate);
		
		double detourDist = ((dToCA + dFromCAToDest) - (umDToCA + umDFromCAToDest));
				
		// Need characteristic walk time to compare this to - use the length of the roads in planning horizon
		double charDist = this.roadLength;
		
		return 1 - (detourDist / charDist);
	}
	
	double caUtility(CrossingAlternative ca) {
		
		double caWT = caRoadsideDetourIndicator(ca);
		double caVE = caVehicleExposureIndicator(ca);
				
		this.sampledCA = ca;
		this.sampledCAve = caVE;
		this.sampledCAwt = caWT;
		
		return this.ped.getAlpha()*caWT + (1-this.ped.getAlpha())*caVE;
	}
	
	
	/*
	 * Method to accumulate activation for crossing alternatives
	 */
	public void accumulateCAActivation() {
		
		// Get index of sampled crossing alternative
		int sampledCAi = sampleCA();
		
		// Loop through activations, decay all by decay factor and increase sampled CA by its utility
		for (int i=0; i<this.caActivations.length; i++) {
			this.caActivations[i] = this.caActivations[i] * this.ped.getGamma();
			
			if (i==sampledCAi) {
				double ui = caUtility(this.cas.get(i));
				this.caActivations[i] += ui;
			}
		}
	}
	
	/*
	 * Check for a dominant activation. If one activation is dominant choose that crossing alternative
	 */
	
	public void chooseCA() {
		
		boolean chooseCA=false;
		
		// Sort the activations
		int nCAs = this.caActivations.length;
		double[] sortedActivations = this.caActivations.clone();
		Arrays.sort(sortedActivations);
		
		// Compare the largest activation to the second largest
		double maxActivation = sortedActivations[nCAs-1];
		
		if ((this.directCrossing | this.endOfRoute) & (this.nAccumulations>ped.getTimeThreshold())) {
			// Force choice of CA with highest activation
			chooseCA=true;
		}
		else if (maxActivation > this.ped.getEpsilon()) {
			// Choose CA if activation exceeds threshold
			chooseCA=true;
		}
		
		// Choose top CA
		if (chooseCA) {
			// initialise variable to record index of chosen crossing alterantive
			Integer choseni = null;
			
			// Choose ca by iterating through activations and finding index of the maximum activation
			for (int i=0; i< nCAs; i++) {
				double a = this.caActivations[i];
				if (a == sortedActivations[nCAs-1]) {
					choseni = i;
				}
			}
			
			// With crossing alternative chosen, set caChosen to true and add the crossing coordinates to the list. Record which crossing the pedestrian chose
			this.setChosenCA(this.cas.get(choseni));
		}
	}
	
	
	/*
	 *  Update the crossing alternative activations and choose a crossing alternative. Also update the ped's ttc with vehicles, only changes when ped is crossing..
	 */
	public void step() {
		// Accumulate activation and chooses a crossing alternative if choice not made yet
		if ((isBlank == false) & (this.caChosen==false) & (this.crossingRequired==true)) {
			this.accumulateCAActivation();
			this.nAccumulations++;
			this.chooseCA();
		}
		
		// Update peds ttc to vehicles when ped is crossing. Once ped stops crossing TTC is set back to null
		if (isCrossing) {
			double[] pLoc= {ped.getLoc().x, ped.getLoc().y};
			HashMap<Vehicle, Double> ttcs = this.chosenCA.vehicleTTCs(pLoc, ped.getV());
			List<Double> values =  ttcs.values().stream().filter(x->x!=null).collect(Collectors.toList());
			if (values.size()>0) {
				double ttc = values.stream().min(Comparator.comparing(Double::valueOf)).get();
				this.ced.setTTC(ttc);
			}
			else {
				this.ced.setTTC(null);
			}
		}
	}
	
	public CrossingAlternative getChosenCA() {
		return this.chosenCA;
	}
	
	public void setChosenCA(CrossingAlternative ca) {
		this.chosenCA = ca;
		this.caChosen = true;
		this.crossingCoordinates.push(this.chosenCA.nearestCoord(this.ped.getLoc()));
		this.crossingCoordinates.push(this.chosenCA.farthestCoord(this.ped.getLoc()));
		this.ped.setYield(false);
	}
	
	public Junction getTargetJunction() {
		return this.targetJunction;
	}
	
	public RepastEdge<Junction> getTargetRouteEdge() {
		return this.targetRouteEdge;
	}
	
	public Junction getDefaultJunction() {
		return this.defaultJunction;
	}
	
	public boolean isBlank() {
		return this.isBlank;
	}
	
	public boolean caChosen() {
		return this.caChosen;
	}
	
	public boolean crossingRequired() {
		return this.crossingRequired;
	}
	
	public List<CrossingAlternative> getCAs() {
		return this.cas;
	}
	
	public LinkedList<Coordinate> getCrossingCoordinates() {
		return this.crossingCoordinates;
	}
	
	public void removeCrossingCoordinate() {
		// if at the start of crossing, initialised a new CrossEventData object
		if (this.isCrossing==false) {
			// Create crossing data collector agent
			CrossEventData ced = new CrossEventData(this.ped.getID(), this.targetRouteEdgeFID(), this.getChosenCA().getType(), this.crossingCoordinates);
			RunState.getInstance().getMasterContext().add(ced);
			this.ced=ced;
		}
		
		// Once crossing coordinates have been updated once, ped has started to cross.
		this.isCrossing = true;
		
		this.crossingCoordinates.removeLast();
		
		// if all crossing coordinates have been passed then agent has finished crossing the road and a crossing is no longer required.
		if(this.crossingCoordinates.size()==0) {
			this.crossingRequired = false;
			this.caChosen = false;
			this.isCrossing = false;
			
			// Remove crossing data collector from context and set to null
			RunState.getInstance().getMasterContext().remove(this.ced);
			this.ced=null;
			
		}
	}
	
	public boolean isCrossing() {
		return this.isCrossing;
	}
	
	public boolean isDirectCrossing() {
		return this.directCrossing;
	}
	
	public String caActivationString() {
		String output = "";
		if (this.isBlank==false) {
			for (int i=0; i< this.caActivations.length; i++) {
				output += Double.toString(this.caActivations[i]) + ";";
			}
		}
		return output;
	}
	
	public String caString() {
		String output = "";
		if (this.isBlank==false) {
			for (int i=0; i< this.cas.size(); i++) {
				CrossingAlternative ca = this.cas.get(i);
				String caString = ca.getType()+Long.toString(ca.getID()) + ";";
				output += caString;
			}
		}
		return output;
	}
	
	public String sampledCAString() {
		String output = "";
		if (this.sampledCA!=null) {
			output += this.sampledCA.getType() + Long.toString(this.sampledCA.getID()) + ";";
			output += Double.toString(this.sampledCAwt) + ";";
			output+= Double.toString(this.sampledCAve) + ";";
		}
		return output;
	}
	
	public String targetRouteEdgeFID() {
		if (this.isBlank) {
			return "";
		}
		else {
			NetworkEdge<Junction> ne = (NetworkEdge<Junction>) this.targetRouteEdge;
			return ne.getRoadLink().getFID();
		}
	}
	
	public void clear() {
		this.ped=null;
		this.cas = null;
		this.targetRouteEdge = null;
		this.targetJunction=null;
		this.defaultJunction=null;
		this.chosenCA = null;
		this.crossingCoordinates = null;
	}

}
