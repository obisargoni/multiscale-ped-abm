/*
 * The AccumulatorRoute class is used to model a pedestrian agent's choice of crossing location, given their origin and destination
 */

package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import com.vividsolutions.jts.geom.Coordinate;

import repastInterSim.agent.Ped;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.Junction;

public class AccumulatorRoute {
	
	private Ped ped;
	
	private TacticalAlternative targetTR = new TacticalAlternative();
	private TacticalAlternative currentTR = new TacticalAlternative();
	
	private double roadLength;
	private List<CrossingAlternative> cas;
	private double[] caActivations;
	
	private Coordinate targetDest;
		
	private boolean caChosen = false;
	
	private boolean isBlank = false;
	
	public AccumulatorRoute() {
		// Blank constructor allows ped agent to be initialised with an AccumulatorRoute object which returns a null initial coordinate
		this.isBlank = true;
	}
	
	public AccumulatorRoute(Ped p, double rL, TacticalAlternative dTR, TacticalAlternative tTR) {
		this.ped = p;
		this.roadLength = rL;
		
		// Agent starts by using the default tactical route but wants to use the target tactical route and must decide where to switch between these
		this.currentTR = dTR;
		this.targetTR = tTR;
		
		this.cas = this.targetTR.getCrossingAlternatives();
		this.caActivations = new double[this.cas.size()];
		for (int i=0;i<this.caActivations.length; i++) {
			this.caActivations[i]=0;
		}
		
		this.targetDest = this.targetTR.getEndJunction().getGeom().getCoordinate();
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
		Double r = new Random().nextDouble();
		
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
	 * Get the crossing exposure indicator value for the input crossing alternative
	 * 
	 * @param CrossingAlternative ca
	 */
	public double caVehicleExposureIndicator(CrossingAlternative ca) {
		
		// Time required for ped to cross the road
		double t_cross = ca.getC1().distance(ca.getC2()) / this.ped.getSpeed();
		
		double vFlow = ca.getvFlow();
		double avVFlow = 1;
		return 1 - (vFlow/avVFlow);
	}
	
	/*
	 * Get the roadside walk time indicator for the input crossing alternative.
	 * 
	 * @param CrossingAlternative ca
	 */
	public double caRoadsideWalkTimeIndicator(CrossingAlternative ca) {
		
		// If the input crossing alternative is unmarked then return 1, as by definition there is no detour using this alternative
		if (ca.getType().contentEquals("unmarked")) {
			return 1;
		}
		
		// Get walk time to crossing alternative and from crossing alternative to destination
		Double dToCA = ca.distanceTo(this.ped.getLoc());
		Double dFromCAToDest = ca.distanceTo(this.targetDest);
				
		// Get the unmarked crossing alternative and calculate the distance to dest using it
		CrossingAlternative umCA = null;
		for (CrossingAlternative ca_: this.cas) {
			if (ca_.getType().contentEquals("unmarked")) {
				umCA = ca_;
			}
		}
		
		Double umDToCA = umCA.distanceTo(this.ped.getLoc());
		Double umDFromCAToDest = umCA.distanceTo(this.targetDest);
		
		double detourWalkTime = ((dToCA + dFromCAToDest) - (umDToCA + umDFromCAToDest)) / this.ped.getSpeed();
				
		// Need characteristic walk time to compare this to - use the length of the roads in planning horizon
		double charWT = this.roadLength / this.ped.getSpeed();
		
		return 1 - (detourWalkTime / charWT);
	}
	
	double caUtility(CrossingAlternative ca) {
		
		double caWT = caRoadsideWalkTimeIndicator(ca);
		double caVE = caVehicleExposureIndicator(ca);
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
		
		// Sort the activations
		int nCAs = this.caActivations.length;
		double[] sortedActivations = this.caActivations.clone();
		Arrays.sort(sortedActivations);
		
		// Compare the largest activation to the second largest
		double diff = sortedActivations[nCAs-1] - sortedActivations[nCAs-2];  
		if (diff > this.ped.getEpsilon()) {
			
			// initialise variable to record index of chosen crossing alterantive
			Integer choseni = null;
			
			// Choose ca by iterating through activations and finding index of the maximum activation
			for (int i=0; i< nCAs; i++) {
				double a = this.caActivations[i];
				if (a == sortedActivations[nCAs-1]) {
					choseni = i;
				}
			}
			
			// With crossing alternative chosen, update the tactical path
			CrossingAlternative chosenCA = this.cas.get(choseni);
			
			// Once crossing choice is made, need to update which tactical route the ped must follow
			// Also need to update the tactical path of the tactical route
			this.targetTR.updatePathToEnd(this.currentTR.getCurrentJunction());
			this.targetTR.addCoordinate(chosenCA.nearestCoord(this.ped.getLoc()));
			this.targetTR.addCoordinate(chosenCA.farthestCoord(this.ped.getLoc()));
			this.currentTR = this.targetTR;
			
			// Update bool indicating whether crossing choice has been made or not
			this.caChosen = true;
			this.ped.setChosenCrossingType(chosenCA.getType());
		}
	}
	
	/*
	 *  Update the crossing alternative activations and choose a crossing alternative
	 */
	public void step() {
		// Accumulate activation and chooses a crossing alternative if choice not made yet
		if ((this.caChosen==false) & (this.cas.size()>0)) {
			this.accumulateCAActivation();
			this.chooseCA();
		}
	}
	
	public Coordinate targetCoordinate() {
		return this.currentTR.getTargetCoordinate();	
	}
	
	public void updateTargetCoordiante() {
		this.currentTR.updateTargetCoordiante();
	}
	
	public Junction getCurrentJunction() {
		return this.currentTR.getCurrentJunction();
	}
	
	public boolean isCrossingChosen() {
		return this.caChosen;
	}
	
	public TacticalAlternative getCurrentTA() {
		return this.currentTR;
	}
	
	public TacticalAlternative getTargetTA() {
		return this.targetTR;
	}
	
	public boolean isBlank() {
		return this.isBlank;
	}

}
