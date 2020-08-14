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

public class AccumulatorRoute {
	
	private Coordinate defaultDestination; // The destination pedestrian walks towards while choosing crossing alternative
	private Ped ped;
	
	private List<CrossingAlternative> cas;
	private double[] caActivations;
	
	private List<Coordinate> routeX; // The list of coordinates that form the pedestrian agent's tactical route
	
	private boolean caChosen = false;
	
	public AccumulatorRoute(Ped p, Coordinate defD) {
		this.defaultDestination = defD;
		this.ped = p;
		this.cas = new ArrayList<CrossingAlternative>();
		this.caActivations = new double[this.cas.size()];
		for (int i=0;i<this.caActivations.length; i++) {
			this.caActivations[i]=0;
		}
		
		routeX = new ArrayList<Coordinate>();
		routeX.add(defaultDestination);
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
			Double d = ca.distanceTo(ped.getLoc());
			
			// Multiply distance by ped attribute lambda to control effect of distance on probability
			distances.add(d*lambda);
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
		
		// Calculate number of vehicles that will pass through ca during that time
		// Need to write methods to do this
		// get speed and position of cars, use to calculate whether they will pass through crossing in allotted time
		// might be tricky because might require checking cars on multiple links, perhaps within the planning horizon. But what about around corners?
		
		// just use constant value for now
		double vFlow = 2;
		
		double avVFlow = 1;
		
		return 1 - (vFlow/avVFlow);
	}
	
	/*
	 * Get the roadside walk time indicator for the input crossing alternative.
	 * 
	 * @param CrossingAlternative ca
	 */
	public double caRoadsideWalkTimeIndicator(CrossingAlternative ca) {
		
		// Get walk time to crossing alternative and from crossing alternative to destination
		Double dToCA = ca.distanceTo(this.ped.getLoc());
		Double dFromCAToDest = ca.distanceTo(ca.getDestination());
		
		double walkTime = (dToCA + dFromCAToDest) / this.ped.getSpeed();
		
		// Need characteristic walk time to compare this to
		double charWT = this.ped.getLoc().distance(this.defaultDestination) / this.ped.getSpeed();
		
		return 1 - (walkTime / charWT);
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
		double[] sortedActivations = this.caActivations;
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
			
			this.routeX = new ArrayList<Coordinate>();
			this.routeX.add(chosenCA.nearestCoord(this.ped.getLoc()));
			this.routeX.add(chosenCA.farthestCoord(this.ped.getLoc()));
			this.routeX.add(chosenCA.getDestination());
			
			// Update bool indicating whether crossing choice has been made or not
			this.caChosen = true;
		}
	}
	
	public List<Coordinate> getRouteX() {
		return routeX;
	}
	
	public boolean isCrossingChosen() {
		return this.caChosen;
	}
	

}
