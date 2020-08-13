/*
 * The AccumulatorRoute class is used to model a pedestrian agent's choice of crossing location, given their origin and destination
 */

package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import repastInterSim.agent.Ped;
import repastInterSim.environment.CrossingAlternative;
import repastInterSim.environment.OD;

public class AccumulatorRoute {
	
	private Ped ped;
	
	private List<CrossingAlternative> cas;
	
	public AccumulatorRoute(Ped p) {
		this.ped = p;
		this.cas = new ArrayList<CrossingAlternative>();
		
	}
	/*
	 * Calculate the probability of sampling each crossing alternative using the softmax function
	 * with the argument given by the distance to the crossing alternative from the pedestrian multiplied by the
	 * pedestrian's lambda factor.
	 */
	private List<Double> caSamplingProbabilities(){
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

}
