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

}
