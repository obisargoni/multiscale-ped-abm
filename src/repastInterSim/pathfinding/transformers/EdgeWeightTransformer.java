package repastInterSim.pathfinding.transformers;

import org.apache.commons.collections15.Transformer;

import repast.simphony.space.graph.RepastEdge;
import repastInterSim.environment.NetworkEdge;

public class EdgeWeightTransformer<T> implements Transformer<RepastEdge<T>,Integer> {
	
	
	public EdgeWeightTransformer() {
		// TODO Auto-generated constructor stub
	}
	
	
	@Override
	public Integer transform(RepastEdge<T> edge) {
		NetworkEdge<T> ne = (NetworkEdge<T>) edge;
		
		return (int) ne.getWeight();
	}
}	