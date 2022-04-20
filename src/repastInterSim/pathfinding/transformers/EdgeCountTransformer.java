package repastInterSim.pathfinding.transformers;

import org.apache.commons.collections15.Transformer;

import repast.simphony.space.graph.RepastEdge;

public class EdgeCountTransformer<T> implements Transformer<RepastEdge<T>,Integer> {

	public EdgeCountTransformer() {
		// TODO Auto-generated constructor stub
	}
	
	@Override
	public Integer transform(RepastEdge<T> edge) {
		return 1;
	}	

}	