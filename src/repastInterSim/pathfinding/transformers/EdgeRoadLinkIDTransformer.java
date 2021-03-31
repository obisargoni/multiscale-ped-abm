package repastInterSim.pathfinding.transformers;

import org.apache.commons.collections15.Transformer;

import repast.simphony.space.graph.RepastEdge;
import repastInterSim.environment.NetworkEdge;

/*
 * This transformer use a network edge has an attribute that indicates the edge crosses a road to determine the weight of the edge.
 */
public class EdgeRoadLinkIDTransformer<T> implements Transformer<RepastEdge<T>,Integer> {
	
	// The road link id to return a different weight value for if the edge crosses them
	private int crossesRoadLinkWeight = 1;
	private int notCrossesRoadLinkWeight = 0;
	
	public EdgeRoadLinkIDTransformer() {
		// TODO Auto-generated constructor stub
	}
	
	public EdgeRoadLinkIDTransformer(int crossesWeight, int notCrossesWeight) {
		this.crossesRoadLinkWeight = crossesWeight;
		this.notCrossesRoadLinkWeight = notCrossesWeight;
	}
	
	@Override
	public Integer transform(RepastEdge<T> edge) {
		NetworkEdge<T> ne = (NetworkEdge<T>) edge;
		
		if (ne.getRoadLink().getPedRLID().contentEquals("")) {
			return this.notCrossesRoadLinkWeight;
		}
		else {
			return this.crossesRoadLinkWeight;
		}
	}
}