package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.collections15.Transformer;

import repast.simphony.space.graph.RepastEdge;
import repastInterSim.environment.NetworkEdge;

public abstract class PathFindingEdgeTransformers {

	
	public class CrossesRoadTransformer<T> implements Transformer<RepastEdge<T>,Double> {
		
		// The road link id to return a different weight value for if the edge crosses them
		private List<String> roadLinkIDs = new ArrayList<String>();;
		private double crossesRoadLinkWeight = Double.MAX_VALUE;
		
		public CrossesRoadTransformer() {
			// TODO Auto-generated constructor stub
		}
		
		public CrossesRoadTransformer(List<String> rlIDs) {
			this.roadLinkIDs = rlIDs;
		}
		
		public CrossesRoadTransformer(List<String> rlIDs, double w) {
			this.roadLinkIDs = rlIDs;
			this.crossesRoadLinkWeight = w;
		}
		
		@Override
		public Double transform(RepastEdge<T> edge) {
			NetworkEdge<T> ne = (NetworkEdge<T>) edge;
			
			boolean intersectsPedRoad = roadLinkIDs.stream().anyMatch(id -> id.contentEquals(ne.getRoadLink().getPedRLID()));
			if (intersectsPedRoad) {
				return this.crossesRoadLinkWeight;
			}
			else {
				return edge.getWeight();
			}
		}
	}
	
	

}