package repastInterSim.pathfinding;

import java.util.List;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import repast.simphony.context.space.graph.ContextJungNetwork;
import repast.simphony.space.graph.JungEdgeTransformer;
import repast.simphony.space.graph.JungNetwork;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repast.simphony.space.projection.ProjectionEvent;
import repast.simphony.space.projection.ProjectionListener;

/*
 * An amended version of the Repast Simphony ShortestPath class that allows bespoke transformers to be used in the shortest path calculation
 */
public class ShortestPath<T> implements ProjectionListener<T> {
	
		private Network<T> net;
		private boolean calc = true;
	  private Transformer<RepastEdge<T>,Double> transformer;
	  private DijkstraShortestPath<T,RepastEdge<T>> dsp;
	  
	  /**
	   * Constructor
	   * 
	   * @param net the Network
	   */
	  
	    public ShortestPath(){
	    	
	    }
		public ShortestPath(Network<T> net){
			init(net);
		}
		
		private void init(Network<T> net){
			this.net = net;
			transformer = new JungEdgeTransformer<T>();
			net.addProjectionListener(this);
		}

		/**
		 * Returns a list of RepastEdges in the shortest path from source to target.
		 * 
		 * @param source
		 * @param target
		 * @return
		 */
		public List<RepastEdge<T>> getPath(T source, T target){
			
			if (calc){
				calcPaths();
				calc = false;
			}
			return dsp.getPath(source, target); 
		}
		
		/**
		 * Gets the path length from the source node to the target node.
		 * 
		 * @param source the node we want to get the path length from
		 * @param target the node we want to get the path length to
		 * @return the path length from the source node to the target node.
		 */
		public double getPathLength(T source, T target){
			if (calc){
				calcPaths();
				calc = false;
			}
			
			Number n = dsp.getDistance(source, target);

			if (n != null)
				return n.doubleValue();
			else
				return Double.POSITIVE_INFINITY;
		}
		
		/**
		 * Creates shortest path info  nodes using the Jung Dijkstra algorithm
		 */
		private void calcPaths(){
			Graph<T, RepastEdge<T>> graph = null;
			
			if (net instanceof JungNetwork)
				graph = ((JungNetwork<T>)net).getGraph();
			else if (net instanceof ContextJungNetwork)
				graph = ((ContextJungNetwork<T>)net).getGraph();
			
			dsp = new DijkstraShortestPath<T,RepastEdge<T>>(graph, transformer);
		}
		
		
		/*
		 * Set the transformer
		 */
		public void setTransformer(Transformer<RepastEdge<T>,Double> t) {
			this.transformer = t;
			this.calc = true;
		}
		
		/*
		 * Set the transformer to be the default transformer the class is initialised with
		 */
		public void setDefaultTransformer() {
			this.transformer = new JungEdgeTransformer<T>();
			this.calc = true;
		}
		
		/**
		 * Called when the network is modified so that this will recalculate the
		 * shortest path info.
		 * 
		 * @param evt
		 */
		public void projectionEventOccurred(ProjectionEvent<T> evt) {
			if (evt.getType() != ProjectionEvent.OBJECT_MOVED) {
				calc = true;
			}
		}
		
		/**
		 * Removes this as a projection listener when this ShortestPath is garbage
		 * collected.
		 */
		public void finalize() {
			if (net != null)
				net.removeProjectionListener(this);
		}
		
		/**
		 * Null the object so that the Garbage Collector recognizes to remove 
		 * the object from the jvm.
		 */
		public static ShortestPath<?> finished(ShortestPath<?> sp){
			sp.finalize();
			sp=null;
			return sp;
		}
		
}
