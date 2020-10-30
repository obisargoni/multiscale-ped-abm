package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
import java.util.function.Predicate;
import java.util.stream.Collectors;

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
public class NetworkPath<T> implements ProjectionListener<T> {
	
		private Network<T> net;
		private boolean calc = true;
		private Transformer<RepastEdge<T>,Double> transformer;
		private DijkstraShortestPath<T,RepastEdge<T>> dsp;
	  
		private Stack<T> connectionPath;
		private List<Stack<T>> connectionPaths;
	  
	  /**
	   * Constructor
	   * 
	   * @param net the Network
	   */
	  
	    public NetworkPath(){
	    	
	    }
		public NetworkPath(Network<T> net){
			init(net);
		}
		
		private void init(Network<T> net){
			this.net = net;
			resetConnectionPaths();
			transformer = new JungEdgeTransformer<T>();
			net.addProjectionListener(this);
		}
		
		/*
		 * Initialise the connection path objects used to record the paths between nodes in the network
		 */
		public void resetConnectionPaths() {
			connectionPath = new Stack<T>();
			connectionPaths = new ArrayList<Stack<T>>();
		}

		/**
		 * Returns a list of RepastEdges in the shortest path from source to target.
		 * 
		 * @param source
		 * @param target
		 * @return
		 */
		public List<RepastEdge<T>> getShortestPath(T source, T target){
			
			if (calc){
				calcShortestPaths();
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
		public double getShortestPathLength(T source, T target){
			if (calc){
				calcShortestPaths();
				calc = false;
			}
			
			Number n = dsp.getDistance(source, target);

			if (n != null)
				return n.doubleValue();
			else
				return Double.POSITIVE_INFINITY;
		}
		
		
		/*
		 * Returns a list of all paths of nodes connecting two nodes in the network. Paths containing duplicated nodes are excluded.
		 * 
		 * @param T node
		 * 		The starting node of the path
		 * @param T targetNode
		 * 		The end node of the path
		 * @param Predicate<? super T> nodeFilter
		 * 		Used to filter which nodes can be including in the paths
		 * 
		 * @return List<Stack<T>>
		 * 		The paths
		 */
		public List<Stack<T>> getSimplePaths(T node, T targetNode, Predicate<? super T> nodeFilter){
			calcSimplePaths(node, targetNode, nodeFilter);
			List<Stack<T>> output = this.connectionPaths;
			resetConnectionPaths(); // Empty the paths
			return output;
		}
		
		/*
		 * Returns a list of all paths of nodes connecting two nodes in the network. Paths containing duplicated nodes are excluded.
		 * 
		 * @param T node
		 * 		The starting node of the path
		 * @param T targetNode
		 * 		The end node of the path
		 * 
		 * @return List<Stack<T>>
		 * 		The paths
		 */
		public List<Stack<T>> getSimplePaths(T node, T targetNode){
			return getSimplePaths(node, targetNode, null);
		}
		
		/*
		 * Given a path of nodes return the corresponding path of edges that connect these nodes
		 */
		public List<RepastEdge<T>> edgePathFromNodes(Stack<T> nodePath) {
			List<RepastEdge<T>> edgePath = new ArrayList<RepastEdge<T>>();
			for (int i = 0; i<nodePath.size()-1;i++) {
				edgePath.add(net.getEdge(nodePath.get(i), nodePath.get(i+1)));
			}
			return edgePath;
		}
		
		/*
		 * Given a list of network edges return the length of this path using the class transformer
		 */
		public double getPathLength(List<RepastEdge<T>> edgePath) {
			double length = 0;
			for (RepastEdge<T> e: edgePath) {
				length += this.transformer.transform(e);
			}
			return length;
		}
		
		/**
		 * Creates shortest path info  nodes using the Jung Dijkstra algorithm
		 */
		private void calcShortestPaths(){
			Graph<T, RepastEdge<T>> graph = null;
			
			if (net instanceof JungNetwork)
				graph = ((JungNetwork<T>)net).getGraph();
			else if (net instanceof ContextJungNetwork)
				graph = ((ContextJungNetwork<T>)net).getGraph();
			
			dsp = new DijkstraShortestPath<T,RepastEdge<T>>(graph, transformer);
		}
		
		/*
		 * Calculate the paths between the node and targetNode that don't contain any duplicate nodes
		 * 
		 * @param T node
		 * 		The starting node of the path
		 * @param T targetNode
		 * 		The end node of the path
		 * @param Predicate<? super T> nodeFilter
		 * 		Used to filter which nodes can be including in the paths
		 */
		private void calcSimplePaths(T node, T targetNode, Predicate<? super T> nodeFilter) {
			List<T> adj = (List<T>)net.getAdjacent(node);
			Iterable<T> validAdj = adj.stream().filter(nodeFilter).collect(Collectors.toList());
		    for (T nextNode : validAdj) {
		       if (nextNode.equals(targetNode)) {
		           Stack<T> temp = new Stack<T>();
		           for (T node1 : connectionPath)
		               temp.add(node1);
		           connectionPaths.add(temp);
		       
		       } else if (!connectionPath.contains(nextNode)) {
		           connectionPath.push(nextNode);
		           calcSimplePaths(nextNode, targetNode, nodeFilter);
		           connectionPath.pop(); // Clears connection path on the way out
		        }
		    }
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
		public static NetworkPath<?> finished(NetworkPath<?> sp){
			sp.finalize();
			sp=null;
			return sp;
		}
		
}
