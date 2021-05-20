package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
import java.util.function.Function;

import org.apache.commons.collections15.Predicate;

import org.apache.commons.collections15.Transformer;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.YenKShortestPath;
import org.jgrapht.alg.shortestpath.YenShortestPathIterator;
import org.jgrapht.graph.AsWeightedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DefaultUndirectedGraph;
import org.jgrapht.graph.builder.GraphTypeBuilder;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.algorithms.filters.FilterUtils;
import edu.uci.ics.jung.algorithms.filters.VertexPredicateFilter;
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
public class NetworkPathFinder<T> implements ProjectionListener<T> {
	
		private Network<T> net;
		private Graph<T, RepastEdge<T>> graph;
		private boolean calc = true;
		private Transformer<RepastEdge<T>,Double> transformer;
		private DijkstraShortestPath<T,RepastEdge<T>> dsp;
	  
		private Stack<T> nodePath;
		private List<Stack<RepastEdge<T>>> edgePaths;
	  
	  /**
	   * Constructor
	   * 
	   * @param net the Network
	   */
	  
	    public NetworkPathFinder(){
	    	
	    }

		public NetworkPathFinder(Network<T> net){
			init(net);
		}
		
		private void init(Network<T> net){
			this.net = net;
			resetConnectionPaths();
			this.setDefaultTransformer();
			this.net.addProjectionListener(this);
			this.graph = this.netToGraph(net);
		}
		
		/*
		 * Initialise the connection path objects used to record the paths between nodes in the network
		 */
		public void resetConnectionPaths() {
			nodePath = new Stack<T>();
			edgePaths = new ArrayList<Stack<RepastEdge<T>>>();
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
		

		public List<Stack<RepastEdge<T>>> getKShortestPaths(T node, T targetNode, int k, Predicate<T> nodeFilter) {
			// Filter graph and convert to JgraphT graph
			this.filterGraph(nodeFilter);
			DefaultUndirectedGraph<T, RepastEdge<T>> jgt = this.getJGraphTGraph();
			
			// Use Yen's K Shortest paths algorithm to get all paths of equally shortest length between node and targetNode
			YenKShortestPath<T, RepastEdge<T>> ksp = new YenKShortestPath<T, RepastEdge<T>>(jgt);
			List<Stack<RepastEdge<T>>> paths = KShortestPaths(ksp, node, targetNode, k);
			return paths;	
		}
		

		public List<Stack<RepastEdge<T>>> getKShortestPaths(T node, T targetNode, int k, Predicate<T> nodeFilter, Transformer<RepastEdge<T>, Integer> transformer) {
			
			// Filter graph and convert to JgraphT graph
			this.filterGraph(nodeFilter);
			DefaultUndirectedGraph<T, RepastEdge<T>> jgt = this.getJGraphTGraph();
			
			// convert to weighted graph with the weights given by the transformer
			Function<RepastEdge<T>, Double> weightFunction = (RepastEdge<T> e)-> {return (double) transformer.transform(e); };
			AsWeightedGraph<T, RepastEdge<T>> wjgt = new AsWeightedGraph<T, RepastEdge<T>>(jgt, weightFunction, false, false);
			
			// Use Yen's K Shortest paths algorithm to get all paths of equally shortest length between node and targetNode
			YenKShortestPath<T, RepastEdge<T>> ksp = new YenKShortestPath<T, RepastEdge<T>>(wjgt);
			List<Stack<RepastEdge<T>>> paths = KShortestPaths(ksp, node, targetNode, k);
			return paths;		
		}
		
		public List<Stack<RepastEdge<T>>> KShortestPaths(YenKShortestPath<T, RepastEdge<T>> ksp, T node, T targetNode, int k) {
			List<GraphPath<T, RepastEdge<T>>> paths = ksp.getPaths(node, targetNode, k);
			
			List<Stack<RepastEdge<T>>> stackPaths = new ArrayList<Stack<RepastEdge<T>>>();
			for (GraphPath<T, RepastEdge<T>> p: paths) {
				Stack<RepastEdge<T>> stackPath = new Stack<RepastEdge<T>>();
				stackPath.addAll(p.getEdgeList());
				stackPaths.add(stackPath);
			}
			return stackPaths;
		}
		
		public List<Stack<RepastEdge<T>>> getAllShortestPaths(T node, Collection<T> targetNodes, Predicate<T> nodeFilter, Transformer<RepastEdge<T>, Integer> transformer) {
			// Filter graph and convert to JgraphT graph
			this.filterGraph(nodeFilter);
			DefaultUndirectedGraph<T, RepastEdge<T>> jgt = this.getJGraphTGraph();
			
			// convert to weighted graph with the weights given by the transformer
			Function<RepastEdge<T>, Double> weightFunction = (RepastEdge<T> e)-> {return (double) transformer.transform(e); };
			AsWeightedGraph<T, RepastEdge<T>> wjgt = new AsWeightedGraph<T, RepastEdge<T>>(jgt, weightFunction, false, false);
			
			// Find shortest path length to each of the target nodes. Then for all those of shortest length, find all shortest paths
			org.jgrapht.alg.shortestpath.DijkstraShortestPath<T, RepastEdge<T>> dsp = new org.jgrapht.alg.shortestpath.DijkstraShortestPath<T, RepastEdge<T>>(wjgt);
			
			List<Double> distances = new ArrayList<Double>();
			List<T> targets = new ArrayList<T>(targetNodes);
			for (T target: targets) {
				Double d = dsp.getPathWeight(node, target);
				distances.add(d);
			}
			
			double minDist = Collections.min(distances); 
			List<T> shortestDistReachable = new ArrayList<T>();
			for (int i=0; i< targets.size(); i++) {
				if (Double.compare(minDist, distances.get(i))==0) {
					shortestDistReachable.add(targets.get(i));
				}
			}
			
			// Now for each of the target nodes that can be reached in within min distance, find all shortest paths to target from source
			List<Stack<RepastEdge<T>>> output = new ArrayList<Stack<RepastEdge<T>>>();
			for (T t: shortestDistReachable) {
				YenShortestPathIterator<T, RepastEdge<T>> iterator = new YenShortestPathIterator<T, RepastEdge<T>>(wjgt, node, t);
				List<Stack<RepastEdge<T>>> paths = allshortestPathsFromYenIterator(iterator);
				for (Stack<RepastEdge<T>> p:paths) {
					output.add(p);
				}
			}
			
			return output;
		}
		
		/*
		 * Methods to get all shortest paths from source to target node. Returns multiple paths if
		 * there are multiple paths with join shortest length. uses JgraphT type graph.
		 */
		public List<Stack<RepastEdge<T>>> getAllShortestPaths(T node, T targetNode, Predicate<T> nodeFilter, Transformer<RepastEdge<T>, Integer> transformer) {
			
			// Filter graph and convert to JgraphT graph
			this.filterGraph(nodeFilter);
			DefaultUndirectedGraph<T, RepastEdge<T>> jgt = this.getJGraphTGraph();
			
			// convert to weighted graph with the weights given by the transformer
			Function<RepastEdge<T>, Double> weightFunction = (RepastEdge<T> e)-> {return (double) transformer.transform(e); };
			AsWeightedGraph<T, RepastEdge<T>> wjgt = new AsWeightedGraph<T, RepastEdge<T>>(jgt, weightFunction, false, false);
			
			// Use Yen's K Shortest paths algorithm to get all paths of equally shortest length between node and targetNode
			YenShortestPathIterator<T, RepastEdge<T>> iterator = new YenShortestPathIterator<T, RepastEdge<T>>(wjgt, node, targetNode);
			List<Stack<RepastEdge<T>>> paths = allshortestPathsFromYenIterator(iterator);
			return paths;		
		}
		
		/*
		 * Methods to get all shortest paths from source to target node. Returns multiple paths if
		 * there are multiple paths with join shortest length. uses JgraphT type graph.
		 */
		public List<Stack<RepastEdge<T>>> getAllShortestPaths(T node, T targetNode, Predicate<T> nodeFilter) {
						
			// Filter graph and convert to JgraphT graph
			this.filterGraph(nodeFilter);
			DefaultUndirectedGraph<T, RepastEdge<T>> jgt = this.getJGraphTGraph();
			
			// Use Yen's K Shortest paths algorithm to get all paths of equally shortest length between node and targetNode
			YenShortestPathIterator<T, RepastEdge<T>> iterator = new YenShortestPathIterator<T, RepastEdge<T>>(jgt, node, targetNode);
			List<Stack<RepastEdge<T>>> paths = allshortestPathsFromYenIterator(iterator);
			return paths;
		}
		
		private List<Stack<RepastEdge<T>>> allshortestPathsFromYenIterator(YenShortestPathIterator<T, RepastEdge<T>> iterator) {
			List<Stack<RepastEdge<T>>> output = new ArrayList<Stack<RepastEdge<T>>>();
			GraphPath<T, RepastEdge<T>> path = iterator.next();
			Stack<RepastEdge<T>> stackPath = new Stack<RepastEdge<T>>();
			stackPath.addAll(path.getEdgeList());
			output.add(stackPath);
			
			Double shortestPathLength = path.getWeight();
			
			boolean notFinished = iterator.hasNext();
			while (notFinished) {
				path = iterator.next();
				Double pathLength = path.getWeight();
				
				// Check if this path length is equal to or less than the first path length (which should be the shortest)
				if( Double.compare(shortestPathLength, pathLength)>=0 ) {
					stackPath = new Stack<RepastEdge<T>>();
					stackPath.addAll(path.getEdgeList());
					notFinished = iterator.hasNext();
				}
				else {
					notFinished = false;
				}
			}
			
			return output;
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
		public List<Stack<RepastEdge<T>>> getSimplePaths(T node, T targetNode, Predicate<T> nodeFilter){
			this.filterGraph(nodeFilter);
			List<Stack<RepastEdge<T>>> paths = getSimplePaths(node, targetNode); 
			return paths;
		}
		
		/*
		 * Returned a list of all paths of nodes connecting source node to any one of the target nodes. paths containing duplicated nodes are excluded.
		 * 
		 * @param T node
		 * 		The starting node of the path
		 * @param Collection<T> targetNodes
		 * 		Paths must terminate at one of these nodes
		 * @param Predicate<? super T> nodeFilter
		 * 		Used to filter which nodes can be including in the paths
		 * 
		 * @return List<Stack<T>>
		 * 		The paths
		 */
		public List<Stack<RepastEdge<T>>> getSimplePaths(T node, Collection<T> targetNodes, Predicate<T> nodeFilter){
			this.filterGraph(nodeFilter);
			List<Stack<RepastEdge<T>>> paths = new ArrayList<Stack<RepastEdge<T>>>();
			for (T tN: targetNodes) {
				getSimplePaths(node, tN).stream().forEach(paths::add);
			}
			return paths;
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
		public List<Stack<RepastEdge<T>>> getSimplePaths(T node, T targetNode){
			calcSimplePaths(node, targetNode);
			List<Stack<RepastEdge<T>>> output = this.edgePaths;
			resetConnectionPaths(); // Empty the paths
			return output;
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
		 * Extract the series of nodes that make up a path of edges
		 */
		public List<T> nodePathFromEdges(List<RepastEdge<T>> edgePath, T startNode) {
			List<T> pathNodes = new LinkedList<T>();
			T prev = startNode;
			pathNodes.add(prev);
			for (RepastEdge<T> e: edgePath) {
				if (e.getSource().equals(prev)) {
					prev = e.getTarget();
				}
				else {
					prev = e.getSource();
				}
				pathNodes.add(prev);
			}
			return pathNodes;
		}
		
		/*
		 * Given a list of network edges return the length of this path using the input transformer
		 */
		public static <T> double getPathLength(List<RepastEdge<T>> edgePath, Transformer<RepastEdge<T>,Double> t) {
			double length = 0;
			for (RepastEdge<T> e: edgePath) {
				length += t.transform(e);
			}
			return length;
		}
		
		/*
		 * Given a list of network edges return the length of this path using the input transformer
		 */
		public static <T> int getIntPathLength(List<RepastEdge<T>> edgePath, Transformer<RepastEdge<T>,Integer> t) {
			int length = 0;
			for (RepastEdge<T> e: edgePath) {
				length += t.transform(e);
			}
			return length;
		}
		
		/*
		 * Given a list of network edges return the length of this path using the class transformer
		 */
		public double getPathLength(List<RepastEdge<T>> edgePath) {
			return getPathLength(edgePath, this.transformer);
		}
		
		/*
		 * Given a list of paths return the shortest. Can include multiple paths if they have equal shortest lengths.
		 */
		public List<Stack<RepastEdge<T>>> getShortestOfMultiplePaths(List<Stack<RepastEdge<T>>> paths, Transformer<RepastEdge<T>, Integer> t) {
			int minLength = Integer.MAX_VALUE;
			List<Stack<RepastEdge<T>>> output = new ArrayList<Stack<RepastEdge<T>>>();
			for (Stack<RepastEdge<T>> path: paths) {
				int pathLength = getIntPathLength(path, t);
				if (pathLength < minLength) {
					output.clear();
					minLength = pathLength;
					output.add(path);
				}
				else if (pathLength == minLength) {
					output.add(path);
				}
			}
			return output;
		}
		
		/**
		 * Creates shortest path info  nodes using the Jung Dijkstra algorithm
		 */
		private void calcShortestPaths(){
			dsp = new DijkstraShortestPath<T,RepastEdge<T>>(this.graph, transformer);
		}
		
		/*
		 * Calculate the paths between the node and targetNode that don't contain any duplicate nodes.
		 * 
		 * Importantly, if the start node and targetNode are the same then this function only returns
		 * a single path containing just the start node
		 * 
		 * @param T node
		 * 		The starting node of the path
		 * @param T targetNode
		 * 		The end node of the path
		 * @param Predicate<? super T> nodeFilter
		 * 		Used to filter which nodes can be including in the paths
		 */
		private void calcSimplePaths(T node, T targetNode) {
	        // First add the node to the path
			nodePath.push(node);
			
			// If target node reached log path and exit
			// Could create edge path here...?
			if (node.equals(targetNode)) {
	           Stack<RepastEdge<T>> temp = new Stack<RepastEdge<T>>();
	           for (int i = 0; i<nodePath.size()-1;i++) {
	        	   T ni = nodePath.get(i);
	        	   T nj = nodePath.get(i+1);
	        	   RepastEdge<T> e = net.getEdge(ni, nj);
	               temp.add(e);
	           }
	           edgePaths.add(temp);
			}
			
			// Otherwise get children nodes and add them to the path
			else {
				// Then find valid children nodes
				Collection<T> adj = graph.getNeighbors(node);
				
				// If node if not contained in the graph (this happens when filter excludes starting node) adj is null
				if (adj != null) {
				    for (T nextNode : adj) {
				    	if (!nodePath.contains(nextNode)) {
				           calcSimplePaths(nextNode, targetNode);
				           nodePath.pop(); // Clears connection path on the way out
					    }
					}
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
		
		public Network<T> getNet() {
			return net;
		}
		
		public Graph<T, RepastEdge<T>> getGraph() {
			return this.graph;
		}
		
		public DefaultUndirectedGraph<T, RepastEdge<T>> getJGraphTGraph() {
			return jgrapgtFromJung(this.graph);
		}
		
		private Graph<T, RepastEdge<T>> netToGraph(Network<T> net) {
			Graph<T, RepastEdge<T>> graph = null;

			if (net instanceof JungNetwork)
				graph = ((JungNetwork<T>)net).getGraph();
			else if (net instanceof ContextJungNetwork)
				graph = ((ContextJungNetwork<T>)net).getGraph();
			
			return graph;
		}
		
	    private DefaultUndirectedGraph<T, RepastEdge<T>> jgrapgtFromJung(Graph<T, RepastEdge<T>> g) {
	    	
	    	DefaultUndirectedGraph<T, RepastEdge<T>> jgtGraph = new DefaultUndirectedGraph<T, RepastEdge<T>>((Class<? extends RepastEdge<T>>) RepastEdge.class);
	    	
	    	// Loop through edges in input graph and add these to the new graph
	    	for (RepastEdge<T> e : g.getEdges()) {
	    		T u = e.getSource();
	    		T v = e.getTarget();
	    		
	    		jgtGraph.addVertex(u);
	    		jgtGraph.addVertex(v);
	    		jgtGraph.addEdge(u, v, e);
	    	}
	    	
	    	return jgtGraph;
	    }
		
		public void filterGraph(Collection<T> vertices) {
			Graph<T, RepastEdge<T>> g = this.netToGraph(net);
			this.graph = FilterUtils.createInducedSubgraph(vertices, g);
			calc = true;
		}
		
		public void filterGraph(Predicate<T> verticesFilter) {
			Graph<T, RepastEdge<T>> g = this.netToGraph(net);
			VertexPredicateFilter<T,RepastEdge<T>> filter = new VertexPredicateFilter<T,RepastEdge<T>>(verticesFilter);
			this.graph = filter.transform(g);
			calc = true;
		}
		
		public void resetgraph() {
			this.graph = this.netToGraph(net);
			calc = true;
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
		public static NetworkPathFinder<?> finished(NetworkPathFinder<?> sp){
			sp.finalize();
			sp=null;
			return sp;
		}		
}
