package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;

public class SimplePath<T> {
	
	Network<T> net;
	private Stack<T> connectionPath;
	private List<Stack<T>> connectionPaths;

	public SimplePath(Network<T> n) {
		net = n;
		init();
	}
	
	public void init() {
		connectionPath = new Stack<T>();
		connectionPaths = new ArrayList<Stack<T>>();
	}
	
	// Push to connectionsPath the object that would be passed as the parameter 'node' into the method below
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
	
	public List<Stack<T>> getSimplePaths(T node, T targetNode, Predicate<? super T> nodeFilter){
		calcSimplePaths(node, targetNode, nodeFilter);
		List<Stack<T>> output = this.connectionPaths;
		init(); // Empty the paths
		return output;
	}
	
	public List<Stack<T>> getSimplePaths(T node, T targetNode){
		return getSimplePaths(node, targetNode, null);
	}

	public List<RepastEdge<T>> edgePathFromNodes(Stack<T> nodePath) {
		List<RepastEdge<T>> edgePath = new ArrayList<RepastEdge<T>>();
		for (int i = 0; i<nodePath.size()-1;i++) {
			edgePath.add(net.getEdge(nodePath.get(i), nodePath.get(i+1)));
		}
		return edgePath;
	}

	public double getPathLength(List<RepastEdge<T>> edgePath) {
		double length = 0;
		for (RepastEdge<T> e: edgePath) {
			length += e.getWeight();
		}
		return length;
	}
}
