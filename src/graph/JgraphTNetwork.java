package graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.stream.Collectors;

import org.jgrapht.Graph;
import org.jgrapht.Graphs;

import repast.simphony.random.RandomHelper;
import repast.simphony.space.graph.DefaultEdgeCreator;
import repast.simphony.space.graph.EdgeCreator;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repast.simphony.space.projection.DefaultProjection;
import repast.simphony.space.projection.ProjectionEvent;
import repast.simphony.space.projection.ProjectionPredicate;

public abstract class JgraphTNetwork<T> extends DefaultProjection<T> implements Network<T> {

	protected Graph<T, RepastEdge<T>> graph;
	private ArrayList<T> tmpRandomList = new ArrayList<T>();
	protected EdgeCreator<? extends RepastEdge<T>, T> creator;

	public JgraphTNetwork(String name) {
		this(name, new DefaultEdgeCreator<T>());
	}

	public JgraphTNetwork(String name, EdgeCreator<? extends RepastEdge<T>, T> creator) {
		super(name);
		this.creator = creator;
	}

	/**
	 * Gets the EdgeCreator used to create edges for
	 * this Network. {@link #addEdge(Object, Object) addEdge} and
	 * {@link #addEdge(Object, Object, double) addEdge} will use
	 * this creator to create edges. Any edge added with
	 * {@link #addEdge(repast.simphony.space.graph.RepastEdge) addEdge} must be of the same
	 * type as that created with this EdgeCreator. By default,
	 * an edge creator that creates RepastEdges is used.
	 * <p/>
	 * The default EdgeCreator will create
	 * RepastEdge
	 *
	 * @return the edge class of this network
	 */
	public EdgeCreator<? extends RepastEdge<T>, T> getEdgeCreator() {
		return creator;
	}

	public void setGraph(Graph<T, RepastEdge<T>> graph) {
		this.graph = graph;
	}

	public Graph<T, RepastEdge<T>> getGraph() {
		return graph;
	}

	public boolean evaluate(ProjectionPredicate predicate) {
		return predicate.evaluate(this);
	}

	public Iterable<T> getAdjacent(T agent) {
		Collection<T> out = new LinkedHashSet<T>();
		out.addAll(Graphs.predecessorListOf(this.graph, agent));
		out.addAll(Graphs.successorListOf(this.graph, agent));
		return out;
	}

	public int getDegree() {
		return graph.edgeSet().size();
	}

	public int getDegree(T agent) {
		return graph.degreeOf(agent);
	}

	public Iterable<RepastEdge<T>> getEdges() {
		return graph.edgeSet();
	}

	public Iterable<RepastEdge<T>> getEdges(T agent) {
		return graph.edgesOf(agent);
	}

	public int getInDegree(T agent) {
		return graph.inDegreeOf(agent);
	}

	public Iterable<RepastEdge<T>> getInEdges(T agent) {
		return graph.incomingEdgesOf(agent);
	}

	public Iterable<T> getNodes() {
		return graph.vertexSet();
	}

	public int getOutDegree(T agent) {
		return graph.outDegreeOf(agent);
	}

	public Iterable<RepastEdge<T>> getOutEdges(T agent) {
		return graph.outgoingEdgesOf(agent);
	}

	public Iterable<T> getPredecessors(T agent) {
		return Graphs.predecessorListOf(graph, agent);
	}

	public T getRandomAdjacent(T agent) {
		Set<T> neighbours = Graphs.neighborSetOf(graph, agent);
		int size = neighbours.size();
		if (size == 0) {
			return null;
		}
		int index = RandomHelper.getUniform().nextIntFromTo(0, size - 1);
		tmpRandomList.clear();
		tmpRandomList.addAll(neighbours);
		return tmpRandomList.get(index);
	}

	public T getRandomPredecessor(T agent) {
		Set<T> preds =  Graphs.predecessorListOf(graph, agent).stream().collect(Collectors.toSet());
		int size = preds.size();
		if (size == 0) {
			return null;
		}
		int index = RandomHelper.getUniform().nextIntFromTo(0, size - 1);
		tmpRandomList.clear();
		tmpRandomList.addAll(preds);
		return tmpRandomList.get(index);
	}

	public T getRandomSuccessor(T agent) {
		Set<T> sucs =  Graphs.successorListOf(graph, agent).stream().collect(Collectors.toSet());
		int size = sucs.size();
		if (size == 0) {
			return null;
		}
		int index = RandomHelper.getUniform().nextIntFromTo(0, size - 1);
		tmpRandomList.clear();
		tmpRandomList.addAll(sucs);
		return tmpRandomList.get(index);
	}

	public Iterable<T> getSuccessors(T agent) {
		return Graphs.successorListOf(graph,  agent);
	}

	public boolean isAdjacent(T first, T second) {
		return graph.containsEdge(first, second) ||
		graph.containsEdge(second, first);
	}

	public boolean isPredecessor(T first, T second) {
		return graph.containsEdge(second, first);
	}

	public boolean isSuccessor(T first, T second) {
		return graph.containsEdge(first, second);
	}

	public int numEdges() {
		return graph.edgeSet().size();
	}

	public void removeEdge(RepastEdge<T> edge) {
	        if (graph.removeEdge(edge)) {
	          fireProjectionEvent(new ProjectionEvent<T>(this, edge,
    				ProjectionEvent.EDGE_REMOVED));
	        }
	}

  /**
   * Returns whether or not this network contains the specified edge.
   *
   * @param edge the edge to check
   * @return true if the network contains the specified edge, otherwise false.
   */
  public boolean containsEdge(RepastEdge<T> edge) {
    return graph.containsEdge(edge);
  }

  public RepastEdge<T> addEdge(RepastEdge<T> edge) {
    RepastEdge<T> oldEdge = graph.getEdge(edge.getSource(), edge.getTarget());
    if (oldEdge != null) removeEdge(oldEdge);    
    graph.addEdge(edge.getSource(), edge.getTarget(), edge);
	fireProjectionEvent(new ProjectionEvent<T>(this, edge, ProjectionEvent.EDGE_ADDED));
	return edge;
  }
  

	public RepastEdge<T> addEdge(T source, T target) {
		return addEdge(source, target, 1);
	}

	public RepastEdge<T> getEdge(T source, T target) {
		return graph.getEdge(source, target);
	}

	public int size() {
		return graph.vertexSet().size();
	}

	public void addVertex(T vertex) {
		graph.addVertex(vertex);
		fireProjectionEvent(new ProjectionEvent<T>(this, vertex,
				ProjectionEvent.OBJECT_ADDED));
	}

	public void removeVertex(T vertex) {
    Collection<RepastEdge<T>> edges = graph.edgesOf(vertex);
    
    ArrayList<RepastEdge<T>> tempEdges = new ArrayList<RepastEdge<T>>();
    for (RepastEdge<T> edge : edges)
    	tempEdges.add(edge);
      
		graph.removeVertex(vertex);
		
		fireProjectionEvent(new ProjectionEvent<T>(this, vertex,ProjectionEvent.OBJECT_REMOVED));
		
		for (RepastEdge<T> edge : tempEdges) 	
			fireProjectionEvent(new ProjectionEvent<T>(this, edge,ProjectionEvent.EDGE_REMOVED));
  }
	
	public void removeEdges() {
		Iterator<RepastEdge<T>> ie = this.getEdges().iterator();
		 while(ie.hasNext()){
			 this.removeEdge(ie.next());
			 ie=this.getEdges().iterator();
		 }
	}
}