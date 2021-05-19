package graph;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Hypergraph;

public class DijkstraMultiDistance<V, E> extends DijkstraDistance<V, E> {

	public DijkstraMultiDistance(Graph<V, E> g) {
		super(g);
		// TODO Auto-generated constructor stub
	}

	public DijkstraMultiDistance(Graph<V, E> g, boolean cached) {
		super(g, cached);
		// TODO Auto-generated constructor stub
	}

	public DijkstraMultiDistance(Hypergraph<V, E> g, Transformer<E, ? extends Number> nev, boolean cached) {
		super(g, nev, cached);
		// TODO Auto-generated constructor stub
	}

	public DijkstraMultiDistance(Hypergraph<V, E> g, Transformer<E, ? extends Number> nev) {
		super(g, nev);
		// TODO Auto-generated constructor stub
	}
	
	

}
