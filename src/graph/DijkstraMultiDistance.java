package graph;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.collections15.Transformer;

import com.jgoodies.common.base.Preconditions;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.util.MapBinaryHeap;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Hypergraph;

public class DijkstraMultiDistance<N, E> extends DijkstraDistance<N, E> {
	
	protected Map<N, SourceData> sourceMap; // a map of source nodes to an instance of SourceData
	private double maxDistance;
	private int maxTargets;
	
	public DijkstraMultiDistance(Hypergraph<N, E> g, Transformer<E, ? extends Number> nev, boolean cached) {
		super(g, nev, cached);
		// TODO Auto-generated constructor stub
	    this.sourceMap = new HashMap<>();
	    this.maxDistance = Double.POSITIVE_INFINITY;
	    this.maxTargets = Integer.MAX_VALUE;
	}


	public DijkstraMultiDistance(Graph<N, E> g, boolean cached) {
		this(g, e -> 1, cached);
	}

	public DijkstraMultiDistance(Hypergraph<N, E> g, Transformer<E, ? extends Number> nev) {
		this(g, nev, true);
	}
	
	public DijkstraMultiDistance(Graph<N, E> g) {
		this(g, e -> 1, true);
	}
	
	
 /**
   * Implements Dijkstra's single-source shortest-path algorithm for weighted graphs. Uses a <code>
   * MapBinaryHeap</code> as the priority queue, which gives this algorithm a time complexity of O(m
   * lg n) (m = # of edges, n = # of nodes). This algorithm will terminate when any of the following
   * have occurred (in order of priority):
   *
   * <ul>
   *   <li>the distance to the specified target (if any) has been found
   *   <li>no more nodes are reachable
   *   <li>the specified # of distances have been found, or the maximum distance desired has been
   *       exceeded
   *   <li>all distances have been found
   * </ul>
   *
   * @param source the node from which distances are to be measured
   * @param numDistances the number of distances to measure
   * @param targets the set of nodes to which distances are to be measured
   * @return a mapping from node to the shortest distance from the source to each target
   */
	@Override
	protected LinkedHashMap<N, Number> singleSourceShortestPath(N source, Collection<N> targets, int numDistances) {
		SourceData sd = getSourceData(source);

		Set<N> toGet = new HashSet<>();
		if (targets != null) {
			toGet.addAll(targets);
			Set<N> existingDistances = sd.distances().keySet();
			for (N o : targets) {
				if (existingDistances.contains(o)) {
					toGet.remove(o);
				}
			}
		}

		// if we've exceeded the max distance or max # of distances we're willing to calculate, or
		// if we already have all the distances we need,
		// terminate
		if (sd.reachedMax()
		    || (targets != null && toGet.isEmpty())
		    || (sd.distances().size() >= numDistances)) {
		  return sd.distances();
		}
		
		while (!sd.unknownNodes().isEmpty() && (sd.distances().size() < numDistances || !toGet.isEmpty())) {
		  Map.Entry<N, Number> p = sd.getNextVertex();
		  N v = p.getKey();
		  double vDist = p.getValue().doubleValue();
		  toGet.remove(v);
		  if (vDist > this.maxDistance) {
		    // we're done; put this node back in so that we're not including
			// a distance beyond what we specified
		    sd.restoreVertex(v, vDist);
		    sd.setReachedMax(true);
		    break;
		  }
		  sd.setDistanceReached(vDist);
		
		  if (sd.distances().size() >= this.maxTargets) {
		    sd.setReachedMax(true);
		    break;
		  }
		
		  for (N w : g.getSuccessors(v)) {
		    for (E e : g.findEdgeSet(v, w)) {
		      if (!sd.distances().containsKey(w)) {
		        double edgeWeight = nev.transform(e).doubleValue();
		        Preconditions.checkArgument(
		            edgeWeight >= 0,
		            "encountered negative edge weight %s for edge %s",
		    nev.transform(e),
		    e);
		        double newDist = vDist + edgeWeight;
		        if (!sd.estimatedDistances().containsKey(w)) {
		        	sd.createRecord(w, e, newDist);
		        } else {
		        	double wDist = (Double) sd.estimatedDistances().get(w);
		        	if (newDist < wDist) { // update tentative distance & path for w
		        		sd.update(w, e, newDist);
		              }
		            }
		          }
		        }
		      }
		    }
		    return sd.distances();
		  }
	@Override
	protected SourceData getSourceData(N source) {
		SourceData sd = sourceMap.get(source);
		if (sd == null) {
			sd = new SourceData(source);
		}
		return sd;
	}
	
	
	protected class SourceData extends DijkstraDistance<N, E>.SourceData {
	    

		protected SourceData(N source) {
			super(source);
		}

		@Override
		protected void createRecord(N w, E e, double new_dist) {
			// TODO Auto-generated method stub
			super.createRecord(w, e, new_dist);
		}

		@Override
		protected Entry<N, Number> getNextVertex() {
			// TODO Auto-generated method stub
			return super.getNextVertex();
		}

		@Override
		protected void restoreVertex(N v, double dist) {
			// TODO Auto-generated method stub
			super.restoreVertex(v, dist);
		}

		@Override
		protected void update(N dest, E tentative_edge, double new_dist) {
			// TODO Auto-generated method stub
			super.update(dest, tentative_edge, new_dist);
		}
		
		
		protected void setDistances(LinkedHashMap<N, Number> d) {
			super.distances = d;
		}
		
		protected void setEstimatedDistances(Map<N, Number> ed) {
			super.estimatedDistances = ed;
		}
		
		
		protected void setUnknownNodes(MapBinaryHeap<N> un) {
			super.unknownVertices = un;
		}
		
		protected void setReachedMax(boolean rm) {
			super.reached_max = rm;
		}
		
		protected void setDistanceReached(double dr) {
			super.dist_reached = dr;
		}
		
		protected LinkedHashMap<N, Number> distances() {
			return super.distances;
		}
		
		protected Map<N, Number> estimatedDistances() {
			return super.estimatedDistances;
		}
		
		
		protected MapBinaryHeap<N> unknownNodes() {
			return super.unknownVertices;
		}
		
		protected boolean reachedMax() {
			return super.reached_max;
		}
		
		protected double distanceReached() {
			return super.dist_reached;
		}
		
	}	
}
