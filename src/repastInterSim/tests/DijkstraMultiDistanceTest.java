package repastInterSim.tests;

import java.util.ArrayList;
import java.util.Map;

import org.junit.jupiter.api.Test;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Point;

import edu.uci.ics.jung.graph.Graph;
import graph.DijkstraMultiDistance;
import repast.simphony.context.Context;
import repast.simphony.context.space.graph.ContextJungNetwork;
import repast.simphony.context.space.graph.NetworkBuilder;
import repast.simphony.space.graph.Network;
import repast.simphony.space.graph.RepastEdge;
import repastInterSim.environment.GISFunctions;
import repastInterSim.environment.Junction;
import repastInterSim.environment.NetworkEdgeCreator;
import repastInterSim.environment.contexts.JunctionContext;

class DijkstraMultiDistanceTest {
	
	Graph<Junction, RepastEdge<Junction>> graph;
	Junction[] nodes;
	
	void initToyGraph() {
		Network<Junction> net;

		Context<Junction> junctionContext = new JunctionContext();
		NetworkBuilder<Junction> builder = new NetworkBuilder<Junction>("net", junctionContext, false);
		builder.setEdgeCreator(new NetworkEdgeCreator<Junction>());
		net = builder.buildNetwork();
		
		// Create nodes of toy graph
		int n=4;
		nodes = new Junction[n];
		for (int i=0;i<n;i++) {
			Coordinate c  = new Coordinate(0,i);
			Point p = GISFunctions.pointGeometryFromCoordinate(c);
			Junction j = new Junction(String.valueOf(i));
			j.setGeom(p);
			junctionContext.add(j);
			nodes[i] = j;
		}
		
		// Add edges to toy graph
		RepastEdge<Junction> e1 = new RepastEdge<Junction>(nodes[0], nodes[1], false, 1);
		RepastEdge<Junction> e2 = new RepastEdge<Junction>(nodes[0], nodes[2], false, 1);
		RepastEdge<Junction> e3 = new RepastEdge<Junction>(nodes[1], nodes[3], false, 1);
		RepastEdge<Junction> e4 = new RepastEdge<Junction>(nodes[2], nodes[3], false, 1);
		
		net.addEdge(e1);
		net.addEdge(e2);
		net.addEdge(e3);
		net.addEdge(e4);
		
		graph = ((ContextJungNetwork<Junction>)net).getGraph();
	}
	
	@Test
	void testShortestMultiPath() {
		initToyGraph();
		
		DijkstraMultiDistance<Junction, RepastEdge<Junction>> dmd = new DijkstraMultiDistance<Junction, RepastEdge<Junction>>(graph);
		
		Junction source = nodes[0];
		ArrayList<Junction> targets = new ArrayList<Junction>();
		targets.add(nodes[3]);
		Map<Junction, Number> dists = dmd.getDistanceMap(source, targets);
		assert dists.entrySet().size()==2;
	}

}
