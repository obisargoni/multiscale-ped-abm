package repastInterSim.pathfinding;

import java.util.ArrayList;
import java.util.List;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.space.graph.RepastEdge;
import repastInterSim.environment.Junction;

public class TacticalRoute {
	
	private Coordinate c;
	private Double costT;
	private Integer parityT;

	private Double costS;
	private Integer parityS;
	
	private List<Junction> routeJunctions;
	private List<RepastEdge<Junction>> routePath;
	private List<RepastEdge<Junction>> remainderPath;
	
	public TacticalRoute(Coordinate c){
		this.c = c;
	}
	
	public TacticalRoute() {
		this.routeJunctions = new ArrayList<Junction>();
		this.routePath = new ArrayList<RepastEdge<Junction>>();
	}
	
	public void addJunction(Junction j) {
		this.routeJunctions.add(j);
	}
	
	public List<Junction> getRouteJunctions(){
		return this.routeJunctions;
	}
	
	public List<RepastEdge<Junction>> getRoutePath() {
		return routePath;
	}

	public void setRoutePath(List<RepastEdge<Junction>> routePath) {
		this.routePath = routePath;
	}

	public Integer getParityT() {
		return parityT;
	}

	public void setParityT(Integer parityT) {
		this.parityT = parityT;
	}

	public Coordinate getC() {
		return c;
	}

	public void setC(Coordinate c) {
		this.c = c;
	}

	public Double getCostT() {
		return costT;
	}

	public void setCostT(Double costT) {
		this.costT = costT;
	}

	public Double getCostS() {
		return costS;
	}

	public void setCostS(Double costS) {
		this.costS = costS;
	}

	public Integer getParityS() {
		return parityS;
	}

	public void setParityS(Integer parityS) {
		this.parityS = parityS;
	}

	public void setRouteRemainderPath(List<RepastEdge<Junction>> path) {
		this.remainderPath = path;
	}
	
	public List<RepastEdge<Junction>> getRouteRemainderPath() {
		return this.remainderPath;
	}
	
	
}