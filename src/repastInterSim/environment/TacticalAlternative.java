package repastInterSim.environment;

import com.vividsolutions.jts.geom.Coordinate;

public class TacticalAlternative {
	
	private Coordinate c;
	private Double costT;
	private Integer parityT;

	private Double costS;
	private Integer parityS;
	
	public TacticalAlternative(Coordinate c){
		this.c = c;
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
	
	
}