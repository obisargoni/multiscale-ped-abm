package repastInterSim.environment;

import com.vividsolutions.jts.geom.Coordinate;

public class TacticalAlternative {
	
	public Coordinate c;
	public Double costT;
	public Integer parityT;

	public Double costS;
	public Integer parityS;
	
	public TacticalAlternative(Coordinate c){
		this.c = c;
	}
	
	public Integer getParityT() {
		return parityT;
	}

	public void setParityT(Integer parityT) {
		this.parityT = parityT;
	}
}