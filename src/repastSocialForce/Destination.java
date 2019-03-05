package repastSocialForce;

import java.util.Random;

import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;

public class Destination {
	private int ext;
	private ContinuousSpace<Object> space;
		
	public Destination(ContinuousSpace<Object> space, int extent) {
		this.space = space;
		this.ext = extent;
	}
	
	public int getExtent() {
		return this.ext;
	}

}
