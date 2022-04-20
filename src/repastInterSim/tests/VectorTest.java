package repastInterSim.tests;

import org.junit.jupiter.api.Test;

import repastInterSim.environment.Vector;

class VectorTest {

	@Test
	void testAngleBetweenTwoVectorsDegree() {
		double[] v1 = {5,0};
		double[] v2 = {1,1};
		
		Double ang = Vector.angleBetweenTwoVectorsDegree(v1, v2);
		assert Math.round(ang) == 45.0;
		
		v2[0] = -1;
		v2[1] = -1;
		
		ang = Vector.angleBetweenTwoVectorsDegree(v1, v2);
		assert Math.round(ang) == 135.0;
		
		v1[0] = -5;
		
		ang = Vector.angleBetweenTwoVectorsDegree(v1, v2);
		assert Math.round(ang) == 45.0;
		
		v1[0] = 5;
		v2[0] = -1;
		v2[1] = 1;
		
		ang = Vector.angleBetweenTwoVectorsDegree(v1, v2);
		assert Math.round(ang) == 135.0;
		
		v2[0] = 2;
		v2[1] = -2;
		
		ang = Vector.angleBetweenTwoVectorsDegree(v1, v2);
		assert Math.round(ang) == 45.0;
		
		v2[0] = 0.1;
		v2[1] = 1;
		
		ang = Vector.angleBetweenTwoVectorsDegree(v1, v2);
		assert Math.round(ang) < 90.0;
		
		v2[0] = -0.1;
		
		ang = Vector.angleBetweenTwoVectorsDegree(v1, v2);
		assert Math.round(ang) > 90.0;
	}
	
	@Test
	void testEdgeTG() {
		// Setup edges and speeds
		double[] e10 = {4,0};
		double[] e11 = {0,0};
		double[] e20 = {4,2};
		double[] e21 = {0,2};
		double[] vS = {1,0};
		
		double[] pLoc = {6,0};
		double[] pV = {0,0};
		
		Double tg = Vector.edgeTG(e10, e11, e20, e21, pLoc, pV, vS);
		assert tg ==-2.0; // -ve since ped is not moving and therefore is second to arrive
		
		// Change ped velocity such that ped would not collide with vehicle
		pV[1] = 3;
		tg = Vector.edgeTG(e10, e11, e20, e21, pLoc, pV, vS);
		assert tg==2.0;
		
		// Again, but make ped move slower from slightly further back. Expect ped to arrive second this time
		pLoc[1] = -1;
		pV[1]=0.1;
		tg = Vector.edgeTG(e10, e11, e20, e21, pLoc, pV, vS);
		assert tg==-10.0;
		
		// Moving in different direction ped now classed as arriving first
		// Vehicle arrives second so time given relates to time that vehicle arrives at collision location
		pV[1]=-0.1;
		tg = Vector.edgeTG(e10, e11, e20, e21, pLoc, pV, vS);
		assert tg==2;
		
		
		// Test movement at angles
		double a = Math.PI/4.0;
		double[] e10_ = {0,0};
		double[] e11_ = {-4*Math.sin(a), 4*Math.cos(a)};
		double[] e20_ = {2*Math.sin(a), 2*Math.cos(a)};
		double[] e21_ = {-2*Math.sin(a), 6*Math.cos(a)};
		double[] vS_ = {6*Math.sin(a), -6*Math.cos(a)};
		
		pLoc[0]=0;
		pLoc[1]=-6;
		pV[0]=1;
		pV[1]=0;
		
		tg = Vector.edgeTG(e10_, e11_, e20_, e21_, pLoc, pV, vS_);
		assert tg == -6*Math.tan(a); // time it takes ped to reach conflict point
		
		// speed up ped so it arrives first
		pV[0] = 5;
		tg = Vector.edgeTG(e10_, e11_, e20_, e21_, pLoc, pV, vS_);
		assert tg == 6/(6*Math.cos(a)); // time is now time it takes vehicle to reach conflict point
		
	}

}
