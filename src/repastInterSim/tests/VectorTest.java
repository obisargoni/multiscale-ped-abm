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
		
		v2[0] = -1;
		v2[1] = 1;
		
		ang = Vector.angleBetweenTwoVectorsDegree(v1, v2);
		assert Math.round(ang) == 135.0;
		
		v2[0] = 2;
		v2[1] = -2;
		
		ang = Vector.angleBetweenTwoVectorsDegree(v1, v2);
		assert Math.round(ang) == 45.0;
		
		
	}

}
