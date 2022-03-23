package repastInterSim.environment;

import org.apache.commons.math3.util.FastMath;

import com.vividsolutions.jts.geom.Coordinate;

// Class for sorting functions that do standard vector operations
public class Vector {
	
	/*
	 * Compute the magnitude of a vector
	 */
	public static double mag(double[] V) {
		double magV = 0;
		for (int i = 0; i<V.length;i++) {
			magV += FastMath.pow(V[i], 2);
		}
		magV = FastMath.sqrt(magV);
		return magV;		
	}
	
	/*
	 * Compute the unit vector for an input vector
	 */
	public static double[] unitV (double[] V) {
		double magV = mag(V);
		
		double[] uV = new double[V.length];
		for (int i = 0; i<V.length;i++) {
			uV[i] = V[i] / magV;
		}
		
		return uV;
	}
	
	/*
	 * Calculate the dot product between two vectors
	 */
	public static double dotProd(double[] U, double[] V) {
		double dP = 0;
		for (int i = 0;i<FastMath.min(U.length, V.length);i++) {
			dP += U[i] * V[i];
		}
		return dP;
	}
	
    /*
     * Generic function to sum two vectors that are the same size
     * 
     * @param a vector of doubles
     * @param b vector of doubles
     * 
     * @return vector of doubles that is the sum of the two input vectors
     */
    public static double[] sumV(double[] a, double[] b) {
    	
    	if (a.length != b.length) {
    		return null;
    	}
    	else {
    		
    	}
        double[] c = new double[a.length];
            for (int i = 0; i < a.length; i++) {
            c[i] = a[i] + b[i];}
        return c;
    }
    
    /*
     * Rotate vector by the given angle
     * 
     * @param double[] u. 2D vector to rotate.
     * @param double theta. The angle in radians to rotate by 
     */
    public static double[] rotate2D(double[] u, double theta) {
    	double[] v = {Math.cos(theta)*u[0] - Math.sin(theta)*u[1],Math.sin(theta)*u[0] + Math.cos(theta)*u[1]};
    	return v;    	
    }
    
    public static double angleBetweenNorthAndUnitVector(double[] u) {
    	double[] n = {0,1};
    	
    	double signX = Math.signum(u[0]);
        double signY = Math.signum(u[1]);
        
        double dp = Vector.dotProd(n, u);
        double a = Math.acos(dp);
        
        // Dot product gives angle between vectors. We need angle clockwise from north
        if (signX == 1) {
        	return a;
        }
        else if (signX == -1) {
        	return 2*Math.PI - a;
        }
        else {
        	return a;
        }
    }
    
    public static double[] difference(double[] u, double[] v) {
    	
    	// Initialise the difference vector
    	double[] diff = new double[Math.min(u.length, v.length)];
    	
    	for (int i = 0 ; i < Math.min(u.length, v.length); i++) {
    		diff[i] = u[i] - v[i];
    	}
    	
    	return diff;
    }
    
    public static double distanceBetweenVectors(double[] u, double[] v) {
    	double[] diff = difference(u,v);
    	
    	double dist = mag(diff);
    	
    	return dist;
    }

	public static double angleBetweenTwoVectors(double[] v1, double[] v2) {
		// Convert to unit vectors
		double cos_ = dotProd(v1,v2)/ mag(v1) / mag(v2);
		double a = Math.acos(cos_);
		
		return a;
	}
	
	public static double angleBetweenTwoVectorsDegree(double[] v1, double[] v2) {
		
		double a = angleBetweenTwoVectors(v1, v2);
		a = Math.toDegrees(a);
		return a;
	}
	
	public static double nonReflexAngleBetweenBearings(double a1, double a2) {
		// Map to range 0-2pi
		a1 = a1%(2*Math.PI);
		a2 = a2%(2*Math.PI);
		
		double r = Math.abs(a1-a2);
	    if (r>Math.PI) {
	        r = 2*Math.PI - r;
	    }
	    return r;
	}
	
	/*
	 * Find time to collision of a point and and edge.
	 */
	public static Double edgeTTC(double[] p, double[] vP, double[] e1, double[] e2, double[] vE) {
		
		double ttc;
		
		double k_denom = 	e2[0] -	e1[0];
		double k_num = 		e2[1] - e1[1];
		
		if (k_denom==0) {
			ttc = (p[0] - e1[0]) / (vP[0] - vE[0]) ; 
		}
		else {
			double k = k_num / k_denom;
			ttc = -1*( (p[1] - e1[1]) - k*(p[0] - e1[0]) ) / ( (vP[1] - vE[1]) - k*(vP[0] - vE[0]) );
		}
		
		// Now check if collision actually occurs by check if position of point at time ttc coincides with edge
		double[] pTTC = {p[0] + vP[0]*ttc, p[1] + vP[1]*ttc};
		double[] e1TTC = {e1[0] + vE[0]*ttc, e1[1] + vE[1]*ttc};
		double[] e2TTC = {e2[0] + vE[0]*ttc, e2[1] + vE[1]*ttc};
		
		boolean xInRange = (Math.min(e1TTC[0], e2TTC[0]) <= pTTC[0]) & (pTTC[0] <= Math.max(e1TTC[0], e2TTC[0]));
		boolean yInRange = (Math.min(e1TTC[1], e2TTC[1]) <= pTTC[1]) & (pTTC[1] <= Math.max(e1TTC[1], e2TTC[1]));
		
		if (xInRange & yInRange) {
			return ttc;
		}
		else {
			return null;
		}
	}
    
	/*
	 * Find the time gap between vehicle edges and pedestrian
	 * 
	 * Time gap is the minimum time required for second arriving agent to reach intersection point
	 * Sign is attached to the value to indicate whether ped or vehicle passes first. If ped passes first then +ve,
	 * if vehicle passes first then -ve.
	 */
	public static Double edgeTG(double[] e10, double[] e11, double[] e20, double[] e21, double[] pLoc, double[] pV, double[] vV) {
		
    	// Find gradients of vehicle edges and ped direction
		// Vehicle edges are parallel so only need to calculate 1 of them
    	Double k1Denom = e11[0] - e10[0];
    	Double kpDenom = pV[0];
    	Double k1Num = e11[1] - e10[1];
    	Double kpNum = pV[1];
    	
    	// Initialise intersection coordinates that we want to find
    	double[] int1 = new double[2];
    	double[] int2 = new double[2];
    	
    	if ( (k1Denom==0) & (kpDenom==0) ) {
    		// Don't have a shared point, return null
    		return null;
    	}
    	else if ( (k1Denom==0) & (kpDenom!=0) ) {
    		// Vehicles x position is constant. Use this to find y pos of intersection
    		Double kp = kpNum / kpDenom;
    		Double cp = pLoc[1] - kp*pLoc[0];
    		
    		double y1p = kp*e10[0] + cp;
    		int1[0] = e10[0];
    		int1[1] = y1p;
    		
    		double y2p = kp*e20[0]+cp;
    		int2[0] = e20[0];
    		int2[1] = y2p;

    	}
    	else if ( (k1Denom!=0) & (kpDenom==0) ) {
    		// Peds x position is constant. Use this to find y pos of intersection
    		Double k1 = k1Num / k1Denom;
    		Double c1 = e10[1] - e10[0]*k1;
    		Double c2 = e20[1] - e20[0]*k1;
    		
    		double y1p = k1*pLoc[0] + c1;
    		int1[0] = pLoc[0];
    		int1[1] = y1p;
    		
    		double y2p = k1*pLoc[0]+c2;
    		int2[0] = pLoc[0];
    		int2[1] = y2p;
    	}
    	else {
    		// Neither gradient is infinite, can calculate using gradients
    		Double k1 = k1Num / k1Denom;
    		Double kp = kpNum / kpDenom;
    		if (k1==kp) {
    			// Moving in same direction assume never going to have shared coordinate, therefore return null
    			return null;
    		}
    		
    		Double c1 = e10[1] - e10[0]*k1;
    		Double c2 = e20[1] - e20[0]*k1;
    		Double cp = pLoc[1] - kp*pLoc[0];
    		
        	// Calculate points of intersection between line 1 and ped
        	double x1p = (c1-cp) / (kp-k1);
        	double y1p = kp*x1p + cp;
        	int1[0]=x1p;
        	int1[1]=y1p;
        	
        	// line 2 and ped
        	double x2p = (c2-cp) / (kp-k1);
        	double y2p = kp*x2p + cp;
        	int2[0]=x2p;
        	int2[1]=y2p;
    	}    	
    	
    	// Now calculate time for ped and vehicle to reach these points
    	double pS = Vector.mag(pV);
    	double vS = Vector.mag(vV);
    	if ( (vS==0) & (pS==0)) {
    		return null;
    	}
    	else if (vS==0) {
    		double t1p = (int1[0] - pLoc[0]) / pV[0];
    		double t1p_ = (int1[1] - pLoc[1]) / pV[1];
    		double t2p = (int2[0] - pLoc[0]) / pV[0];
    		double t2p_ = (int2[1] - pLoc[1]) / pV[1];
    		return FastMath.min(t1p, t2p);
    	}
    	else if (pS==0) {
    		double t1v = (int1[0] - e10[0]) / vV[0];
    		double t1v_ = (int1[1] - e10[1]) / vV[1];
    		double t2v = (int2[0] - e20[0]) / vV[0];
    		double t2v_ = (int2[1] - e20[1]) / vV[0];
    		return -1*FastMath.min(t1v, t2v);
    	}
    	else {
    		double t1p;
    		double t2p;
    		if (pV[0]!=0) {
    			t1p = (int1[0] - pLoc[0]) / pV[0];
    			t2p = (int2[0] - pLoc[0]) / pV[0];
    		}
    		else {
    			t1p = (int1[1] - pLoc[1]) / pV[1];
    			t2p = (int2[1] - pLoc[1]) / pV[1];
    		}
    		
    		double t1v;
    		double t2v;
    		if (vV[0]!=0) {
    			t1v = (int1[0] - e10[0]) / vV[0];
    			t2v = (int2[0] - e20[0]) / vV[0];
    		}
    		else {
    			t1v = (int1[1] - e10[1]) / vV[1];
    			t2v = (int2[1] - e20[1]) / vV[0];
    		}
        	
        	// Time gap is the minimum time required for second arriving agent to reach intersection point
        	// Calculated using min-max
        	double tg = FastMath.min(FastMath.max(t1p,  t1v), FastMath.max(t2p, t2v));
        	
        	// Return -1*tg if vehicle arrives first, which occurs if t1p or t2p is greater than the corresponding vehicle times
        	if ( (Double.compare(tg, t1p)==0) | (Double.compare(tg, t2p)==0) ) {
        		return -1*tg;
        	}
        	else {
        		return tg;
        	}
    	}
	}
}
