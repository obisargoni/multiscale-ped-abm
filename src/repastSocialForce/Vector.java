package repastSocialForce;

import org.apache.commons.math3.util.FastMath;

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

}
