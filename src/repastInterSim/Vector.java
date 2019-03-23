package repastInterSim;

import org.apache.commons.math3.util.FastMath;
import org.geotools.factory.Hints;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import com.vividsolutions.jts.geom.Geometry;

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
        
        // Need to get the correct angle
        if (signX == 1 & signY == 1) {
        	return a;
        }
        else if (signY == -1) {
        	return a + Math.PI;
        }
        else if (signX == -1 & signY == 1) {
        	return 2*Math.PI - a;
        }
        else {
        	return a;
        }
    }
    
}
