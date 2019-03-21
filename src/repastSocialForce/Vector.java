package repastSocialForce;

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
     * Transform the coordinate reference system of the geometry
     * 
     * @param Geometry G. The geometry object to transform
     * @param String epsgIn. The CRS of the input geometry in the format "EPSG:xxxx".
     * @param String epsgOut. THe CRS to transform the input geometry to in the format "EPSG:xxxx"
     * 
     * @return Geometry. The transformed geometry object.
     */
    public static Geometry geomCRSTrans(Geometry G, String epsgIn, String epsgOut) {

		CoordinateReferenceSystem sourceCRS = null;
		try {
			sourceCRS = CRS.decode(epsgIn);
		} catch (FactoryException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		CoordinateReferenceSystem targetCRS = null;
		try {
			targetCRS = CRS.decode(epsgOut);
		} catch (FactoryException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Not sure what this line does and whether it is required
		Hints.putSystemDefault(Hints.FORCE_LONGITUDE_FIRST_AXIS_ORDER, Boolean.TRUE);
		
		MathTransform transform = null;
		try {
			transform = CRS.findMathTransform(sourceCRS, targetCRS);
		} catch (FactoryException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Geometry targetGeometry = null;
		try {
			targetGeometry = JTS.transform( G, transform);
		} catch (MismatchedDimensionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (TransformException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return targetGeometry;
    }
}
