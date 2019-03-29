package styles;

import java.awt.Color;

import gov.nasa.worldwind.render.SurfacePolygon;
import gov.nasa.worldwind.render.SurfaceShape;
import repast.simphony.visualization.gis3D.style.SurfaceShapeStyle;
import repastInterSim.Road;

public class RoadStyler implements SurfaceShapeStyle<Road> {

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.SurfaceShapeStyle#getSurfaceShape(java.lang.Object, gov.nasa.worldwind.render.SurfaceShape)
	 */
	@Override
	public SurfaceShape getSurfaceShape(Road object, SurfaceShape shape) {
		// TODO Auto-generated method stub
		return new SurfacePolygon();
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.SurfaceShapeStyle#getFillColor(java.lang.Object)
	 */
	@Override
	public Color getFillColor(Road obj) {
		// TODO Auto-generated method stub
		return Color.magenta;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.SurfaceShapeStyle#getFillOpacity(java.lang.Object)
	 */
	@Override
	public double getFillOpacity(Road obj) {
		// TODO Auto-generated method stub
		return 0.25;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.SurfaceShapeStyle#getLineColor(java.lang.Object)
	 */
	@Override
	public Color getLineColor(Road obj) {
		// TODO Auto-generated method stub
		return Color.BLACK;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.SurfaceShapeStyle#getLineOpacity(java.lang.Object)
	 */
	@Override
	public double getLineOpacity(Road obj) {
		// TODO Auto-generated method stub
		return 1.0;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.SurfaceShapeStyle#getLineWidth(java.lang.Object)
	 */
	@Override
	public double getLineWidth(Road obj) {
		// TODO Auto-generated method stub
		return 3;
	}

}
