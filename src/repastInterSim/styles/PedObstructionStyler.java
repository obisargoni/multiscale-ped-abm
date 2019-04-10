package repastInterSim.styles;

import java.awt.Color;

import gov.nasa.worldwind.render.SurfacePolyline;
import gov.nasa.worldwind.render.SurfaceShape;
import repast.simphony.visualization.gis3D.style.SurfaceShapeStyle;
import repastInterSim.environment.PedObstruction;

public class PedObstructionStyler implements SurfaceShapeStyle<PedObstruction> {
	
	@Override
	public SurfaceShape getSurfaceShape(PedObstruction object, SurfaceShape shape) {
	  return new SurfacePolyline();
	}

	@Override
	public Color getFillColor(PedObstruction obj) {
		return null;
	}

	@Override
	public double getFillOpacity(PedObstruction obj) {
		return 0;
	}

	@Override
	public Color getLineColor(PedObstruction obj) {
		return Color.RED;
	}

	@Override
	public double getLineOpacity(PedObstruction obj) {
		return 1.0;
	}

	@Override
	public double getLineWidth(PedObstruction obj) {
		return 2;
	}
}
