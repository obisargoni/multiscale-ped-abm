package repastSocialForce;

import java.awt.Color;

import repast.simphony.visualizationOGL2D.DefaultStyleOGL2D;
import saf.v3d.ShapeFactory2D;

public class PedStyler extends DefaultStyleOGL2D {

	/* (non-Javadoc)
	 * @see repast.simphony.visualizationOGL2D.DefaultStyleOGL2D#init(saf.v3d.ShapeFactory2D)
	 */
	@Override
	public void init(ShapeFactory2D factory) {
		// TODO Auto-generated method stub
		super.init(factory);
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualizationOGL2D.DefaultStyleOGL2D#getColor(java.lang.Object)
	 */
	@Override
	public Color getColor(Object agent) {
		if (agent instanceof Ped) {
			Color pedCol = ((Ped)agent).getColor();
			return pedCol;
		}
		// TODO Auto-generated method stub
		return super.getColor(agent);
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualizationOGL2D.DefaultStyleOGL2D#getScale(java.lang.Object)
	 */
	@Override
	public float getScale(Object object) {
		if (object instanceof Ped) {
			int Scale = 5;
			return Scale;
		}
		// TODO Auto-generated method stub
		return super.getScale(object);
	}

}
