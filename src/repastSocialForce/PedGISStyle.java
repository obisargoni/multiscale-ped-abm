package repastSocialForce;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.image.BufferedImage;
import java.util.HashMap;
import java.util.Map;

import org.geotools.renderer.style.MarkStyle2D;

import gov.nasa.worldwind.WorldWind;
import gov.nasa.worldwind.avlist.AVKey;
import gov.nasa.worldwind.render.BasicWWTexture;
import gov.nasa.worldwind.render.Material;
import gov.nasa.worldwind.render.Offset;
import gov.nasa.worldwind.render.PatternFactory;
import gov.nasa.worldwind.render.WWTexture;
import repast.simphony.visualization.gis3D.PlaceMark;
import repast.simphony.visualization.gis3D.style.MarkStyle;

public class PedGISStyle implements MarkStyle<Ped> {
	
	private Map<String, WWTexture> textureMap;
	private Offset labelOffset;
	
	public PedGISStyle(){
		
		labelOffset = new Offset(1.2d, 0.6d, AVKey.FRACTION, AVKey.FRACTION);

		
		/**
		 * Use of a map to store textures significantly reduces CPU and memory use
		 * since the same texture can be reused.  Textures can be created for different
		 * agent states and re-used when needed.
		 */		
		textureMap = new HashMap<String, WWTexture>();
		
		BufferedImage image = PatternFactory.createPattern(PatternFactory.PATTERN_CIRCLE, 
				new Dimension(50, 50), 0.7f,  Color.BLUE);
		
		textureMap.put("blue circle", new BasicWWTexture(image));
		
		image = PatternFactory.createPattern(PatternFactory.PATTERN_CIRCLE, 
				new Dimension(10, 10), 0.7f,  Color.RED);
		
		textureMap.put("red circle", new BasicWWTexture(image));
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getTexture(java.lang.Object, gov.nasa.worldwind.render.WWTexture)
	 */
	@Override
	public WWTexture getTexture(Ped object, WWTexture texture) {
		/*
		if (((Ped)object).getColor() == Color.RED){
			return textureMap.get("red circle");
		}
		else{
			return textureMap.get("blue circle");
		}
		*/
		
		return textureMap.get("blue circle");
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getElevation(java.lang.Object)
	 */
	@Override
	public double getElevation(Ped obj) {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getPlaceMark(java.lang.Object, repast.simphony.visualization.gis3D.PlaceMark)
	 */
	@Override
	public PlaceMark getPlaceMark(Ped object, PlaceMark mark) {
		// TODO Auto-generated method stub
		// PlaceMark is null on first call.
		if (mark == null)
			mark = new PlaceMark();
		
		mark.setAltitudeMode(WorldWind.RELATIVE_TO_GROUND);
		mark.setLineEnabled(false);
		
		return mark;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getIconOffset(java.lang.Object)
	 */
	@Override
	public Offset getIconOffset(Ped obj) {
		// TODO Auto-generated method stub
		return Offset.CENTER;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getScale(java.lang.Object)
	 */
	@Override
	public double getScale(Ped obj) {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getHeading(java.lang.Object)
	 */
	@Override
	public double getHeading(Ped obj) {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getLabel(java.lang.Object)
	 */
	@Override
	public String getLabel(Ped obj) {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getLabelColor(java.lang.Object)
	 */
	@Override
	public Color getLabelColor(Ped obj) {
		// TODO Auto-generated method stub
		return ((Ped)obj).getColor();
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getLabelFont(java.lang.Object)
	 */
	@Override
	public Font getLabelFont(Ped obj) {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getLabelOffset(java.lang.Object)
	 */
	@Override
	public Offset getLabelOffset(Ped obj) {
		// TODO Auto-generated method stub
		return labelOffset;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getLineWidth(java.lang.Object)
	 */
	@Override
	public double getLineWidth(Ped obj) {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see repast.simphony.visualization.gis3D.style.MarkStyle#getLineMaterial(java.lang.Object, gov.nasa.worldwind.render.Material)
	 */
	@Override
	public Material getLineMaterial(Ped obj, Material lineMaterial) {
		// TODO Auto-generated method stub
		if (lineMaterial == null){
			lineMaterial = new Material(Color.RED);
		}
		
		return lineMaterial;
	}
	
	

}
