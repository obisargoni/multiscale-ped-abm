package repastInterSim.environment.contexts;

import org.geotools.coverage.grid.GridEnvelope2D;

import repast.simphony.context.DefaultContext;
import repastInterSim.main.GlobalVars;

public class GridEnvelopeContext extends DefaultContext<GridEnvelope2D> {
	
	public GridEnvelopeContext() {
		super(GlobalVars.CONTEXT_NAMES.GRID_ENVELOPE_CONTEXT);
	}

}