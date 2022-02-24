package repastInterSim.environment.contexts;

import repast.simphony.context.DefaultContext;
import repastInterSim.datasources.DefaultDataRecorder;
import repastInterSim.main.GlobalVars;

public class DataRecorderContext extends DefaultContext<DefaultDataRecorder> {

	public DataRecorderContext() {
		super(GlobalVars.CONTEXT_NAMES.DATA_RECORDER_CONTEXT);
	}

}
