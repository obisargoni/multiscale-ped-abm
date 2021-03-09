package repastInterSim.environment;

import java.util.HashMap;

import repast.simphony.engine.schedule.ScheduledMethod;

public class Signal {
	
	protected String signalID; // The id of the road link this road object should be linked with
	
	// The state of the signal will indicate which road links cars can progress on
	protected HashMap<String, Boolean> state = new HashMap<String, Boolean>();;
	
	protected HashMap<String, Boolean>[] phases;
	protected int[] phaseDurations; 
	protected int phaseIndex;
	
	protected int ticksSincePhaseChange = 0;
	
	/*
	 * Blank instance method.
	 */
	public Signal() {
	}
	
	public Signal(HashMap<String, Boolean>[] p, int[] pD, int pI) {
		init(p, pD, pI);
	}
	
	public Signal(HashMap<String, Boolean>[] p, int[] pD) {
		init(p, pD, 0);
	}
	
	private void init(HashMap<String, Boolean>[] p, int[] pD, int pI) {
		this.phases = p;
		this.phaseDurations = pD;
		this.phaseIndex = pI;
		setState(p[this.phaseIndex]);
	}
		
	/*
	 * Get the state of the signal for a given road link
	 * 
	 * @param String rlID
	 * 		The road link ID to get the signal state for
	 * 
	 * @return Boolean. The state of the signal
	 */
	public boolean getState(String rlID) {
		return this.state.get(rlID);
	}
	
	/*
	 * Set the state of the Signal.
	 * 
	 * @param HashMap<String, Boolean> s
	 * 		The state to set the signal to. Contains lookup from road link IDs to indicator of whether vehicles on this road links can progress 
	 * when the reach the signal.
	 */
	public void setState(HashMap<String, Boolean> s) {
		this.state = s;
	}
	
	/*
	 * Reverse the state of the signal.
	 */
	public void switchState() {
		this.state.replaceAll((k,v) -> !v);
	}
	
	public void switchState(String rlID) {
		Boolean switched = !this.state.get(rlID);
		this.state.replace(rlID, switched);
	}
	
	public void switchState(String[] rlIDs) {
		for(int i=0; i< rlIDs.length; i++) {
			this.switchState(rlIDs[i]);
		}
	}
	
	
	/*
	 * Change the current state of the signal to be the next phase when the current phases duration has ended.
	 */
	@ScheduledMethod(start = 0, interval = 1)
	public void changePhase() {
		// Check if it's time to change phase
		if(this.ticksSincePhaseChange >= this.phaseDurations[this.phaseIndex]) {
			this.ticksSincePhaseChange = 0;
			this.phaseIndex++;
			if(this.phaseIndex == this.phases.length) {
				this.phaseIndex = 0;
			}
			
			this.state = this.phases[this.phaseIndex];
		}
		else {
			this.ticksSincePhaseChange+=1;
		}
	}
	
	public String getSignalID() {
		return this.signalID;
	}
	
	public void setSignalID(String id) {
		this.signalID = id;
	}
}
