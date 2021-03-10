package repastInterSim.environment;

import com.vividsolutions.jts.geom.Coordinate;

import repast.simphony.engine.schedule.ScheduledMethod;

public class Signal {
		
	// The itn road links this signal controls
	protected String[] itnLinkIDs;
	
	// Nested array, each elements contains a phase state
	protected char[][] phases;
	protected int[] phaseDurations;
	
	// Used to access the state of the signal, which indicates which road links cars can progress on
	protected int phaseIndex;
	
	protected int ticksSincePhaseChange = 0;
	
	protected Coordinate signalLoc;
	
	/*
	 * Blank instance method.
	 */
	public Signal() {
	}
	
	public Signal(String[] rlIDs, char[][] phs, int[] pD, int pI) {
		init(rlIDs, phs, pD, pI);
	}
	
	public Signal(String[] rlIDs, char[][] phs, int[] pD) {
		init(rlIDs, phs, pD, 0);
	}
	
	private void init(String[] rlIDs, char[][] phs, int[] pD, int pI) {
		this.itnLinkIDs = rlIDs;
		this.phases = phs;
		this.phaseDurations = pD;
		this.phaseIndex = pI;
	}
		
	/*
	 * Get the state of the signal for a given road link
	 * 
	 * @param String rlID
	 * 		The road link ID to get the signal state for
	 * 
	 * @return char. The state of the signal
	 */
	public char getState(String rlID) {
		// get index this road link appears at and return state value at this index
		for(int i=0; i<this.itnLinkIDs.length; i++) {
			if(this.itnLinkIDs[i].contentEquals(rlID)) {
				return this.phases[this.phaseIndex][i];
			}
		}
		// return u for unknown if can't match road link id
		return 'u';
	}
	
	/*
	 * Get the start across all itn road links this signal controls
	 * 
	 * @return char[]
	 * 		Array of characters indicating state for each road link this signal controls
	 */
	public char[] getState() {
		return this.phases[this.phaseIndex];
	}
	
	
	/*
	 * Set the state of the Signal.
	 * 
	 * @param char[[ s
	 * 		The state to set the signal to. Array of characters, each character refers to the state the signal assigns a road link. 
	 * when the reach the signal.
	 */
	public void setPhaseIndex(int i) {
		this.phaseIndex = i;
	}
	
	/*
	 * Change the current state of the signal to be the next phase when the current phases duration has ended.
	 */
	@ScheduledMethod(start = 0, interval = 1)
	public void step() {
		// Check if it's time to change phase
		if(this.ticksSincePhaseChange >= this.phaseDurations[this.phaseIndex]) {
			this.ticksSincePhaseChange = 0;
			this.phaseIndex++;
			if(this.phaseIndex == this.phases.length) {
				this.phaseIndex = 0;
			}
		}
		else {
			this.ticksSincePhaseChange+=1;
		}
	}
	
	public String[] getITNRoadLinkIDs() {
		return this.itnLinkIDs;
	}
	
	public char[][] getPhases() {
		return this.phases;
	}
	
	public int[] getPhaseDurations() {
		return this.phaseDurations;
	}
	
	public Coordinate getSignalLoc() {
		return this.signalLoc;
	}
}
