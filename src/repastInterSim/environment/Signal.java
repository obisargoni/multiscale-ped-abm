package repastInterSim.environment;

import repast.simphony.engine.schedule.ScheduledMethod;

public class Signal {
		
	// The itn road links this signal controls
	protected String[] itnLinkIDs;
	
	// The state of the signal will indicate which road links cars can progress on
	protected char[] state;
	
	// Nested array, each elements contains a phase state
	protected char[][] phases;
	protected int[] phaseDurations; 
	protected int phaseIndex;
	
	protected int ticksSincePhaseChange = 0;
	
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
		setState(phs[this.phaseIndex]);
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
				return this.state[i];
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
		return this.state;
	}
	
	
	/*
	 * Set the state of the Signal.
	 * 
	 * @param char[[ s
	 * 		The state to set the signal to. Array of characters, each character refers to the state the signal assigns a road link. 
	 * when the reach the signal.
	 */
	public void setState(char[] s) {
		this.state = s;
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
	
	public String[] getITNRoadLinkIDs() {
		return this.itnLinkIDs;
	}
	
	public char[][] getPhases() {
		return this.phases;
	}
	
	public int[] getPhaseDurations() {
		return this.phaseDurations;
	}
}
