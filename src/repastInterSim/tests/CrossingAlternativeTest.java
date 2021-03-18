package repastInterSim.tests;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import org.junit.jupiter.api.Test;

import repastInterSim.environment.CrossingAlternative;
import repastInterSim.main.SpaceBuilder;

class CrossingAlternativeTest {

	@Test
	void testLoadCrossingAlternatives1() {
		try {
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		assert SpaceBuilder.caGeography != null;
		assert SpaceBuilder.caGeography.size() == 5;
		
		CrossingAlternative ca1 = null;
		for (CrossingAlternative ca_: SpaceBuilder.caGeography.getAllObjects()) {
			if(ca_.getID() == 3) {
				ca1 = ca_;
				break;
			}
		}
		
		// Test ca1 setup as expected
		assert ca1.getRoadLinkID().contentEquals("CCBB6F91-5D00-46AA-8835-CC06A8D97B91_0");
		
		String[] expectedITNLinks = {"osgb4000000030238817","osgb4000000030238818"};
		for (int i=0; i<expectedITNLinks.length; i++) {
			assert ca1.getITNRoadLinkIDs()[i].contentEquals(expectedITNLinks[i]);
		}
		
		char[] phase1 = {'r','r'};
		char[] phase2 = {'g','g'};
		char[][] expectedPhases = {phase1, phase2};
		for(int i=0; i< expectedPhases.length; i++) {
			for (int j=0; j<expectedPhases[i].length; j++) {
				assert ca1.getPhases()[i][j] == expectedPhases[i][j];
			}
		}
		
		int[] expectedDurations = {10,30};
		for (int i=0; i<expectedDurations.length; i++) {
			assert ca1.getPhaseDurations()[i] == expectedDurations[i];
		}
		
		// Finally check initial state of signal
		assert ca1.getState("osgb4000000030238817") == 'r';
		assert ca1.getState("osgb4000000030238818") == 'r';
		assert ca1.getState("wrong_link_id") == 'u';
	}
	
	@Test
	void testSignalPhaseChange() {
		try {
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		assert SpaceBuilder.caGeography != null;
		assert SpaceBuilder.caGeography.size() == 5;
		
		CrossingAlternative ca1030 = null;
		CrossingAlternative ca530a = null;
		CrossingAlternative ca530b = null;
		for (CrossingAlternative ca_: SpaceBuilder.caGeography.getAllObjects()) {
			if(ca_.getID() == 3) {
				ca1030 = ca_;
			}
			else if (ca_.getID() == 4) {
				ca530a = ca_;
			}
			else if (ca_.getID() == 5) {
				ca530b = ca_;
			}
		}
		
		// Check initial states of the signals then check that state changes as expected
		assert ca1030.getState("osgb4000000030238817") == 'r';
		assert ca1030.getState("osgb4000000030238818") == 'r';
		
		assert ca530a.getState("osgb4000000030238839") == 'r';
		assert ca530b.getState("osgb4000000030418774") == 'r';
		
		// After 5 ie on the 6th tick) ticks two of the signals should change to green
		for (int i=0; i<6; i++) {
			ca1030.step();
			ca530a.step();
			ca530b.step();
		}
		
		assert ca1030.getState("osgb4000000030238817") == 'r';
		assert ca1030.getState("osgb4000000030238818") == 'r';
		
		assert ca530a.getState("osgb4000000030238839") == 'g';
		assert ca530b.getState("osgb4000000030418774") == 'g';
		
		// After 5 more ticks the other signal should also be green
		for (int i=0; i<6; i++) {
			ca1030.step();
			ca530a.step();
			ca530b.step();
		}

		assert ca1030.getState("osgb4000000030238817") == 'g';
		assert ca1030.getState("osgb4000000030238818") == 'g';
		
		assert ca530a.getState("osgb4000000030238839") == 'g';
		assert ca530b.getState("osgb4000000030418774") == 'g';
		
		// After a further 25 steps 530 signals should revert to red
		for (int i=0; i<26; i++) {
			ca1030.step();
			ca530a.step();
			ca530b.step();
		}
		
		assert ca1030.getState("osgb4000000030238817") == 'g';
		assert ca1030.getState("osgb4000000030238818") == 'g';
		
		assert ca530a.getState("osgb4000000030238839") == 'r';
		assert ca530b.getState("osgb4000000030418774") == 'r';
		
		// After 5 more ticks 1030 signal should be on red whilst the others are back to green
		for (int i=0; i<6; i++) {
			ca1030.step();
			ca530a.step();
			ca530b.step();
		}
		assert ca1030.getState("osgb4000000030238817") == 'r';
		assert ca1030.getState("osgb4000000030238818") == 'r';
		
		assert ca530a.getState("osgb4000000030238839") == 'g';
		assert ca530b.getState("osgb4000000030418774") == 'g';
		
	}
	
	@Test
	void testUnsignalisedCrossings() {
		try {
			EnvironmentSetup.setUpCrossingAlternatives("CrossingAlternatives.shp");
		} catch (MalformedURLException | FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		assert SpaceBuilder.caGeography != null;
		assert SpaceBuilder.caGeography.size() == 5;
		
		CrossingAlternative us1 = null;
		CrossingAlternative us2 = null;
		for (CrossingAlternative ca_: SpaceBuilder.caGeography.getAllObjects()) {
			if(ca_.getID() == 1) {
				us1 = ca_;
			}
			else if (ca_.getID() == 2) {
				us2 = ca_;
			}
		}
		
		// Test that these crossing alternatives return 'u' when trying to get their state
		assert us1.getState("osgb4000000030238839") == 'u';
		assert us2.getState("osgb4000000030343781") == 'u';
		
	}
}
