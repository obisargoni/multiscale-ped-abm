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
		assert SpaceBuilder.caGeography.size() == 3;
		
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
	}

}
