package repastInterSim.tests;

import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import org.junit.jupiter.api.Test;

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
	}

}
