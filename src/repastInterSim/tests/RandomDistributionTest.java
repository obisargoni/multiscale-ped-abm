package repastInterSim.tests;

import org.junit.jupiter.api.Test;

import cern.jet.random.Uniform;
import cern.jet.random.engine.RandomEngine;
import repast.simphony.random.RandomHelper;

class RandomDistributionTest {

	@Test
	void testUniformGenerator() {
		RandomEngine engCASample = RandomHelper.registerGenerator("caSampleDistribution", 1);
		Uniform caSampleUniform = new Uniform(0, 1, engCASample);
		RandomHelper.registerDistribution("caSampleDistribution", caSampleUniform);

		
		RandomEngine engCASample2 = RandomHelper.registerGenerator("caSampleDistribution2", 1);
		Uniform caSampleUniform2 = new Uniform(0, 1, engCASample2);
		RandomHelper.registerDistribution("caSampleDistribution2", caSampleUniform2);
		
		// Compare values between distributions
		for (int i=0; i<50; i++) {
			double u1 = RandomHelper.getDistribution("caSampleDistribution").nextDouble();
			double u2 = RandomHelper.getDistribution("caSampleDistribution2").nextDouble();
			
			assert Double.compare(u1, u2)==0;
		}
		
		// Create anoth uniform generator witha  different seed and check that all values are not the same
		RandomEngine engCASample3 = RandomHelper.registerGenerator("caSampleDistribution3", 2);
		Uniform caSampleUniform3 = new Uniform(0, 1, engCASample3);
		RandomHelper.registerDistribution("caSampleDistribution3", caSampleUniform3);
		
		RandomEngine engCASample4 = RandomHelper.registerGenerator("caSampleDistribution4", 1);
		Uniform caSampleUniform4 = new Uniform(0, 1, engCASample4);
		RandomHelper.registerDistribution("caSampleDistribution4", caSampleUniform4);
		
		for (int i=0; i<50; i++) {
			double u4 = RandomHelper.getDistribution("caSampleDistribution4").nextDouble();
			double u3 = RandomHelper.getDistribution("caSampleDistribution3").nextDouble();
			
			assert Double.compare(u3, u4)!=0;
		}
		
	}

}

