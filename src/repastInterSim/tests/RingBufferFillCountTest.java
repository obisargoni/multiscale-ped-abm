package repastInterSim.tests;

import org.junit.jupiter.api.Test;

import repastInterSim.util.RingBufferFillCount;

class RingBufferFillCountTest {
	
	String[] alphabet = {"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"};

	RingBufferFillCount<String> initStringBuffer(int length) {
		String[] a = new String[length];
		for(int i=0;i<length;i++) {
			a[i] = alphabet[i];
		}
		RingBufferFillCount<String> rbfc = new RingBufferFillCount<String>(a);
		return rbfc;
	}
	
	@Test
	void testInit() {
		
		// Initialise a Ring Buffer
		String[] a = new String[5];
		RingBufferFillCount<String> rb5 = new RingBufferFillCount<String>(a);
		
		assert rb5.readPos() == 0;
		assert rb5.writePos() == 0;
		assert rb5.count() == 0;
	}
	
	@Test
	void testAddElemsReset() {
		
		// Initialise a Ring Buffer
		int length = 20;
		String[] a = new String[length];
		RingBufferFillCount<String> rb = new RingBufferFillCount<String>(a);
		
		// Add some elements
		int nElems = length / 2;
		for (int i=0; i<nElems; i++) {
			rb.put(alphabet[i]);
		}
		
		assert rb.readPos() == 0;
		assert rb.writePos() == nElems;
		assert rb.count() == nElems;
		
		rb.reset();
		
		assert rb.readPos() == 0;
		assert rb.writePos() == 0;
		assert rb.count() == 0;
	}
	
	@Test
	void testAddRemoveElems() {
		
		// Initialise a Ring Buffer
		int length = 20;
		String[] a = new String[length];
		RingBufferFillCount<String> rb = new RingBufferFillCount<String>(a);
		
		// Add some elements
		int nElems = length / 2;
		for (int i=0; i<nElems; i++) {
			rb.put(alphabet[i]);
		}
		
		rb.take();
		
		assert rb.readPos() == 1;
		assert rb.writePos() == nElems;
		assert rb.count() == nElems-1;		
	}
	
	@Test
	void testAddElemsOverwrite() {
		
		// Initialise a Ring Buffer
		int length = 20;
		String[] a = new String[length];
		RingBufferFillCount<String> rb = new RingBufferFillCount<String>(a);
		
		// Add elements up to capacity
		int nElems = length;
		for (int i=0; i<nElems; i++) {
			assert rb.put(alphabet[i]);
		}
		
		// Try adding another, confirm that not possible to exceed capacity
		assert rb.put(alphabet[length+1]) == false;
		
		// Remove some elems then add some
		rb.take();
		rb.take();
		assert rb.put(alphabet[length+1]);
		assert rb.put(alphabet[length+2]);
		
		// Now test read and write pos	
		assert rb.readPos() == 2;
		assert rb.writePos() == 2;
		assert rb.count() == length;	
	}
	
	/*
	 * Tests getting element ahead when buffer has not looped round
	 */
	@Test
	void testGetElementAhead1() {
		
		// Initialise a Ring Buffer
		int length = 20;
		String[] a = new String[length];
		RingBufferFillCount<String> rb = new RingBufferFillCount<String>(a);
		
		// Add some elements
		int nElems = length / 2;
		for (int i=0; i<nElems; i++) {
			rb.put(alphabet[i]);
		}
		
		// Try getting the element ahead of 'c'
		String ea = rb.getElementAhead(2);
		
		assert ea.contentEquals("b");
		
		// Try getting the element ahead of the read position, should be null since nothing is ahead of read pos
		ea = rb.getElementAhead(rb.readPos());
		assert ea == null;
		
		// What about getting element ahead of position behind write pos
		ea = rb.getElementAhead(rb.writePos()+1);
		assert ea == null;
	}
	
	/*
	 * Tests getting element ahead when buffer has looped round
	 */
	@Test
	void testGetElementAhead2() {
		
		// Initialise a Ring Buffer
		int length = 20;
		String[] a = new String[length];
		RingBufferFillCount<String> rb = new RingBufferFillCount<String>(a);
		
		// Add some elements
		int nElems = length;
		for (int i=0; i<nElems; i++) {
			assert rb.put(alphabet[i]);
		}
		
		// take two and add two elems
		rb.take();
		rb.take();
		assert rb.put(alphabet[length+1]);
		assert rb.put(alphabet[length+2]);
		
		// Try getting the element ahead of 'l'
		String ea = rb.getElementAhead(11);
		assert ea.contentEquals("k");
		
		// Element ahead of c should now be null since c should be at head
		ea = rb.getElementAhead(2);
		assert ea == null;
		
		// Element ahead of 0 should be 20th letter of the alphabet
		ea = rb.getElementAhead(0);
		assert ea.contentEquals("t");
	}
	
	/*
	 * Tests getting element at end of buffer
	 */
	@Test
	void testGetEndElement() {
		
		// Initialise a Ring Buffer
		int length = 20;
		String[] a = new String[length];
		RingBufferFillCount<String> rb = new RingBufferFillCount<String>(a);
		
		// Add some elements
		int nElems = length-1;
		for (int i=0; i<nElems; i++) {
			assert rb.put(alphabet[i]);
		}
		
		// get element at end, 19ths character of alphabet
		assert rb.getEndElement().contentEquals("s");
		
		// take one and add two elems
		rb.take();
		assert rb.put(alphabet[nElems]);
		
		assert rb.writePos() == 0;
		
		// Get element at end, should now have looped around
		String ea = rb.getEndElement();
		assert ea.contentEquals("t");
		
		assert rb.put(alphabet[nElems+1]);
		
		ea = rb.getEndElement();
		assert ea.contentEquals("u");

	}
	
	
}
