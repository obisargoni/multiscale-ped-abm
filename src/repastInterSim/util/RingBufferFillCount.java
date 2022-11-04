package repastInterSim.util;

public class RingBufferFillCount<T> {
    
	public T[] elements = null;

    private int capacity  = 0;
    private int writePos  = 0;
    private int count = 0;

    public RingBufferFillCount(T[] array) {
        this.elements = array;
        this.capacity = array.length;        
    }

    public void reset() {
        this.writePos = 0;
        this.count = 0;
    }

    public int capacity() { return this.capacity; }
    public int count(){ return this.count; }

    public int remainingCapacity() {
        return this.capacity - this.count;
    }
    
    public int writePos() {
    	return this.writePos;
    }
    
    public int readPos() {
		int nextSlot = writePos - count;
		if(nextSlot < 0){
		    nextSlot += capacity;
		}
		return nextSlot; 
    }

    public Integer put(T element){
    	Integer elemPos = null;
        if(count < capacity){
            if(writePos >= capacity){
                writePos = 0;
            }
            elements[writePos] = element;
            elemPos = writePos;
            writePos++;
            count++;
        }

        return elemPos;
    }

    public T take() {
        if(count == 0){
            return null;
        }
        int nextSlot = this.readPos();
        T nextObj = elements[nextSlot];
        elements[nextSlot] = null; // Need to set to null so that vehicle object does not have references pointing to it
        count--;
        return nextObj;
    }

	public T getElementAhead(int elemPos) {
		if (this.readPos() == elemPos) {
			return null;
		}
		else {
			int inFront = elemPos - 1;
			if (inFront < 0) {
				inFront = capacity-1;
			}
			return elements[inFront];
		}
	}
	
	private int endPos() {
		int end = this.writePos - 1;
		if (end<0) {
			end+=capacity;
		}
		return end;
	}
	
	public T getEndElement() {
		if (this.count==0) {
			return null;
		}
		int end = this.endPos();
		return this.elements[end];
	}
	
	public T getStartElement() {
		if (this.count==0) {
			return null;
		}
		int start = this.readPos();
		return this.elements[start];
	}
	
	public boolean hasCapacity() {
		return this.capacity > this.count;
	}
}