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

    public boolean put(T element){

        if(count < capacity){
            if(writePos >= capacity){
                writePos = 0;
            }
            elements[writePos] = element;
            writePos++;
            count++;
            return true;
        }

        return false;
    }

    public T take() {
        if(count == 0){
            return null;
        }
        int nextSlot = this.readPos();
        T nextObj = elements[nextSlot];
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
}