package util;

import java.util.concurrent.CountDownLatch;

public class MultiThreads<T extends Runnable> implements Runnable {
	private T[] contains;
	private CountDownLatch latch;
	
	public MultiThreads(T[] contains){
		this.setContains(contains);
	}
	
	public void setRun(CountDownLatch latch){
		this.latch = latch;		
	}
	
	public void run() {
		for(T t: this.contains){
			t.run();
		}
		latch.countDown();
	}

	public T[] getContains() {
		return contains;
	}

	public void setContains(T[] contains) {
		this.contains = contains;
	}

}
