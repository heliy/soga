package util;

public class AllotThreads {
	public int[] allot(int total, int threads){
		int[] nums = new int[threads];
		for(int i = 0; i < total; i++){
			nums[i%threads]++;
		}
		return nums;
	}

}
