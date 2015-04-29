package parameter;

public class GenoType {
	int A;
	int B;
	int N = 0;
	
	public GenoType(int ANo, int BNo){
		A = ANo;
		B = BNo;
	}
	
	public int getNo(int first, int second){
		if(second == N){ //NN
			return 0;
		}else if(first == A){ //A?
			if(second == A){  //AA
				return 1;
			}else{
				return 2;     //AB
			}
		}else{
			return 3;         //BB
		}
	}
	
	public int[] getType(int type){
		int bases[] = new int[2];
		if(type == 1){
			bases[0] = A;
			bases[1] = A;
		}else if(type == 2){
			bases[0] = A;
			bases[1] = B;
		}else if(type == 3){
			bases[0] = B;
			bases[1] = B;
		}
		return bases;
	}

}
