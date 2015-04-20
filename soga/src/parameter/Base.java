package parameter;

public class Base {
	char[] bases = {'N', 'A', 'T', 'G', 'C'};

	public int N(){
		return 0;
	}
	
	public int A(){
		return 1;
	}
	public int T(){
		return 2;
	}
	public int G(){
		return 3;
	}
	public int C(){
		return 4;
	}

	
	public int baseNo(char L){
		int i;
		for(i = 0; i < 5; i++){
			if(L == bases[i]){
				return i;
			}
		}
		return -1;
	}
	
	public char getBase(int i){
		return bases[i];
	}
};
