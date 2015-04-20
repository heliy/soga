package dna;

public class Scope implements Comparable<Scope> {

	private int begin;
	private int end;
	private int len;
	
	public Scope(int b,int e){
		begin = b;
		end = e;
		len = e - b + 1;
	}
	
	public int compareTo(Scope other) { // 降序
		if(len < other.len){
			return 1;
		}
		if(len > other.len){
			return -1;
		}
		return 0;
	}
	
	public String toString(){
		return len + " : " + begin + ", " + end;
	}
	
	public int getBegin() {
		return begin;
	}

	public int getEnd() {
		return end;
	}
}
