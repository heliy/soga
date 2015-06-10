package population;

public class Sample {
	
	private int no;	
	private String name;
	
	private boolean pass;
	private double ratio = 0;
	
	private int statu = -1;
	
	public Sample(int no){
		this.no = no;
		this.name = "sample_"+no;
	}
	
	public Sample(int no, String name){
		this.no = no;
		this.name = name;
	}

	public int getNo() {
		return no;
	}

	public void setNo(int no) {
		this.no = no;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public boolean isPass() {
		return pass;
	}

	public void setPass(boolean pass) {
		this.pass = pass;
	}

	public double getRatio() {
		return ratio;
	}

	public void setRatio(double ratio) {
		this.ratio = ratio;
	}

	public int getStatu() {
		return statu;
	}

	public void setStatu(int statu) {
		this.statu = statu;
	}
	
	

}
