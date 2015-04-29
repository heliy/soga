package parameter;

public class Setting {
	public int HEAD = 4;
	public double NNRatio = 0.05;
	public double MAF = 0.01;
	public double HWE = 0.05;
	public double MAXN = 0.8;
	public int SAMPLES;
	public long SNPS;
	
	public int getHEAD() {
		return HEAD;
	}

	public void setHEAD(int hEAD) {
		HEAD = hEAD;
	}

	public double getMAF() {
		return MAF;
	}

	public void setMAF(double mAF) {
		MAF = mAF;
	}

	public double getHWE() {
		return HWE;
	}

	public void setHWE(double hWE) {
		HWE = hWE;
	}

	public double getMAXN() {
		return MAXN;
	}

	public void setMAXN(double mAXN) {
		MAXN = mAXN;
	}

	public int[] getValids() {
		return valids;
	}

	public void setValids(int[] valids) {
		this.valids = valids;
	}

	public void setNNRatio(double nNRatio) {
		NNRatio = nNRatio;
	}

	public long getSNPS() {
		return SNPS;
	}

	public void setSNPS(long sNPS) {
		SNPS = sNPS;
	}

	public int[] valids;
	public int valid;

	public int getSAMPLES() {
		return SAMPLES;
	}

	public void setSAMPLES(int sAMPLES) {
		SAMPLES = sAMPLES;
	}

	public String getFileSplit() {
		// TODO Auto-generated method stub
		return " ";
	}

	public double getNNRatio() {
		return NNRatio;
	}
	

}
