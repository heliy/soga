package dna;

public class LD {
	private double rsq;
	private double dprime;
	private double ciLow;
	private double ciHigh;
	private double lod;
	private int type;

	public LD(double crsq, double cdprime, double cciLow,
			double cciHigh, double clod, int ctype) {
		rsq = crsq;
		dprime = cdprime;
		ciLow = cciLow;
		ciHigh = cciHigh;
		lod = clod;
		type = ctype;
	}

	public double getRsq() {
		return rsq;
	}

	public double getDprime() {
		return dprime;
	}

	public double getCiLow() {
		return ciLow;
	}

	public double getCiHigh() {
		return ciHigh;
	}

	public double getLod() {
		return lod;
	}

	public int getType() {
		return type;
	}

	public String toString() {
		String s = dprime + "\t" + rsq + "\t" + ciLow + "\t" + ciHigh + "\t" + lod;
		if(type == 0){
			return s + "\totherwise\n";
		}else if(type == 1){
			return s + "\tLD\n";
		}else{
			return s + "\tEHR\n";
		}
	}

}
