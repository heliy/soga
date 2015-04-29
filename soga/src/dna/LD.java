package dna;

public class LD {
	double rsq;
	double dprime;
	double ciLow;
	double ciHigh;
	double lod;	
	
	public LD(double crsq, double cdprime, double cciLow, double cciHigh, double clod){
		rsq = crsq;
		dprime = cdprime;
		ciLow = cciLow;
		ciHigh = cciHigh;
		lod = clod;
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
	
}
