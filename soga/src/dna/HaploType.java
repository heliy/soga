package dna;

import parameter.Base;

public class HaploType implements Comparable<HaploType> {
	private int[] alleles;
	private int caseNum;
	private int controlNum;
	private int unknownNum;
	private double pvalue = -1;
	private int sampleNo;
	private double or;
	private double orLow;
	private double orHigh;

	public HaploType(int[] lleles, int sampleno) {
		alleles = lleles;
		caseNum = 0;
		controlNum = 0;
		unknownNum = 0;
		sampleNo = sampleno;
	}

	public HaploType(int sampleno) {
		alleles = new int[0];
		caseNum = 0;
		controlNum = 0;
		unknownNum = 0;
		this.sampleNo = sampleno;
	}

	public HaploType clone() {
		HaploType ht = new HaploType(this.alleles, this.sampleNo);
		ht.setPvalue(pvalue, or, orLow, orHigh);
		ht.caseNum = this.caseNum;
		ht.controlNum = this.controlNum;
		ht.unknownNum = this.unknownNum;
		return ht;
	}

	public int getNum() {
		return caseNum + controlNum + unknownNum;
	}

	public double getPvalue() {
		return pvalue;
	}

	public double getOR() {
		return or;
	}

	public double getOrLow() {
		return orLow;
	}

	public double getOrHigh() {
		return orHigh;
	}

	public void setPvalue(double pvalue, double or, double orLow, double orHigh) {
		this.pvalue = pvalue;
		this.or = or;
		this.orLow = orLow;
		this.orHigh = orHigh;
	}

	public int getCaseNum() {
		return caseNum;
	}

	public int getControlNum() {
		return controlNum;
	}

	public int getUnknownNum() {
		return unknownNum;
	}

	public int[] getAlleles() {
		return alleles;
	}

	public String toString(int total) {
		Base base = new Base();
		StringBuilder sb = new StringBuilder();
		if (pvalue < 0) {
			sb.append("----\t----");
		} else {
			sb.append(pvalue + "\t" + caseNum + "/" + controlNum);
		}
		sb.append("\t" + (double) getNum() / total);
		for (int c : alleles) {
			sb.append("\t" + base.getBase(c));
		}
		if (pvalue >= 0) {
			sb.append("\t" + or + "\t" + orLow + "\t" + orHigh);
		}
		sb.append("\n");
		return sb.toString();
	}

	public boolean equals(HaploType obj) {
		int i, len = alleles.length;
		int[] another = obj.getAlleles();
		if (len != another.length) {
			return false;
		}
		for (i = 0; i < len; i++) {
			if (alleles[i] != another[i]) {
				return false;
			}
		}
		return true;
	}

	public void add(int statu) {
		switch (statu) {
		case -1:
			addUnknown();
			return;
		case 0:
			addControl();
			return;
		case 1:
			addCase();
			return;
		default:
			return;
		}
	}

	public void addUnknown() {
		unknownNum++;
	}

	public void addCase() {
		caseNum++;
	}

	public void addControl() {
		controlNum++;
	}

	public int compareTo(HaploType other) { // 数目降序 数目相同时按碱基类型排列
		int d = other.getNum() - this.getNum(), i;
		if (d != 0) {
			return d;
		}
		int l = alleles.length;
		int[] o = other.getAlleles();
		for (i = 0; i < l; i++) {
			if (alleles[i] != o[i]) {
				return alleles[i] - o[i];
			}
		}
		return 0;
	}

	public int getSampleNo() {
		return sampleNo;
	}

	public void extend(HaploType next) {
		if(next.getSampleNo() != this.sampleNo){
			return;
		}
		
		int[] nextAlleles = next.getAlleles();
		int[] newAlleles = new int[nextAlleles.length+this.alleles.length];
		int i, j, n;
		for(i = 0, n = this.alleles.length; i < n; i++){
			newAlleles[i] = this.alleles[i];
		}
		for(j = 0, n = nextAlleles.length; j < n; j++){
			newAlleles[i+j] = nextAlleles[j];
		}
		this.alleles = newAlleles;		
	}

}
