package dna;

import dna.LD;
import exceptions.GenoTypeException;
import parameter.Setting;
import parameter.Base;

import java.util.Random;
import java.util.WeakHashMap;

import calculate.FitTest;
import calculate.OR;
import calculate.PairwiseCal;

public class Snp {
	private static int NN = 0;
	private static int AA = 1;
	private static int AB = 2;
	private static int BB = 3;
	private static int N = 0;
	
	private String rs; // rs + number
	private int A;
	private int B;
	private String chr; // chromesome
	private long position; // position in chromesome
	private Base base;
	private int genos[]; // genotypes data
	private int types[][]; // 与AB对应的基因型 用来计算LD
	private WeakHashMap<Snp, LD> pws; // pairwise to other snps
	private int alleleCount[]; // count of alleles
	private double alleleFreq[]; // frequency of alleles
	private double maf; // min allele frequency
	private int genoCount[]; // count of genotypes of genos
	private double genoFreq[]; // frequency of genotypes
	private double hw; // hardy-weinberg p value
	private double aPValue; // allele A p value of case/control test
	private double bPValue; // allele b p value of case/control test
	private double or; // ORrate
	private double orLow;  // OR置信区间
	private double orHigh; // OR置信区间
	private boolean isBad; // is bad snp

	public Snp(String rsNo, int a, int b, String Chr, String sposition,
			String[] parts, Setting rc) throws GenoTypeException {
		// rs, alleles, chr, position, genos....
		rs = rsNo;
		position = Integer.parseInt(sposition); // 所在位置
		chr = Chr;
		base = new Base();
		A = a;
		B = b;
		typeCount(parts);
		maf = Math.min(alleleFreq[A], alleleFreq[B]); // 最小等位频率
		hw = calHW(); // HW平衡p值
		isBad = judgeBad(rc); // 是否通过QC
		pws = new WeakHashMap<Snp, LD>(); // 与其他snp的LD值
		or = -1;
		if(rc.doCC()){
			setCC(rc.getStatus(), rc.getFreqs());
		}
	}
	
	public int hashCode(){
		return (int)position;
	}
	

	public String toString() {
		String s = rs + "\t" + base.getBase(A) + "/" + base.getBase(B) + "\t"
				+ chr + "\t" + position + "\t" + base.getBase(A) + "\t"
				+ alleleCount[A] + "\t" + alleleFreq[A] + "\t"
				+ base.getBase(B) + "\t" + alleleCount[B] + "\t"
				+ alleleFreq[B] + "\t" + alleleCount[N] + "\t"
				+ alleleFreq[N] + "\t" + genoCount[AA] + "\t"
				+ genoFreq[AA] + "\t" + genoCount[AB] + "\t" + genoFreq[AB] + "\t"
				+ genoCount[BB] + "\t" + genoFreq[BB] + "\t" + genoCount[NN]
				+ "\t" + genoFreq[NN] + "\t" + maf + "\t" + hw + "\t";
		if(or >= 0){
			s += aPValue+"\t"+bPValue+"\t"+or+"\t"+orLow+"\t"+orHigh+"\t";
		}
		if (isBad) {
			s += "NO";
		} else {
			s += "YES";
		}
		s += "\n";
		return s;
	}
	
	// 0/1等位
	public int getDis(int s){
		if(s == A){
			return 0;
		}else if(s == B){
			return 1;
		}else{
			return -1;
		}
	}
	
	public int[][] getTypes() {
		return types;
	}
		
	public int[] getAlleleCount() {
		// return [N, A, B]
		int[] alls = new int[3];
		alls[N] = alleleCount[N];
		alls[1] = alleleCount[A];
		alls[2] = alleleCount[B];
		return alls;
	}

	public double getOR(){
		return or;
	}
	
	public double getOrLow() {
		return orLow;
	}

	public double getOrHigh() {
		return orHigh;
	}

	public String getRs() {
		return rs;
	}

	public long getPosition() {
		return position;
	}

	public double getMaf() {
		return maf;
	}

	public int[] getGenoCount() {
		return genoCount;
	}

	public double[] getGenoFreq() {
		return genoFreq;
	}

	public double getHWE() {
		return hw;
	}

	public boolean isBad() {
		return isBad;
	}
	
	public int getA() {
		return A;
	}

	public int getB() {
		return B;
	}

	public String getChr() {
		return chr;
	}
	
	public boolean hasLD(Snp another){
		if(pws.containsKey(another)){
			return true;
		}else{
			return false;
		}
		
	}
	
	public LD getLD(PairwiseCal cal, Snp another) {
		if(pws.containsKey(another)){
			return pws.get(another);
		}else{
			LD ld = cal.pw(this, another);
			this.setLD(another, ld);
			another.setLD(this, ld);
			return ld;
		}
	}
	
	private void setLD(Snp another, LD ld){
		pws.put(another, ld);
	}

	private void typeCount(String[] parts) throws GenoTypeException {
		String s;
		int i, len = parts.length;
		alleleCount = new int[5];
		genoCount = new int[4];
		genos = new int[len];
		types = new int[len][2];
		for (i = 0; i < len; i++) {
			s = parts[i];
			int first = base.baseNo(s.charAt(0));
			int second = base.baseNo(s.charAt(1));
//			System.out.println(A+", "+B);
			if((first != A && first != B && first != N) || (second != A && second != B && second != N)){
				throw new GenoTypeException(s);
			}
			
			if(first == N || second == N){
				// 缺失的
				genoCount[NN]++;
				alleleCount[N] += 2;
				types[i][0] = types[i][1] = -1;
				genos[i] = -1;
				continue;
			}
			if (first > second) { // 保证 first 的序号小于 second
				int tmp = second;
				second = first;
				first = tmp;
			}
			
			alleleCount[first]++;
			alleleCount[second]++;
			types[i][0] = (first == A) ? 0 : 1;
			types[i][1] = (second == A) ? 0 : 1;
			genos[i] = types[i][0] + types[i][1] + 1; 
			genoCount[genos[i]]++;
		}
		double n = (double) alleleCount[N];
		alleleFreq = new double[5];
		alleleFreq[N] = n / len * 2; // 缺失率
		alleleFreq[A] = alleleCount[A] / (len * 2 - n); // A 的比例
		alleleFreq[B] = 1 - alleleFreq[A];

		genoFreq = new double[4];
		n = (double) genoCount[NN];
		genoFreq[0] = n / len; // 缺失率
		double sum = len - n;
		for (i = 1; i < 4; i++) {
			genoFreq[i] = genoCount[i] / sum;
		}
		estimateNN();
	}
	
	private void estimateNN(){
		// TODO: 更高档的算法
		Random r = new Random();
		double p = this.alleleFreq[A];
		for(int ts[]: types){
			if(ts[0] == -1){
				ts[0] = r.nextDouble() < p ? (0) : (1);
				ts[1] = r.nextDouble() < p ? (0) : (1);
			}
		}
	}

	private double calHW() {
		int homc = Math.max(genoCount[AA],
				genoCount[BB]); // AA/BB中大的一个
		int homr = Math.min(genoCount[AA],
				genoCount[BB]); // AA/BB中小的一个
		int hets = genoCount[AB]; // AB
		// TODO: 转移
		return hwe(hets, homc, homr);
	}

	private double hwe(int hets, int homc, int homr) {
		int rares = 2 * homr + hets; // 稀有基因的出现次数
		int diplo = hets + homc + homr; // 所有

		double[] hetProbs = new double[rares + 1];
		int mid = (rares * (2 * diplo - rares) / (2 * diplo));

		if (((rares & 1) ^ (mid & 1)) != 0) {
			mid++;
		}

		int currHets = mid;
		int currHomr = (rares - mid) / 2;
		int currHomc = diplo - currHets - currHomr;

		hetProbs[mid] = 1.0;
		double sum = 1.0;
		for (currHets = mid; currHets > 1; currHets -= 2) {
			hetProbs[currHets - 2] = hetProbs[currHets] * currHets
					* (currHets - 1.0)
					/ (4.0 * (currHomr + 1.0) * (currHomc + 1.0));
			sum += hetProbs[currHets - 2];
			currHomr++;
			currHomc++;
		}

		currHets = mid;
		currHomr = (rares - mid) / 2;
		currHomc = diplo - currHets - currHomr;

		for (currHets = mid; currHets <= rares - 2; currHets += 2) {
			hetProbs[currHets + 2] = hetProbs[currHets] * 4.0 * currHomr
					* currHomc / ((currHets + 2.0) * (currHets + 1.0));
			sum += hetProbs[currHets + 2];
			currHomr--;
			currHomc--;
		}

		int i;
		double pvalue = 0.0;
		for (i = 0; i <= rares; i++) {
			hetProbs[i] /= sum;
		}
		for (i = 0; i <= rares; i++) {
			if (hetProbs[i] > hetProbs[hets])
				continue;
			pvalue += hetProbs[i];
		}

		if (pvalue > 1.0) {
			return 1.0;
		} else {
			return pvalue;
		}

	}

	private boolean judgeBad(Setting rc) {
		if (genoFreq[NN] > rc.getNNRatio() || maf < rc.getMAF()
				|| hw < rc.getHWE()) {
			return true;
		} else {
			return false;
		}
	}
	
	private void setCC(int[] status, double[] excepted){
		int a, b, c, d, i, l = status.length;
		a = b = c = d = 0;
		for(i = 0; i < l; i++){
			if(status[i] == -1 || genos[i] == 0){
				continue;
			}
			if(status[i] == 0){          // control
				switch(genos[i]){
				case 1: c += 2; break;   // AA
				case 2: c++; d++; break; // AB
				case 3: d += 2; break;
				default: break;
				}
			}else{                       // case
				switch(genos[i]){
				case 1: a += 2; break;   // AA
				case 2: a++; b++; break; // AB
				case 3: b += 2; break;
				default: break;
				}
			}
		}
		/*
		 * -----------------
		 *        | A  | B
		 * -----------------
		 *  case  | a  | b
		 * -----------------
		 * control| c  | d
		 * -----------------
		 */

		OR orv = new OR();
		double[] ors = orv.getOR_CI(a, b, c, d);
		or = ors[0];
		orLow = ors[1];
		orHigh = ors[2];
		
		FitTest fit = new FitTest(excepted);
		double[] observed = new double[2];
		observed[0] = a;
		observed[1] = c;
		this.aPValue = fit.test(observed);
		observed[0] = b;
		observed[1] = d;
		this.bPValue = fit.test(observed);
	}

}
