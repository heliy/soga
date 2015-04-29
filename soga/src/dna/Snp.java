package dna;

import dna.LD;
import parameter.Setting;
import parameter.Base;
import parameter.GenoType;

import java.util.HashMap;

public class Snp {
	String rs; // rs + number
	int A;
	int B;
	String chr; // chromesome
	long position; // position in chromesome
	Base base;
	GenoType genotype;
	int genos[]; // genotypes data
	int types[][]; // 与AB对应的基因型 用来计算LD
	@SuppressWarnings("rawtypes")
	HashMap pws; // pairwise to other snps
	int alleleCount[]; // count of alleles
	double alleleFreq[]; // frequency of alleles
	double maf; // min allele frequency
	int genoCount[]; // count of genotypes of genos
	double genoFreq[]; // frequency of genotypes
	double hw; // hardy-weinberg p value
	boolean isBad; // is bad snp

	public Snp(String rsNo, String Chr, String sposition, String alleles,
			String[] parts, Setting rc) {
		// rs, alleles, chr, position, genos....
		rs = rsNo;
		parseAlleles(alleles); // 等位的碱基 [A, B, N]
		position = Long.parseLong(sposition); // 所在位置
		chr = Chr;
		base = new Base();
		genotype = new GenoType(A, B);
		typeCount(parts);
		maf = Math.min(alleleFreq[A], alleleFreq[B]); // 最小等位频率
		hw = calHW(); // HW平衡p值
		isBad = judgeBad(rc); // 是否通过QC
		pws = new HashMap<Snp, LD>(); // 与其他snp的LD值

	}

	public void display() {
		// TODO
		System.out.println(rs + "\t" + chr + "\t" + position + "\t"
				+ base.getBase(A) + "\t" + base.getBase(B) + "\t["
				+ alleleCount[A] + ", " + alleleCount[B] + ", "
				+ alleleCount[base.N()] + "]\n");

		System.out.println("frequency of alleles: " + alleleFreq[0] + ", "
				+ alleleFreq[3] + ", " + alleleFreq[4]);
		System.out.println("counts of genotypes: " + genoCount[0] + ", "
				+ genoCount[1] + ", " + genoCount[2] + ", " + genoCount[3]);
		System.out.println("frequency of genotypes: " + genoFreq[0] + ", "
				+ genoFreq[1] + ", " + genoFreq[2] + ", " + genoFreq[3]);
		System.out.println("MAF:" + maf);
		System.out.println("H-W:" + hw);
		if (isBad) {
			System.out.println("Bad");
		} else {
			System.out.println("OK");
		}
	}
	
	public int[][] getTypes(){
		return types;
	}

	private void parseAlleles(String ts) {
		// parse string to two alleles
		// A/C -> A C
		// TODO: error control and insert deletion
		Base base = new Base();
		A = base.baseNo(ts.charAt(0));
		B = base.baseNo(ts.charAt(2));
		if (A > B) { // 保证 A 的序号小于 B
			int tmp = B;
			B = A;
			A = tmp;
		}
	}

	private void typeCount(String[] parts) {
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
			if (first > second) { // 保证 first 的序号小于 second
				int tmp = second;
				second = first;
				first = tmp;
			}
			alleleCount[first]++;
			alleleCount[second]++;
			int type = genotype.getNo(first, second);
			genos[i] = type;
			types[i] = new int[2];
			if(first == A){
				types[i][0] = 0;
			}else if(first == B){
				types[i][0] = 1;
			}else{
				types[i][0] = -1;
			}
			if(second == A){
				types[i][1] = 0;
			}else if (second == B){
				types[i][1] = 1;
			}else{
				types[i][1] = -1;
			}
			
			genoCount[type]++;
		}
		double n = (double) alleleCount[base.N()] * 2;
		alleleFreq = new double[5];
		alleleFreq[base.N()] = n / len * 2; // 缺失率
		alleleFreq[A] = alleleCount[A] / (len * 2 - n); // A 的比例
		alleleFreq[B] = 1 - alleleFreq[A];

		genoFreq = new double[4];
		n = (double) genoCount[0];
		genoFreq[0] = n / len; // 缺失率
		double sum = len - n;
		for (i = 1; i < 4; i++) {
			genoFreq[i] = genoCount[i] / sum;
		}
	}

	private double calHW() {
		int homc = Math.max(genoCount[genotype.getNo(A, A)],
				genoCount[genotype.getNo(B, B)]); // AA/BB中大的一个
		int homr = Math.min(genoCount[genotype.getNo(A, A)],
				genoCount[genotype.getNo(B, B)]); // AA/BB中小的一个
		int hets = genoCount[genotype.getNo(A, B)]; // AB
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
		if (alleleFreq[base.N()] > rc.MAXN || maf < rc.MAF || hw < rc.HWE) {
			return true;
		} else {
			return false;
		}
	}

}
