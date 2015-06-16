package util;

import java.util.Random;

import parameter.Setting;
import dna.Snp;

public class SelectH {
//	private ArrayList<Snp> snplist;
	private int K;
	private int N;
	private int L;
//	private ArrayList<int[]> H;
//	private ArrayList<int[]> genotypes;
//	private ArrayList<Integer> mafs;
	
	private int[][] genotypes;
	private double[] mafs;
	private int[][] H;
	
	// 本来想着添加SNP然后快速算出候选的单体型分布的
	
	// 现在是简化版
	// 没有添加这一选项
	// 构造 -> 得背景分布
	
	public SelectH(Setting rc, Snp[] snps){
		K = rc.getK();
		N = rc.getSAMPLES();
		L = snps.length;
		genotypes = new int[N][L];
		mafs = new double[L];
		
		H = new int[K][L];
		
		int i, j;
		for(i = 0; i < L; i++){
			mafs[i] = snps[i].getMaf();
		}
		for (i = 0; i < L; i++) {
			int[][] types = snps[i].getTypes();
			for(j = 0; j < N; j++){
    			genotypes[j][i] = types[j][0]+types[j][1];
			}
		}
//		for(i = 0;i < N; i++){
//			for(j = 0; j < L; j++){
//    			System.out.print(genotypes[i][j]+", ");				
//			}
//			System.out.println();
//		}
//		System.out.println();
		int[] pows, hetes;
		pows = new int[N]; 
		hetes = new int[N];
		
		for (i = 0; i < N; i++) {
			hetes[i] = 0;
			for (j = 0; j < L; j++) {
				hetes[i] += (genotypes[i][j] == 1) ? (1) : (0);
			}
			// System.out.println(i+": "+hetes[i]);
			pows[i] = ((int) Math.pow(2, hetes[i])) % K;
		}

		H = new int[K][L];

		Random r = new Random();
		int hNum = 0;
		for (i = 0; i < N; i++) {
			double selectp = Math.max(1.0 / K, K / (N * pows[i]));
			for (j = 0; j < pows[i]; j++) {
				if (hNum == K) {
					break;
				}
				if (r.nextDouble() > selectp) {
					continue;
				}
				getCandidate(i, hNum);
				hNum++;
			}
		}
		while (hNum < K) {
			i = r.nextInt(N);
			getCandidate(i, hNum);
			hNum++;
		}	
	}

	private void getCandidate(int index, int h) {
		int[] genotype = genotypes[index];
		int i, L = genotype.length;
		Random r = new Random();
		for (i = 0; i < L; i++) {
			if (genotype[i] == 1) {
				if (r.nextDouble() < this.mafs[i]) {
					H[h][i] = 1;
				} else {
					H[h][i] = 0;
				}
			} else {
				H[h][i] = genotype[i] / 2;
			}
		}
	}
	
	public int[][] getH(){
		return H;
	}
	
//	public SelectH(Setting rc){
//		snplist = new ArrayList<Snp>();
//		K = rc.getK();
//		N = rc.getSAMPLES();
//		H = new ArrayList<int[]>();
//		genotypes = new ArrayList<int[]>();
//		mafs = new ArrayList<Integer>();
//	}
	
//	public int[][] select(int[][] genotypes){
//		if(this.H == null){
//			
//		}
//	}
//	
//	private void setH(){
//		int i, j, sum, hNum;
//		int L = this.snplist.size();
//		int[] pows, hetes;
//		pows = new int[N]; 
//		hetes = new int[N];
//		double selectp;
//		
//		sum = 0;
//		for (i = 0; i < N; i++) {
//			hetes[i] = 0;
//			for (j = 0; j < L; j++) {
//				hetes[i] += (genotypes[i][j] == 1) ? (1) : (0);
//			}
//			// System.out.println(i+": "+hetes[i]);
//			pows[i] = ((int) Math.pow(2, hetes[i])) % K;
//			sum += pows[i];
//		}
//		K = (sum < K) ? (sum) : (K);
//
//		H = new int[K][L];
//
//		Random r = new Random();
//		hNum = 0;
//		for (i = 0; i < N; i++) {
//			selectp = Math.max(1.0 / K, K / (N * pows[i]));
//			for (j = 0; j < pows[i]; j++) {
//				if (hNum == K) {
//					break;
//				}
//				if (r.nextDouble() > selectp) {
//					continue;
//				}
//				getCandidate(i, hNum);
//				hNum++;
//			}
//		}
//		while (hNum < K) {
//			i = r.nextInt(N);
//			getCandidate(i, hNum);
//			hNum++;
//		}		
//	}
//
//	private void getCandidate(int index, int h) {
//		int[] genotype = genotypes.get(index);
//		int i, L = genotype.length;
//		double p = this.mafs.get(index);
//		Random r = new Random();
//		for (i = 0; i < L; i++) {
//			if (genotype[i] == 1) {
//				if (r.nextDouble() < p) {
//					H[h][i] = 1;
//				} else {
//					H[h][i] = 0;
//				}
//			} else {
//				H[h][i] = genotype[i] / 2;
//			}
//		}
//	}
	

}
