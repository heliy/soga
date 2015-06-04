package util;

import java.util.ArrayList;
import java.util.Random;

import parameter.Setting;
import dna.Snp;

public class SelectH {
	private ArrayList<Snp> snplist;
	private int K;
	private int N;
	private ArrayList<int[]> H;
	private ArrayList<int[]> genotypes;
	private ArrayList<Integer> mafs;
	
	public SelectH(Setting rc){
		snplist = new ArrayList<Snp>();
		K = rc.getK();
		N = rc.getSAMPLES();
		H = new ArrayList<int[]>();
		genotypes = new ArrayList<int[]>();
		mafs = new ArrayList<Integer>();
	}
	
//	public int[][] select(Snp[] snps){
//		if(this.H == null){
//			
//		}
//	}
	
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
//	
	private void clear(){
		this.snplist.clear();
		this.H.clear();
		this.genotypes.clear();
		this.mafs.clear();
	}

}
