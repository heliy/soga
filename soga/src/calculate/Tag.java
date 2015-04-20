package calculate;

import java.util.ArrayList;
import java.util.Iterator;

import dna.LD;
import dna.Snp;
import parameter.Setting;


public class Tag {
	private static double threshold = 0.9;
	private int maxNum;
//	private double weight;
	private double R[][];
	private boolean N[][];
	private PairwiseCal cal;
	
	public Tag(Setting rc){
		maxNum = rc.getSIZE() + 2;
		cal = new PairwiseCal(rc);
		R = new double[maxNum][maxNum];
		N = new boolean[maxNum][maxNum];
//		weight = 0.05;
	}
	
	public boolean[] getTags(Snp[] snps){
		int num = snps.length;
//		double threshold = 1 - weight*Math.exp(num*weight);
//		System.out.println(num+", "+threshold);
		boolean[] isTags = new boolean[num];   // all false
		
		ArrayList<Integer> tests = new ArrayList<Integer>();
		ArrayList<Integer> remains = new ArrayList<Integer>();
		ArrayList<Double> mris = new ArrayList<Double>();
		
		int i, j, m, n;
		double min, tmp;
		for(i = 0; i < num; i++){
			remains.add(i);
			for(j = 0; j < num; j++){
				if(i == j){
					N[i][j] = true;
					R[i][j] = 1;
					continue;
				}
				LD ld = snps[i].getLD(cal, snps[j]);
				R[i][j] = ld.getRsq();
				if(R[i][j] > threshold){
					N[i][j] = true;
				}else{
					N[i][j] = false;
				}
			}
		}
				
		while(remains.size() != 0){
			tests.clear();
			tests.add(remains.get(0));
			remains.remove(0);
			for(i = 0; i < num; i++){
				for(j = 0; j < num; j++){
					m = index(tests, i);
					n = index(remains, j);
					if((m > -1)&&(n > -1)&&(N[i][j])){
						tests.add(j);
						remains.remove(n);
					}
				}
			}
			
			min = num;
			m = 0;
			Iterator<Integer> iter = tests.iterator();
			while(iter.hasNext()){
				i = iter.next();
				tmp = sum(R[i], num);
				mris.add(tmp);
				if(tmp < min){
					min = tmp;
					m = i;
				}
			}
			isTags[m] = true;
		}
		
		return isTags;
	}
	
	private Double sum(double[] ds, int num) {
		double s = 0;
		int i;
		for(i = 0; i < num; i++){
			s += ds[i];
		}
		return s;
	}

	private int index(ArrayList<Integer> al, int n){
		Iterator<Integer> iter = al.iterator();
		int i = 0;
		while(iter.hasNext()){
			if(iter.next().equals(n)){
				return i;
			}
			i++;
		}
		return -1;
	}

}
