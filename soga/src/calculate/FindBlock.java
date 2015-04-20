package calculate;

import java.util.ArrayList;
import java.util.Arrays;
//import java.util.concurrent.CountDownLatch;

import dna.Block;
import dna.LD;
import dna.Scope;
import dna.Snp;
import parameter.Setting;

public class FindBlock implements Runnable {
	private double TOLERENCE = 1e-8;
	private int SAMPLES;

	private int win;  // migpp算法的小窗口的大小
	private double d;
	private int last;
	
	private double[] w;
	private int TARGET;
	
	private PairwiseCal cal;
	
	private Snp[] snps;
	private boolean isLast;
	private Block[] blocks = null;
	private boolean isBlock;
	
//	private CountDownLatch latch;
	
	public FindBlock(Setting rc, boolean isblock){
		
		isBlock = isblock;
		SAMPLES = rc.getSAMPLES();
		d = rc.getD();
		cal = new PairwiseCal(rc);
		last = 1;
		w = new double[3];      
		if(isBlock){                // 查找单体型块
			w[0] = 0;               // otherwise
			w[1] = 1 - d;           // strong LD
			w[2] = -d;              // EHR
			TARGET = 1;
		}else{                      // 查找重组块
			w[0] = 0;               // otherwise
			w[1] = -d;              // strong LD
			w[2] = 1 - d;           // EHR
			TARGET = 2;
		}
	}

	public Block[] mig(Snp[] snps, boolean isLast) {
	
		int type, i, j, n = snps.length;
		double s;
		double[] W = new double[n];
		ArrayList<Scope> scopes = new ArrayList<Scope>();

		for (j = 1; j < n; j++) {
			s = 0;
			for (i = j - 1; i >= 0; i--) {
				LD l = snps[j].getLD(cal, snps[i]);
				type = l.getType();
				s += w[type];
				W[i] += s;
				if ((type == TARGET) && W[i] >= 0) {
					scopes.add(new Scope(i,j));
				}
			}
		}
		return getBlocks(scopes, snps);
	}

	public Block[] migpp(Snp[] snpS, boolean islast) {
		snps = snpS;
		isLast = islast;
		
		int n, newWin, calculations, b, newB, type;
		int i, j;
		int[] T, B;
		double wMax, tmp;
		double[] W, S;
		ArrayList<Scope> scopes = new ArrayList<Scope>();

		n = snps.length;
		newWin = 0;
//		win = ceil((n-1)*(1-d)/2)
		tmp = (n - 1) * (1 - d) / 2;
		if (Math.abs((int) tmp - tmp) <= TOLERENCE) {
			win = (int) tmp;
		} else {
			win = (int) tmp + 1;
		}
//System.out.println(win);
		// W S T B H 初始化
		W = new double[n];
		S = new double[n];
//		wmaxs = new double[n];
		T = new int[n];
		B = new int[n];
		for (i = 0; i < n; i++) {
			T[i] = i; // 伪代码中T/ B从2至n 并且后续均为T[j-1] / B[j-1] <- 下标从1开始
			B[i] = i; // 改为从1到n 后续为T[j] / B[j]
//			System.out.println(i+":" + T[i]+", "+B[i]);
		}
		
//		wmaxs[0] = 0;

		calculations = 1;
		while (calculations > 0) {
			newWin += win;
			calculations = 0;
			b = 0;
			newB = 0;
			for (j = last; j < n; j++) { // 这里从上次计算的结尾开始
				if (newB == B[j]) {
					B[j] = b;
					b = T[j];
					newB = b;
					continue;
				}
				if (j - newB > newWin) {
					B[j] = j - newWin;
					b = B[j];
				} else {
					B[j] = b;
					b = newB;
				}

				newB = T[j];
//				System.out.println(newB+", "+b);
				for (i = T[j]-1; i >= b; i--) {
					type = snps[i].getLD(cal, snps[j]).getType();
//					System.out.println(w[type]);
					S[j] += w[type];
					W[i] += S[j];
//					wmaxs[j] = wmaxs[i]+W[i];
					if ((type == TARGET) && (W[i] >= 0)) {
						
//    					LD ld = snps[i].getLD(cal, snps[j]);
//    					System.out.println(W[i]+", "+ld);
//						if(ld.getCiLow() < 0.1 ||)
						scopes.add(new Scope(i, j));
//						int snpsn = j-i+1;
//						System.out.println(i+", "+j+", "+TARGET);
//						for(int a = i; a <= j; a++){
//							for(int c = i; c <= j; c++){
//									System.out.print(snps[a].getLD(cal, snps[c]).getType()+" ");
//							}
//							System.out.println();
//						}
//						System.out.println();
					}
					wMax = W[i]+w[1]*((j-i+1)+(n-i))*(n-j)/2;
//					if(isBlock){
//					wMax = W[i] + 0.5*d*(n-i)*(n+i-2*j);
//					}else{
////					System.out.println(i+": "+W[i]);
//					    wMax = W[i]+w[snps[i].getLD(cal, snps[j]).getType()];
//					    int k = j+1;
//					    while (wMax < 0 && k < n) {
////						System.out.println(i+", "+k+": "+snps[i].hasLD(snps[k]));
//						type = snps[i].getLD(cal, snps[k]).getType();
//						wMax += w[type];
//						k++;
//					    }
//					}
					if(wMax >= 0){
						newB = i;
					}
					calculations++;
				}
				T[j] = b;
			}
		}
		return getBlocks(scopes, snps);
	}

	private Block[] getBlocks(ArrayList<Scope> scopelist, Snp[] snps) {
		int[] H = new int[snps.length];
		int i;
		ArrayList<Block> blocklist = new ArrayList<Block>();
		
		Scope[] scopes = new Scope[scopelist.size()];
		scopelist.toArray(scopes);
		Arrays.sort(scopes);
		
		for(Scope s:scopes){
			if((H[s.getBegin()] == 0) && (H[s.getEnd()] == 0)){
				if(!isLast && s.getEnd() == (snps.length - 1)){
					continue;
				}
				blocklist.add(new Block(s.getBegin(), s.getEnd(), snps));
				for(i = s.getBegin(); i <= s.getEnd(); i++){
					H[i] = 1;
				}
			}
		}
		
		Block[] blocks = new Block[blocklist.size()];
		blocklist.toArray(blocks);
		Arrays.sort(blocks);
		return blocks;
	}

	public void setLast(int last) {
		this.last = last;
	}
	
	public int getWin(){
		return win;
	}
	
	public void setRun(Snp[] snpS, int begin, int end, boolean islast){
		snps = new Snp[end-begin];
		int i;
		for(i = begin; i < end; i++){
			snps[i-begin] = snpS[i];
		}
		isLast = islast;
	}

	@Override
	public void run() {
//		if(isBlock){
		blocks = migpp(snps, isLast);
//		}else{
//			blocks = mig(snps, isLast);
//		}
//		latch.countDown();
	}

	public Block[] getBlocks(){
		return blocks;
	}
}
