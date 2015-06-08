package calculate;

import java.util.Random;
import java.util.concurrent.CountDownLatch;

import parameter.Setting;
import dna.HaploType;
import dna.Snp;

public class Phase implements Runnable {
	private HaploType[] hts; // 所有单体型

	static private int B; // B个一段
	static private int T; // T = 2**B
	static private double LIMIT = 0.0001;// 最小概率
	static private int Ne = 15000; // 群体个数

	private int L; // snp个数
	private int N; // 样本个数
	private int K; // 考虑的备选单体型数
	private static int bits[][] = null; // 段内单体型的0/1
	private int hete; // 样本的AB型个数
	private boolean ishete[]; // 内是否是AB型

	private int genotypes[][]; // 基因型
	private int H[][] = null; // 备选单体型

	private double lamb; // lambda
	private double thetas[]; // e^(-ro/K)
	private double unthetas[]; // (1-e^(-ro/K))/K
	private double mafs[]; // 最小等位频率

	private double PZ1; // P(Z1 = u) = 1/K

	CountDownLatch latch;
	private int endType;
	private Snp[] endSnps;
	private Snp[] snps;
	private Snp[] next;
	
	private int sample;
//	private int maxDistance;

	public void setRun(Snp[] snps, CountDownLatch latch) {
		this.next = snps;
		this.latch = latch;
	}
	
	public void clearHistory(){
		this.endType = -1;
		this.endSnps = new Snp[0];
		this.H = null;
	}
	
	public Phase(int sample, int[][] H, Setting rc){
		this.construct(sample, rc);
		this.H = H;
	}
	
	private void construct(int sample, Setting rc) {
		B = rc.getB();
		T = (int) Math.pow(2, B);
		N = rc.getSAMPLES();
		K = rc.getK();
		
		this.sample = sample;		
		this.clearHistory();
	}

	public Phase(int sample, Setting rc) {
		this.construct(sample, rc);
	}
	
	public boolean hasCleared(){
		return this.endSnps.length == 0;
	}
	
	private void addSnps(Snp[] snps){
		int i, j, n;
		this.snps = new Snp[snps.length+this.endSnps.length];
		L = this.snps.length;
		for(i = 0, n = this.endSnps.length; i < n; i++){
			this.snps[i] = this.endSnps[i];
		}
		for(j = 0, n = snps.length; j < n; j++){
			this.snps[i+j] = snps[j];
		}		
	}
	
	public void setH(int[][] H){
		// TODO
		this.H = H;
	}
	
	public HaploType[] phase(Snp[] snps2, int[][] h){
		this.H = h;
		return this.phase(snps2);
	}
	
	public HaploType[] phase(Snp[] snps2){

		HaploType[] hts = new HaploType[2];
		int i, j, types[][];
		hete = i = 0;
		this.addSnps(snps2);
		this.genotypes = new int[N][L];
		this.ishete = new boolean[L];
		
//		System.out.println(this.sample+", "+this.snps.length);
		// 基因型
		for (i = 0; i < L; i++) {
			types = this.snps[i].getTypes();
			for(j = 0; j < N; j++){
    			genotypes[j][i] = types[sample][0]+types[sample][1];
			}
			ishete[i] = (genotypes[sample][i] == 1);
			hete += (ishete[i]?1:0);
		}

		if(hete < 2){
			j = snps.length;
			int[][] diplos = new int[2][j];
			int k = endSnps.length;
			for (i = 0; i < j; i++) {
				if (genotypes[sample][i+k] == 1) {
					diplos[0][i] = snps[i].getA();
					diplos[1][i] = snps[i].getB();
				} else if (genotypes[sample][i+k] == 0) {
					diplos[1][i] = diplos[0][i] = snps[i].getA();
				} else {
					diplos[1][i] = diplos[0][i] = snps[i].getB();
				}
			}
			hts[0] = new HaploType(diplos[0], sample);
			hts[1] = new HaploType(diplos[1], sample);

			this.clearHistory();
			return hts;			
		}
		
		
		int sum, hNum;
		double selectp, distance;
		lamb = getLambda();
		thetas = new double[L];
		unthetas = new double[L];
		mafs = new double[L];

		// System.out.println(L + ", " + N + ", " + K);

		GeneticDistance gd = new GeneticDistance();
		mafs[0] = this.snps[0].getMaf();
		for (i = 1; i < L; i++) {
			distance = gd.ditance(this.snps[i - 1], this.snps[i]);
			thetas[i] = Math.exp((-4) * Ne * distance / (double) K);
			unthetas[i] = (1 - thetas[i]) / K;
			mafs[i] = this.snps[i].getMaf();
		}

		if(H == null){
			int[] pows, hetes;
			pows = new int[N]; 
			hetes = new int[N];
			
			sum = 0;
			for (i = 0; i < N; i++) {
				hetes[i] = 0;
				for (j = 0; j < L; j++) {
					hetes[i] += (genotypes[i][j] == 1) ? (1) : (0);
				}
				// System.out.println(i+": "+hetes[i]);
				pows[i] = ((int) Math.pow(2, hetes[i])) % K;
				sum += pows[i];
			}
	
			H = new int[K][L];
	
			Random r = new Random();
			hNum = 0;
			for (i = 0; i < N; i++) {
				selectp = Math.max(1.0 / K, K / (N * pows[i]));
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

		if(bits == null){
			bits = new int[T][B];
			for (i = 0; i < T; i++) {
				sum = i;
				// i二进制化
				// eg i=5时为 [1, 0, 1]
				for (j = 0; j < B; j++) {
					bits[i][j] = sum % 2;
					sum /= 2;
				}
			}
		}
		PZ1 = 1.0 / K;
		
		return this.shapeit();
	}

	public void run() {
		this.hts = phase(this.next);
		this.latch.countDown();
		return;
	}

	public HaploType[] getHts() {
		return hts;
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


	private HaploType[] shapeit() {
		// 计算用
		int S[] = new int[L]; // 分段
		int A[][] = new int[K][L]; // 断内可能单体型
		double alpha[][][] = new double[L][T][K]; // 前向概率
		double beta[][][] = new double[L][T][K]; // 后向概率
		double sumsOneAll[] = new double[T]; // 行和
		double sumsAllOne[] = new double[K]; // 列和
		double inturn[][] = new double[T][T]; // 内部转移矩阵
		double outturn[][] = new double[T][T]; // 外部转移矩阵
		double tmpturn[][] = new double[T][T]; // 临时转移矩阵

		int[] genotype = genotypes[this.sample];
		int i, j, k, u, u1, u2;
		int inNo, segNo, m, total, maxi, maei, maej;
		double nextalpha, nextbeta, pmax, tmp, tmpsum;

		// 分段
		inNo = segNo = 0;
		for (i = 0; i < L; i++) {
			if (genotype[i] == 1) {
				inNo++;
			}
			S[i] = segNo;
			// System.out.println(i + ": " + S[i]);
			if (inNo == B) {
				segNo++;
				inNo = 0;
			}
		}

		// 断内可能单体型
		for (i = 0; i < T; i++) {
			for (k = j = 0; j < L; j++) {
				if (genotype[j] == 1) {
					A[i][j] = bits[i][k % B];
					k++;
				} else {
					A[i][j] = genotype[j] / 2;
				}
			}
		}

		// 前向
		m = total = 0;
		for (i = 0; i < T; i++) {
			for (u = 0; u < K; u++) {
				alpha[m][i][u] = limit(PZ1 * pXlZl(A, m, i, u));
				total += alpha[m][i][u];
				sumsOneAll[i] += alpha[m][i][u];
				sumsAllOne[u] += alpha[m][i][u];
			}
		}
		for (m = 1; m < L; m++) {
			for (i = 0; i < T; i++) {
				for (u = 0; u < K; u++) {
					nextalpha = pXlZl(A, m, i, u);
					if (S[m - 1] == S[m]) {
						nextalpha *= (sumsOneAll[i] * unthetas[m] + thetas[m]
								* alpha[m - 1][i][u]);
					} else {
						nextalpha *= (total * unthetas[m] + thetas[m]
								* sumsAllOne[u]);
					}
					alpha[m][i][u] = limit(nextalpha);
				}
			}
			clear(sumsOneAll, T);
			clear(sumsAllOne, K);
			total = 0;
			for (i = 0; i < T; i++) {
				for (u = 0; u < K; u++) {
					total += alpha[m][i][u];
					sumsOneAll[i] += alpha[m][i][u];
					sumsAllOne[u] += alpha[m][i][u];
				}
			}
		}

		// 后向
		clear(sumsOneAll, T);
		clear(sumsAllOne, K);
		total = 0;
		m = L - 1;
		for (i = 0; i < T; i++) {
			for (u = 0; u < K; u++) {
				beta[m][i][u] = 1;
				nextbeta = pXlZl(A, m, i, u);
				sumsOneAll[i] += nextbeta;
				sumsAllOne[u] += nextbeta;
				total += nextbeta;
			}
		}
		for (m = L - 2; m >= 0; m--) {
			for (i = 0; i < T; i++) {
				for (u = 0; u < K; u++) {
					if (S[m] == S[m + 1]) {
						nextbeta = unthetas[m + 1] * sumsOneAll[i]
								+ thetas[m + 1] * beta[m + 1][i][u]
								* pXlZl(A, m + 1, i, u);
					} else {
						nextbeta = unthetas[m + 1] * total + thetas[m + 1]
								* sumsAllOne[u];
					}
					beta[m][i][u] = limit(nextbeta);
				}
			}
			clear(sumsOneAll, T);
			clear(sumsAllOne, K);
			total = 0;
			for (i = 0; i < T; i++) {
				for (u = 0; u < K; u++) {
					beta[m][i][u] = 1;
					nextbeta = pXlZl(A, m, i, u) * beta[m][i][u];
					sumsOneAll[i] += nextbeta;
					sumsAllOne[u] += nextbeta;
					total += nextbeta;
				}
			}
		}

		// 概率
		int[] X1 = new int[segNo + 1];
		int[] X2 = new int[segNo + 1];

		m = maxi = 0;
		pmax = 0.0;
		if (this.endType < 0) {
			for (i = 0; i < T / 2; i++) {
				j = getCp(i);
				for (tmp = u = 0; u < K; u++) {
					tmp += alpha[m][i][u] * beta[m][i][u];
				}
				nextalpha = tmp;
				for (tmp = u = 0; u < K; u++) {
					tmp += alpha[m][j][u] * beta[m][i][u];
				}
				tmp *= nextalpha;
				if (tmp > pmax) {
					pmax = tmp;
					maxi = i;
				}
			}
			X1[0] = maxi;
			X2[0] = getCp(maxi);
		}else{
			X1[0] = this.endType;
			X2[0] = getCp(this.endType);
		}
		// System.out.println(S[L-1]);
		if (S[L - 1] > 0) {
			while (true) {
				m++;
				// System.out.println(m);
				// 计算前一个位点到这个位点的转移概率
				for (i = 0; i < T; i++) {
					for (j = 0; j < T; j++) {
						tmpsum = 0;
						for (u1 = 0; u1 < K; u1++) {
							for (u2 = 0; u2 < K; u2++) {
								tmpsum += (alpha[m - 1][i][u1]
										* pZlZlm1(m, u2, u1)
										* pXlZl(A, m, j, u2) * beta[m][j][u2]);
							}
						}
						inturn[i][j] = tmpsum;
					}
				}
				// 内部转移与外部转移合并
				for (i = 0; i < T; i++) {
					tmpsum = 0;
					for (j = 0; j < T; j++) {
						tmp = 0;
						for (k = 0; k < T; k++) {
							tmp += outturn[i][k] * inturn[k][j];
						}
						tmpturn[i][j] = tmp;
						tmpsum += tmp;
					}
					for (j = 0; j < T; j++) {
						outturn[i][j] = tmpturn[i][j] / tmpsum;
					}
				}
				// 是否跨段
				if (S[m] != S[m - 1]) {
					// 跨段 算单体型概率 最大加入
					maei = X1[S[m - 1]];
					maej = X2[S[m - 1]];
					pmax = 0;
					for (i = 0; i < T; i++) {
						j = getCp(i);
						tmp = outturn[maei][i] * outturn[maej][j];
						if (tmp > pmax) {
							pmax = tmp;
							maxi = i;
						}
					}
					X1[S[m]] = maxi;
					X2[S[m - 1]] = getCp(maxi);
					if (S[m] == S[L - 1]) {
						break;
					}
					for (i = 0; i < T; i++) {
						for (j = 0; j < T; j++) {
							outturn[i][j] = 0;
						}
						outturn[i][i] = 1;
					}
				}
			}
		}

		HaploType[] hts = new HaploType[2];
		int e = this.endSnps.length;
		int[] ht1 = new int[L-e];
		int[] ht2 = new int[L-e];
		int a, b;
		for (i = e; i < L; i++) {
			a = snps[i].getA();
			b = snps[i].getB();
			ht1[i-e] = (A[X1[S[i]]][i] == 0) ? (a) : (b);
			ht2[i-e] = (A[X2[S[i]]][i] == 0) ? (a) : (b);
		}
		
		this.setEnd(ht1, ht2);
		
		hts[0] = new HaploType(ht1, this.sample);
		hts[1] = new HaploType(ht2, this.sample);

		return hts;
		
	}
	
	private void setEnd(int[] ht1, int[] ht2){
		int i, j, e = this.endSnps.length;
	    int[] hetes = new int[B];
	    this.endType = 0;
	    for(i = L-1, j = 0; i >= e && j < B; i--){
	    	if(ht1[i-e] != ht2[i-e]){
	    		hetes[j] = i;
	    		if(this.snps[i].getA() == ht1[i-e]){
    	    		this.endType += (int)Math.pow(2, j);
	    		}
	    		j++;
	    	}
	    }
	    if(i < e){
	    	i = e;
	    }
		endSnps = new Snp[L-i];
		
//		for(i = 0, a = ht1.length; i < a; i++){
//			System.out.println(ht1[i]+" | " +ht2[i]);
//		}
//		System.out.println();
		
		for(j = i; i < L; i++){
//			System.out.println(e+": "+i);
			endSnps[i-j] = this.snps[i];
		}		
	}

	private int getCp(int i) {
		return T - i - 1;
	}

	private void clear(double[] array, int n) {
		int i;
		for (i = 0; i < n; i++) {
			array[i] = 0;
		}
	}

	// P(Zl = u| Zl-1 = v)
	// P(Zsnp = ht1| Zsnp-1 = ht2)
	private double pZlZlm1(int snpIndex, int ht1, int ht2) {
		// System.out.println(snpIndex);
		double theta = thetas[snpIndex];
		double untheta = unthetas[snpIndex];
		if (ht1 == ht2) {
			return theta + untheta;
		} else {
			return untheta;
		}
	}

	// P(Xl = i| Zl = u,H)
	// P(Xsnp = ht1| Zsnp = ht2, H)
	private double pXlZl(int[][] A, int snp, int ht1, int ht2) {
		return (A[ht1][snp] == H[ht2][snp]) ? (1 - lamb) : (lamb);
	}

	// p(Xl = i, Zl = u| Xl-1 = j, Zl-1 = v, H)
	// p(Xsnp = nexthtx, Zsnp = nexthtz| Xsnp-1 = prevhtx, Zsnp-1 = prevhtz)
	// private double pXzlXzlm1(int snp, int nexthtx, int nexthtz, int prevhtx,
	// int prevhtz){
	// if(S[snp] == S[snp-1] && nexthtx != prevhtx){
	// return 0;
	// }else{
	// return limit(pXlZl(snp, nexthtx, nexthtz)*pZlZlm1(snp, nexthtz,
	// prevhtz));
	// }
	// }

	private double limit(double p) {
		return Math.max(p, LIMIT);
	}

	private double getLambda() {
		lamb = 0;
		int i;
		for (i = 1; i < K; i++) {
			lamb += 1.0 / i;
		}
		lamb = 1 / lamb;
		lamb = lamb / (2 * (lamb + K));
		return lamb;
	}

	public Snp[] getBound(Snp[] snps2) {
		int i, t, n, m = this.endSnps.length;
		for(i = t = 0, n = snps2.length; i < n && t < B; i++){
			int[] type = snps2[i].getTypes()[this.sample];
			if(type[0] != type[1]){
				t++;
			}
		}
		Snp[] bound = new Snp[i+m];
		for(n = i, i = 0; i < m; i++){
			bound[i] = this.endSnps[i];
		}
		for(i = 0; i < n; i++){
			bound[m+i] = snps2[i];
		}
		return bound;
	}
	
	public void setHt(Snp[] snps, HaploType[] bound, HaploType[] haploTypes) {
		this.addSnps(snps);
		this.hts = new HaploType[2];
		if(bound == null){
			this.setEnd(haploTypes[0].getAlleles(), haploTypes[1].getAlleles());
			this.hts[0] = new HaploType(haploTypes[0].getAlleles(), this.sample);
			this.hts[1] = new HaploType(haploTypes[1].getAlleles(), this.sample);
			return;
		}else{
			int[] ht1alleles = bound[0].getAlleles();
			int[] ht2alleles = bound[0].getAlleles();
			int e = this.endSnps.length;
			for(int i = 0; i+e < ht1alleles.length; i++){
				if(ht1alleles[i+e] != ht2alleles[i+e]){
					if(haploTypes[0].getAlleles()[i] == ht1alleles[i]){
						this.setEnd(haploTypes[0].getAlleles(), haploTypes[1].getAlleles());
						this.hts[0] = new HaploType(haploTypes[0].getAlleles(), this.sample);
						this.hts[1] = new HaploType(haploTypes[1].getAlleles(), this.sample);						
					}else{
						this.setEnd(haploTypes[1].getAlleles(), haploTypes[0].getAlleles());
						this.hts[0] = new HaploType(haploTypes[1].getAlleles(), this.sample);
						this.hts[1] = new HaploType(haploTypes[0].getAlleles(), this.sample);												
					}
					return;
				}
			}
		}
		this.setEnd(haploTypes[0].getAlleles(), haploTypes[1].getAlleles());
		this.hts[0] = new HaploType(haploTypes[0].getAlleles(), this.sample);
		this.hts[1] = new HaploType(haploTypes[1].getAlleles(), this.sample);								
	}

}
