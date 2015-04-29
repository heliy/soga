package calculate;

import java.util.Random;
import java.util.concurrent.CountDownLatch;

import parameter.Setting;
import dna.HaploType;
import dna.Snp;

public class Phase implements Runnable {
	private HaploType[] htall; // 所有单体型
	private int htbegin;
	private int htend;
	
	static private int B; // B个一段
	static private int T; // T = 2**B
	static private double LIMIT = 0.0001;// 最小概率
	static private int Ne = 15000; // 群体个数

	private int L; // snp个数
	private int N; // 样本个数
	private int K; // 考虑的备选单体型数
	private int bits[][]; // 段内单体型的0/1
	private int hetes[]; // 每个样本的AB型个数

	private int genotypes[][]; // 基因型
	private int H[][]; // 备选单体型
	private int getH[][]; // 表示H[j]是从第i个基因型中选出的
	private int genoHs[]; // 从第i个基因型中选出几个

	private Snp snps[];
	private double lamb; // lambda
	private double thetas[]; // e^(-ro/K)
	private double unthetas[]; // (1-e^(-ro/K))/K
	private double mafs[]; // 最小等位频率

	private double PZ1; // P(Z1 = u) = 1/K

	CountDownLatch latch;
	
	public Phase(Snp[] Snps, int htBegin, int htEnd, Setting rc, CountDownLatch latch) {
		this.latch = latch;
		
		int i, j, sum, hNum;
		double selectp, distance;
		int types[][], pows[];
		
		htbegin = htBegin;
		htend = htEnd; 
		B = rc.getB();
		T = (int) Math.pow(2, B);
		N = rc.getSAMPLES();
		K = rc.getK();

		snps = Snps;
		L = snps.length;
		lamb = getLambda();
		thetas = new double[L];
		unthetas = new double[L];
		mafs = new double[L];
		genotypes = new int[N][L];
		hetes = new int[N];
		pows = new int[N];

//		System.out.println(L + ", " + N + ", " + K);

		GeneticDistance gd = new GeneticDistance();
		mafs[0] = snps[0].getMaf();
		for (i = 1; i < L; i++) {
			distance = gd.ditance(snps[i - 1], snps[i]);
			thetas[i] = Math.exp((-4) * Ne * distance / (double) K);
			unthetas[i] = (1 - thetas[i]) / K;
			mafs[i] = snps[i].getMaf();
		}

		// 基因型
		for (i = 0; i < L; i++) {
			types = snps[i].getTypes();
			for (j = 0; j < N; j++) {
				genotypes[j][i] = types[j][0] + types[j][1];
			}
		}

		sum = 0;
		for (i = 0; i < N; i++) {
			hetes[i] = 0;
			for (j = 0; j < L; j++) {
				hetes[i] += (genotypes[i][j] == 1) ? (1) : (0);
			}
//			System.out.println(i+": "+hetes[i]);
			pows[i] = ((int) Math.pow(2, hetes[i])) % K;
			sum += pows[i];
		}
		K = (sum < K) ? (sum) : (K);

		H = new int[K][L];
		getH = new int[N][K];
		genoHs = new int[N];

		Random r = new Random();
		hNum = 0;
		for (i = 0; i < N; i++) {
			selectp = Math.max(1.0 / K, K / (N * pows[i]));
			genoHs[i] = 0;
			for (j = 0; j < pows[i]; j++) {
				if (r.nextDouble() > selectp) {
					continue;
				}
				getCandidate(i, hNum);
				getH[i][genoHs[i]++] = hNum;
				hNum++;
			}
			if (hNum == K) {
				break;
			}
		}
		while (hNum < K) {
			i = r.nextInt(N);
			getCandidate(i, hNum);
			getH[i][genoHs[i]++] = hNum;
			hNum++;
		}

		bits = new int[T][B];
		for (i = 0; i < T; i++) {
			sum = i;
			// i二进制化
			for (j = 0; j < B; j++) {
				bits[i][j] = sum % 2;
				sum /= 2;
			}
		}
		PZ1 = 1.0 / K;
	}
	
	public void run(){
		int i, n = 0;
		htall = new HaploType[(htend-htbegin)*2];
//		System.out.println("----"+htend+", "+htbegin+": "+htall.length);
		HaploType[] hts;
		for(i = htbegin; i < htend; i++){
			hts = phase(i);
//			System.out.println(hts[]);
			htall[n++] = hts[0];
			htall[n++] = hts[1];
		}
		latch.countDown();
	}
	
	public HaploType[] getHts(){
		return htall;
	}

	public HaploType[] phase(int index) {
//		System.out.println(index + ": " + hetes[index]);
		if (hetes[index] < 2) {
			HaploType[] hts = new HaploType[2];
			int[][] diplos = new int[2][L];
			int i;
			for (i = 0; i < L; i++) {
				if (genotypes[index][i] == 1) {
					diplos[0][i] = snps[i].getA();
					diplos[1][i] = snps[i].getB();
				} else if (genotypes[index][i] == 0) {
					diplos[1][i] = diplos[0][i] = snps[i].getA();
				} else {
					diplos[1][i] = diplos[0][i] = snps[i].getB();
				}
			}
			hts[0] = new HaploType(diplos[0], index);
			hts[1] = new HaploType(diplos[1], index);
			return hts;
		} else {
			return shapeit(index);
		}
	}

	private HaploType[] shapeit(int genoIndex) {
		// 计算用
	    int S[] = new int[L]; // 分段
		int A[][] = new int[K][L] ; // 断内可能单体型
		double alpha[][][] = new double[L][T][K]; // 前向概率
		double beta[][][] = new double[L][T][K]; // 后向概率
		double sumsOneAll[] = new double[T]; // 行和
		double sumsAllOne[] = new double[K]; // 列和
		double inturn[][] = new double[T][T]; // 内部转移矩阵
		double outturn[][] = new double[T][T]; // 外部转移矩阵
		double tmpturn[][] = new double[T][T]; // 临时转移矩阵

		int[] genotype = genotypes[genoIndex];
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
//			System.out.println(i + ": " + S[i]);
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
			X1[0] = maxi;
			X2[0] = getCp(maxi);
		}
 //       System.out.println(S[L-1]);		
		if (S[L-1] > 0) {
			while (true) {
				m++;
//                System.out.println(m);
				// 计算前一个位点到这个位点的转移概率
				for (i = 0; i < T; i++) {
					for (j = 0; j < T; j++) {
						tmpsum = 0;
						for (u1 = 0; u1 < K; u1++) {
							for (u2 = 0; u2 < K; u2++) {
								tmpsum += (alpha[m - 1][i][u1]
										* pZlZlm1(m, u2, u1) * pXlZl(A, m, j, u2) * beta[m][j][u2]);
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
					if (S[m] == S[L-1]) {
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
		int[] ht1 = new int[L];
		int[] ht2 = new int[L];
		int a, b;
		for (i = 0; i < L; i++) {
			a = snps[i].getA();
			b = snps[i].getB();
			ht1[i] = (A[X1[S[i]]][i] == 0) ? (a) : (b);
			ht2[i] = (A[X2[S[i]]][i] == 0) ? (a) : (b);
		}
		hts[0] = new HaploType(ht1, genoIndex);
		hts[1] = new HaploType(ht2, genoIndex);
		
		update(genoIndex, ht1, ht2);
		
		return hts;
	}
	
	private void update(int geno, int[] ht1, int[] ht2){
        int i, j, k, n = genoHs[geno], ht[];
        Random r = new Random();
        for(i = 0; i < n; i++){
        	if(r.nextBoolean()){
        		ht = ht1;
        	}else{
        		ht = ht2;
        	}
        	k = getH[geno][i];
        	for(j = 0; j < L; j++){
        		H[k][j] = ht[j];
        	}
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
//		System.out.println(snpIndex);
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

	private void getCandidate(int index, int h) {
		int[] genotype = genotypes[index];
		int i;
		Random r = new Random();
		for (i = 0; i < L; i++) {
			if (genotype[i] == 1) {
				if (r.nextDouble() < mafs[i]) {
					H[h][i] = 1;
				} else {
					H[h][i] = 0;
				}
			} else {
				H[h][i] = genotype[i] / 2;
			}
		}
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

}
