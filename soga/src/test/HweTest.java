package test;

public class HweTest {
	private static double hwe(int hets, int homc, int homr) {
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
		System.out.println(mid);

		hetProbs[mid] = 1.0;
		double sum = 1.0;
		for (currHets = mid; currHets > 1; currHets -= 2) {
			hetProbs[currHets - 2] = hetProbs[currHets] * currHets
					* (currHets - 1.0)
					/ (4.0 * (currHomr + 1.0) * (currHomc + 1.0));
			sum += hetProbs[currHets - 2];
//			System.out.println(sum);
			
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
//			System.out.println(pvalue);
		}

		if (pvalue > 1.0) {
			return 1.0;
		} else {
			return pvalue;
		}

	}

/*	public static void main(String[] args) {
		// TODO Auto-generated method stub
      System.out.println(hwe(102,156,11));
	}
*/
}
