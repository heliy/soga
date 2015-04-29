package calculate;

import parameter.Setting;
import dna.Snp;
import dna.LD;

public class PairwiseCal {
	private double TOLERENCE = 1e-8;
	private double LN10 = Math.log(10);
	private double PROBMIN = 1e-10;
	private int COUNTMAX = 1000;
	private double LOGINIT = (-(1e9) - 1);
	private int LSURLEN = 100;
	
//	private int[] lpCount;
//	private double[] lpFreqs;

	private int SAMPLES;
	private int N = -1;
	private int AA = 0;
	private int AB = 1;
	private int BA = 2;
	private int BB = 3;
	
	private double LDCL;
	private double LDCU;
	private double EHRCU;
	
	private static double z = 1.75; // 0.975
	
	public PairwiseCal(Setting rc) {
		SAMPLES = rc.getSAMPLES();
		LDCL = rc.getLDCL();
		LDCU = rc.getLDCU();
		EHRCU = rc.getEHRCU();
//		lps = new int
	}

	public LD pw(Snp s1, Snp s2) {
		int unknown = 0;
		int[] known = new int[4];
		int[][] onetypes = s1.getTypes();
		int[][] twotypes = s2.getTypes();

		int i;
		for (i = 0; i < SAMPLES; i++) {
			int[] typeOne = onetypes[i];
			int[] typeTwo = twotypes[i];
			if (typeOne[0] == N || typeOne[1] == N || typeTwo[0] == N
					|| typeTwo[1] == N) {
				continue;
			}
			if (typeOne[0] != typeOne[1] && typeTwo[0] != typeTwo[1]) {
				unknown++;
			} else {
				known[typeOne[0] * 2 + typeTwo[0]]++;
				known[typeOne[1] * 2 + typeTwo[1]]++;
			}
		}
//		System.out.println(testcal(known, unknown).toString());
//		System.out.println(cal(known, unknown).toString());

		return cal(known, unknown);
	}

	private LD testcal(int[] known, int unknown){
		double total_chroms, p_A1, p_A2, p_B1, p_B2, p_multi, const_prob;
		double prob_haps[] = { 0.1, 0.1, 0.1, 0.1 };
		double num_haps[] = new double[4], temp_array[] = new double[4];
		double loglike, old_loglike, aa_bb, ab_ba, D, Dmax;
		double aa, bb, ab, ba, vd, vdprime, f1, f2, f3, sum_prob;
		int type, count, i;
		double rsq, dprime, ciLow, ciHigh, lod;

		total_chroms = known[0] + known[1] + known[2] + known[3] + 2 * unknown;
		p_A1 = (known[AA] + known[AB] + unknown) / total_chroms;
		p_A2 = (known[AA] + known[BA] + unknown) / total_chroms;
		p_B1 = 1 - p_A1;
		p_B2 = 1 - p_A2;
		p_multi = p_A1 * p_A2 * p_B1 * p_B2;

		const_prob = 0.1;
		count_haps(num_haps, known, prob_haps, unknown, 0);
		estimate_p(prob_haps, num_haps, const_prob);

		const_prob = 0.0;
		loglike = LOGINIT;
		for (count = 1; count < COUNTMAX; count++) {
			old_loglike = loglike;
			count_haps(num_haps, known, prob_haps, unknown, count);
			aa_bb = prob_haps[AA] * prob_haps[BB];
			ab_ba = prob_haps[AB] * prob_haps[BA];

			sum_prob = 0.0;
			for (i = 0; i < 4; i++) {
				temp_array[i] = known[i] * Math.log(prob_haps[i]);
				sum_prob += temp_array[i];
			}
			loglike = (sum_prob + unknown * Math.log(aa_bb + ab_ba))
					/ LN10;
			if (Math.abs(loglike - old_loglike) <= TOLERENCE) {
				break;
			}
			estimate_p(prob_haps, num_haps, const_prob);
		}

		lod = loglike
				- (known[AA] * Math.log(p_A1 * p_A2) + known[AB]
						* Math.log(p_A1 * p_B2) + known[BA]
						* Math.log(p_B1 * p_A2) + known[BB]
						* Math.log(p_B1 * p_B2) + unknown
						* Math.log(2 * p_multi)) / LN10;
		
		aa = prob_haps[AA] + prob_haps[AB];            // fu1
		ab = prob_haps[BA] + prob_haps[BB];            // fu2
		ba = prob_haps[AA] + prob_haps[BA];            // fv1
		bb = prob_haps[AB] + prob_haps[BB];            // fv2
		
		// 按定义来 怎么不对。。。
//		System.out.println(aa+","+ab+","+ba+","+bb);
//		D = prob_haps[AA] - aa*ab;
//		Dmax = (D < 0) ? Math.min(aa*ba, ab*bb) : Math.min(aa*bb, ab*ba); 
//		dprime = D/Dmax;
//		rsq = D*D/(aa*ab*ba*bb);
//		System.out.println(D+","+Dmax);
				
		aa_bb = prob_haps[AA] * prob_haps[BB];
		ab_ba = prob_haps[AB] * prob_haps[BA];
		D = aa_bb - ab_ba;
//		System.out.println(D);

		if (D < 0) {
			swap(prob_haps);
			double temp = p_A2;
			p_A2 = p_B2;
			p_B2 = temp;
			swap(num_haps);
			aa_bb = prob_haps[AA] * prob_haps[BB];
			ab_ba = prob_haps[AB] * prob_haps[BA];
			D = aa_bb - ab_ba;
			swap(known);
		}
		aa = (prob_haps[AA] + prob_haps[AB]);
		bb = (prob_haps[BA] + prob_haps[BB]);
		ab = (prob_haps[AA] + prob_haps[AB]);
		ba = (prob_haps[AB] + prob_haps[BB]);
		aa_bb = aa*bb;
		ab_ba = ab*ba;
		Dmax = (aa_bb < ab_ba) ? aa_bb : ab_ba;
		dprime = D / Dmax;
		rsq = D * D / p_multi;

		f1 = (dprime > 0) ? ba : bb;
		f2 = ba + bb - f1;
		if(Dmax == aa*ba){
			f3 = prob_haps[AA];
		}else if(Dmax == ab*bb){
			f3 = prob_haps[BB];
		}else if(Dmax == aa*bb){
			f3 = prob_haps[AB];
		}else{
			f3 = prob_haps[BA];
		}

        vd = (aa*ab*ba*bb + D*(ab-aa)*(bb-ba) - D*D)/SAMPLES;
        vdprime = ((1-Math.abs(dprime))*(SAMPLES*vd-Math.abs(dprime)*Dmax*(aa*f1+ab*f2-2*Math.abs(D)))
        		   +Math.abs(dprime)*f3*(1-f3))/(Dmax*Dmax*SAMPLES);
        vdprime = z*Math.sqrt(vdprime);
        ciLow = dprime - vdprime;
        ciHigh = dprime + vdprime;
        aa = Math.min(Math.abs(ciLow), Math.abs(ciHigh));
        ab = Math.max(Math.abs(ciLow), Math.abs(ciHigh));
        ciLow = aa;
        ciHigh = ab;
		if(ciLow >= LDCL && ciHigh >= LDCU){
			type = 1;  // strong LD
		}else if(ciHigh < EHRCU){
			type = 2;  // strong EHR
		}else{
			type = 0;  // otherwise
		}

		LD ld = new LD(rsq, dprime, ciLow, ciHigh, lod, type);
		return ld;
	}
	private LD cal(int[] known, int unknown) {
		double total_chroms, p_A1, p_A2, p_B1, p_B2, p_multi, const_prob;
		double prob_haps[] = { 0.1, 0.1, 0.1, 0.1 };
		double num_haps[] = new double[4], temp_array[] = new double[4], lsurface[] = new double[LSURLEN + 1];
		double loglike, old_loglike, aa_bb, ab_ba, num, temp;
		double aa, bb, ab, ba;
		int type, count, i, j, low_i = 0, high_i = 0;
		double denom, t_array[] = new double[4], total_prob, sum_prob;
		double rsq, dprime, ciLow, ciHigh, lod;

		total_chroms = known[0] + known[1] + known[2] + known[3] + 2 * unknown;
		p_A1 = (known[AA] + known[AB] + unknown) / total_chroms;
		p_A2 = (known[AA] + known[BA] + unknown) / total_chroms;
		p_B1 = 1 - p_A1;
		p_B2 = 1 - p_A2;
		p_multi = p_A1 * p_A2 * p_B1 * p_B2;

		const_prob = 0.1;
		count_haps(num_haps, known, prob_haps, unknown, 0);
		estimate_p(prob_haps, num_haps, const_prob);

		const_prob = 0.0;
		loglike = LOGINIT;
		for (count = 1; count < COUNTMAX; count++) {
			old_loglike = loglike;
			count_haps(num_haps, known, prob_haps, unknown, count);
			aa_bb = prob_haps[AA] * prob_haps[BB];
			ab_ba = prob_haps[AB] * prob_haps[BA];

			sum_prob = 0.0;
			for (i = 0; i < 4; i++) {
				temp_array[i] = known[i] * Math.log(prob_haps[i]);
				sum_prob += temp_array[i];
			}
			loglike = (sum_prob + unknown * Math.log(aa_bb + ab_ba))
					/ LN10;
			if (Math.abs(loglike - old_loglike) <= TOLERENCE) {
				break;
			}
			estimate_p(prob_haps, num_haps, const_prob);
		}

		lod = loglike
				- (known[AA] * Math.log(p_A1 * p_A2) + known[AB]
						* Math.log(p_A1 * p_B2) + known[BA]
						* Math.log(p_B1 * p_A2) + known[BB]
						* Math.log(p_B1 * p_B2) + unknown
						* Math.log(2 * p_multi)) / LN10;
				
		aa_bb = prob_haps[AA] * prob_haps[BB];
		ab_ba = prob_haps[AB] * prob_haps[BA];
		num = aa_bb - ab_ba;
	//	System.out.println(num);

		if (num < 0) {
			swap(prob_haps);
			temp = p_A2;
			p_A2 = p_B2;
			p_B2 = temp;
			swap(num_haps);
			aa_bb = prob_haps[AA] * prob_haps[BB];
			ab_ba = prob_haps[AB] * prob_haps[BA];
			num = aa_bb - ab_ba;
			swap(known);
		}
		aa = (prob_haps[AA] + prob_haps[AB]);
		bb = (prob_haps[BA] + prob_haps[BB]);
		ab = (prob_haps[AA] + prob_haps[AB]);
		ba = (prob_haps[AB] + prob_haps[BB]);
		aa_bb = aa*bb;
		ab_ba = ab*ba;
		denom = (aa_bb < ab_ba) ? aa_bb : ab_ba;
		dprime = num / denom;
		rsq = num * num / p_multi;
//		System.out.println(num+","+denom);

		total_prob = 0;
		for (i = 0; i < LSURLEN + 1; i++) {
			temp = i * 0.01 * denom + p_A1 * p_A2;
			t_array[AA] = temp > PROBMIN ? temp : PROBMIN;
			temp = p_A1 - t_array[AA];
			t_array[AB] = temp > PROBMIN ? temp : PROBMIN;
			temp = p_A2 - t_array[AA];
			t_array[BA] = temp > PROBMIN ? temp : PROBMIN;
			temp = p_B1 - t_array[BA];
			t_array[BB] = temp > PROBMIN ? temp : PROBMIN;
			
            sum_prob = 0.0;
			for (j = 0; j < 4; j++) {
				temp_array[j] = known[j] * Math.log(t_array[j]);
				sum_prob += temp_array[j];
			}
			lsurface[i] = (sum_prob + unknown
					* Math.log(t_array[AA] * t_array[BB] + t_array[AB]
							* t_array[BA]))
					/ LN10;
			lsurface[i] = Math.pow(10, lsurface[i] - loglike);
			total_prob += lsurface[i];
		}

		sum_prob = 0.0;
		for (i = 0; i <= LSURLEN; i++) {
			sum_prob += lsurface[i];
			if (sum_prob > 0.05 * total_prob
					&& sum_prob - lsurface[i] < 0.05 * total_prob) {
				low_i = i - 1;
				break;
			}
		}
		sum_prob = 0.0;
		for (i = LSURLEN; i >= 0; i--) {
			sum_prob += lsurface[i];
			if (sum_prob > 0.05 * total_prob
					&& sum_prob - lsurface[i] < 0.05 * total_prob) {
				high_i = i + 1;
				break;
			}
		}

		if (high_i > 100) {
			high_i = 100;
		}
		ciLow = ((float) low_i) / LSURLEN;
		ciHigh = ((float) high_i) / LSURLEN;
	
//		if(ciLow > 0.3){
		if(ciLow >= LDCL && ciHigh >= LDCU){
			type = 1;  // strong LD
//		}else if(ciHigh < 0.3){
		}else if(ciHigh < EHRCU){
			type = 2;  // strong EHR
		}else{
			type = 0;  // otherwise
		}

		LD ld = new LD(rsq, dprime, ciLow, ciHigh, lod, type);
		return ld;
	}

	private void swap(int[] array) {
		int temp;
		temp = array[0];
		array[0] = array[1];
		array[1] = temp;
		temp = array[2];
		array[2] = array[3];
		array[3] = temp;
	}

	private void swap(double[] array) {
		double temp;
		temp = array[0];
		array[0] = array[1];
		array[1] = temp;
		temp = array[2];
		array[2] = array[3];
		array[3] = temp;
	}

	private void estimate_p(double[] probs, double[] nums, double const_prob) {
		double total;
		int i;
		total = nums[0] + nums[1] + nums[2] + nums[3] + 4.0 * const_prob;
		for (i = 0; i < 4; i++) {
			probs[i] = (nums[i] + const_prob) / total;
			probs[i] = (probs[i] < PROBMIN) ? PROBMIN : probs[i];
		}
	}

	private void count_haps(double[] nums, int[] known, double[] haps,
			int unknown, int em_record) {
		  int i;
		  double aa_bb,ab_ba;
		  double probs[] = new double[4];

		  if(em_record > 0){
		    aa_bb = haps[AA] * haps[BB];
		    ab_ba = haps[AB] * haps[BA];
		    probs[0] = aa_bb / (aa_bb + ab_ba);
		    probs[1] = 1 - probs[0];
		    probs[2] = probs[1];
		    probs[3] = probs[0];

		    for(i = 0; i < 4; i++){
		      nums[i] = known[i] + probs[i] * unknown;
		    }
		  }else{
		    for(i = 0; i < 4; i++){
		      nums[i] = known[i];
		    }
		  }

	}

}
