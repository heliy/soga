package dna;

import parameter.Setting;
import parameter.Base;
import dna.Snp;
import dna.LD;

public class PairwiseCal {
	double TOLERENCE = 1e-8;
	double LN10 = Math.log(10);
	double PROBMIN = 1e-10;
	int COUNTMAX = 1000;
	double LOGINIT = (-(1e9) - 1);
	int LSURLEN = 100;

	int SAMPLES;
	Base base;
	int N = -1;
	int AA = 0;
	int AB = 1;
	int BA = 2;
	int BB = 4;

	public PairwiseCal(Setting rc) {
		SAMPLES = rc.getSAMPLES();
		base = new Base();
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

		return cal(known, unknown);
	}

	private LD cal(int[] known, int unknown) {
		double total_chroms, p_A1, p_A2, p_B1, p_B2, p_multi, const_prob;
		double prob_haps[] = { 0.1, 0.1, 0.1, 0.1 };
		double num_haps[] = new double[4], temp_array[] = new double[4], lsurface[] = new double[LSURLEN + 1];
		double loglike, old_loglike, aa_bb, ab_ba, num, temp;
		int count, i, j, low_i = 0, high_i = 0;
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
		aa_bb = (prob_haps[AA] + prob_haps[BA])
				* (prob_haps[BA] + prob_haps[BB]);
		ab_ba = (prob_haps[AA] + prob_haps[AB])
				* (prob_haps[AB] + prob_haps[BB]);
		denom = (aa_bb < ab_ba) ? aa_bb : ab_ba;
		dprime = num / denom;
		rsq = num * num / p_multi;

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

		LD ld = new LD(rsq, dprime, ciLow, ciHigh, lod);
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
