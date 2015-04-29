package calculate;

public class FitTest {
	private double[] distrib;
	private static double TRESH = 1e-8;
	private static double MIN = 1e-14;
	private static double RECIP = 1/Math.E;
	private static double TWOPI = Math.PI*2;
	
	public FitTest(double[] dis){
		distrib = dis;
	}

	public double test(double[] observed){
		int i, len = observed.length;
		double sum = 0, chi2 = 0, expected;
		for(double d : observed){
			sum += d;
		}
		for(i = 0; i < len; i++){
			expected = distrib[i] * sum;
			chi2 += (observed[i] - expected)*(observed[i] - expected)/expected;
		}
		return pvalue(chi2, len - 1);
	}
	
	private double pvalue(double chi2, int dof){
		if(chi2 < 0 || dof < 1){
			return 0.0;
		}

		double k, v;
		k = dof * 0.5;
		v = chi2 * 0.5;
		if(dof == 2){
			return Math.exp((-1)*v);
		}
		
		double incomplete = log_igf(k,v);
		double expIncomplete = Math.exp(incomplete);
        if(expIncomplete <= TRESH || Double.isInfinite(expIncomplete) || Double.isNaN(expIncomplete)){
        	return 1.0;
        }
        
        double gamma = Math.log(approx(k));
        incomplete -= gamma;
		expIncomplete = Math.exp(incomplete);
		if(expIncomplete > 1){
			return MIN;
		}
		return 1 - expIncomplete;
	}

	private double approx(double k) {
		double d = 1/(10*k);
		d = 1/(12 * k - d);
		d = (d + k) * RECIP;
		d = Math.pow(d, k);
		d *= Math.sqrt(TWOPI/k);
		return d;
	}

	private double log_igf(double s, double z) {
		if(z < 0.0){
			return 0.0;
		}
		double sc, k;
		sc = (Math.log(z)*s) - z - Math.log(s);
		k = KM(s, z);
		return Math.log(k) + sc;
	}

	private double KM(double s, double z) {
		double sum, nom, denom, logNom, logDenom, logS, logZ;
		sum = nom = denom = 1.0;
		logNom =  Math.log(nom);
		logDenom = Math.log(denom);
		logS = Math.log(s);
		logZ = Math.log(z);
		
		int i;
		for(i = 0; i < 1000; ++i){
			logNom += logZ;
			logS  = Math.log(++s);
			logDenom += logS;
			sum += Math.exp(logNom - logDenom);
		}
		return sum;
	}
}
