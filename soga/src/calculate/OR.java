package calculate;

public class OR {
        private double LIMIT1 = 0.00000001;
        private double LIMIT2 = 999999999.0;
	public double[] getOR_CI(double counts, double d, double counts2, double e){
		double[] ors = new double[3];
		ors[0] = limit((counts*e+0.01)/(d*counts2+0.01));
		double t = 0;
		if(counts != 0){ t += 1.0/counts; }
		if(d != 0){ t += 1.0/d; }
		if(counts2 != 0){ t += 1.0/counts2; }
		if(e != 0){ t += 1.0/e; }
		t = Math.sqrt(t);
		
		ors[1] = limit(Math.exp(Math.log(ors[0])-1.96*t));
		ors[2] = limit(Math.exp(Math.log(ors[0])+1.96*t));
		return ors;
	}
        private double limit(double v){
	    if(v < LIMIT1){
		return LIMIT1;
	    }else if(v > LIMIT2){
		return LIMIT2;
	    }else{
		return v;
	    }
	}

}
