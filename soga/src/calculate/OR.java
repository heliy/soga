package calculate;

public class OR {
	public double[] getOR_CI(double counts, double d, double counts2, double e){
		double[] ors = new double[3];
		ors[0] = (counts*e+0.01)/(d*counts2+0.01);
		double t = 0;
		if(counts != 0){ t += 1.0/counts; }
		if(d != 0){ t += 1.0/d; }
		if(counts2 != 0){ t += 1.0/counts2; }
		if(e != 0){ t += 1.0/e; }
		t = Math.sqrt(t);
		
		ors[1] = Math.exp(Math.log(ors[0])-1.96*t);
		ors[2] = Math.exp(Math.log(ors[0])+1.96*t);
		return ors;
	}

}
