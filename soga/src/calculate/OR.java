package calculate;

public class OR {
        private double LIMIT1 = 0.00000001;
        private double LIMIT2 = 999999999.0;
		/*
		 * -----------------
		 *        | A  | B
		 * -----------------
		 *  case  | a  | b
		 * -----------------
		 * control| c  | d
		 * -----------------
		 */

	public double[] getOR_CI(double a, double b, double c, double d){
		double[] ors = new double[3];
		ors[0] = limit((a*d+0.01)/(b*c+0.01));
		double t = 0;
		if(a != 0){ t += 1.0/a; }
		if(b != 0){ t += 1.0/b; }
		if(c != 0){ t += 1.0/c; }
		if(d != 0){ t += 1.0/d; }
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
