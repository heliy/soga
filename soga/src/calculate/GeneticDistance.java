package calculate;

import dna.Snp;

public class GeneticDistance {
    private double DIS = 357042;  // Giraud et al, 2011
    
    public double ditance(Snp snp1, Snp snp2){
    	if(snp1.getChr().equals(snp2.getChr())){
        	return (snp1.getPosition()-snp2.getPosition())/DIS;
    	}else{
    		return -1;
    	}
    }
}
