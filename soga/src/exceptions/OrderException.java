package exceptions;

import dna.Snp;

public class OrderException extends Exception {
	private static final long serialVersionUID = 4L;
		
	public OrderException(Snp mae, Snp ato){
		System.out.println("Have Unvalided SNP order in input file:");
		System.out.println("SNP1: "+mae.getRs()+", "+mae.getChr()+", "+mae.getPosition());
		System.out.println("SNP2: "+ato.getRs()+", "+ato.getChr()+", "+ato.getPosition());
		if(mae.getChr().equals(ato.getChr())){
			System.out.println(mae.getPosition()+" -> "+ato.getPosition()+", in same chromosome, SNPs MUST in ASCENDING order.");
		}else{
			System.out.println(ato.getChr()+", SNPs in this chromosome have been appeared before.");
		}
	}
	
	public OrderException(){
		// Nothing ...
	}
	
}
