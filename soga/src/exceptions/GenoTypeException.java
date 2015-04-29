package exceptions;

public class GenoTypeException extends Exception{
	
	private static final long serialVersionUID = 2L;

	private String genotype;
	
	public GenoTypeException(String geno){
		genotype = geno;
	}
	
	public String getGeno(){
		return genotype;
	}
	
	public GenoTypeException(String file, String line, String geno){
		System.out.println("Have Unvalided Genotype:");
		System.out.println("Snp: "+line);
		System.out.println("In File: "+file);
		System.out.println("At: "+geno);
	}

}
