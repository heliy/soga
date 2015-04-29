package exceptions;

public class AlleleException extends Exception {
	private static final long serialVersionUID = 1L;
	
	public AlleleException(){
	}
	
	public AlleleException(String file, String line, String alleles){
		System.out.println("Have Unvalided Alleles:");
		System.out.println("Snp: "+line);
		System.out.println("In File: "+file);
		System.out.println("At: "+alleles);
	}
}
