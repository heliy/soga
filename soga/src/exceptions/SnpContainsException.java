package exceptions;

public class SnpContainsException extends Exception {

	private static final long serialVersionUID = 3L;
	
	public SnpContainsException(String file, String line, int need){
		System.out.println("Have Unvalided Snp Contains:");
		System.out.println("Snp: "+line);
		System.out.println("In File: "+file);
		System.out.println("Need: "+need);		
	}

}
