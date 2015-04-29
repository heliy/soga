package exceptions;

public class ChromosomeSegmentException extends Exception {
	private static final long serialVersionUID = 4L;
	
	private String chr;
	
	public ChromosomeSegmentException(String c){
		chr = c;
	}
	
	public String getChr(){
		return chr;
	}
	
	public ChromosomeSegmentException(String file, String line, String ch){
		System.out.println("Have Unvalided Multi-Segment Chromosome:");
		System.out.println("Snp: "+line);
		System.out.println("In File: "+file);
		System.out.println("At: "+ch);
	}

}
