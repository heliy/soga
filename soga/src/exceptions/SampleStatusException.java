package exceptions;

public class SampleStatusException extends Exception {

	private static final long serialVersionUID = 6L;
	
	public SampleStatusException(int index, int num){
		System.out.println("Statu of Sample No."+index+" is Invalid.");
		System.out.println("Statu Num: "+num);
		System.out.println("Should be -1/0/1");
		System.out.println("More help: `--help`");
	}

}
