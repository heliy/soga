package exceptions;

public class SampleException extends Exception {

	private static final long serialVersionUID = 7L;
	
	public SampleException(String line){
		System.out.println("Ooooooooops! You have unvalided information about samples in input sample file!");
		System.out.println(line);
		System.out.println("You MUST have NAME and STAUS in every line (split by '\\t')");
		System.out.println("More help: `--help`");	
	}

}
