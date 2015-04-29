package exceptions;

public class ArgsException extends Exception {

	private static final long serialVersionUID = 5L;
	
	public ArgsException(String info){
		System.out.println("Ooooops! You have some wrong in your parameters.");
		System.out.println(info);
		System.out.println("More help: `--help`");
	}

}
