package control;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import exceptions.AlleleException;
import exceptions.ArgsException;
import exceptions.ChromosomeSegmentException;
import exceptions.GenoTypeException;
import exceptions.SampleStatusException;
import exceptions.SnpContainsException;
import parameter.Setting;
import parameter.Summary;

public class Main {

	public static void main(String[] args) throws FileNotFoundException, ArgsException, SampleStatusException, AlleleException, GenoTypeException, SnpContainsException, ChromosomeSegmentException, InterruptedException {
		Setting rc = new Setting();
		if(rc.parseArgs(args)){
			System.out.println("==============================================");
			Take run = new Take(rc);
			Summary summary = run.run();
			System.out.println("==============================================");
			System.out.println(summary);
			PrintWriter writer = new PrintWriter(rc.getOutput()+".SUMMARY");
			writer.write(summary.toString());
			writer.close();
		}else{
			System.out.println("HAVE NOTHING TO RUN.");
		}
	}

}
