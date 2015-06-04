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
			Summary summary = new Summary(rc);
			if(rc.getPhasedFile() != null){              // 输入的文件是已分相的
				Take run = new Take(true, rc, summary);
				summary = run.run();
			}else{
     			Take run = new Take(false, rc, summary);         
	    		summary = run.run();				 
		    	if(rc.needRescanPhasedFile()){
		    		rc.setPhasedFile(summary.getPhasedFile());
		    		rc.cancelFULL();
		    		rc.cancelCHECK();
		    		rc.cancelFILTER();
		    		rc.cancelBADDATA();
			    	Take run2 = new Take(true, rc, summary);
				    summary = run2.run();
		    	}
			}
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
