package control;

import io.SampleQCFp;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

import exceptions.AlleleException;
import exceptions.ArgsException;
import exceptions.GenoTypeException;
import exceptions.OrderException;
import exceptions.SampleException;
import exceptions.SampleStatusException;
import exceptions.SnpContainsException;
import parameter.Setting;
import parameter.Summary;
import population.Sample;

public class Main {

	public static void main(String[] args) throws FileNotFoundException, ArgsException, SampleStatusException, AlleleException, GenoTypeException, SnpContainsException, OrderException, InterruptedException, SampleException {
		Setting rc = new Setting();
		if(!rc.parseArgs(args)){
			return;
		}
		if(rc.needCheck()){
			Check check = new Check(rc);
			Sample[] samples = check.run();
			if(rc.doSAMPLEQC()){
				rc.setSamples(samples);
				SampleQCFp writer = new SampleQCFp(rc);
				for(Sample sample: rc.getSamples()){
					writer.write(sample);
				}
				writer.close();
			}else{
				rc.passAllSample();
			}
		}else{
			rc.passAllSample();
		}
		
		System.out.println("==============================================");
		if(rc.display(false)){
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
