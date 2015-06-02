package control;

import java.io.FileNotFoundException;

import dna.Snp;
import exceptions.AlleleException;
import exceptions.ChromosomeSegmentException;
import exceptions.GenoTypeException;
import exceptions.SnpContainsException;
import parameter.Setting;
import parameter.Summary;
import io.BadFp;
import io.CheckFp;
import io.ExtractSnp;
import io.FilterFp;


//  总的相关处理的类
//  run()

public class Take {
	private Setting rc;
	private ExtractSnp input = null;
	
	private PhasePro full = null;
	
	private LDPro ld = null;
	private HaploPro block = null;
	private HaploPro recom = null;

	private CheckFp check = null;
	private FilterFp filter = null;
	private BadFp bad = null;
	
	private String chr;
	private boolean renew;
	
	private Summary summary;
	
	public Take(Setting Rc) throws FileNotFoundException{
		rc = Rc;
		summary = new Summary(rc);

		input = new ExtractSnp(rc);
		if(rc.doCHECK()){
			check = new CheckFp(rc);
			summary.add(check);
		}
    	if(rc.doFILTER()){
	    	filter = new FilterFp(rc);
	    	summary.add(filter);
		}
		if(rc.doBAD()){
			bad = new BadFp(rc);
		    summary.add(bad);
		}
		
		if(rc.doFULL()){
			full = new PhasePro(rc, summary);
		}else{
			if(rc.doLD()){
				ld = new LDPro(rc, summary);
			}
			if(rc.doBLOCK()){
				block = new HaploPro(rc, summary, true);
			}
			if(rc.doRECOM()){
				recom = new HaploPro(rc, summary, false);
			}
		}
		chr = null;
		renew = false;
	}

	public Summary run() throws FileNotFoundException, AlleleException, GenoTypeException, SnpContainsException, ChromosomeSegmentException, InterruptedException{
		String line;
   		Snp snp = input.nextSnp();
		while(snp != null){
//			System.out.print(chr+", "+snp.getChr());
			if(snp.getChr().equals(chr)){
				summary.add(snp, chr, true);				
			}else{
				renew = true;
				summary.add(snp, snp.getChr(), false);
			}
//			System.out.println(renew);
			chr = snp.getChr();
			line = input.getThisline() + "\n";
			if(check != null){
				check.write(snp);
			}
			if(bad != null && snp.isBad()){
				bad.write(line);
			}
			if(filter != null && !snp.isBad()){
				filter.write(line);
			}
			if(renew){
				    if(full != null){
				    	full.restart();
				    }else{
				    	if(ld != null){
				    		ld.restart();
				    	}
				    	if(block != null){
				    		block.restart(summary);
				    	}
				    	if(recom != null){
				    		recom.restart(summary);
				    	}
				    }
			}
			if(!snp.isBad()){
				if(full != null){
					full.add(snp); 
				}else{
					if(ld != null){
						ld.add(snp, summary);
					}
					if(block != null){
						block.add(snp, summary);
					}
					if(recom != null){
						recom.add(snp, summary);
					}
				}
			}
			if(renew){
				System.out.println("New Chromesome: "+chr);					
				renew = false;				
			}
			snp = input.nextSnp();
		}
		if(full != null){
			full.close(true);
		}
		if(check != null){
			check.close();
		}
		if(filter != null){
			filter.close();
		}
		if(bad != null){
			bad.close();
		}
		if(ld != null){
			ld.close();
		}
		if(block != null){
			block.close(summary);
		}
		if(recom != null){
			recom.close(summary);
		}
		
		summary.end();
		return summary;
	}	

}
