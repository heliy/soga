package control;

import java.io.FileNotFoundException;

import dna.Snp;
import exceptions.AlleleException;
import exceptions.GenoTypeException;
import exceptions.OrderException;
import exceptions.SnpContainsException;
import parameter.Setting;
import parameter.Summary;
import io.BadFp;
import io.ExtractSnp;
import io.FilterFp;
import io.SNPFp;


//  总的相关处理的类
//  run()

public class Take {
	private boolean isPhased;
	
	private Setting rc;
	private ExtractSnp input = null;
	
	private PhasePro full = null;
	
	private LDPro ld = null;
	private HaploPro block = null;
	private HaploPro recom = null;

	private SNPFp snpFp = null;
	private FilterFp filter = null;
	private BadFp bad = null;
	
	private String chr;
	private Snp mae = null;
	private boolean renew;
	
	private Summary summary;
	
	
	public void take(boolean isPhased, Setting rc, Summary summary) throws FileNotFoundException{
		this.isPhased = isPhased;
		this.summary = summary;
		this.rc = rc;
		this.input = new ExtractSnp(this.rc, this.isPhased);
		if(rc.doSNPQC()){
			this.snpFp = new SNPFp(this.rc);
			this.summary.add(this.snpFp);
		}
    	if(rc.doFILTER()){
    		this.filter = new FilterFp(this.rc);
    		this.summary.add(this.filter);
		}
		if(rc.doBAD()){
			this.bad = new BadFp(this.rc);
			this.summary.add(this.bad);
		}
		
		if(rc.doFULL()){
			this.full = new PhasePro(this.rc, this.summary);
		}else{
			if(rc.doLD()){
				this.ld = new LDPro(this.rc, this.summary);
			}
			if(rc.doBLOCK()){
				this.block = new HaploPro(this.rc, this.summary, true);
			}
			if(rc.doRECOM()){
				this.recom = new HaploPro(this.rc, this.summary, false);
			}
		}
		this.chr = null;
		this.renew = false;
		if(this.isPhased){
			this.full = null;
		}
	}
	public Take(boolean isPhased, Setting Rc) throws FileNotFoundException{
		summary = new Summary(Rc);
		take(isPhased, rc, summary);
	}

	public Take(boolean isPhased, Setting rc2, Summary summary2) throws FileNotFoundException {
		this.take(isPhased, rc2, summary2);		
	}
	
	public Summary run() throws FileNotFoundException, AlleleException, GenoTypeException, SnpContainsException, OrderException, InterruptedException{
		String line;
   		Snp snp = input.nextSnp();
		while(snp != null){
//			System.out.print(chr+", "+snp.getChr());
			try{
				if(snp.getChr().equals(chr)){
					summary.add(snp, chr, true);				
				}else{
					renew = true;
					summary.add(snp, snp.getChr(), false);
				}
				if(this.mae != null){
					if(!renew && this.mae.getPosition() >= snp.getPosition()){
						throw new OrderException();
					}
				}
			}catch (OrderException e){
				throw new OrderException(this.mae, snp);
			}
			this.mae = snp;
//			System.out.println(renew);h
			chr = snp.getChr();
			line = input.getThisline() + "\n";
			if(this.snpFp != null){
				this.snpFp.write(snp);
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
		if(snpFp != null){
			snpFp.close();
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
