package parameter;

import java.math.BigInteger;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import io.FileOutput;
import dna.Block;
import dna.LD;
import dna.Snp;
import exceptions.ChromosomeSegmentException;

public class Summary {
	private Date begin;
	private Date end;
	
	private int samples;
	private int cases;
	private int controls;
	private int unknowns;
	
	private long snps;
	private long oksnps;
	private HashMap<String, Long[]> snpPerChr;
	private BigInteger lds;
	private long blocks;
	private long hotspots;
	private HashMap<String, Long> blockPerChr;
	private HashMap<String, Long> hotspotPerChr;
	
	private String snpfile;
	private String samplefile;
	private String[] outputs;
	private int outnum;
	
	public Summary(Setting rc){
		begin = new Date();
		samples = rc.getSAMPLES();
		cases = rc.getCases();
		controls = rc.getControls();
		unknowns = rc.getUnknowns();
		
		snpfile = rc.getSnpFile();
		samplefile = rc.getSampleFile();
		outputs = new String[100];
		outnum = 0;
		
		oksnps = snps = blocks = hotspots = 0;
		
		lds = BigInteger.valueOf(0);
		snpPerChr = new HashMap<String, Long[]>();
		blockPerChr = new HashMap<String, Long>();
		hotspotPerChr = new HashMap<String, Long>();
	}
	
	public void add(FileOutput output){
		outputs[outnum++] = output.getFileName();
	}

	public void add(Snp snp, String chr, boolean chrIsIn) throws ChromosomeSegmentException {
		Long[] nums;
		snps++;
		if(!snp.isBad()){
			oksnps++;
		}
		if(chrIsIn){
			nums = snpPerChr.get(chr);
			nums[0]++;
			if(!snp.isBad()){
				nums[1]++;
			}
			snpPerChr.put(chr, nums);
		}else{
			if(snpPerChr.containsKey(chr)){
				throw new ChromosomeSegmentException(chr);
			}
			nums = new Long[2];
			nums[0]= (long) 1;
			if(!snp.isBad()){
				nums[1] = (long) 1;
			}else{
				nums[1] = (long) 0;
			}
			snpPerChr.put(chr, nums);
		}
	}
	
	public void add(Block block, boolean isBlock){
		String chr = block.getSnps()[0].getChr();
		Long num;
		if(isBlock){
    		blocks++;
    		if(blockPerChr.containsKey(chr)){
    			num = blockPerChr.get(chr);
    			num++;
    		}else{
    			num = (long)1;
    		}
    		blockPerChr.put(chr, num);
		}else{
			hotspots++;
    		if(hotspotPerChr.containsKey(chr)){
    			num = hotspotPerChr.get(chr);
    			num++;
    		}else{
    			num = (long)1;
    		}
    		hotspotPerChr.put(chr, num);
		}
	}
	
	public void end(){
		end = new Date();
	}

	public void add(LD ld) {
		lds.and(BigInteger.valueOf(0));
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append("BEGIN: "+begin + "\n");
		sb.append("END: "+end+"\n");
		sb.append("Snp Files: "+snpfile+"\n");
		sb.append("Sample Info File: "+samplefile+"\n");
		sb.append("OUTPUT File:\n");
		for(int i = 0; i < outnum; i++){
			sb.append(outputs[i]+"\n");
		}
		sb.append("Samples:\n"+"\tTotal: "+samples+"\n\tCases: "+cases+"\n\tControls: "+controls+"\n\tUnknowns: "+unknowns+"\n");
		sb.append("SNPS: \tCount \tpassQC\n");
		sb.append("total: \t"+snps+"\t"+oksnps+"\n");
		Set<String> chrs = snpPerChr.keySet();
		String chr;
		Long[] nums;
		Iterator<String> iter = chrs.iterator();
		while(iter.hasNext()){
			chr = iter.next();
			nums = snpPerChr.get(chr);
			if(nums == null){
				nums = new Long[2];
				nums[0] = nums[1] = (long)0;
			}
			sb.append(chr+": \t"+nums[0]+"\t"+nums[1]+"\n");
		}
		if(blocks != 0){
    		sb.append("BLOCKS: \tCount\t\n");
	    	sb.append("total: \t"+blocks+"\n");
		    iter = chrs.iterator();
		    while(iter.hasNext()){
			    chr = iter.next();
			    Long num = blockPerChr.get(chr);
			    if(num == null){
				    num = (long) 0;
			    }
			    sb.append(chr+": \t"+num+"\n");
		    }
			if(hotspots != 0){
	    		sb.append("HOTSPOTS: \tCount\t\n");
		    	sb.append("total: \t"+hotspots+"\n");
			    iter = chrs.iterator();
			    while(iter.hasNext()){
				    chr = iter.next();
				    Long num = hotspotPerChr.get(chr);
				    if(num == null){
					    num = (long) 0;
				    }
				    sb.append(chr+": \t"+num+"\n");
			    }
			}
		}
		return sb.toString();
	}

}
