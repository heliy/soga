package control;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;

import calculate.FindBlock;
import calculate.Tag;
import dna.Block;
import dna.PhasedRange;
import dna.Snp;
import io.BlockFp;
import io.FileOutput;
import io.PhasedFp;
import io.RecomFp;
import parameter.Setting;
import parameter.Summary;

public class HaploPro {
	private FileOutput writer;
	private ArrayList<Snp> snplist;
	private boolean isBlock;

	private static int SIZE;
	private static int threads;
	private static int MAXSIZE;
	private Tag tag = null;
	private boolean isCC;
	private Setting rc;
	
	private int blockNo;
	private PhasedFp phased = null;

	private PhasedRange range = null;

	public HaploPro(Setting Rc, Summary summary, boolean isblock) throws FileNotFoundException {
		snplist = new ArrayList<Snp>();
		rc = Rc;
		isBlock = isblock;
		SIZE = rc.getSIZE(isblock);
		blockNo = 0;
		threads = rc.getThreads();
		MAXSIZE = SIZE*threads;
		if(rc.doFULL()){
			this.range = new PhasedRange(rc);
		}
		if(rc.doTAG()){
			tag = new Tag(rc);
		}
		if(isBlock){
			writer = new BlockFp(rc);
			if(rc.doPHASE()){
				phased = new PhasedFp(rc);
				summary.add(phased);
			}
			isCC = rc.doCC();
		}else{
			writer = new RecomFp(rc);
		}
		summary.add(writer);
	}

	private void pro(Summary summary, boolean isLast) throws InterruptedException, FileNotFoundException{
		if(snplist.size() == 0){
			return;
		}
		int THREADS;
		if(snplist.size()%SIZE == 0){
			THREADS = snplist.size()/SIZE;
		}else{
			THREADS = snplist.size()/SIZE+1;
		}
//		System.out.println(snplist.size());
		Snp[][] snps = new Snp[THREADS][SIZE*2];
		Iterator<Snp> iter = snplist.iterator();
		int[] ends = new int[THREADS];
		int i, j = 0;
		for(i = 0; i < THREADS; i++){
			for(j = 0; j < SIZE; j++){
				if(iter.hasNext()){
					snps[i][j] = iter.next();
				}else{
					break;
				}
			}
			ends[i] = j;
		}
	    ends[THREADS-1] = j;
		FindBlock[] finds = new FindBlock[THREADS];
		Thread[] runs = new Thread[THREADS];
		for(i = 0; i < THREADS; i++){
			finds[i] = new FindBlock(rc, isBlock);
//			System.out.println(snps[i][0].getPosition()+", "+snps[i][ends[i]-1].getPosition());
			finds[i].setRun(snps[i], 0, ends[i], isLast);
			runs[i] = new Thread(finds[i]);
		}
		for(i = 0; i < THREADS; i++){
			runs[i].start();
		}
    	try {
    		for(i = 0; i < THREADS; i++){
                runs[i].join();
    		}
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    	ArrayList<Block> blocklist = new ArrayList<Block>();
    	Block[] blocks;
    	int[] already = new int[THREADS];
    	int to;
//		System.out.println(THREADS);
    	for(i = 0; i < THREADS; i++){
    		blocks = finds[i].getBlocks();
    		if(blocks == null || blocks.length == 0){
    			already[i] = 0;
    			ends[i] = 0;
    			continue;
    		}
    		already[i] = blocks.length;
    		if(i > 0){
    		    to = blocks[0].getUp();
    		    for(j = 0; j < to; j++){
//    				System.out.println(i+", "+ends[i-1]+", "+j);
    		    	snps[i-1][ends[i-1]++] = snps[i][j]; 
    		    }
    		}    		
    		for(Block block: blocks){
        		blocklist.add(block);
    		}
    		ends[i] = 0;
    		if(i < THREADS-1){
    			for(j = blocks[blocks.length-1].getDown(); j < SIZE; j++){
//    				System.out.println(i+", "+ends[i]+", "+j);
    				snps[i][ends[i]++] = snps[i][j];
    			}
    		}
    	}
    	for(i = 0; i < THREADS-1; i++){
    		finds[i] = new FindBlock(rc, isBlock);
//			System.out.println(snps[i][0].getPosition()+", "+snps[i][ends[i]-1].getPosition());
    		finds[i].setRun(snps[i], 0, ends[i], isLast);
    		runs[i] = new Thread(finds[i]);
    	}
		for(i = 0; i < THREADS-1; i++){
			runs[i].start();
		}
    	try {
    		for(i = 0; i < THREADS-1; i++){
                runs[i].join();
    		}
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    	j = 0;
    	for(i = 0; i < THREADS-1; i++){
    		j += already[i];
    		blocks = finds[i].getBlocks();
    		if(blocks == null){
    			continue;
    		}
    		for(Block block: blocks){
    			blocklist.add(j++, block);
    		}
    	}
    	blocks = new Block[blocklist.size()];
    	blocklist.toArray(blocks);
    	if(isBlock){
    		System.out.print("Find new Haplotype blocks: ");
    	}else{
    		System.out.print("Find new Recombination Hotspots:");
    	}
		System.out.println(blocks.length);
		if (blocks.length != 0) {
			for(Block b : blocks){
				summary.add(b, isBlock);
				writer.writeline(">> "+blockNo+"\n");
				if(phased != null){
					phased.writeline(">> "+blockNo+"\n");
					b.phase(rc, phased);
				}
				blockNo++;
				if(tag != null){
					b.setTags(tag.getTags(b.getSnps()));
				}
				if(isCC){
					b.ccTest();
				}
				writer.write(b);
			}
		}
		if(!isLast){
			int limit;
			if(blocks.length == 0){
				limit = (SIZE-1)/20;
			}else{
				limit = blocks[blocks.length - 1].getDown();
			}
			limit += SIZE*(THREADS-1);
			for(i = 0; i < limit; i++){
				snplist.remove(0);
			}
//			System.out.println(snplist.size());
		}
	}
	
	public void add(Snp snp, Summary summary) throws InterruptedException, FileNotFoundException {
		snplist.add(snp);
		if (snplist.size() == MAXSIZE) {
			pro(summary, false);
		}
	}

	public void restart(Summary summary) throws InterruptedException, FileNotFoundException{
//        System.out.println(snplist.size());
		pro(summary, true);
		snplist.clear();
	}
	
	public void close(Summary summary) throws InterruptedException, FileNotFoundException {
        pro(summary, true);
		writer.close();
		if(phased != null){
			phased.close();
		}
	}

	public void add(PhasedRange range, Summary summary) throws FileNotFoundException, InterruptedException {
		this.range.extend(range);
		
	}
	
}
