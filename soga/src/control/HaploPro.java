package control;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import calculate.FindBlock;
import calculate.Tag;
import dna.Block;
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

	private FindBlock find;

	private static int SIZE;
	private static int threads;
	private int MAXSIZE;
	private Tag tag = null;
	private boolean isCC;
	private Setting rc;
	
	private static int DISCARD = 4;

	private int blockNo;
	private PhasedFp phased = null;

	public HaploPro(Setting Rc, Summary summary, boolean isblock)
			throws FileNotFoundException {
		snplist = new ArrayList<Snp>();
		rc = Rc;
		isBlock = isblock;
		SIZE = rc.getSIZE(isblock);
		blockNo = 1;
		threads = rc.getThreads();
		MAXSIZE = SIZE * threads;
		this.MAXSIZE = SIZE;
		if (rc.doTAG()) {
			tag = new Tag(rc);
		}
		find = new FindBlock(rc, isBlock);
		if (isBlock) {
			writer = new BlockFp(rc);
			if (rc.doPHASE()) {
				phased = new PhasedFp(rc);
				summary.add(phased);
			}
			isCC = rc.doCC();
		} else {
			writer = new RecomFp(rc);
		}
		summary.add(writer);
	}

	// 非并行 弹性窗口版
	private void pro(Summary summary, boolean isLast)
			throws FileNotFoundException, InterruptedException {
		if (this.snplist.size() == 0) {
			return;
		}
		Snp[] snps = new Snp[this.snplist.size()];
		this.snplist.toArray(snps);

		Block[] blocks = find.migpp(snps, isLast);
		int limit = 0;
		if(isLast){
			for(Block b: blocks){
				summary.add(b, isBlock);
				this.useBlock(b);
			}
		}else{
			int st = 0;
			for(Block b: blocks){
				st += (b.getUp()-b.getDown());
			}
			if(st >= 0.8*this.MAXSIZE){
				this.MAXSIZE *= 2;
//				System.out.println(this.MAXSIZE);
				return;
			}
			int n = blocks.length;
			if(n < DISCARD){
				this.MAXSIZE *= 2;
			}else{
				n -= DISCARD/2;
				if (isBlock) {
					System.out.print("Find new Haplotype blocks: ");
				} else {
					System.out.print("Find new Recombination Hotspots:");
				}
				System.out.println(blocks.length);
				for(int i = 0; i < n; i++){
					summary.add(blocks[i], isBlock);
					this.useBlock(blocks[i]);
					limit = blocks[i].getDown();
				}
				for(int i = 0; i <= limit; i++){
					this.snplist.remove(0);
				}
				this.MAXSIZE = SIZE;
			}
		}		
//		System.out.println(this.MAXSIZE);
	}

	private void useBlock(Block b) throws FileNotFoundException, InterruptedException {
		writer.writeline("# BLOCK " + blockNo + "\n");
//		System.out.println("("+this.MAXSIZE+") "+blockNo+": "+b.getUp()+", "+b.getDown());
		if (phased != null) {
			phased.writeline("# BLOCK " + blockNo + "\n");
			b.phase(rc, phased);
		}
		blockNo++;
		if (tag != null) {
			b.setTags(tag.getTags(b.getSnps()));
		}
		if (isCC) {
			b.ccTest();
		}
		writer.write(b);		
	}

	public void add(Snp snp, Summary summary) throws InterruptedException,
			FileNotFoundException {
		snplist.add(snp);
		if (snplist.size() >= MAXSIZE) {
			pro(summary, false);
		}
	}

	public void restart(Summary summary) throws InterruptedException,
			FileNotFoundException {
		// System.out.println(snplist.size());
		pro(summary, true);
		snplist.clear();
	}

	public void close(Summary summary) throws InterruptedException,
			FileNotFoundException {
		pro(summary, true);
		writer.close();
		if (phased != null) {
			phased.close();
		}
	}

	// 多线程版
	private void proMulti(Summary summary, boolean isLast)
			throws InterruptedException, FileNotFoundException {
		return;
//		if (snplist.size() == 0) {
//			return;
//		}
//		int THREADS;
//		if (snplist.size() % SIZE == 0) {
//			THREADS = snplist.size() / SIZE;
//		} else {
//			THREADS = snplist.size() / SIZE + 1;
//		}
//		// System.out.println(snplist.size());
//		Snp[][] snps = new Snp[THREADS][SIZE * 2];
//		Iterator<Snp> iter = snplist.iterator();
//		int[] ends = new int[THREADS];
//		int i, j = 0;
//		for (i = 0; i < THREADS; i++) {
//			for (j = 0; j < SIZE; j++) {
//				if (iter.hasNext()) {
//					snps[i][j] = iter.next();
//				} else {
//					break;
//				}
//			}
//			ends[i] = j;
//		}
//		ends[THREADS - 1] = j;
//		FindBlock[] finds = new FindBlock[THREADS];
//		Thread[] runs = new Thread[THREADS];
//		for (i = 0; i < THREADS; i++) {
//			finds[i] = new FindBlock(rc, isBlock);
//			// System.out.println(snps[i][0].getPosition()+", "+snps[i][ends[i]-1].getPosition());
//			finds[i].setRun(snps[i], 0, ends[i], isLast);
//			runs[i] = new Thread(finds[i]);
//		}
//		for (i = 0; i < THREADS; i++) {
//			runs[i].start();
//		}
//		try {
//			for (i = 0; i < THREADS; i++) {
//				runs[i].join();
//			}
//		} catch (InterruptedException e) {
//			e.printStackTrace();
//		}
//		ArrayList<Block> blocklist = new ArrayList<Block>();
//		Block[] blocks;
//		int[] already = new int[THREADS];
//		int to;
//		// System.out.println(THREADS);
//		for (i = 0; i < THREADS; i++) {
//			blocks = finds[i].getBlocks();
//			if (blocks == null || blocks.length == 0) {
//				already[i] = 0;
//				ends[i] = 0;
//				continue;
//			}
//			already[i] = blocks.length;
//			if (i > 0) {
//				to = blocks[0].getUp();
//				for (j = 0; j < to; j++) {
//					// System.out.println(i+", "+ends[i-1]+", "+j);
//					snps[i - 1][ends[i - 1]++] = snps[i][j];
//				}
//			}
//			for (Block block : blocks) {
//				blocklist.add(block);
//			}
//			ends[i] = 0;
//			if (i < THREADS - 1) {
//				for (j = blocks[blocks.length - 1].getDown(); j < SIZE; j++) {
//					// System.out.println(i+", "+ends[i]+", "+j);
//					snps[i][ends[i]++] = snps[i][j];
//				}
//			}
//		}
//		for (i = 0; i < THREADS - 1; i++) {
//			finds[i] = new FindBlock(rc, isBlock);
//			// System.out.println(snps[i][0].getPosition()+", "+snps[i][ends[i]-1].getPosition());
//			finds[i].setRun(snps[i], 0, ends[i], isLast);
//			runs[i] = new Thread(finds[i]);
//		}
//		for (i = 0; i < THREADS - 1; i++) {
//			runs[i].start();
//		}
//		try {
//			for (i = 0; i < THREADS - 1; i++) {
//				runs[i].join();
//			}
//		} catch (InterruptedException e) {
//			e.printStackTrace();
//		}
//		j = 0;
//		for (i = 0; i < THREADS - 1; i++) {
//			j += already[i];
//			blocks = finds[i].getBlocks();
//			if (blocks == null) {
//				continue;
//			}
//			for (Block block : blocks) {
//				blocklist.add(j++, block);
//			}
//		}
//		blocks = new Block[blocklist.size()];
//		blocklist.toArray(blocks);
//		if (isBlock) {
//			System.out.print("Find new Haplotype blocks: ");
//		} else {
//			System.out.print("Find new Recombination Hotspots:");
//		}
//		System.out.println(blocks.length);
//		if (blocks.length != 0) {
//			for (Block b : blocks) {
//				summary.add(b, isBlock);
//				writer.writeline(">> " + blockNo + "\n");
//				if (phased != null) {
//					phased.writeline(">> " + blockNo + "\n");
//					b.phase(rc, phased);
//				}
//				blockNo++;
//				if (tag != null) {
//					b.setTags(tag.getTags(b.getSnps()));
//				}
//				if (isCC) {
//					b.ccTest();
//				}
//				writer.write(b);
//			}
//		}
//		if (!isLast) {
//			int limit;
//			if (blocks.length == 0) {
//				limit = (SIZE - 1) / 20;
//			} else {
//				limit = blocks[blocks.length - 1].getDown();
//			}
//			limit += SIZE * (THREADS - 1);
//			for (i = 0; i < limit; i++) {
//				snplist.remove(0);
//			}
//			// System.out.println(snplist.size());
//		}
	}

}
