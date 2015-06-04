package control;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.concurrent.CountDownLatch;

import calculate.Phase;
import dna.HaploType;
import dna.PhasedRange;
import dna.Snp;
import io.PhasedFp;
import parameter.Setting;
import parameter.Summary;
import util.MultiThreads;

public class PhasePro {
	private int THREADS;
	private PhasedFp writer;
	private Phase[] phases;
	private int window;
	private MultiThreads<Phase>[] runs;
	private ArrayList<Snp> snplist;
	
	private Summary summary;

	public PhasePro(Setting rc, Summary summary) throws FileNotFoundException {
		this.writer = new PhasedFp(rc);
		this.summary = summary;
		this.summary.add(this.writer);
		this.summary.setPhasedFile(this.writer.getFileName());
		this.alloc(rc);
	}

	public PhasePro(Setting rc, PhasedFp writer) throws FileNotFoundException {
		this.writer = writer;
		this.alloc(rc);
	}
	
	public PhasedRange add(Snp snp) throws InterruptedException, FileNotFoundException {
		snplist.add(snp);
		if (this.needCal()) {
			PhasedRange range = this.run();
			writer.write(range);
			this.snplist.clear();
			return range;
		} else {
			return null;
		}
	}

	public PhasedRange restart() throws FileNotFoundException, InterruptedException {
		PhasedRange range = this.run();
		this.snplist.clear();
		for (Phase phase: this.phases) {
			phase.clearHistory();
		}
		return range;
	}

	public PhasedRange close(boolean writerNeed) throws InterruptedException, FileNotFoundException {
		PhasedRange range = this.restart();
		if (writerNeed) {
			writer.close();
		}
		return range;
	}

	@SuppressWarnings("unchecked")
	private void alloc(Setting rc) throws FileNotFoundException {
		int i, j, r, n = rc.getSAMPLES();
		this.THREADS = rc.getThreads();
		this.phases = new Phase[n];
		this.window = rc.getPhaseWINDOW();
		this.snplist = new ArrayList<Snp>();
		int[] threadsNum = new int[this.THREADS];

		for (i = 0; i < n; i++) {
			this.phases[i] = new Phase(i, rc);
		}
		
		this.runs = new MultiThreads[this.THREADS];
		for(i = 0; i < n; i++){
			threadsNum[i%this.THREADS]++;
		}
		j = 0;
		for(r = 0; r < this.THREADS; r++){
//			System.out.println(threadsNum[r]);
			Phase[] phs = new Phase[threadsNum[r]];
			for(i = 0; i < phs.length; i++){
				phs[i] = this.phases[j+i];
			}
			j += phs.length;
			this.runs[r] = new MultiThreads<Phase>(phs);
		}
	}

	private boolean needCal() {
		return (this.snplist.get(this.snplist.size() - 1).getPosition()
				- this.snplist.get(0).getPosition() > window);
	}

	private PhasedRange run() throws InterruptedException, FileNotFoundException {
		if(this.snplist.size() == 0){
			return null;
		}
		System.out.println(this.snplist.size());
		Snp[] snps = new Snp[snplist.size()];
		snplist.toArray(snps);
		
		int i, n = this.phases.length;
		HaploType[][] hts = new HaploType[n][2];
		
		for(i = 0; i < n; i++){
			this.phases[i].setRun(snps, new CountDownLatch(1));
		}
		
		CountDownLatch latch = new CountDownLatch(this.THREADS);
		Thread[] threads = new Thread[this.THREADS];
		for(i = 0; i < this.THREADS; i++){
			this.runs[i].setRun(latch);
		}
		for(i = 0; i < this.THREADS; i++){
			threads[i] = new Thread(this.runs[i]);
		}
		for(i = 0; i < this.THREADS; i++){
			threads[i].start();
		}
		latch.await();
		
		for(i = 0; i < n; i++){
			hts[i] = this.phases[i].getHts();
		}
		
		PhasedRange range = new PhasedRange(snps, hts);
		this.writer.write(range);
		return range;

	}

}
