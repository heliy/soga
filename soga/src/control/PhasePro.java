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

public class PhasePro {
	private int THREADS;
	private PhasedFp writer;
	private Phase[] phases;
	private int window;
	private int[] rounds;
	private int[] threadsPerRounds;
	private Thread[] runs;
	private ArrayList<Snp> snplist;
	
	private Summary summary;

	public PhasePro(Setting rc, Summary summary) throws FileNotFoundException {
		this.writer = new PhasedFp(rc);
		this.summary = summary;
		this.summary.add(this.writer);
		this.alloc(rc);
	}

	public PhasePro(Setting rc, PhasedFp writer) {
		this.writer = writer;
		this.alloc(rc);
	}
	
	public PhasedRange add(Snp snp) throws InterruptedException {
		if(this.snplist.size() != 0 && snp.getChr() != this.snplist.get(0).getChr()){
			PhasedRange range = this.close(false);
			this.snplist.add(snp);
			return range;
		}
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

	public void restart() {
		this.snplist.clear();
		for (Phase phase: this.phases) {
			phase.clearHistory();
		}
	}

	public PhasedRange close(boolean writerNeed) throws InterruptedException {
		PhasedRange range = this.run();
		this.restart();
		if (writerNeed) {
			writer.close();
		}
		return range;
	}

	private void alloc(Setting rc) {
		int i, j, r, n = rc.getSAMPLES();
		this.THREADS = rc.getThreads();
		this.phases = new Phase[n];
		this.rounds = new int[n];
		this.window = rc.getPhaseWindow();
		this.snplist = new ArrayList<Snp>();

		for (i = 0; i < n; i++) {
			this.phases[i] = new Phase(i, rc);
		}
		
		for(i = 0, r = 0, j = 0; i < n; i++){
			this.rounds[i] = r;
			j++;
			if(j == THREADS){
				r++;
				j = 0;
			}
		}
		this.threadsPerRounds = new int[r+1];
		for(i = 0; i <= r; i++){
			this.threadsPerRounds[i] = this.THREADS;
		}
		if(j != this.THREADS && j != 0){
			this.threadsPerRounds[r] = j;
		}
		
//		for(i = 0; i < n; i++){
//			System.out.print(this.rounds[i]+", ");
//		}
//		System.out.println();
//		for(int s : this.threadsPerRounds){
//			System.out.print(s+", ");
//		}
//		System.out.println();
		
	}

	private boolean needCal() {
		return (this.snplist.get(this.snplist.size() - 1).getPosition()
				- this.snplist.get(0).getPosition() > window);
	}

	private PhasedRange run() throws InterruptedException {
		Snp[] snps = new Snp[snplist.size()];
		snplist.toArray(snps);
		
		int i, j, k = 0, r = 0, n = this.phases.length;
		int maxR = this.rounds[n-1];
		HaploType[][] hts = new HaploType[n][2];
		for(i = 0; i <= maxR; i++){
			r = this.threadsPerRounds[i];
			CountDownLatch latch = new CountDownLatch(r);
			k = i*this.THREADS;
			this.runs = new Thread[r];
			
			for(j = 0; j < r; j++){
				this.phases[k+j].setRun(snps, latch);
				this.runs[j] = new Thread(this.phases[k+j]);
			}
			
			for(j = 0; j < r; j++){
				runs[j].start();
			}
			latch.await();
			
			for(j = 0; j < r; j++){
				hts[k+j] = this.phases[k+j].getHts();
			}
		}
		PhasedRange range = new PhasedRange(snps, hts);
		this.writer.write(range);
		return range;

	}

}
