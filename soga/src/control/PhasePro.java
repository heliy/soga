package control;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.CountDownLatch;

import com.sun.org.apache.bcel.internal.generic.Select;

import calculate.Phase;
import dna.GenoType;
import dna.HaploType;
import dna.PhasedRange;
import dna.Snp;
import io.PhasedFp;
import parameter.Setting;
import parameter.Summary;
import util.AllotThreads;
import util.MultiThreads;
import util.SelectH;

public class PhasePro {
	private int THREADS;
	private PhasedFp writer;
	private Phase[] phases;
	private int window;
	private MultiThreads<Phase>[] runs;
	private ArrayList<Snp> snplist;
	private AllotThreads allot;
	
	private Setting rc;
	private Summary summary;

	public PhasePro(Setting rc, Summary summary) throws FileNotFoundException {
		this.rc = rc;
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
		this.allot = new AllotThreads();
		this.THREADS = rc.getThreads();
		this.phases = new Phase[n];
		this.window = rc.getPhaseWINDOW();
		this.snplist = new ArrayList<Snp>();
		this.runs = new MultiThreads[this.THREADS];
		int[] threadsNum = this.allot.allot(n, this.THREADS);
		this.rc = rc;

		for (i = 0; i < n; i++) {
			this.phases[i] = new Phase(i, rc);
		}
		for(j = 0, r = 0; r < this.THREADS; r++){
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
				- this.snplist.get(0).getPosition() >= window);
	}
	
	private PhasedRange run() throws InterruptedException, FileNotFoundException {
		if(this.snplist.size() == 0){
			return null;
		}
//		System.out.println(this.snplist.size());
		Snp[] snps = new Snp[snplist.size()];
		snplist.toArray(snps);
		
		int i, j, r, n = this.phases.length, m = snps.length;
		HaploType[][] hts = new HaploType[n][2];
		int[][] alleles = new int[n][m];

		// 找出这段的基因型的所有类型
		// 因为大体是连锁不平衡快内的不同样本的基因型会有重复
		for(i = 0; i < m; i++){
			int[][] types = snps[i].getTypes();
			for(j = 0; j < n; j++){
				alleles[j][i] = types[j][0]+types[j][1];
			}
		}
		GenoType[] genotypes = new GenoType[n];
		for(i = 0; i < n; i++){
			genotypes[i] = new GenoType(alleles[i], i);
		}
		Arrays.sort(genotypes);  // 相同的基因型会相进排列
		boolean[] needCal = new boolean[n];
		int[] refs = new int[n];  // 第i个样本参考 standalone[refs[i]]
		int news = 0;
		needCal[0] = true;
		refs[genotypes[0].getSample()] = news++;
		for(i = 1; i < n; i++){
			if(genotypes[i].equals(genotypes[i-1])){
				refs[genotypes[i].getSample()] = refs[genotypes[i-1].getSample()];  // 第i个样本参考 standalone[refs[i]]
				needCal[i] = false;
			}else{
				refs[genotypes[i].getSample()] = news++;
				needCal[i] = true;
			}
		}
//		System.out.println(news);
		SelectH select = new SelectH(rc, snps);
		Phase[] standalone = new Phase[news];
		for(i = 0, j = 0; i < n; i++){
			if(needCal[i]){
//				System.out.println(j);
				standalone[j++] = new Phase(i, select.getH(), rc); 
			}
		}
		
		// 对这些“独特”的基因型分相
		int[] multis = this.allot.allot(news, this.THREADS);
		for(i = 0; i < news; i++){
			standalone[i].setRun(snps, new CountDownLatch(1));
		}
		for(r = 0, j = 0; r < this.THREADS; r++){
			Phase[] phs = new Phase[multis[r]];
			for(i = 0; i < multis[r]; i++){
				phs[i] = standalone[j+i];
			}
			j += multis[r];
			this.runs[r] = new MultiThreads<Phase>(phs);
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
		
		// 剩下的找到有相同的基因型的样品的分相情况 并接上头
		for(i = news = 0; i < n; i++){
			if(this.phases[i].hasCleared()){
				needCal[i] = false; // 没有头 直接加尾
				this.phases[i].setHt(snps, null, standalone[refs[i]].getHts());
			}else{
				needCal[i] = true;
				news++;
			}
		}
		if(news > 0){
			Phase[] boundPhase = new Phase[news];
			multis = this.allot.allot(news, this.THREADS);
			for(i = 0, j = 0; i < n; i++){
				if(needCal[i]){
					Snp[] bound = this.phases[i].getBound(snps); // 能确定怎么接的SNP区域
					boundPhase[j] = new Phase(i, rc);
					boundPhase[j++].setRun(bound, new CountDownLatch(1));
				}
			}
			for(r = 0, j = 0; r < this.THREADS; r++){
				Phase[] phs = new Phase[multis[r]];
				for(i = 0; i < multis[r]; i++){
					phs[i] = boundPhase[j+i];
				}
				j += multis[r];
				this.runs[r] = new MultiThreads<Phase>(phs);
			}
			latch = new CountDownLatch(this.THREADS);
			threads = new Thread[this.THREADS];
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
			// 加尾
			for(i = j = 0; i < n; i++){
				if(needCal[i]){
					this.phases[i].setHt(snps, boundPhase[j++].getHts(), standalone[refs[i]].getHts());
				}
			}
		}
		
		for(i = 0; i < n; i++){
			hts[i] = this.phases[i].getHts();
		}
		
		PhasedRange range = new PhasedRange(snps, hts);
		this.writer.write(range);
		return range;
	}
	
}
