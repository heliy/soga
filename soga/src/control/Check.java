package control;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;

import io.ExtractSnp;
import dna.Snp;
import exceptions.AlleleException;
import exceptions.GenoTypeException;
import exceptions.OrderException;
import exceptions.SnpContainsException;
import parameter.Setting;
import population.Sample;

public class Check {
	// 检查文件 和 样本的质量控制
	private boolean[] sampleQC;
	private int SAMPLE;
	private double ratio;

	private String chr = null;
	private ArrayList<String> chrlist = null;
	private Snp mae = null;

	private long total;
	private long[] nns;
	private Sample[] samples;

	private ExtractSnp input;

	public Check(Setting rc) throws FileNotFoundException {
		rc.passAllSample();
		this.SAMPLE = rc.getSAMPLES();
		this.sampleQC = new boolean[this.SAMPLE];
		this.total = 0;
		this.nns = new long[this.SAMPLE];
		this.ratio = rc.getSampleRatio();
		this.input = new ExtractSnp(rc, rc.isPhased());
		this.chrlist = new ArrayList<String>();
		this.samples = rc.getSamples();
	}

	public Sample[] run() throws FileNotFoundException, AlleleException,
			GenoTypeException, SnpContainsException, OrderException {
		Snp snp = input.nextSnp();
		while (snp != null) {
			if(chr != null){
				if (snp.getChr().equals(chr)) {
					if (snp.getPosition() <= this.mae.getPosition()) {
						throw new OrderException(this.mae, snp);
					}
				} else {
					Iterator<String> iter = this.chrlist.iterator();
					while (iter.hasNext()) {
						if (iter.next().equals(snp.getChr())) {
							throw new OrderException(this.mae, snp);
						}
					}
					this.chrlist.add(snp.getChr());
					System.out.println("New Chromosome: " + snp.getChr());
				}
			}else{
				System.out.println("New Chromosome: " + snp.getChr());				
			}
			this.chr = snp.getChr();
			this.mae = snp;
			this.total++;
			int[] genos = snp.getGenos();
			for (int i = 0; i < this.SAMPLE; i++) {
				if (genos[i] == -1) {
					this.nns[i]++;
				}
			}
			snp = input.nextSnp();
		}
		long limit = (long) (this.total * this.ratio);
		for (int i = 0; i < this.SAMPLE; i++) {
			this.sampleQC[i] = (this.nns[i] < limit);
			this.samples[i].setPass(this.sampleQC[i]);
			this.samples[i].setRatio((float) this.nns[i] / this.total);
		}
		return this.samples;
	}


}
