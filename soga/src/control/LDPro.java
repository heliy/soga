package control;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;

import calculate.PairwiseCal;
import dna.LD;
import dna.Snp;
import parameter.Setting;
import parameter.Summary;
import io.LDFp;

public class LDPro {
	private int maxDistance;
	private long begin; // 最先的snp的位置

	private ArrayList<Snp> snps;
	private LDFp writer;
	private PairwiseCal cal;

	public LDPro(Setting rc, Summary summary) throws FileNotFoundException {
		maxDistance = rc.getDISTANCE();
		begin = 0;
		snps = new ArrayList<Snp>();
		writer = new LDFp(rc);
		cal = new PairwiseCal(rc);
		summary.add(writer);
	}
	
	

	public void add(Snp snp, Summary summary) {
		long ima = snp.getPosition();
		while (ima - begin > maxDistance) {  // 保证维护的snps内的距离在范围内
			if (snps.isEmpty()) {
				begin = ima;
				break;
			} else {
				snps.remove(0);
			}
			if(snps.isEmpty()){
				begin = ima;
				break;
			}else{
			    begin = snps.get(0).getPosition();
			}
		}
		Iterator<Snp> iter = snps.iterator();
		while (iter.hasNext()) {
			Snp mae = iter.next();
			String s = mae.getChr() + "\t" + mae.getRs() + "\t" + mae.getPosition() + "\t"
					+ snp.getRs() + "\t" + snp.getPosition() + "\t";
			LD ld = mae.getLD(cal, snp);
			summary.add(ld);
			writer.write(s + ld.toString());
		}
		snps.add(snp);
	}
	
	public void restart(){
		begin = 0;
		snps.clear();
	}

	public void close() {
		writer.close();
	}

}
