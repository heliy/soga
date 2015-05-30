package dna;

import parameter.Base;

public class PhasedRange {
	private Snp[] snps;
	private HaploType[][] hts;
	
	public PhasedRange(Snp[] snps, HaploType[][] hts){
		this.snps = snps;
		this.hts = hts;
	}
	
	public void merge(PhasedRange nextRange){
		Snp[] nextSnps, newSnps;
		HaploType[][] nextHts = nextRange.getHaploTypes();
		nextSnps = nextRange.getSnps();
		newSnps = new Snp[this.snps.length+nextSnps.length];
		
		int i, j, n;
		for(i = 0, n = this.snps.length; i < n; i++){
			newSnps[i] = this.snps[i];
		}
		for(j = 0, n = nextSnps.length; j < n; j++){
			newSnps[i+j] = nextSnps[j];
		}
		this.snps = newSnps;
		
		for(i = 0, n = this.hts.length; i < n; i++){
			hts[i][0].extend(nextHts[i][0]);
			hts[i][1].extend(nextHts[i][1]);
		}
	}
	
	public Snp[] getSnps() {
		return snps;
	}

	public HaploType[][] getHaploTypes() {
		return hts;
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		int i, j, num = snps.length, samples = hts.length;
		Base base = new Base();
		for(i = 0; i < num; i++){
			sb.append(snps[i].getChr()+"\t"+snps[i].getPosition());
			for(j = 0; j < samples; j++){
//				System.out.println(haplotypes.length);
				sb.append("\t"+base.getBase(hts[j][0].getAlleles()[i])+"|"+base.getBase(hts[j][1].getAlleles()[i]));
			}
			sb.append("\n");
		}
		return sb.toString();
		
	}

}
