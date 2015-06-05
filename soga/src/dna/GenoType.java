package dna;

public class GenoType implements Comparable<GenoType>{
	private int[] alleles;
	private int sample;
	
	public GenoType(int[] alleles, int sampleNo){
		this.setAlleles(alleles);
		this.setSample(sampleNo);
	}

	@Override
	public int compareTo(GenoType arg0) {
		int[] another = arg0.getAlleles();
		if(this.alleles.length > another.length){
			return 1;
		}else if(this.alleles.length < another.length){
			return -1;
		}else{
			int i, n = this.alleles.length;
			for(i = 0; i < n; i++){
				if(this.alleles[i] != another[i]){
					return(this.alleles[i]-another[i]);
				}
			}
			return 0;
		}
	}
	
	public boolean equals(GenoType another){
		int[] a = another.getAlleles();
		if(a.length != this.alleles.length){
			return false;
		}
		int i, n = this.alleles.length;
		for(i = 0; i < n; i++){
			if(this.alleles[i] != a[i]){
				return false;
			}
		}
		return true;
	}

	public int[] getAlleles() {
		return alleles;
	}

	public void setAlleles(int[] alleles) {
		this.alleles = alleles;
	}

	public int getSample() {
		return sample;
	}

	public void setSample(int sample) {
		this.sample = sample;
	}

}
