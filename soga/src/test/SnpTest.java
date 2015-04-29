package test;

import dna.Snp;
import parameter.Setting;

public class SnpTest {

	public static void main(String[] args) {
		String rs = "rs62224610";
		String alleles = "G/C";
		String chr = "chr22";
		String position = "16051347";
		String line = "CG GC GG GC GC GG GG GG GG GC GC GG GG GG GG GG GG GC GG GC GG GG GC CC GC GG CC GC GG GC GC GG GC GG GC GG GC GC GG GG GC GC GC GG GG GG GG GG CC GC GC GC GG GC GG GC GG GG GG GC GG GG GG GG CG GC GG GG GC GG GG GG GC GC GC GG GG GC GC GC GG GC GG GG GG CG CG GG GG GG CG CG GC GG GG CG CG CC GG GG GG GG GG GG GG GG GC GC GC GC GC GG GC GG GC GG GG GG CG GC GC GC GC GC GG GG GG CC GG GG GC GC GG GG GG GG GG GG GG GG GC GG GC GG GG GC GG GG GG GC GG GG GG GG GG CG GC GG GG GC GG GG GC GG GC GG GG GG GC GC GC CG GG GC GG GC GG GC GC GG GG GG CC CC GG GG GG GG GC CC GC GG GC GC GG GC GC CG GG GG GG GG GC GC CG GG GC GG GC GC GG GG CG GG GG GG GC GG CG GG GC GG GG GG GC GG GC GC GC GC GG GG GC GG GG GC GG GC GG GG GG GG CG GG GG GG GG GG GG GG GG GG GC GG GG GC GG GC GG GG GC CC GC GC GC GC GG GC GC GC GG CC GC GG GG GC GC GC CG GG CC GG GC GG GG GC";
        Setting rc = new Setting();
        String []parts = line.split(" ");
        Snp snp = new Snp(rs, chr, position, alleles, parts, rc);
        snp.display();
	}

}
