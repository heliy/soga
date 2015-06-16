package io;

import java.io.FileNotFoundException;

import parameter.Setting;

public class SNPFp extends FileOutput {
	static String head = "# rs#\talleles\tchr\tposition\t"
			+ "A\tA_count\tA_freq\tB\tB_count\tB_freq\t"
			+ "N_count\tDeletion\tAA_count\tAA_freq\t"
			+ "AB_count\tAB_freq\tBB_count\tBB_freq\t" 
			+ "NN_count\tDeletion\tMAF\tH-W\t";

	public SNPFp(Setting rc) throws FileNotFoundException {
		super(rc, "SNP", head);
		if(rc.doCC()){
			super.writeline("A_pvalue\tB_pvalue\tOR\tOR_CI_Low\tOR_CI_High\tPASS\n");
		}else{
			super.writeline("PASS\n");
		}
	}

}
