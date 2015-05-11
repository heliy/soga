package io;

import java.io.FileNotFoundException;

import parameter.Setting;

public class LDFp extends FileOutput {

	static String head = "# chr\tSNP1_rs#\tSNP1_position"
			+ "\tSNP2_rs#\tSNP2_position\t"
			+ "dprime\tr2\tci_low\tci_high\tlod\ttype\n";

	public LDFp(Setting rc) throws FileNotFoundException {
		super(rc, "LD", head);
	}

}
