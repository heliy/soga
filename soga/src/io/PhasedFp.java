package io;

import java.io.FileNotFoundException;

import parameter.Setting;

public class PhasedFp extends FileOutput {

	public PhasedFp(Setting rc) throws FileNotFoundException {
		super(rc, "PHASED");
	}

}
