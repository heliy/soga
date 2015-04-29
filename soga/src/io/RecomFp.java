package io;

import java.io.FileNotFoundException;

import parameter.Setting;

public class RecomFp extends FileOutput {

	public RecomFp(Setting rc) throws FileNotFoundException {
		super(rc, "RECOM");
	}

}
