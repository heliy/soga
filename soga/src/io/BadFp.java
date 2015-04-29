package io;

import java.io.FileNotFoundException;

import parameter.Setting;

public class BadFp extends FileOutput {

	public BadFp(Setting rc) throws FileNotFoundException {
		super(rc, "BADDATA");
	}

}
