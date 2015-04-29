package io;

import java.io.FileNotFoundException;

import parameter.Setting;

public class FilterFp extends FileOutput {

	public FilterFp(Setting rc) throws FileNotFoundException {
		super(rc, "FILTER");
	}

}
