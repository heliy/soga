package io;

import java.io.FileNotFoundException;

import parameter.Setting;

public class BlockFp extends FileOutput {
	public BlockFp(Setting rc)
			throws FileNotFoundException {
		super(rc, "BLOCK");
	}

}
