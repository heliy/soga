package io;

import java.io.FileNotFoundException;

import parameter.Setting;

public class SampleQCFp extends FileOutput {
	private static String label = "SAMPLE";
	private static String head = "# NAME\tSTATU\tDeletion\tPASS\n";
	public SampleQCFp(Setting rc) throws FileNotFoundException{
		super(rc, label, head);
	}
}
