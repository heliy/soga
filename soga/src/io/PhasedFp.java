package io;

import java.io.FileNotFoundException;

import parameter.Setting;
import population.Sample;

public class PhasedFp extends FileOutput {

	public PhasedFp(Setting rc) throws FileNotFoundException {
		super(rc, "PHASED");
		this.addhead(rc.getSamples());
	}
	
	public void addhead(Sample[] samples){
		this.write("# rs\talleles\tchromosome\tposition");
		for(Sample s: samples){
			if(s.isPass()){
				this.write("\t"+s.getName());
			}
		}
		this.write("\n");
	}

}
