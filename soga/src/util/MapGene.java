package util;

import io.FileOutput;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import dna.Block;
import dna.Snp;
import parameter.Setting;
import parameter.Summary;

public class MapGene {
	private String current = null;
	private String chr = null;
	
	private File f;
	private Scanner in;
	
    private FileOutput snps = null;
    private FileOutput blocks = null;
    
	public MapGene(Setting rc, Summary summary) throws FileNotFoundException{
		this.f = new File(rc.getGeneFile());
		this.in = new Scanner(f);
		this.snps = new FileOutput(rc, "SNPMAP");
		summary.add(snps);
		if(rc.doPHASE()){
			this.blocks = new FileOutput(rc, "BLPCKMAP");
			summary.add(blocks);
		}
	}
	
	public void add(Snp snp){
		// TODO
	}
	
	public void add(Block b){
		// TODO
	}
	
	public void close(){
		in.close();
		this.snps.close();
		if(this.blocks != null){
			this.blocks.close();
		}
	}
	
	
}
