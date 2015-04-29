package io;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import dna.Snp;
import parameter.Setting;

public class FileIn {
	String name;
	File f;
	Scanner in;
	String splitFlag;
	QCSample qc;
	int valid;
	int[] valids;
	int head;
	int len;
	
	public FileIn(String filename, Setting rc) throws FileNotFoundException{
		name = filename;
		f = new File(name);
		in = new Scanner(f);
		splitFlag = rc.getFileSplit();
		qc = new QCSample(rc);
		head = rc.getHEAD();
		len = rc.getSAMPLES();
	}

	private void setQCSamples(Setting rc) throws FileNotFoundException{
		String line;
		int i;
		String[] types = new String[len];
		int head = rc.getHEAD();
		
		while(in.hasNextLine()){
			line = in.nextLine();
			String[] parts = line.split(rc.getFileSplit());
			for(i = 0; i < len; i++){
				types[i] = parts[i + head];
			}
			qc.update(types);
		}
		rc.setValids(qc.valids());
		valids = rc.getValids();
	    valid = rc.valids.length;
		rc.setSNPS(qc.totalSnp());
		in.close();
		in = new Scanner(f);
	}

	
	public Snp nextSnp(Setting rc){
		if(in.hasNext()){
			String line = in.nextLine();
			if(line.charAt(0) == '#'){
				return this.nextSnp(rc);
			}else{
				String[] parts = line.split(splitFlag);
				int i;
				String types[] = new String[valid];
				for(i = 0; i < valid; i++){
					types[i] = parts[valids[i] + head];
				}
				return new Snp(parts[0], parts[1], parts[2], parts[3], types, rc);
			}
		}else{
			return null;
		}
	}

}