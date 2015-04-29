package io;

import java.io.FileNotFoundException;
import java.util.Scanner;

import parameter.Setting;

public class sampleQC {
	double nnRatio;
	int[] samples;
	int len;
	long total;
	
    public sampleQC(Setting rc){
    	nnRatio = rc.getNNRatio();
    	len = rc.getSAMPLES();
    	samples = new int[len];
    	total = 0;
    }
    
    public void update(String[] types){
    	int i;
    	for(i = 0; i < len; i++){
    		if(isNN(types[i])){
    			samples[i]++;
    		}
    	}
    	total++;
    }
    
    public int[] valids(){
    	int stand = (int)(nnRatio * len);
    	int valid = 0, i;
    	for(i = 0; i < len; i++){
    		if(samples[i] >= stand){
    			samples[i] = 0;
    		}else{
    			samples[i] = 1;
    			valid++;
    		}
    	}
    	int[] vas = new int[valid];
    	for(i = 0, valid = 0; i < len; i++){
    		if(samples[i] == 1){
    			vas[valid++] = i;
    		}
    	}
    	return vas;
    }

    public long totalSnp(){
    	return total;
    }
    
	private boolean isNN(String string) {
		if(string.toUpperCase().equals("NN")){
			return true;
		}else{
			return false;
		}
	}
    
//	private void setQCSamples(Setting rc) throws FileNotFoundException{
//		String line;
//		int i;
//		String[] types = new String[len];
//		int head = rc.getHEAD();
//		
//		while(in.hasNextLine()){
//			line = in.nextLine();
//			String[] parts = line.split(rc.getFileSplit());
//			for(i = 0; i < len; i++){
//				types[i] = parts[i + head];
//			}
//			qc.update(types);
//		}
//		rc.setValids(qc.valids());
//		valids = rc.getValids();
//	    valid = rc.getValids().length;
//		rc.setSNPS(qc.totalSnp());
//		in.close();
//		in = new Scanner(f);
//	}


}
