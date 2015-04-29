package io;

import parameter.Setting;

public class QCSample {
	double nnRatio;
	int[] samples;
	int len;
	long total;
	
    public QCSample(Setting rc){
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
    
    
}
