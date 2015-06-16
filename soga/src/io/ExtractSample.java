package io;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import exceptions.SampleException;
import population.Sample;

public class ExtractSample {
	public Sample[] getSample(int n){
		Sample[] samples = new Sample[n];
		for(int i = 1; i <= n; i++){
			samples[i-1] = new Sample(i);
		}
		return samples;
	}
	
	@SuppressWarnings("resource")
	public Sample[] getSample(String filename) throws FileNotFoundException, SampleException{
		File f = new File(filename);
		Scanner in = new Scanner(f);
		int total = 0;
		while(in.hasNext()){
			String line = in.nextLine();
			if(line.charAt(0) == '#'){
				continue;
			}
			total++;
		}
		in.close();
		Sample[] samples = new Sample[total];
		total = 0;
		in = new Scanner(f);
		while(in.hasNext()){
			String line = in.nextLine();
			if(line.charAt(0) == '#'){
				continue;
			}
			total++;
			String[] parts = line.split("\t");
			if(parts.length != 2){
				throw new SampleException(line);
			}
			samples[total-1] = new Sample(total, parts[0]);
			samples[total-1].setStatu(Integer.parseInt(parts[1]));
		}
		in.close();
		return samples;
	}

}
