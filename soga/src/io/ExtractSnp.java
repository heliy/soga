package io;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import dna.Snp;
import exceptions.AlleleException;
import exceptions.GenoTypeException;
import exceptions.SnpContainsException;
import parameter.Base;
import parameter.Setting;

public class ExtractSnp {
	private String snpfile;
//	private int fileNo;
	private File f;
	private Scanner in;
	private String splitFlag;
	private int head;
	private int samples;
	private String thisline;
	private Setting rc;
	private boolean ignoreGenotypeException;
	
	public ExtractSnp(Setting rrc) throws FileNotFoundException{
		rc = rrc;
		snpfile = rc.getSnpFile();
//		fileNo = 0;
		f = new File(snpfile);
		in = new Scanner(f);
		splitFlag = rc.getFileSplit();
		head = rc.getHEAD();
		samples = rc.getSAMPLES();
	}

	public Snp nextSnp() throws FileNotFoundException, AlleleException, GenoTypeException, SnpContainsException{
		if(in.hasNext()){
			String line = in.nextLine();
			if(line.charAt(0) == '#'){
				return this.nextSnp();
			}else{
				thisline = line;
				String[] parts = line.split(splitFlag);
//				System.out.println(samples);
				if(parts.length != (samples+head)){
					throw new SnpContainsException(snpfile, thisline, samples+head);
				}
				int i;
				String types[] = new String[samples];
				for(i = 0; i < samples; i++){
					types[i] = parts[i + head];
				}
				int[] alleles;
				try {
					alleles = parseAlleles(parts[1]);
				} catch (AlleleException e1) {
					throw new AlleleException(snpfile, thisline, parts[1]);
				}
				try {
					return new Snp(parts[0], alleles[0], alleles[1], parts[2], parts[3], types, rc);
				} catch (GenoTypeException e) {
					if(ignoreGenotypeException){
						return this.nextSnp();
					}else{
						System.out.println("This line have invalide genotype: ");
						System.out.println(thisline);
						System.out.print("Continue?[Y/N]:");
						@SuppressWarnings("resource")
						Scanner in = new Scanner(System.in);
						if(in.next().toLowerCase().charAt(0) == 'y'){
							return this.nextSnp();
						}else{
							String geno = e.getGeno();
							throw new GenoTypeException(snpfile, thisline, geno);							
						}
					}
				}
			}
		}else{
			in.close();
//			fileNo++;
//			if(fileNo == snpfiles.length){
//				return null;
//			}else{
//				f = new File(snpfiles[fileNo]);
//				in = new Scanner(f);
//				System.out.println("New File: "+snpfiles[fileNo]);
//				return this.nextSnp();
//			}
			return null;
		}
	}
	
	public String getThisline() {
		return thisline;
	}
		
	private int[] parseAlleles(String ts) throws AlleleException {
		// parse string to two alleles
		// A/C -> A C
		// TODO: error control and insert deletion
		Base base = new Base();
		int[] alleles = new int[2];
		if(ts.length() != 3 || ts.charAt(1) != '/'){
			throw new AlleleException();			
		}
		alleles[0] = base.baseNo(ts.charAt(0));
		alleles[1] = base.baseNo(ts.charAt(2));
		if(alleles[0] < 0 || alleles[1] < 0){
			throw new AlleleException();
		}
		if (alleles[0] > alleles[1]) { // 保证 A 的序号小于 B
			int tmp = alleles[0];
			alleles[0] = alleles[1];
			alleles[1] = tmp;
		}
		return alleles;
	}

}