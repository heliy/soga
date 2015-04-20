package test;
//
//import io.ExtractSnp;
//
//import java.io.FileNotFoundException;
//
//import calculate.FindBlock;
//import calculate.Phase;
//import dna.Block;
//import dna.HaploType;
//import dna.Snp;
//import parameter.Setting;
//
public class shapeitTest {
//	static String lit = "/home/hly/pro/resoff/ld/test-data/lit";
//	static String chr22 = "/home/hly/pro/resoff/ld/test-data/chr22_merged_Asia.txt";
//	static int litSample = 286;
//
//	public static void main(String[] args) throws FileNotFoundException {
//		    Setting rc = new Setting();
//			rc.setInputfile(chr22);
//			rc.setSAMPLES(litSample);
//			rc.setBLOCK();
//			rc.setPHASE();
//			ExtractSnp next = new ExtractSnp(rc);
//			Snp[] snps = new Snp[200];
//			int i, j, k, h;
//			FindBlock find = new FindBlock(rc, true);
//			for(i = 0; i < 5; i++){
//  			    for(j = 0; j < 200; j++){
//				   snps[j] = next.nextSnp();
//			    }
//			}
//			for(i = 5; i < 200; i++){
//  			    for(j = 0; j < 200; j++){
//				   snps[j] = next.nextSnp();
//			    }
//    			Block[] blocks = find.migpp(snps, false);
//			    System.out.println(i+", "+j);
//	    		for(k = 0; k < blocks.length; k++){
////			    k = 2;
//				    System.out.println(i+", "+j+", "+k);
//		   		    Block block = blocks[k];
//		   		    for(h = 0; h < litSample; h++){
////		   		    h = 1;
//		   		    	if(h%50 == 0){
//   					    System.out.println(i+", "+j+", "+k+", "+h);
//		   		    	}
//					    Phase phase = new Phase(block.getSnps(), rc);
//					    phase.phase(h);
//		   		    }
////				    System.out.println(hts[0].toString());
////				    System.out.println(hts[1].toString());				
//    			}
//			}
//	}
//
}
