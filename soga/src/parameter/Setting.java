package parameter;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Scanner;

import exceptions.ArgsException;
import exceptions.SampleStatusException;

public class Setting {
	private int threads = 2; // 线程数

	// 文件数据
	private int HEAD = 4;
	private int SAMPLES = 0;
//	private int[] valids;
	
	private String fileSplit = "\t";

	// QC
	private double NNRatio = 0.05;
	private double MAF = 0.01;
	private double HWE = 0.001;
	private double MAXN = 0.8;

	//LD
	private int DISTANCE = 25000;
	// 这三个暂不作为用户参数设定
	private double LDCL = 0.7;
	private double LDCU = 0.98;
	private double EHRCU = 0.9;
	
	// Block
	private int HTWINDOW = 4000;
	private int RHWINDOW = 200;
//	private int MAXSIZE = 2000;
//	private int WIN = -1;
	private double D = 0.95; //暂不作为用户参数设定
	
	// Phase 暂不作为用户参数设定
	private int B = 3;
	private int K = 200;
	// 窗口大小， 参数
	private int phaseWINDOW = 2000000;
	
	// case-control
	private int[] status;  // -1 未知 0 正常 1对照 
	private int cases;
	private int controls;
	private int unknowns;
	private double[] freqs;
	
	// 全SNP分相
	
	// IO
	private String snpFile = null;
	private String sampleFile = null;
	private String output = null;
	private String phasedFile = null;
	private double htminratio = 0.05;
	private String geneFile = null;
	
	// 处理
	private boolean FULL = false;
	
	private boolean CHECK = false;
	private boolean FILTER = false;
	private boolean BADDATA = false;

	private boolean LD = false;
	private boolean BLOCK = false;
	private boolean RECOM = false;
	private boolean PHASE = false;
	private boolean TAG = false;
	private boolean CC = false;
	
	private boolean MAP = false;
	
	// 运行
	private boolean ignoreGenotypeException = false;
	private boolean silence = false;

	public boolean parseArgs(String[] args) throws ArgsException, FileNotFoundException, SampleStatusException{
		int i, l = args.length;
		if(l < 2){
			return display(true);
		}
		for(i = 0; i < l; i++){
			String arg = args[i].toLowerCase();
			if(arg.equals("--help") || arg.equals("-h")){
				return display(true);
			}else if(arg.equals("-snpfile")){
				i++;
				if(i == l){
					throw new ArgsException("No Input Snp File.");					
				}
				snpFile = args[i];
				if((i == args.length) || snpFile.charAt(0) == '-'){
					throw new ArgsException("No Input Snp File.");
				}
			}else if(arg.equals("-phasedfile")){
				i++;
				if(i == l){
					throw new ArgsException("No Input Sample Info File.");
				}
				this.phasedFile = args[i];
				if((i == args.length) || phasedFile.charAt(0) == '-'){
					throw new ArgsException("No Input Phased SNP File.");
				}
				parseStatus();
			}else if(arg.equals("-sampleinfo")){
				i++;
				if(i == l){
					throw new ArgsException("No Input Sample Info File.");
				}
				sampleFile = args[i];
				if((i == args.length) || sampleFile.charAt(0) == '-'){
					throw new ArgsException("No Input Sample Info File.");
				}
				parseStatus();
			}else if(arg.equals("-geneposition")){
				i++;
				if(i == l){
					throw new ArgsException("No Input Snp File.");					
				}
				geneFile = args[i];
				if((i == args.length) || snpFile.charAt(0) == '-'){
					throw new ArgsException("No Gene Position File.");
				}
			}else if(arg.equals("-filesplit")){
				i++;
				if(i == l){
					throw new ArgsException("No File Split Tag.");
				}
				fileSplit = args[i];
				if((i == args.length) || fileSplit.charAt(0) == '-'){
					throw new ArgsException("No File Split Tag.");
				}
			}else if(arg.equals("-output")){
				i++;
				if(i == l){
					throw new ArgsException("No Output File Name.");
				}
				output = args[i];
				if((i == args.length) || output.charAt(0) == '-'){
					throw new ArgsException("No Output File Name.");
				}
			}else if(arg.equals("-nnratio")){
				i++;
				if((i == args.length) || args[i].charAt(0) == '-'){
					throw new ArgsException("No Limit Ratio of NN.");
				}
				NNRatio = Double.parseDouble(args[i]);
				if(NNRatio < 0 || NNRatio > 1){
					throw new ArgsException("Unvalid Ratio of NN.");					
				}
			}else if(arg.equals("-maf")){
				i++;
				if((i == args.length) || args[i].charAt(0) == '-'){
					throw new ArgsException("No Limit Value of MAF.");
				}
				MAF = Double.parseDouble(args[i]);
				if(MAF < 0 || MAF > 1){
					throw new ArgsException("Unvalid value of MAF.");					
				}
			}else if(arg.equals("-hwe")){
				i++;
				if((i == args.length) || args[i].charAt(0) == '-'){
					throw new ArgsException("No Limit Value of Hardy–Weinberg equilibrium.");
				}
				HWE = Double.parseDouble(args[i]);
				if(HWE < 0 || HWE > 1){
					throw new ArgsException("Unvalid Ratio of Hardy–Weinberg equilibrium.");					
				}
//			}else if(arg.equals("-maxn")){
//				i++;
//				if((i == args.length) || args[i].charAt(0) == '-'){
//					throw new ArgsException("No Limit Ratio of N in one Snp.");
//				}
//				MAXN = Double.parseDouble(args[i]);
//				if(MAXN < 0 || MAXN > 1){
//					throw new ArgsException("Unvalid Limit Ratio of N in one Snp.");					
//				}
			}else if(arg.equals("-minht")){
				i++;
				if((i == args.length) || args[i].charAt(0) == '-'){
					throw new ArgsException("No Limit Ratio of Haplotype in Block.");
				}
				htminratio = Double.parseDouble(args[i]);
				if(htminratio < 0 || htminratio > 1){
					throw new ArgsException("Unvalid Limit Ratio of Haplotype in Block.");					
				}
			}else if(arg.equals("-ld-distance")){
				i++;
				if((i == args.length) || args[i].charAt(0) == '-'){
					throw new ArgsException("No Max Distance Limitation in the Process of Caculation of LD.");
				}
				DISTANCE = Integer.parseInt(args[i]);
				if(DISTANCE < 0){
					throw new ArgsException("Unvalid Max Distance Limitation.");					
				}
			}else if(arg.equals("-ht-window")){
				i++;
				if((i == args.length) || args[i].charAt(0) == '-'){
					throw new ArgsException("No window size while using sliding window approach to find haplotype block.");
				}
				this.HTWINDOW = this.parseWindow(args[i]);
				if(this.HTWINDOW < 0){
					throw new ArgsException("Unvalid Window Size.");					
				}
			}else if(arg.equals("-phase-window")){
				i++;
				if((i == args.length) || args[i].charAt(0) == '-'){
					throw new ArgsException("No window size while using sliding window approach to phase.");
				}
				this.phaseWINDOW = this.parseWindow(args[i]);
				if(this.phaseWINDOW < 0){
					throw new ArgsException("Unvalid Window Size.");					
				}
			}else if(arg.equals("-threads")){
				i++;
				if((i == args.length) || args[i].charAt(0) == '-'){
					throw new ArgsException("No Number of Threads Input.");
				}
				threads = Integer.parseInt(args[i]);
				if(threads < 0){
					throw new ArgsException("Unvalid Number of Threads.");					
				}
			}else if(arg.equals("--check")){
				CHECK = true;
			}else if(arg.equals("--ld")){
				LD = true;
			}else if(arg.equals("--filter")){
				FILTER = true;
			}else if(arg.equals("--badddata")){
				BADDATA = true;
			}else if(arg.equals("--full")){
				FULL = true;
			}else if(arg.equals("--gwhas")){
				BLOCK = true;
				PHASE = true;
				CC = true;
				TAG = true;
				CHECK = true;
			}else if(arg.equals("--block")){
				BLOCK = true;
//			}else if(arg.equals("--recom")){
//				RECOM = true;
			}else if(arg.equals("--phase")){
				PHASE = true;
			}else if(arg.equals("--cc")){
				CC = true;
			}else if(arg.equals("--tag")){
				TAG = true;
//			}else if(arg.equals("--map")){
//				MAP = true;
			}else if(arg.equals("--ignoregenotypeexception")){
				ignoreGenotypeException = true;
			}else if(arg.equals("--silence")){
				silence = true;
			}else{
				throw new ArgsException("UNKNOWN Parameter: "+args[i]);
			}
		}
		if(snpFile == null){
			throw new ArgsException("MUST have snp input file!");
		}
		if(sampleFile == null && CC){
			throw new ArgsException("If you want to take case/control test, you MUST have sample info file!");
		}
//		if(!CC && MAP){
//			throw new ArgsException("If you want to Map SNP and Haplotype block to gene(s), you MUST take case/control test!");						
//		}
//		if(MAP && this.geneFile == null){
//			throw new ArgsException("If you want to Map SNP and Haplotype block to gene(s), you MUST have Gene position file!");			
//		}
		if(output == null){
			output = snpFile;
		}
		if(this.FULL && this.BLOCK){
			this.PHASE = true;
		}
		if(!BLOCK){
			if(PHASE){
				throw new ArgsException("If you try to PHASE, you have to process BLOCK!");
			}
		}
		return display(false);
	}
	
	private int parseWindow(String argument){
		return Integer.parseInt(argument);
	}
	
	@SuppressWarnings("resource")
	private void parseStatus() throws FileNotFoundException, SampleStatusException{
		Scanner in = new Scanner(new File(sampleFile));
		ArrayList<Integer> ss = new ArrayList<Integer>();
		int i = 1;
		while(in.hasNext()){
			int n = in.nextInt();
			if(n < -1 || n > 1){
				throw new SampleStatusException(i, n);
			}
			ss.add(n);
			i++;
		}
		i--;
		SAMPLES = i;
		status = new int[i];
		Iterator<Integer> iter = ss.iterator();
		i = 0;
		while(iter.hasNext()){
			status[i++] = iter.next();
		}
		setStatus(status);
	}

	private Boolean display(boolean needHelp){
		if(needHelp){
			System.out.println(new Info());
			return false;
		}
		System.out.println("SNP data file: "+snpFile);
		System.out.println("Sample Info file: "+sampleFile);
		System.out.println("OUTPUT: "+output);
		System.out.println("cases: "+cases);
		System.out.println("controls: "+controls);
		System.out.println("unknowns: "+unknowns);
		if(CC){
			System.out.println("excepted distribution: "+freqs[0]+", "+freqs[1]);
		}
		if(BLOCK){
			System.out.print("find block, ");
		}
		if(PHASE || FULL){
			System.out.print("phase, ");
		}
		if(CC){
			System.out.print("case-control test, ");
		}
		if(TAG){
			System.out.print("tag, ");
		}
		if(LD){
			System.out.print("ld, ");
		}
		if(RECOM){
			System.out.print("recombination hotspots, ");
		}
		if(MAP){
			System.out.print("map gene, ");
		}
		System.out.println();
		if(!CHECK && !FILTER && !BADDATA && !LD && !BLOCK && !TAG && !FULL && !RECOM){
			return false;
		}
		if(this.silence){
			return true;
		}
		System.out.print("caculate? [Y/N]:");	
		@SuppressWarnings("resource")
		Scanner in = new Scanner(System.in);
		String answer = in.next();
		if(answer == null || answer.length() == 0 || (answer.charAt(0) != 'y' && answer.charAt(0) != 'Y')){
			return false;
		}else{
			return true;
		}
	}

	public void cancelBADDATA(){
		this.BADDATA = false;
	}
	
	public void cancelFILTER(){
		this.FILTER = false;
	}
	
	public void cancelFULL(){
		this.FULL = false;
	}
	
	public void cancelCHECK(){
		this.CHECK = false;
	}

	public boolean doFULL() {
		return FULL;
	}

	public void setFULL() {
		FULL = true;
	}
	
	public boolean doCC() {
		return CC;
	}

	public void setCC() {
		CC = true;
	}

	public boolean doTAG() {
		return TAG;
	}

	public void setTAG() {
		TAG = true;
	}

	public boolean doFILTER() {
		return FILTER;
	}

	public void setFILTER() {
		FILTER = true;
	}

	public boolean doCHECK() {
		return CHECK;
	}

	public void setCHECK() {
		CHECK = true;
	}

	public boolean doBLOCK() {
		return BLOCK;
	}

	public void setBLOCK() {
		BLOCK = true;
	}

	public boolean doPHASE() {
		return PHASE;
	}

	public void setPHASE() {
		PHASE = true;
	}

	public boolean doLD() {
		return LD;
	}

	public void setLD() {
		LD = true;
	}

	public void setFileSplit(String fileSplit) {
		this.fileSplit = fileSplit;
	}

	public double getD() {
		return D;	
	}

	public void setD(float d) {
		D = d;
	}

//	public int getWIN() {
//		return WIN;
//	}
//
//	public void setWIN(int wIN) {
//		WIN = wIN;
//	}
//
	public int getHEAD() {
		return HEAD;
	}

	public void setHEAD(int hEAD) {
		HEAD = hEAD;
	}

	public double getMAF() {
		return MAF;
	}

	public void setMAF(double mAF) {
		MAF = mAF;
	}

	public double getHWE() {
		return HWE;
	}

	public void setHWE(double hWE) {
		HWE = hWE;
	}

	public double getMAXN() {
		return MAXN;
	}

	public void setMAXN(double mAXN) {
		MAXN = mAXN;
	}

/*	public int[] getValids() {
		return valids;
	}

	public void setValids(int[] valids) {
		this.valids = valids;
	}
	
	 默认时全体均设为合格的SNP 
	public void setValids(){
		int i;
		valids = new int[SAMPLES];
		for(i = 0; i < SAMPLES; i++){
			valids[i] = i;
		}
	}
*/
	public void setNNRatio(double nNRatio) {
		NNRatio = nNRatio;
	}

	public int getSAMPLES() {
		return SAMPLES;
	}

	public void setSAMPLES(int sAMPLES) {
		SAMPLES = sAMPLES;
	}

	public String getFileSplit() {
		return fileSplit;
	}

	public double getNNRatio() {
		return NNRatio;
	}

	public String getOutput() {
		return output;
	}

	public void setOutput(String output) {
		this.output = output;
	}

	public boolean doBAD() {
		return BADDATA;
	}

	public void setBAD() {
		BADDATA = true;
	}

	public int getDISTANCE() {
		return DISTANCE;
	}

	public void setDISTANCE(int dISTANCE) {
		DISTANCE = dISTANCE;
	}

	public double getLDCL() {
		return LDCL;
	}

	public void setLDCL(double lDCL) {
		LDCL = lDCL;
	}

	public double getLDCU() {
		return LDCU;
	}

	public void setLDCU(double lDCU) {
		LDCU = lDCU;
	}

	public double getEHRCU() {
		return EHRCU;
	}

	public void setEHRCU(double eHRCU) {
		EHRCU = eHRCU;
	}

//	public int getMAXSIZE() {
//		return MAXSIZE;
//	}
//
//	public void setMAXSIZE(int mAXSIZE) {
//		MAXSIZE = mAXSIZE;
//	}
//
	public boolean doRECOM() {
		return RECOM;
	}

	public void setRECOM() {
		RECOM = true;
	}

	public double[] getFreqs(){
		return freqs;
	}
	
	public int[] getStatus() {
		return status;
	}

	public void setStatus(int[] status) {
		this.status = status;
		cases = controls = unknowns = 0;
		for(int i: status){
			if(i == 0){
				controls++;
			}else if(i == 1){
				cases++;
			}else{
				unknowns++;
			}
		}
		if(cases+controls != 0){
			freqs = new double[2];
			freqs[0] = (double)(cases)/(cases+controls);
			freqs[1] = (double)(controls)/(cases+controls);
		}else{
			freqs = null;			
		}
	}

	public String getSnpFile() {
		return snpFile;
	}

	public String getSampleFile() {
		return sampleFile;
	}

	public int getCases() {
		return cases;
	}

	public int getControls() {
		return controls;
	}

	public int getUnknowns() {
		return unknowns;
	}

	public int getB() {
		return B;
	}

	public void setB(int b) {
		B = b;
	}

	public int getK() {
		return K;
	}

	public void setK(int k) {
		K = k;
	}

	public int getThreads() {
		return threads;
	}

	public void setThreads(int threads) {
		this.threads = threads;
	}

	public boolean isIgnoreGenotypeException() {
		return ignoreGenotypeException;
	}

	public void setIgnoreGenotypeException(boolean ignoreGenotypeException) {
		this.ignoreGenotypeException = ignoreGenotypeException;
	}

	public double getHtminratio() {
		return htminratio;
	}

	public void setHtminratio(double htminratio) {
		this.htminratio = htminratio;
	}

	public boolean isSilence() {
		return silence;
	}

	public void setSilence(boolean silence) {
		this.silence = silence;
	}

	public int getPhaseWINDOW() {
		return phaseWINDOW;
	}

	public void setPhaseWINDOW(int phaseWindow) {
		this.phaseWINDOW = phaseWindow;
	}

	public int getHTWINDOW() {
		return HTWINDOW;
	}

	public void setHTWINDOW(int hTWINDOW) {
		HTWINDOW = hTWINDOW;
	}

	public int getRHWINDOW() {
		return RHWINDOW;
	}

	public void setRHWINDOW(int rHWINDOW) {
		RHWINDOW = rHWINDOW;
	}
	
	public int getSIZE(boolean isBlock){
		if(isBlock){
			return this.HTWINDOW;
		}else{
			return this.RHWINDOW;
		}
	}

	public String getPhasedFile() {
		return phasedFile;
	}

	public void setPhasedFile(String phasedFile) {
		this.phasedFile = phasedFile;
	}

	public boolean needRescanPhasedFile() {
		if(this.FULL && (this.LD || this.BLOCK)){
			return true;
		}else{
			return false;
		}
	}

	public String getGeneFile() {
		return geneFile;
	}

	public void setGeneFile(String geneFile) {
		this.geneFile = geneFile;
	}

	public boolean isMAP() {
		return MAP;
	}

	public void setMAP(boolean mAP) {
		MAP = mAP;
	}

}
