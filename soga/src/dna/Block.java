package dna;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.CountDownLatch;

import parameter.Base;
import parameter.Setting;
import calculate.FitTest;
import calculate.OR;
import calculate.Phase;

public class Block implements Comparable<Block> {
	private static int threads;
	private static int[][] genosPerThreads;
	
	private int up;
	private int down;
	private long begin;
	private long end;
	private int num;
	private Snp[] snps;
	private int samples;

	private static double[] freqs;
	private static int caseHts;
	private static int controlHts;
	private HaploType[] haplotypes = null;
	private HaploType[] htall;
	private boolean[] isTag = null;
	private int minht;

	public Block(int Up, int Down, Snp[] snplist) {
		up = Up;
		down = Down;
		begin = snplist[up].getPosition();
		end = snplist[down].getPosition();
		num = down - up + 1;
		snps = new Snp[num];
		int i;
		for (i = up; i <= down; i++) {
			snps[i - up] = snplist[i];
		}
	}
	
	public void phase(Setting rc) throws InterruptedException {
		samples = rc.getSAMPLES();
		freqs = rc.getFreqs();
		caseHts = rc.getCases()*2;
		controlHts = rc.getControls()*2;
		
		minht = (int)(rc.getHtminratio()*this.samples*2);
//		System.out.println(minht);
		multithreads(rc);
		int i, n = htall.length;
		HaploType[] hts = new HaploType[n];
		for(i = 0; i < n; i++){
			hts[i] = htall[i];
    	}
		
		int[] status = rc.getStatus();
		Arrays.sort(hts);
		ArrayList<HaploType> htlist = new ArrayList<HaploType>();
		HaploType mae;
		i = 0; 
		mae = new HaploType(hts[0].getAlleles(), 0);
		mae.add(status[mae.getSampleNo()]);
		for(i = 1; i < n; i++){
//			System.out.println(hts[i]);
			if(hts[i].equals(mae)){
				mae.add(status[hts[i].getSampleNo()]);
			}else{
				htlist.add(mae.clone());
				mae = new HaploType(hts[i].getAlleles(), 0);
				mae.add(status[mae.getSampleNo()]);
			}
		}
	    htlist.add(mae);
		haplotypes = new HaploType[htlist.size()];
		htlist.toArray(haplotypes);
		Arrays.sort(haplotypes);
		n = 0;
		for(HaploType ht: haplotypes){
			n += ht.getNum();
		}
//		System.out.println(i+", "+n);
	}
	
	public void setTags(boolean[] bs){
		isTag = bs;
	}

	public void ccTest() {
		if(haplotypes == null){
			return;
		}
		FitTest ft = new FitTest(freqs);
		OR orc = new OR();
		double[] ors;
		double[] counts = new double[2];
		for (HaploType ht : haplotypes) {
			counts[0] = ht.getCaseNum();
			counts[1] = ht.getControlNum();
			if(counts[0]+counts[1] > 0){
				ors = orc.getOR_CI(counts[0], caseHts-counts[0],
						counts[1], controlHts-counts[1]);
    			ht.setPvalue(ft.test(counts), ors[0], ors[1], ors[2]);
			}
		}
	}

	

	public String toString() {
		int i;
		StringBuffer s = new StringBuffer();
		s.append(snps[0].getChr() + "\t"+ begin + "\t" + end + "\t" + num);
		for (i = 0; i < num; i++) {
			s.append("\t" + snps[i].getRs());
		}
		s.append("\n");
		
		if(haplotypes != null){
			i = 1;
			for(HaploType ht: haplotypes){
				if(ht.getNum() < minht){
					break;
				}
				s.append(i+"\t"+ht.toString(samples*2));
				i++;
			}
		}
		if(isTag != null){
			s.append("----\t----\t----\t");
			int tagnum = 0;
			for(boolean istag: isTag){
				if(istag){
					tagnum++;
				}
			}
			s.append(tagnum);
			for(i = 0; i < num; i++){
				if(isTag[i]){
					s.append("\tYES");
				}else{
					s.append("\tNO");
				}
			}
			s.append("\n\n");
		}
		return s.toString();
	}
	

	// 按起始位置升续
	public int compareTo(Block arg) {
		if (begin < arg.getBegin()) {
			return -1;
		}
		if (begin > arg.getBegin()) {
			return 1;
		}
		return 0;
	}

	public Snp[] getSnps() {
		return snps;
	}

	public int getUp() {
		return up;
	}

	public int getDown() {
		return down;
	}

	public long getBegin() {
		return begin;
	}

	public long getEnd() {
		return end;
	}

	public int getNum() {
		return num;
	}
	
	public boolean[] getIsTag(){
		return isTag;
	}

	public HaploType[] getHaplotypes() {
		return haplotypes;
	}

	
//	public void outspread(boolean[] status) {
//	HaploType[] haplos;
//	ArrayList<HaploType> haplolist = new ArrayList<HaploType>();
//	Snp snp;
//	long[] positions = new long[num];
//
//	int i, j;
//	int[][] G = new int[samples][num];
//
//	for (i = 0; i < num; i++) {
//		snp = snps[i];
//		int[][] types = snp.getTypes();
//		positions[i] = snp.getPosition();
//		for (j = 0; j < samples; j++) {
//			G[j][i] = types[j][0]+types[j][1];
//		}
//	}
//
//	for (j = 0; j < samples; j++) {
//		// G[j] = {0, 1, 2, 1, 2, ...} 0为AA，1为AB，2为BB
//		GenoType genotype = new GenoType(G[j], positions);
//		haplos = genotype.phase();
//		insert(haplolist, haplos[0], status[j]);
//		insert(haplolist, haplos[1], status[j]);
//	}
//
//	haplotypes = new HaploType[haplolist.size()];
//	haplolist.toArray(haplotypes);
//	Arrays.sort(haplotypes);
//
//	isTag = new boolean[num]; // 默认初值为false
//	int[] partion = new int[haploNum];
//
//	ArrayList<Integer> NRS = new ArrayList<Integer>(); // 单个非冗余的SNP的编号
//	ArrayList<int[]> singles = new ArrayList<int[]>(); // 对应的分类
//
//	// SNP分类 得到非冗余的SNPs（NRS）
//	for (i = 0; i < num; i++) {
//		snp = snps[i];
//		for (j = 0; j < haploNum; j++) {
//			partion[j] = snp.getDis(haplotypes[j].getAlleles()[i]) + 1;
//		}
//		j = 0;
//		Iterator<int[]> iter = singles.iterator();
//		while (iter.hasNext()) {
//			int[] al = iter.next();
//			if (Arrays.equals(al, partion)) {
//				j = 1;
//				break;
//			} else {
//				j++;
//			}
//		}
//		if (j == 0) {
//			singles.add(partion.clone());
//			NRS.add(i);
//		}
//	}
//
//	if (NRS.size() < 3) {
//		// 两个以内 不用选
//		isTag[NRS.get(0)] = true;
//		if (NRS.size() == 2) {
//			isTag[NRS.get(1)] = true;
//		}
//		return;
//	}
//
//	ArrayList<ArrayList<Integer>> schemes = new ArrayList<ArrayList<Integer>>(); // 选定的模式们
//	ArrayList<int[]> partions = new ArrayList<int[]>(); // 模式的分类情况
//	ArrayList<ArrayList<Integer>> backlogs = new ArrayList<ArrayList<Integer>>(); // 模式的可能加入的其他SNP
//	
//	int NRSs = NRS.size();
//	
//	ArrayList<Integer> template = new ArrayList<Integer>();
//	for(i = 0; i < NRSs; i++){
//		template.add(i);
//	}
//	
//	// 按单个SNP的方案补完
//	for(i = 0; i < NRSs; i++){
//		ArrayList<Integer> sc = new ArrayList<Integer>(i);
//		schemes.add(sc);
//		partions.add(singles.get(i));
//		template.remove(i);
//		backlogs.add(template);
//		template.add(i, i);
//	}
//
//	while (true) {
//		// 新一个
//		ArrayList<Integer> scheme = schemes.get(0);
//		Iterator<int[]> partionsIter = partions.iterator();
//		partion = partionsIter.next();
//		Iterator<ArrayList<Integer>> backlogsIter = backlogs.iterator();
//		ArrayList<Integer> backlog = backlogsIter.next();
//			
//		// 每个可能新加入的SNP
//		Iterator<Integer> backlogIter = backlog.iterator();
//		while(backlogIter.hasNext()){
//			i = backlogIter.next();
//			int[] sin = singles.get(i);
//			int[] L = new int[haploNum];
//			int[] newPartion = new int[haploNum];
//			ArrayList<Integer> flags = new ArrayList<Integer>();
//				
//			// 合并
//			for(j = 0; j < haploNum; j++){
//				L[j] = sin[j] + partion[j] * 2;
//				if(!flags.contains(L[j])){
//					flags.add(L[j]);
//				}
//			}
//			for(j = 0; j < haploNum; j++){
//				newPartion[j] = flags.indexOf(L[j]) + 1;
//			}
//				
//			// 满足完全分型要求
//			if(flags.size() == haploNum){
//				isTag[i] = true;
//				backlogIter = scheme.iterator(); // 借用一下
//				while(backlogIter.hasNext()){
//					isTag[backlogIter.next()] = true;
//				}
//				return;
//			}
//				
//			j = 0;
//			while(partionsIter.hasNext()){
//				j++;
//				int[] anotherPartion = partionsIter.next();
//				// 同型 将两者中的可能加入的SNP合并（交集）
//				if(Arrays.equals(newPartion, anotherPartion)){						
//					ArrayList<Integer> bl = backlogsIter.next();
//					ArrayList<Integer> n = new ArrayList<Integer>();						
//					Iterator<Integer> mae = bl.iterator();
//					while(mae.hasNext()){
//						int k = mae.next();
//						if(backlog.contains(k)){
//							n.add(k);
//						}
//					}
//					n.remove(i);
//					backlogs.set(j, n);
//				}else{
//					backlogsIter.next();
//				}
//			}
//			
//			// 新类
//			if(j == partions.size()){
//				ArrayList<Integer> ns = new ArrayList<Integer>(i);
//				ns.addAll(scheme);
//				schemes.add(ns);
//				partions.add(newPartion);
//				ArrayList<Integer> nb = new ArrayList<Integer>();
//				nb.addAll(backlog);
//				nb.remove(i);
//				backlogs.add(nb);
//			}
//		}
//		
//		// 该分类完结
//		schemes.remove(0);
//		partions.remove(0);
//		backlogs.remove(0);
//	}
//
//}

//private void insert(ArrayList<HaploType> haplolist, HaploType ht,
//		boolean isCase) {
//	Iterator<HaploType> iter = haplolist.iterator();
//	while (iter.hasNext()) {
//		HaploType in = iter.next();
//		if (in.equals(ht)) {
//			in.add(isCase);
//			return;
//		}
//	}
//	ht.add(isCase);
//	haplolist.add(ht);
//	haploNum++;
//}

}
