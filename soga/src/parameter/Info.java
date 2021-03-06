package parameter;

public class Info {
	private static String name;
	private static String discription;
	private static String[][] ioparameters;
	private static String[][] qcparameters;
	private static String[][] proparameters;
	private static String[][] inputrequire;
//	private static String version;
//	private static String[] authors;
//	private static String license;
	private static String[] alsosee;
//	private static String[] todos;
	
	public Info(){
		name = "**SOGA** - SOftware of Genome-wide Analysis for SNP Haplotype";
		
		discription = "A toolkit to caculate Linkage Disequilibrium and find haplotype block in SNP data.";
		
		ioparameters = new String[5][2];
		ioparameters[0][0] = "**-snpfile**";
		ioparameters[0][1] = "specifies the file of SNP data";
		ioparameters[1][0] = "**-phasedfile";
		ioparameters[1][1] = "specifies the PHASED file of SNP data";
		ioparameters[2][0] = "**-sampleinfo**";
		ioparameters[2][1] = "specifies the file of sample infomation";
		ioparameters[3][0] = "**-filesplit**";
		ioparameters[3][1] = "use delimiter instead of **defualt TAB** for process SNP file.";
		ioparameters[4][0] = "**-output**";
		ioparameters[4][1] = "specifies the files where to write the results.";
		
		qcparameters = new String[8][2];		
		qcparameters[0][0] = "**-sample-nn**";
		qcparameters[0][1] = "upper bound of ratio of nn genotype in one SAMPLE, **default 0.20**.";
		qcparameters[1][0] = "**-snp-nn**";
		qcparameters[1][1] = "upper bound of ratio of nn genotype in one SNP, **default 0.05**.";
		qcparameters[2][0] = "**-maf**";
		qcparameters[2][1] = "lower bound of MAF(Minor Allele Frequency), **default 0.01.**";
		qcparameters[3][0] = "**-hwe**";
		qcparameters[3][1] = "lower bound of HWE(Hardy–Weinberg equilibrium), **default 0.001.**";
		qcparameters[4][0] = "**--ignoreGenotypeException**";
		qcparameters[4][1] = "the SNP which have invalid genotype will be discarded silently.";
		qcparameters[5][0] = "**-minHt**";
		qcparameters[5][1] = "lower bound of ratio of haplotype when write haplotypes in __OUTPUT__.BLOCK file, **default 0.05**.";
		qcparameters[6][0] = "**--cancel-check**";
		qcparameters[6][1] = "do NOT check order in dataset.";
		qcparameters[7][0] = "**--cancel-sample-qc**";
		qcparameters[7][1] = "do NOT take quality control on samples";
   
		proparameters = new String[16][2];
		proparameters[0][0] = "**-threads**";
		proparameters[0][1] = "the number of threads that SOGA uses in PHASE, **defualt 2**.";
		proparameters[1][0] = "**--full**";
		proparameters[1][1] = "to phase all SNPs.";
		proparameters[2][0] = "**--block**";
		proparameters[2][1] = "to find Haplotype Blocks, output as __OUTPUT__.BLOCK";
		proparameters[3][0] = "**--phase**";
		proparameters[3][1] = "to phase samples in blocks, output as __OUTPUT__.PHASED, MUST also have `--block`.";
		proparameters[4][0] = "**--cc**";
		proparameters[4][1] = "to do case-control test(include OR and  chi-square) in 1)SNP, in __OUTPUT__.CHECK, and, 2)haplotypes, "
				+ " in __OUTPUT__.BLOCK; if you have `--block` and `--phase` in parameters, SOGA will do 1) and 2),"
				+ " otherwise, it just do 1).";
		proparameters[5][0] = "**--tag**";
		proparameters[5][1] = "to find tag SNP in blocks, output in __OUTPUT__.BLOCK, MUST also have `--block`";
		proparameters[6][0] = "**--GWHAS**";
		proparameters[6][1] = "do GWAS and GHWA";
		proparameters[7][0] = "**--snp-qc**";
		proparameters[7][1] = "QC results of SNPs, output as __OUTPUT__.SNPQC.";
		proparameters[8][0] = "**--filter**";
		proparameters[8][1] = "SNPs which PASS QC, output as __OUTPUT__.FILTER.";
		proparameters[9][0] = "**--baddata**";
		proparameters[9][1] = "SNPs which NOT PASS QC, output as __OUTPUT__.BADDATA.";
		proparameters[10][0] = "**--ld**";
		proparameters[10][1] = "to do LD caculation between two SNPs which distance is less than DISTANCE, output as __OUTPUT__.LD.";
		proparameters[11][0] = "**-ld-distance**";
		proparameters[11][1] = "MAX distance in calculate LDs in adjacent SNPs, **DEFAULT 1000**(bp).";
		proparameters[12][0] = "**-ht-window**";
		proparameters[12][1] = "the size of snps in one window while recongnizing haplotype blocks, **DEFAULT 200.**";
		proparameters[13][0] = "**-phase-window**";
		proparameters[13][1] = "the size of snps in one windo while phasing, **DEFAULT 2,000,000.**";
		proparameters[14][0] = "**--silence**";
		proparameters[14][1] = "running without asking confirmation";
		proparameters[15][0] = "**--help**";
		proparameters[15][1] = "print THIS help.";
		
		inputrequire = new String[2][2];
		inputrequire[0][0] = "**SNP DATA**";
		inputrequire[0][1] = "One line is a SNP, eg. `rs3013006 C/G chr21 16054667 CG GC GC GG GC CC GC CC CC........` is "
				+ "`#rs allele1/allele2 chromosome position genotype1 genotype2 genotyp3 ...`. "
				+ "You can have many different chromosomes, but MUST in one continuous segment,"
				+ " NOT have many segments in many places.";
		inputrequire[1][0] = "**Sample INFO**";
		inputrequire[1][1] = "Info file must have `0`(normal)/`1`(case)/`-1`(unknown)s, "
				+ "represents status in samples order by genotypes in SNP data. NO more contents."
				+ "Number of digits MUST as SAME as genotypes.";
		
//		version = "0.8.0";
//		
//		authors = new String[1];
//		authors[0] = "Heliy<heliy7#163.com>";
//		
//		license = "GPLv3";
		
		alsosee = new String[5];
		alsosee[0] = "The algorithm of LD caculation is same as [Haploview](http://www.ncbi.nlm.nih.gov/pubmed/15297300).";
		alsosee[1] = "The algorithm of haplotype block finding is [Taliun et al, 2014](http://www.ncbi.nlm.nih.gov/pubmed/24423111).";
		alsosee[2] = "The algorithm of phase is [SHAPEIT2](http://www.ncbi.nlm.nih.gov/pubmed/22138821/).";
		alsosee[3] = "The algorithm of tag SNP is [Yan et al, 2015](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4373116/).";
		alsosee[4] = "This is the work of Heliy's undergraduate graduation project, [article, chinese]().";
		
		
//		todos = new String[5];
//		todos[0] = "Find Recombination Hotspots.";
//		todos[1] = "Caculation iHS.(for Positive Selection)";
//		todos[2] = "Input Genetic Map to caculate recombination rates between SNPs.";
//		todos[3] = "Process in X/Y Chromosome and Family Dataset.";
//		todos[4] = "Input File from PLINK format.";
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append("NAME\n");
		sb.append("  "+name+"\n");
		sb.append("DISCRIPTION\n");
		sb.append("  "+discription+"\n");

		sb.append("PARAMETERS\n");
		sb.append(paras(proparameters));
		sb.append(paras(ioparameters));
		sb.append(paras(qcparameters));
		
		sb.append("\nMore information, see: https://www.bioapp.org/soga\n");
//		
//		sb.append("\n**INPUT REQUIREMENT**\n");
//		sb.append(paras(inputrequire));
//		
//		sb.append("\n**VERSION**\t"+version+"\n");
//		sb.append("\n**LICENSE**\t"+license+"\n");
//		sb.append("\n**AUTHORS**\n");
//		sb.append(list(authors));
//		
//		sb.append("\n**ALSO SEE**\n");
//		sb.append(list(alsosee));
//		sb.append("\n**TODO**\n");
//		sb.append(list(todos));
		return sb.toString();
	}

	private StringBuffer paras(String[][] parameters){
		StringBuffer sb = new StringBuffer();
		
		int i, n;
		n = parameters.length;
		for(i = 0; i < n; i++){
			sb.append("\n    ");
			sb.append(parameters[i][0]);
			sb.append(",  ");
			sb.append(parameters[i][1]);
			sb.append("\n");
		}
		return sb;
	}
	
	private StringBuffer list(String[] l){
		int i, n = l.length;
		StringBuffer sb = new StringBuffer();
		
		for(i = 0; i < n; i++){
			sb.append("\n  "+(i+1)+". ");
			sb.append(l[i]);
			sb.append("\n");
		}
		return sb;
	}
}
