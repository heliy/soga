package parameter;

public class Info {
	private static String name;
	private static String discription;
	private static String[][] ioparameters;
	private static String[][] qcparameters;
	private static String[][] proparameters;
	private static String[][] inputrequire;
	private static String version;
	private static String[] authors;
	private static String license;
	private static String[] alsosee;
	private static String[] todos;
	
	public Info(){
		name = "**SOGA** - SOftware of Genome-wide Analysis for caculating LD and identifying haplptype block";
		
		discription = "A toolkit to caculate Linkage Disequilibrium and find haplotype block in SNP data.";
		
		ioparameters = new String[4][2];
		ioparameters[0][0] = "**-snpfile**";
		ioparameters[0][1] = "specifies the file of SNP data, see input reqirements.";
		ioparameters[1][0] = "**-sampleinfo**";
		ioparameters[1][1] = "specifies the file of sample infomation, see input requirements.";
		ioparameters[2][0] = "**-filesplit**";
		ioparameters[2][1] = "use delimiter instead of **defualt TAB** for process SNP file.";
		ioparameters[3][0] = "**-output**";
		ioparameters[3][1] = "specifies the files where to write the results.";
		
		qcparameters = new String[6][2];		
		qcparameters[0][0] = "**-maxn**";
		qcparameters[0][1] = "uppper bound of ratio of n in one SNP, **default 0.8**.";
		qcparameters[1][0] = "**-nnratio**";
		qcparameters[1][1] = "upper bound of ratio of nn genotype in one SNP, **default 0.05**.";
		qcparameters[2][0] = "**-maf**";
		qcparameters[2][1] = "lower bound of MAF(Minor Allele Frequency), **default 0.01.**";
		qcparameters[3][0] = "**-hwe**";
		qcparameters[3][1] = "lower bound of HWE(Hardy–Weinberg equilibrium), **default 0.001.**";
		qcparameters[4][0] = "**--ignoreGenotypeException**";
		qcparameters[4][1] = "the SNP which have invalid genotype will be discarded silently.";
		qcparameters[5][0] = "**-minHt**";
		qcparameters[5][1] = "lower bound of ratio of haplotype when write haplotypes in __OUTPUT__.BLOCK file, **default 0.05**.";
   
		proparameters = new String[12][2];
		proparameters[0][0] = "**-threads**";
		proparameters[0][1] = "the number of threads that SOGA uses, **defualt 2**.";
		proparameters[1][0] = "**--block**";
		proparameters[1][1] = "to find Haplotype Blocks, output as __OUTPUT__.BLOCK";
		proparameters[2][0] = "**--phase**";
		proparameters[2][1] = "to phase samples in blocks, output as __OUTPUT__.PHASED, MUST also have `--block`.";
		proparameters[3][0] = "**--cc**";
		proparameters[3][1] = "to do case-control test in 1)SNP, have OR value in __OUTPUT__.CHECK, or, 2)haplotypes, output in __OUTPUT__.BLOCK;"
				+ "if you have `--block` and `--phase` in parameters, it will do 1) and 2), otherwise, soga just do 1).";
		proparameters[4][0] = "**--tag**";
		proparameters[4][1] = "to find tag SNP in blocks, output in __OUTPUT__.BLOCK, MUST also have `--block`";
		proparameters[5][0] = "**--check**";
		proparameters[5][1] = "QC results of SNPs, output as __OUTPUT__.CHECK.";
		proparameters[6][0] = "**--filter**";
		proparameters[6][1] = "SNPs which PASS QC, output as __OUTPUT__.FILTER.";
		proparameters[7][0] = "**--baddata**";
		proparameters[7][1] = "SNPs which NOT PASS QC, output as __OUTPUT__.BADDATA.";
		proparameters[8][0] = "**--ld**";
		proparameters[8][1] = "to do LD caculation between two SNPs which distance is less than DISTANCE, output as __OUTPUT__.LD.";
		proparameters[9][0] = "**-ld-distance**";
		proparameters[9][1] = "MAX distance in calculate LDs in adjacent SNPs, **DEFAULT 1000**(bp).";
		proparameters[10][0] = "**-window**";
		proparameters[10][1] = "the size of snps in one window, **DEFAULT 1000.**";
		proparameters[11][0] = "**--help**";
		proparameters[11][1] = "print THIS help.";
		
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
		
		version = "0.8.0";
		
		authors = new String[1];
		authors[0] = "Heliy<heliy7#163.com>";
		
		license = "GPLv3";
		
		alsosee = new String[5];
		alsosee[0] = "The algorithm of LD caculation is same as [Haploview](http://www.ncbi.nlm.nih.gov/pubmed/15297300).";
		alsosee[1] = "The algorithm of haplotype block finding is [Taliun et al, 2014](http://www.ncbi.nlm.nih.gov/pubmed/24423111).";
		alsosee[2] = "The algorithm of phase is [SHAPEIT2](http://www.ncbi.nlm.nih.gov/pubmed/22138821/).";
		alsosee[3] = "The algorithm of tag SNP is [Yan et al, 2015](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4373116/).";
		alsosee[4] = "This is the work of Heliy's undergraduate graduation project, [article, chinese]().";
		
		
		todos = new String[5];
		todos[0] = "Find Recombination Hotspots.";
		todos[1] = "Caculation iHS.(for Positive Selection)";
		todos[2] = "Input Genetic Map to caculate recombination rates between SNPs.";
		todos[3] = "Process in X/Y Chromosome and Family Dataset.";
		todos[4] = "Input File from PLINK format.";
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append("**NAME**\n\t");
		sb.append(name+"\n");
		sb.append("**DISCRIPTION**\n\t");
		sb.append(discription+"\n");

		sb.append("\n**PARAMETERS\n\t**CACULATION**\n");
		sb.append(paras(proparameters));
		sb.append("\t**Input/Output**\n");
		sb.append(paras(ioparameters));
		sb.append("\t**Quality Control**\n");
		sb.append(paras(qcparameters));
		
		sb.append("\n**INPUT REQUIREMENT\n");
		sb.append(paras(inputrequire));
		
		sb.append("\n**VERSION**\t"+version+"\n");
		sb.append("\n**LICENSE**\t"+license+"\n");
		sb.append("\n**AUTHORS**\n");
		sb.append(list(authors));
		
		sb.append("\n**ALSO SEE**\n");
		sb.append(list(alsosee));
		sb.append("\n**TODO**\n");
		sb.append(list(todos));
		return sb.toString();
	}

	private StringBuffer paras(String[][] parameters){
		StringBuffer sb = new StringBuffer();
		
		int i, n;
		n = parameters.length;
		for(i = 0; i < n; i++){
			sb.append("\n\t\t");
			sb.append(parameters[i][0]);
			sb.append(",\t");
			sb.append(parameters[i][1]);
			sb.append("\n");
		}
		return sb;
	}
	
	private StringBuffer list(String[] l){
		int i, n = l.length;
		StringBuffer sb = new StringBuffer();
		
		for(i = 0; i < n; i++){
			sb.append("\n\t"+(i+1)+". ");
			sb.append(l[i]);
			sb.append("\n");
		}
		return sb;
	}
}
