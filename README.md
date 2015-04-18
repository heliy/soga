**NAME**
	**SOGA** - SOftware of Genome-wide Analysis for caculating LD and identifying haplptype block
**DISCRIPTION**
	A toolkit to caculate Linkage Disequilibrium and find haplotype block in SNP data.

**PARAMETERS
	**CACULATION**
		**-threads**,	the number of threads that SOGA uses, **defualt 2**.
		**--block**,	to find Haplotype Blocks, output as __OUTPUT__.BLOCK
		**--phase**,	to phase samples in blocks, output as __OUTPUT__.PHASED, MUST also have `--block`.
		**--cc**,	to do case-control test in haplotypes, output in __OUTPUT__.BLOCK, MUST also have `--block` AND `--phase`.
		**--tag**,	to find tag SNP in blocks, output in __OUTPUT__.BLOCK, MUST also have `--block`
		**--check**,	QC results of SNPs, output as __OUTPUT__.CHECK.
		**--filter**,	SNPs which PASS QC, output as __OUTPUT__.FILTER.
		**--baddata**,	SNPs which NOT PASS QC, output as __OUTPUT__.BADDATA.
		**--ld**,	to do LD caculation between two SNPs which distance is less than DISTANCE, output as __OUTPUT__.LD.
		**-ld-distance**,	MAX distance in calculate LDs in adjacent SNPs, **DEFAULT 1000**(bp).
		**-window**,	the size of snps in one window, **DEFAULT 1000.**
		**--help**,	print THIS help.
	**Input/Output**
		**-snpfile**,	specifies the file of SNP data, see input reqirements.
		**-sampleinfo**,	specifies the file of sample infomation, see input requirements.
		**-filesplit**,	use delimiter instead of **defualt TAB** for process SNP file.
		**-output**,	specifies the files where to write the results.
	**Quality Control**
		**--maxn**,	uppper bound of ratio of n in one SNP, **default 0.8**.
		**--nnratio**,	upper bound of ratio of nn genotype in one SNP, **default 0.05**.
		**--maf**,	lower bound of MAF(Minor Allele Frequency), **default 0.01.**
		**--hwe**,	lower bound of HWE(Hardy–Weinberg equilibrium), **default 0.001.**

**INPUT REQUIREMENT
		**SNP DATA**,	One line is a SNP, eg. `rs3013006 C/G chr21 16054667 CG GC GC GG GC CC GC CC CC........` is `#rs allele1/allele2 chromosome position genotype1 genotype2 genotyp3 ...`. You can have many different chromosomes, but MUST in one continuous segment, NOT have many segments in many places.
		**Sample INFO**,	Info file must have `0`(normal)/`1`(case)/`-1`(unknown)s, represents status in samples order by genotypes in SNP data. NO more contents.Number of digits MUST as SAME as genotypes.

**VERSION**	0.8.0

**LICENSE**	GPLv3

**AUTHORS**
	1. Heliy<heliy7#163.com>

**ALSO SEE**
	1. The algorithem of LD caculation is same as [Haploview](http://www.ncbi.nlm.nih.gov/pubmed/15297300).
	2. The algorithem of haplotype block finding is [Taliun et al, 2014](http://www.ncbi.nlm.nih.gov/pubmed/24423111).
	3. The algorithem of phase is [SHAPEIT2](http://www.ncbi.nlm.nih.gov/pubmed/22138821/).
	4. The algorithem of tag SNP is [Yan et al, 2015](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4373116/).
	5. This is the work of Heliy's undergraduate graduation project, [article, chinese]().

**TODO**
	1. Find Recombination Hotspots.
	2. Caculation iHS.(for Positive Selection)
	3. Input Genetic Map to caculate recombination rates between SNPs.
	4. Process in X/Y Chromosome and Family Dataset.
	5. Input File from PLINK format.

