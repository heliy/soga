package dna;

public class Block {
	int up;
	int down;
	long begin;
	long end;
	Snp[] snps;

	public Block(int cup, int cdown, Snp[] snplist){
		up = cup;
		down = cdown;
		begin = snplist[up].position;
		end = snplist[down].position;
		snps = new Snp[down - up + 1];
		int i;
		for(i = up; i <= down; i++){
			snps[i - up] = snplist[i];
		}
	}

}
