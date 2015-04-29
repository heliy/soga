import sys

def tag(tagline):
    cons = tagline.split("\t")[4:]
    tags = []
    for term in cons:
        if term[0] == 'Y':
            tags.append(1)
        else:
            tags.append(0)
    return tags

class Haplotype:
    def __init__(self, line):
        cons = line.split("\t")
        self.pvalue = cons[1][0] != '-' and float(cons[1]) or None
        self.nums = cons[2][0]!= '-' and cons[2] or None
        self.ratio = float(cons[3])
        self.alleles = cons[1][0] != '-' and cons[4:-3] or cons[4:]
        self.OR = cons[1][0] != '-' and float(cons[-3]) or None
        self.low = cons[1][0] != '-' and float(cons[-2]) or None
        self.high = cons[1][0] != '-' and float(cons[-1]) or None

class Block:
    def __init__(self, head, haplolines=None, tagline=None, text = None):
        self.text = text
        cons = head.split("\t")
        self.chr = cons[0]
        self.begin = int(cons[1])
        self.end = int(cons[2])
        self.nums = int(cons[3])
        self.snps = cons[4:]
        self.haplotypes = haplolines and [Haplotype(line) for line in haplolines] or []
        if tagline:
            self.tags = tag(tagline)
            self.tag_num = sum(self.tags)
        else:
            self.tags = None
            self.tag_num = 0
        if self.haplotypes and self.haplotypes[0].pvalue:
            self.pvalue = min(ht.pvalue for ht in self.haplotypes)
        else:
            self.pvalue = None

    def length(self):
        return self.end - self.begin

def getblocks(filename):
    cons = open(filename).read().split('>>')
    if len(cons) > 1:
        blocks = []
        for text in cons[1:]:
            con = [ line for line in text.splitlines() if line]
            head = con[1]
            tagline = con[-1][0] == '-' and con[-1] or None
            if tagline:
                con.remove(tagline)
            haplolines = con[2][0] != '-' and con[2:] or None
            blocks.append(Block(head, haplolines, tagline, text))
        return blocks
    else:
        return None

def getsummary(filename):
    snps = {}
    lines = open(filename).read().splitlines()
    need = False
    for line in lines:
        if not line:
            continue
        if line.startswith("SNPS:"):
            need = True
        elif line.startswith("BLOCKS:"):
            break
        elif need:
            chro = line.split(":")[0]
            num = int(line.split("\t")[-1])
            snps[chro] = num
    return snps

class Chromosome:
    def __init__(self, chro, snps):
        self.chro = chro
        self.snps = snps
        self.blocks = []

    def addblock(self, block):
        self.blocks.append(block)

    def summary(self):
        self.blocks_num = len(self.blocks)
        lens = [block.length() for block in self.blocks]
        self.max_len = max(lens)
        self.ave_len = sum(lens)/(self.blocks_num)
        snp_nums = [block.nums for block in self.blocks]
        self.total_snps = sum(snp_nums)
        self.max_num = max(snp_nums)
        self.ave_num = self.total_snps/self.blocks_num
        self.snp_ratio = self.total_snps/float(self.snps)
        self.tags = sum(block.tag_num for block in self.blocks)
        self.tag_ratio = self.tags/float(self.total_snps)
        return [self.blocks_num, self.max_len, self.ave_len, self.max_num, self.ave_num, self.total_snps, self.tags, self.snp_ratio, self.tag_ratio]
        


def get_chromosomes(blocks, snps):
    chrs = {chro:Chromosome(chro, snps[chro]) for chro in snps}
    for block in blocks:
        chrs[block.chr].addblock(block)
    return chrs

def block_key(block):
    return block.pvalue

def sign_blocks(blocks, pvalue = 0.05, OR_value = 1.25, low_value = 1):
    if len(blocks) == 0 or not blocks[0].pvalue:
        return []
    signs = []
    for block in blocks:
        for haplotype in block.haplotypes:
            if haplotype.pvalue < pvalue and haplotype.OR > OR_value and haplotype.low >= low_value:
                signs.append(block)
    signs.sort(key = block_key)
    return signs

if __name__ == "__main__":
    blockfile = sys.argv[1]
    summaryfile = sys.argv[2]
    pvalue = len(sys.argv) == 3 and 0.0001 or float(sys.argv[4])
    ORValue = len(sys.argv) == 3 and 1.25 or float(sys.argv[5])
    lowValue = len(sys.argv) == 3 and 1 or float(sys.argv[5])
    blocks = getblocks(blockfile)
    snps = getsummary(summaryfile)
    chrs = get_chromosomes(blocks, snps)
    print "Chromosomes\t#blocks\tmax_len\tave_len\tmax_snps\tave_snps\t#snps_in_block\t#tags\t%snps_in_chromosomes\t%tags_in_blocks"
    cc = False
    for chro in sorted(snps.keys()):
        if len(chrs[chro].blocks):
            if chrs[chro].blocks[0].pvalue:
                cc = True
            print chro, " ".join([ str(i) for i in chrs[chro].summary()])

    if cc:
        signs = sign_blocks(blocks, pvalue)
        print
        print "Blocks which have haplotypes is significant:", len(signs)
        for block in signs:
            print ">>"+block.text,
