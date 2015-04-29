import sys
from math import sqrt, exp, log

def OR(cos):
    cos = [float(i) for i in cos]
    OR = (cos[0]*cos[3]+0.01)/(cos[2]*cos[3]+0.01)
    t = sqrt(sum(i and 1/i or 0 for i in cos))
    low = exp(log(OR)-1.96*t)
    high = exp(log(OR)+1.96*t)
    return [OR, low, high]

blockfile = sys.argv[1]
cases = int(sys.argv[2])*2
controls = int(sys.argv[3])*2

need = -1

for line in open(blockfile).readlines():
    if need == 1:
        if line[0] == '-':
            print line
            need = -1
        else:
            [a, c] = [int(i) for i in line.split("\t")[2].split("/")]
            b = cases - a
            d = controls - c
            print line[:-1]+"\t"+"\t".join(str(i) for i in OR([a, b, c, d]))
    else:
        print line
        if line[0] == '>':
            need = 0
        else:
            if not need:
                need = 1
                
