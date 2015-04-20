# coding:UTF-8
#!/bin/python

import sys

filename = sys.argv[1]

base = {
    "A" : "1",
    "C" : "2",
    "G" : "3",
    "T" : "4",
    }

f = open(filename)
ped = open(filename + ".ped", "wb+")
info = open(filename + ".info", "wb+")

line = f.readline()
parts = line.split()
sample = len(parts) - 4
peds = [[]] * sample
    
for i in range(sample):
    peds[i] = peds[i] + [parts[i + 4][0]]
    peds[i] = peds[i] + [parts[i + 4][1]]

info.write(parts[0] + "\t" + parts[3] + "\n")
        
for line in f.readlines():
    if line[0] == "#":
        continue    
    parts = line.split()
    info.write(parts[0] + "\t" + parts[3] + "\n")
    for i in range(sample):
        peds[i] = peds[i] + [parts[i + 4][0]]
        peds[i] = peds[i] + [parts[i + 4][1]]
print peds[1]
print len(peds[1])
for i in range(sample):
    ped.write( str(i+1) + "\t" + str(i+1) + "\t0\t0\t0\t0\t" + "\t".join(peds[i]) + "\n")

f.close()
ped.close()
info.close()
