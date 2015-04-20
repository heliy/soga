import sys

filename = sys.argv[1]
before = sys.argv[2]
splitflag = sys.argv[3]

isignore = True
f = open(filename)

for line in f.readlines():
    if not isignore:
        print line
        continue
    if before == line.split()[0]:
        isignore = False

f.close()
