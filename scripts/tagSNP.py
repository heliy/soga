import random

R = [
      [1.0000000000000004, 0.9266923107650029,],
      [0.9266923107650029, 0.9999999999999994,],
    ]

num = 2
r = 0.99

N, Q, B, S = [], [], [], []

for i in xrange(num):
    N.append([])
    for j in xrange(num):
        if R[i][j] > r:
            N[i].append(1)
        else:
            N[i].append(0)

B = range(num)
while len(B):
    Q = random.sample(B, 1)
    B.remove(Q[0])
    for m in xrange(num):
        for n in xrange(i, num):
            if m in Q and n in B and N[m][n]:
                Q.append(n)
                B.remove(n)
    mRis = []
    for i in Q:
        mRis.append( sum(R[i]))
    minRi = min(mRis)
    S.append(Q[mRis.index(minRi)])
    print B

print S
        
