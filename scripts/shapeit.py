#coding:UTF-8

import math
import random

genotype = [1, 1, 0, 1, 1, 0, 1, 1, 1]
genotypes = [genotype]
# [0, 1, 6][7, 6, 1]
L = len(genotype)
N = 20
hetes = [genotype.count(1)]

for i in xrange(1, N):
    genotypes.append([])
    c = 0
    for j in xrange(L):
        t = random.random()
        if t < 0.3:
            genotypes[i].append(0)
        elif t < 0.8:
            genotypes[i].append(1)
            c += 1
        else:
            genotypes[i].append(2)
    hetes.append(c)

from bak import genotypes

hetes = [g.count(1) for g in genotypes]

poss = [-3.12, -4.53, -2.70, -3.77, -3.170, -3, -4.0, -2, -1,75]
mafs = [0.45, 0.3, 0.25, 0.23, 0.4, 0.2, 0.35, 0.22, 0.16]
B = 3

T2 = 2**(B)
K = min(sum(2**(n) for n in hetes), 200)
# print K

H = []     # H[K][L], H[i][j] 为第j个SLP在第i个单体型上是0/1
getH = []  # getH[K][...] H[i][j] 表示H[j]是从第i个基因型中选出的
for i in xrange(0, N):
    # print i, hetes[i], float(K)/(N*(2**hetes[i]))
    selectp = max(float(K)/(N*(2**hetes[i])), 1.0/K)
    # print i, selectp
    genotype = genotypes[i]
    getH.append([])
    for j in xrange(2**(hetes[i])):
        haplotype = []
        if random.random() > selectp:
            continue
        for (index, loc) in enumerate(genotype):
            if loc == 1:
                if random.random() < mafs[index]:
                    haplotype.append(1)
                else:
                    haplotype.append(0)
            else:
                haplotype.append(loc/2)
        H.append(haplotype)
        getH[i].append(len(H)-1)
        # print haplotype
    # print getH[i]
K = len(H)
# print K

def update(i, hs):
    for j in getH[i]:
        if random.random() < 0.5:
            H[j] = hs[0]
        else:
            H[j] = hs[1]

thetas, unthetas = [0], [0] # e^(-ro/K), (1-e^(-ro/K))/K
##################### TODO: 遗传距离
for pos in poss[1:]:
    thetas.append(math.exp(pos/K))
    unthetas.append((1-thetas[-1])/K)

lamb = 0
for i in xrange(1, K):
    lamb += 1.0/i
lamb = 1/lamb
lamb = lamb/(2*(lamb+K))

def limit(p):
    return max(p, 0.001)

# 以上与基因型无关 放入构造中
#######################################

def shapeit(n):
    genotype = genotypes[n]
    # print genotype

    # 分段
    S = []        # 每个位点在第几段上
    inLo = 0
    segLo = 0
    for index in xrange(L):
        if genotype[index] == 1:
            inLo += 1
        S.append(segLo)
        if inLo == B:
            segLo += 1
            inLo = 0
    #     print segLo, inLo, S[-1]
    # print S
    # return

    # hete = hetes[0]
    A = []   # A[T2][L], A[i][j] 为第j个SLP在第i个单体型上是0/1
    for i in xrange(T2):
        m = i     # 保存为m号单体型
        dits = [] # 0/1 化
        for j in xrange(B-1, -1, -1):
            dits.append(m%2)
            m /= 2
        A.append([])
        # print dits
        k = 0
        for j in xrange(L): # 考入HA中
            if genotype[j] == 1:
                A[i].append(dits[k%B])
                k += 1
            else:
                A[i].append(genotype[j]/2)
        # print A[i]

    # P(Z1 = u)
    # P(Z1 = ht)
    def p_z1(ht):
        return 1.0/K

    # P(Zl = u| Zl-1 = v)
    # P(Zsnp = ht1| Zsnp-1 = ht2)
    def p_zl_zlm1(snp, ht1, ht2):
        theta, untheta = thetas[snp], unthetas[snp]
        if ht1 == ht2:
            return theta + untheta
        else:
            return untheta
        
    # P(Xl = i| Zl = u,H)
    # P(Xsnp = ht1| Zsnp = ht2, H)
    def p_xl_zl(snp, ht1, ht2):
        if A[ht1][snp] == H[ht2][snp]:
            return 1-lamb
        else:
            return lamb

    # p(Xl = i, Zl = u| Xl-1 = j, Zl-1 = v, H)
    # p(Xsnp = nexthtx, Zsnp = nexthtz| Xsnp-1 = prevhtx, Zsnp-1 = prevhtz)
    def p_xzl_xzlm1(snp, nexthtx, nexthtz, prevhtx, prevhtz):
        if S[snp] == S[snp-1] and nexthtx != prevhtx:
            return 0
        return limit(p_xl_zl(snp, nexthtx, nexthtz)*p_zl_zlm1(snp, nexthtz, prevhtz))

    sums_one_all, sums_all_one = [0]*T2, [0]*K 

    ###########
    ## 前向
    ###########
    alpha = [] # alpha[SNP个数][单体型总个数][单体型对个数]

    # 初始
    alpha.append([])
    m, total = 0, 0
    for i in  xrange(T2):
        alpha[m].append([])
        for u in xrange(K):
            # a[0](i, u) = P(Xl = i| Zl = u,H)*P(Z1 = u)
            alpha[m][i].append(limit(p_z1(u)*p_xl_zl(m, i, u)))
            total += alpha[m][i][u]
            sums_one_all[i] += alpha[m][i][u]
            sums_all_one[u] += alpha[m][i][u]
            # print alpha[m][i][u]
        
    # 递推
    for m in xrange(0, L-1):
        alpha.append([])
        for i in xrange(T2):
            alpha[m+1].append([])
            for u in xrange(K):
                nextalpha = p_xl_zl(m+1, i, u)
                if S[m] == S[m+1]:
                    am_all = sums_one_all[i]
                    nextalpha *= (am_all*unthetas[m+1] + thetas[m+1]*alpha[m][i][u])
                else:
                    am_all, am_all_all = sums_all_one[u], total
                    nextalpha *= (am_all_all*unthetas[m+1]+ thetas[m+1]*am_all)
                alpha[m+1][i].append(limit(nextalpha))
            # print alpha[m+1][i]
        sums_one_all, sums_all_one, total = [0]*T2, [0]*K, 0
        for i in xrange(T2):
            for u in xrange(K):
                sums_one_all[i] += alpha[m][i][u]
                sums_all_one[u] += alpha[m][i][u]
                total += alpha[m][i][u]
        
    ##################
    ### 后退
    ##################
    beta = []  # beta[hete个数][单体型总数][单体型对个数]

    # 初始
    beta.append([])
    m = L-1
    sums_one_all, sums_all_one, total = [0]*T2, [0]*K, 0
    for i in xrange(T2):
        beta[0].append([])
        for u in xrange(K):
            beta[0][i].append(1)
            p = p_xl_zl(m, i, u)
            sums_one_all[i] += p*beta[0][i][u]
            sums_all_one[u] += p*beta[0][i][u]
            total += p*beta[0][i][u]

    # 递推
    for m in xrange(L-2, -1, -1):
        beta.insert(0, [])
        for i in xrange(T2):
            beta[0].append([])
            for u in xrange(K):
                if S[m] == S[m+1]:
                    bm_all = sums_one_all[i]
                    # print 0, m, i, u, bm_all, p_xl_zl(m+1, i, u),
                    nextbeta = unthetas[m+1]*bm_all + thetas[m+1]*beta[1][i][u]*p_xl_zl(m+1, i, u)
                else:
                    bm_all, bm_all_all = sums_all_one[u], total
                    # print 1, m, i, u, bm_all, bm_all_all, unthetas[m+1], thetas[m+1],
                    nextbeta = unthetas[m+1]*bm_all_all + thetas[m+1]*bm_all
                # print nextbeta
                beta[0][i].append(limit(nextbeta))
        sums_one_all, sums_all_one, total = [0]*T2, [0]*K, 0
        for i in xrange(T2):
            for u in xrange(K):
                p = p_xl_zl(m, i, u)
                sums_one_all[i] += p*beta[0][i][u]
                sums_all_one[u] += p*beta[0][i][u]
                total += p*beta[0][i][u]
                
    ###################
    ## 概率
    ###################
    X1, X2 = [], []
    tmpturn, inturn, outturn = [], [], []  # 临时、内、外部转移矩阵

    for i in xrange(T2):
        inturn.append([0]*T2)
        outturn.append([0]*T2)
        tmpturn.append([0]*T2)
        outturn[i][i] = 1

    # 初始

    def getCp(i):
        return T2-i-1

    def formatlize(l, m):
        total = sum(l)
        for i in xrange(m):
            l[i] = l[i]/total
        
    m, pmax, maxi = 0, 0, 0
    for i in xrange(T2/2):
        j = getCp(i)
        tmp = sum(alpha[m][i][u] * beta[m][i][u] for u in xrange(K))
        tmp *= sum(alpha[m][j][u] * beta[m][j][u] for u in xrange(K))
        if(tmp > pmax):
            pmax = tmp
            maxi = i
        # print A[i], tmp
    # print maxi, pmax
    X1.append(maxi)
    X2.append(getCp(maxi))
    
    while True:
        m += 1
        # 计算前一个位点到这个位点的转移概率
        for i in xrange(T2):
            for j in xrange(T2):
                tmpsum = 0
                for u1 in xrange(K):
                    for u2 in xrange(K):
                        tmp = alpha[m-1][i][u1]
                        tmp *= p_zl_zlm1(m, u2, u1)
                        tmp *= p_xl_zl(m, j, u2)
                        tmp *= beta[m][j][u2]
                        tmpsum += tmp
                inturn[i][j] = tmpsum
            # formatlize(inturn[i], T2)
        # 内部转移与外部转移合并
        # 外部*内部
        for i in xrange(T2):
            for j in xrange(T2):
                tmpturn[i][j] = sum(outturn[i][k]*inturn[k][j] for k in xrange(T2))
            formatlize(tmpturn[i], T2)
        for i in xrange(T2):
            for j in xrange(T2):
                outturn[i][j] = tmpturn[i][j]
        if S[m] != S[m-1]:
            # 跨断 算单体型概率 最大加入
            maei, maej = X1[-1], X2[-1]
            maxp = 0
            for i in xrange(T2):
                j = getCp(i);
                p = outturn[maei][i]*outturn[maej][j]
                if p > maxp:
                    maxi = i
                    maxp = p
            X1.append(maxi)
            X2.append(getCp(maxi))
            if len(X1) > S[-1]:
                break
            for i in xrange(T2):
                outturn[i] = [0]*T2
                outturn[i][i] = 1

    haplotype1 = []
    haplotype2 = []
    # inNo = 0
    # print S
    for i in xrange(L):
        haplotype1.append(A[X1[S[i]]][i])
        haplotype2.append(A[X2[S[i]]][i])
    # print S
    # print X1, X2
    print n, haplotype1, haplotype2

    update(n, [haplotype1, haplotype2])
    return [haplotype1, haplotype2]

def tonum(l):
    p = 1
    s = 0
    for n in l:
        if n:
            s += (p*2)
        p *= 2
    return s

haplotypes = []
for i in xrange(N):
    haplotypes.extend([tonum(h) for h in shapeit(i)])

d = {}
for h in haplotypes:
    if h in d:
        d[h] += 1
    else:
        d[h] = 1
        
for key in d:
    print key, d[key]

