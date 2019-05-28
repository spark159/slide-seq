import sys
import random
import math
import copy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def ID_cmp(a, b):
    win1, loc1 = a.split('-')
    win2, loc2 = b.split('-')
    loc1, loc2 = int(loc1), int(loc2)
    if win1 < win2:
        return -1
    elif win1 == win2:
        if loc1 < loc2:
            return -1
        else:
            return 1
    else:
        return 1

def acf (samples):
    mu, std = np.mean(samples), np.std(samples)
    result = []
    for step in range(int(len(samples)/10)):
        corr = 0.0
        count = 0
        for i in range(len(samples) - step):
            corr += (samples[i]-mu)*(samples[i+step]-mu)
            count += 1
        corr = float(corr)/count
        corr = float(corr)/std
        result.append(corr)
    return result

def read_file(filenames, ref_length=225, dyad_offset=52, choice=['valid']):    
    cutmaps, dyadmaps = {}, {}
    for i in range(len(filenames)):
        filename = filenames[i]
        for line in open(filename, 'r'):
            if line.strip():
                read_id, type, mapped_id, cutloc, seq = line.strip().split()
                if type not in choice:
                    continue            
                cols=cutloc.split(':')
                if int(cols[1]) < 0 or int(cols[1]) >= ref_length:
                    continue
                side, cut = cols[0], int(cols[1])
                if side == 'L':
                    offset = -dyad_offset
                    #nt = id_seq[id][cut+1]
                    #weight = weight_dict[side][nt]
                else:
                    assert side == 'R'
                    offset = dyad_offset
                    #nt = id_seq[id][cut-1]
                    #weight = weight_dict[side][nt]
                if cut+offset < 0  or cut+offset >= ref_length:
                    continue
                if mapped_id not in dyadmaps:
                    dyadmaps[mapped_id] = [0.0]*ref_length
                dyadmaps[mapped_id][cut+offset] += 1.0 #/ weight
                if mapped_id not in cutmaps:
                    cutmaps[mapped_id] = {}
                    cutmaps[mapped_id]['L'], cutmaps[mapped_id]['R'] = [0.0]*ref_length, [0.0]*ref_length
                cutmaps[mapped_id][side][cut] += 1
    dcutmaps = copy.deepcopy(cutmaps)
    mcutmaps = {}
    for mapped_id in cutmaps:
        if mapped_id not in mcutmaps:
            mcutmaps[mapped_id] = {}
            mcutmaps[mapped_id]['L'], mcutmaps[mapped_id]['R'] = [np.NaN]*ref_length, [np.NaN]*ref_length
        win, loc = mapped_id.split('-')
        size, loc = len(win), int(loc)           
        if size <= 0:
            continue
        st, ed = loc, loc + size
        for i in range(st, ref_length):
            try:
                counts = cutmaps[mapped_id]['L'][i + 2*dyad_offset]
                dyadmaps[mapped_id][i+dyad_offset] += counts
            except:
                counts = 0.0
            mcutmaps[mapped_id]['R'][i] = counts
            dcutmaps[mapped_id]['R'][i] = np.NaN
        for i in range(0, ed):
            try:
                counts = cutmaps[mapped_id]['R'][i - 2*dyad_offset]
                dyadmaps[mapped_id][i-dyad_offset] += counts
            except:
                counts = 0.0
            mcutmaps[mapped_id]['L'][i] = counts
            dcutmaps[mapped_id]['L'][i] = np.NaN
    print >> sys.stderr, "Reading is done"
    return dcutmaps, mcutmaps, dyadmaps

"""
def MetroHasting (alpha, beta, A, B_list, C_list, D_list, N=1000, sigma=100, mode=1):
    def ratio_func1 (x, y): # beta-like group
        if x<=0 or x>=1:
            return 0
        output = -A*(x-y) + (alpha-1)*(np.log(x)-np.log(y)) + (beta-1)*(np.log(1-x)-np.log(1-y))
        for i in range(len(B_list)):
            if D_list[i] == 0:
                continue
            output += D_list[i]*(np.log(B_list[i]*x + C_list[i]) - np.log(B_list[i]*y + C_list[i]))
        return np.exp(output)
    def ratio_func2 (x, y): # gamma-like group
        if x<=0:
            return 0
        output = -beta*(x-y) - A*(x-y) + (alpha-1)*(np.log(x)-np.log(y))
        for i in range(len(B_list)):
            if D_list[i] == 0:
                continue
            output += D_list[i]*(np.log(B_list[i]*x + C_list[i]) - np.log(B_list[i]*y + C_list[i]))
        return np.exp(output)
    def update (mu, sigma):
        return random.gauss(mu=mu, sigma=sigma)
    if mode == 1:
        prev = 0.01
        sigma = 0.1
        N = 1000
    elif mode == 2:
        if beta + A >0:
            prev = float(alpha)/(beta + A)
        else:
            prev = 10
        sigma = 100
        N = 1000
    samples = []
    prev = 0.01
    while len(samples) < N:
        current = update(mu=prev, sigma=sigma)
        #A = min(1, float(prob_func(current))/prob_func(prev))
        if mode == 1:
            threshold = min(1, ratio_func1(current, prev))
        elif mode == 2:
            threshold = min(1, ratio_func2(current, prev))
        if random.uniform(0,1) <= threshold:
            samples.append(current)
            prev = copy.deepcopy(current)
        else:
            samples.append(prev)
    return samples
"""

# beta-like group
def beta_ABCD (i, dyad, free, efficiency, cleavage, mcleavage, dyad_offset, mode=1):
    A = 0
    Blist, Clist, Dlist = [], [], []
    for j in range(len(dyad)):
        if i+dyad_offset < 0 or i+dyad_offset >= len(dyad[j]):
            dyadsig = 0.0
        else:
            dyadsig = dyad[j][i+dyad_offset]
        A += dyadsig
        Blist.append(dyadsig)
        total = 0.0
        for k in range(len(dyad[j])):
            if k == i + dyad_offset:
                continue
            total += dyad[j][k]
        total += free[j][0]
        Clist.append(total)
        Cut = cleavage[j][i]
        if np.isnan(Cut):
            Cut = mcleavage[j][i]
        #print j, i
        assert not np.isnan(Cut)
        Dlist.append(Cut)
    if mode == 1:
        Clist = [Clist[u]*efficiency[i] for u in range(len(Clist))]
        return A, Blist, Clist, Dlist
    else:
        Blist = [Blist[u]*efficiency[i] for u in range(len(Blist))]
        return sum(Clist), Clist, Blist, Dlist

#gamma-like group
def gamma_ABCD (j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset): 
    if i - dyad_offset < 0:
        t, ep, T = 0.0, 0.0, 0.0
    else:
        t = tlist[i-dyad_offset]
        ep = eplist[i-dyad_offset]
        T = top[j][i-dyad_offset]
    if i + dyad_offset >= len(dyad[j]):
        b, de, B = 0.0, 0.0, 0.0
    else:
        b = blist[i+dyad_offset]
        de = delist[i+dyad_offset]
        B = bott[j][i+dyad_offset]
    if np.isnan(T):
        T = mtop[j][i-dyad_offset]
    if np.isnan(B):
        B = mbott[j][i+dyad_offset]
    A = t + b
    Blist = [t, b]
    total = 0.0
    for u in range(len(dyad[j])):
        if u == i:
            continue
        total += dyad[j][u]
    total += free[j][0]
    Clist = [total*ep, total*de]
    Dlist = [T, B]
    return A, Blist, Clist, Dlist

# beta-like group
def accr_beta (x, y, A, B_list, C_list, D_list, alpha, beta, a=0, b=1):
    if x<=a or x>=b:
        return 0
    output = -A*(x-y) + (alpha-1)*(np.log(x-a)-np.log(y-a)) + (beta-1)*(np.log(b-x)-np.log(b-y))
    for i in range(len(B_list)):
        if D_list[i] == 0:
            continue
        output += D_list[i]*(np.log(B_list[i]*x + C_list[i]) - np.log(B_list[i]*y + C_list[i]))
    return np.exp(output)

# gamma-like group
def accr_gamma (x, y, A, B_list, C_list, D_list, alpha, beta):
    if x<=0:
        return 0
    output = -beta*(x-y) - A*(x-y) + (alpha-1)*(np.log(x)-np.log(y))
    for i in range(len(B_list)):
        if D_list[i] == 0:
            continue
        output += D_list[i]*(np.log(B_list[i]*x + C_list[i]) - np.log(B_list[i]*y + C_list[i]))
    return np.exp(output)

def MetroHasting (N, initial_pt, sigma, accr_func, *args, **kwargs):
    def update (mu, sigma):
        return random.gauss(mu=mu, sigma=sigma)
    samples = []
    prev = initial_pt
    #success = 0
    while len(samples) < N:
        current = update(mu=prev, sigma=sigma)
        threshold = min(1, accr_func(current, prev, *args, **kwargs))
        if random.uniform(0,1) <= threshold:
            #success += 1
            samples.append(current)
            prev = copy.deepcopy(current)
        else:
            samples.append(prev)
    #print float(success)/N
    return samples

def t_conditional(i, dyad, free, eplist, top, mtop, dyad_offset, alpha, beta, a=0, b=1):
    #return 1.0
    A, Blist, Clist, Dlist = beta_ABCD(i, dyad, free, eplist, top, mtop, dyad_offset)
    #newt = beta_conditional(alpha, beta, A, Blist, Clist, Dlist, N=1)[0]
    newt = MetroHasting (500, b/100.0, b/10.0, accr_beta, A, Blist, Clist, Dlist, alpha, beta, a=a, b=b)
    #print "t"
    #plt.plot(acf(newt))
    #plt.show()
    #print acf(newt)[:3]
    #print acf(newt)[-3:]
    #newt = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, mode=1)
    #print acf(newt)[-1]
    newt = newt[-1]
    return newt

def b_conditional(i, dyad, free, delist, bott, mbott, dyad_offset, alpha, beta, a=0, b=1):
    #return 0.5
    A, Blist, Clist, Dlist = beta_ABCD(i, dyad, free, delist, bott, mbott, -dyad_offset)
    #newb = beta_conditional(alpha, beta, A, Blist, Clist, Dlist, N=1)[0]
    newb = MetroHasting (500, b/100.0, b/10.0, accr_beta, A, Blist, Clist, Dlist, alpha, beta, a=a, b=b)
    #print "b"
    #print acf(newb)[:3]
    #print acf(newb)[-3:]
    #newb = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, mode=1)
    #print acf(newb)[-1]
    newb = newb[-1]
    return newb

def ep_conditional(i, dyad, free, tlist, top, mtop, dyad_offset, alpha, beta, a=0, b=0.01):
    #return 0.0
    A, Blist, Clist, Dlist = beta_ABCD (i, dyad, free, tlist, top, mtop, dyad_offset, mode=2)
    #newep = beta_conditional(alpha, beta, A, Blist, Clist, Dlist, N=1)[0]
    newep = MetroHasting (500, b/100.0, b/10.0, accr_beta, A, Blist, Clist, Dlist, alpha, beta, a=a, b=b)
    #newep = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, mode=1)
    #plt.plot(acf(newep))
    #plt.show()
    #print acf(newep)[-1]
    #print 
    newep = newep[-1]
    return newep

def de_conditional(i, dyad, free, blist, bott, mbott, dyad_offset, alpha, beta, a=0, b=0.01):
    #return 0.0
    A, Blist, Clist, Dlist = beta_ABCD (i, dyad, free, blist, bott, mbott, -dyad_offset, mode=2)
    #newde = beta_conditional(alpha, beta, A, Blist, Clist, Dlist, N=1)[0]
    newde = MetroHasting (500, b/100.0, b/10.0, accr_beta, A, Blist, Clist, Dlist, alpha, beta, a=a, b=b)
    #newde = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, mode=1)
    #plt.plot(acf(newde))
    #plt.show()
    #print acf(newde)[-1]
    newde = newde[-1]
    return newde

def D_conditional(j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset, alpha, beta):
    A, Blist, Clist, Dlist = gamma_ABCD (j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset)
    #newD = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, mode=2)
    if beta + A > 0:
        initial_pt = float(alpha)/(beta + A)
    else:
        initial_pt = 10
    newD = MetroHasting (1000, initial_pt, 100, accr_gamma, A, Blist, Clist, Dlist, alpha, beta)
    #print "D"
    #plt.plot(acf(newD))
    #plt.show()
    #print acf(newD)[3:]
    #print acf(newD)[-3:]
    #print acf(newD)[-1]
    newD = newD[-1]
    return newD

def T_conditional(j, i, dyad, free, tlist, eplist, dyad_offset):
    if i + dyad_offset >= len(dyad[j]):
        dyadsig = 0.0
    else:
        dyadsig = dyad[j][i+dyad_offset]
    lamb = dyadsig*tlist[i]
    total = 0.0
    for u in range(len(dyad[j])):
        if u == i + dyad_offset:
            continue
        total += dyad[j][u]
    total += free[j][0]
    total *= eplist[i]
    lamb += total
    return np.random.poisson(lamb,1)[0]

def B_conditional(j, i, dyad, free, blist, delist, dyad_offset):
    if i - dyad_offset < 0:
        dyadsig = 0.0
    else:
        dyadsig = dyad[j][i-dyad_offset]
    lamb = dyadsig*blist[i]
    total = 0.0
    for u in range(len(dyad[j])):
        if u == i - dyad_offset:
            continue
        total += dyad[j][u]
    total += free[j][0]
    total *= delist[i]
    lamb += total
    return np.random.poisson(lamb,1)[0]

def Gibbs (N, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset):
    t = open("tlist.txt", 'w')
    b = open("blist.txt", 'w')
    ep = open("eplist.txt", 'w')
    de = open("delist.txt", 'w')
    d = open("dyad.txt", 'w')
    mT = open("mtop.txt", 'w')
    mB = open("mbott.txt", 'w')    
    for k in range(N):
        print >> t, '>' + str(k)
        print >> b, '>' + str(k)
        print >> ep, '>' + str(k)
        print >> de, '>' + str(k)
        print >> d, '>' + str(k)
        print >> mT, '>' + str(k)
        print >> mB, '>' + str(k)
        s = ""
        for i in range(len(tlist)):
            if k > 0 and i >= st - dyad_offset and i < ed - dyad_offset:
                tlist[i] = t_conditional(i, dyad, free, eplist, top, mtop, dyad_offset, 3, 3)
            s += str(tlist[i]) + ','
            #tlist_List.append(copy.deepcopy(tlist))
            #print 't', i
        print >> t, s[:-1]
        s = ""
        for i in range(len(blist)):
            if k > 0 and i >= st + dyad_offset and i < ed + dyad_offset:
                blist[i] = b_conditional(i, dyad, free, delist, bott, mbott, dyad_offset, 3, 3)
            s += str(blist[i]) + ','
            #blist_List.append(copy.deepcopy(blist))
            #print 'b', i
        print >> b, s[:-1]
        s = ""
        for i in range(len(eplist)):
            if k > 0:
                eplist[i] = ep_conditional(i, dyad, free, tlist, top, mtop, dyad_offset, 1, 10)
            s += str(eplist[i]) + ','
            #eplist_List.append(copy.deepcopy(eplist))
            #print 'ep', i
        print >> ep, s[:-1]
        s = ""
        for i in range(len(delist)):
            if k > 0:
                delist[i] = de_conditional(i, dyad, free, blist, bott, mbott, dyad_offset, 1, 10)
            s += str(delist[i]) + ','
            #delist_List.append(copy.deepcopy(delist))
            #print 'de', i
        print >> de, s[:-1]
        for j in range(len(dyad)):
            s = ""
            for i in range(len(dyad[j])):
                if k > 0 and i >= st and i < ed:
                    dyad[j][i] = D_conditional(j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset, alphas[j][i], betas[j][i])
                    #dyad[j][i] = D_conditional(j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset, 1, 0)
                s += str(dyad[j][i]) + ','
                #dyad_List.append(copy.deepcopy(dyad))
                #print 'dyad', j, i
            print >> d, s[:-1]        
        for j in range(len(mtop)):
            s1, s2 = "", ""
            for i in range(ref_length):
                if not np.isnan(mtop[j][i]):
                    if k > 0:
                        mtop[j][i] = T_conditional(j, i, dyad, free, tlist, eplist, dyad_offset)
                    s1 += str(mtop[j][i]) + ','
                    #mtop_List.append(copy.deepcopy(mtop))
                    #print 'mT', j,i
                if not np.isnan(mbott[j][i]):
                    if k > 0 :
                        mbott[j][i] = B_conditional(j, i, dyad, free, blist, delist, dyad_offset)
                    s2 += str(mbott[j][i]) + ','
                    #mbott_List.append(copy.deepcopy(mbott))
                    #print 'mB', j,i
            if s1:
                print >> mT, s1[:-1]
            if s2:
                print >> mB, s2[:-1]
        
        print >> sys.stderr, k
        #print "dyad", dyad[0][225/2]
        #print "mtop", mtop[0][225/2 - 52]
        #print "bott", bott[0][225/2 + 52]
        #print 
    t.close()
    b.close()
    ep.close()
    de.close()
    d.close()
    mT.close()
    mB.close()
    return

ref_length = 225
dyad_offset = 52
#left_bound, right_bound = 225/2, 225/2
left_bound, right_bound = 147/2, 147/2
#left_bound, right_bound = 0, 0
st, ed = left_bound, ref_length - right_bound
#cutmaps, mcutmaps, dyadmaps = read_file(['/home/spark159/../../media/spark159/sw/Ascan-5min_S1_L001_R.sort'])
cutmaps, mcutmaps, dyadmaps = read_file(['batest15.sort'])
top, bott = [], []
mtop, mbott = [], []
dyad = []
free = []
for key in sorted(cutmaps.keys(), cmp=ID_cmp):
    #if not (key.startswith('AAA-46') or key.startswith("AAAAAAAAAAAA-157")):
    #    continue
    top.append(cutmaps[key]['R'])
    bott.append(cutmaps[key]['L'])
    mtop.append(mcutmaps[key]['R'])
    mbott.append(mcutmaps[key]['L'])
    dyad.append([0.0]*left_bound + dyadmaps[key][st:ed] + [0.0]*right_bound)
    free.append([0])
tlist, blist, eplist, delist = [0.5]*ref_length, [0.5]*ref_length, [0.0]*ref_length, [0.0]*ref_length

# get prior for dyad
alphas, betas = [], []
for j in range(len(dyad)):
    var = np.std(dyad[j])**2
    temp1, temp2 = [], []
    for i in range(len(dyad[j])):
        mean = dyad[j][i]
        temp1.append(float(mean**2)/var)
        temp2.append(float(mean)/var)
    alphas.append(temp1)
    betas.append(temp2)

#top = np.nan_to_num(np.asarray(top)) + np.nan_to_num(np.asarray(mtop))
#bott = np.nan_to_num(np.asarray(bott)) + np.nan_to_num(np.asarray(mbott))
#top = np.nan_to_num(np.asarray(top))
#bott = np.nan_to_num(np.asarray(bott))

f = open("top.txt", 'w')
print >> f, ">0"
for i in range(len(top)):
    #print >> f, '>' + str(i)
    s = ""
    for j in range(len(top[i])):
        if np.isnan(top[i][j]):
            continue
        s += str(top[i][j]) + ','
    print >> f, s[:-1]
f.close()

f = open("bott.txt", 'w')
print >> f, ">0"
for i in range(len(bott)):
    #print >> f, '>' + str(i)
    s = ""
    for j in range(len(bott[i])):
        if np.isnan(bott[i][j]):
            continue
        s += str(bott[i][j]) + ','
    print >> f, s[:-1]
f.close() 
Gibbs (1000, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset)
