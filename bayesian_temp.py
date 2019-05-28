import sys
import random
import math
import copy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#from scipy.special import gamma as gamma_func
#from scipy.stats import beta as beta_dist
#from scipy.stats import gamma as gamma_dist

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

def beta_conditional (alpha, beta, A, B_list, C_list, D_list, N=1):
    assert len(B_list) == len(C_list) == len(D_list)
    def prob_func (x):
        output = (x**(alpha-1))*((1-x)**(beta-1))*math.exp(-A*x)
        for i in range(len(B_list)):
            output *= (B_list[i]*x + C_list[i])**D_list[i]
        return output
    def env_func (x):
        output = (x**(alpha-1))*((1-x)**(beta-1))
        for i in range(len(B_list)):
            output *= (B_list[i] + C_list[i])**D_list[i]
        return output
    def ratio_func (x):
        output = -A*x
        for i in range(len(B_list)):
            output += D_list[i]*(math.log(B_list[i]*x + C_list[i]) - math.log(B_list[i] + C_list[i]))
        return np.exp(output)
    count = 0
    samples = []
    while len(samples) < N :
        count += 1
        t = beta_dist.rvs(alpha, beta)
        #cutoff = prob_func(t) / env_func(t)
        cutoff = ratio_func(t)
        assert cutoff <= 1
        if random.uniform(0,1) < cutoff:
            samples.append(t)
    print count
    return samples

def gamma_conditional (alpha, beta, A, B_list, C_list, D_list, N=1):
    assert len(B_list) == len(C_list) == len(D_list)
    def product(x):
        output = 1.0
        for i in range(len(B_list)):
            output *= (B_list[i]*x + C_list[i])**D_list[i]
        return output
    def prob_func (x):
        output = (x**(alpha-1))*math.exp(-beta*x)*math.exp(-A*x)*product(x)
        return output
    def env_func (x):
        prev = 0.001
        prev_value = math.exp(-A*prev)*product(prev)
        while True:
            current = prev*1.5
            current_value = math.exp(-A*current)*product(current)
            if current_value < prev_value:
                break
            prev = copy.deepcopy(current)
            prev_value = copy.deepcopy(current_value)
        print math.exp(-A*prev)*product(current)
        output = (x**(alpha-1))*math.exp(-beta*x)*math.exp(-A*prev)*product(current)
        return output
    count = 0
    samples = []
    while len(samples) < N :
        count +=1
        t = gamma_dist.rvs(a=alpha, scale=1.0/beta)
        cutoff = prob_func(t) / env_func(t)
        assert cutoff <= 1
        if random.uniform(0,1) < cutoff:
            samples.append(t)
    print count
    return samples

def MetroHasting (alpha, beta, A, B_list, C_list, D_list, N=1000, sigma=100, mode=1):
    def prob_func (x):
        if x <= 0:
            return 0
        def product(x):
            output = 1.0
            for i in range(len(B_list)):
                output *= (B_list[i]*x + C_list[i])**D_list[i]
            return output
        output = (x**(alpha-1))*math.exp(-beta*x)*math.exp(-A*x)*product(x)
        return output
    def ratio_func1 (x, y):
        if x<=0:
            return 0
        output = -beta*(x-y) - A*(x-y) + (alpha-1)*(math.log(x)-math.log(y))
        for i in range(len(B_list)):
            if D_list[i] == 0:
                continue
            output += D_list[i]*(math.log(B_list[i]*x + C_list[i]) - math.log(B_list[i]*y + C_list[i]))
        return np.exp(output)
    def ratio_func2 (x, y):
        if x<=0 or x>=1:
            return 0
        output = -A*(x-y) + (alpha-1)*(math.log(x)-math.log(y)) + (beta-1)*(math.log(1-x)-math.log(1-y))
        for i in range(len(B_list)):
            if D_list[i] == 0:
                continue
            output += D_list[i]*(math.log(B_list[i]*x + C_list[i]) - math.log(B_list[i]*y + C_list[i]))
        return np.exp(output)
    def update (mu, sigma):
        return random.gauss(mu=mu, sigma=sigma)
    samples = []
    #prev = float(alpha)/(beta + A)
    prev = 0.01
    while len(samples) < N:
        current = update(mu=prev, sigma=sigma)
        #A = min(1, float(prob_func(current))/prob_func(prev))
        if mode == 1:
            A = min(1, ratio_func1(current, prev))
        else:
            A = min(1, ratio_func2(current, prev))
        if random.uniform(0,1) <= A:
            samples.append(current)
            prev = copy.deepcopy(current)
        else:
            samples.append(prev)
    return samples

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
        Lrange = range(min(st + 2*dyad_offset, ref_length), ref_length)
        Rrange = range(0, max(ed -1 - 2*dyad_offset, -1) + 1)
        for i in Lrange:
            counts = cutmaps[mapped_id]['L'][i]
            cutmaps[mapped_id]['R'][i - 2*dyad_offset] = np.NaN
            mcutmaps[mapped_id]['R'][i - 2*dyad_offset] = counts
            dyadmaps[mapped_id][i - dyad_offset] += counts
        for i in Rrange:
            counts = cutmaps[mapped_id]['R'][i]
            cutmaps[mapped_id]['L'][i + 2*dyad_offset] = np.NaN
            mcutmaps[mapped_id]['L'][i + 2*dyad_offset] = counts
            dyadmaps[mapped_id][i + dyad_offset] += counts

    return cutmaps, mcutmaps, dyadmaps

def get_ABCD (i, dyad, free, efficiency, cleavage, mcleavage, dyad_offset, mode=1):
    A = 0
    Blist, Clist, Dlist = [], [], []
    for j in range(len(dyad)):
        try:
            dyadsig = dyad[j][i+dyad_offset]
        except:
            dyadsig = 0            
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

def get_ABCD_dyad (j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset):
    try:
        t = tlist[i-dyad_offset]
        ep = eplist[i-dyad_offset]
        T = top[j][i-dyad_offset]
    except:
        t = 0.0
        ep = 0.0
        T = 0.0
    try:
        b = blist[i+dyad_offset]
        de = delist[i+dyad_offset]
        B = bott[j][i+dyad_offset]
    except:
        b = 0.0
        de = 0.0
        B = 0.0
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

def t_conditional(i, dyad, free, eplist, top, mtop, dyad_offset, alpha, beta):
    A, Blist, Clist, Dlist = get_ABCD (i, dyad, free, eplist, top, mtop, dyad_offset)
    #newt = beta_conditional(alpha, beta, A, Blist, Clist, Dlist, N=1)[0]
    newt = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, N=1000, sigma=0.1, mode=2)
    #print acf(newt)[-1]
    newt = newt[-1]
    return newt

def b_conditional(i, dyad, free, delist, bott, mbott, dyad_offset, alpha, beta):
    A, Blist, Clist, Dlist = get_ABCD (i, dyad, free, delist, bott, mbott, -dyad_offset)
    #newb = beta_conditional(alpha, beta, A, Blist, Clist, Dlist, N=1)[0]
    newb = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, N=1000, sigma=0.1, mode=2)
    #print acf(newb)[-1]
    newb = newb[-1]
    return newb

def ep_conditional(i, dyad, free, tlist, top, mtop, dyad_offset, alpha, beta):
    A, Blist, Clist, Dlist = get_ABCD (i, dyad, free, tlist, top, mtop, dyad_offset, mode=2)
    #newep = beta_conditional(alpha, beta, A, Blist, Clist, Dlist, N=1)[0]
    newep = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, N=1000, sigma=0.1, mode=2)
    #print acf(newep)[-1]
    newep = newep[-1]
    return newep

def de_conditional(i, dyad, free, blist, bott, mbott, dyad_offset, alpha, beta):
    A, Blist, Clist, Dlist = get_ABCD (i, dyad, free, blist, bott, mbott, -dyad_offset, mode=2)
    #newde = beta_conditional(alpha, beta, A, Blist, Clist, Dlist, N=1)[0]
    newde = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, N=1000, sigma=0.1, mode=2)
    #print acf(newde)[-1]
    newde = newde[-1]
    return newde

def D_conditional(j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset, alpha, beta):
    A, Blist, Clist, Dlist = get_ABCD_dyad (j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset)
    newD = MetroHasting (alpha, beta, A, Blist, Clist, Dlist, N=1000, sigma=100)
    #print acf(newD)[-1]
    newD = newD[-1]
    return newD

def T_conditional(j, i, dyad, free, tlist, eplist, dyad_offset):
    try:
        dyadsig = dyad[j][i+dyad_offset]
    except:
        dyadsig = 0.0
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
    try:
        dyadsig = dyad[j][i-dyad_offset]
    except:
        dyadsig = 0.0
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
    dyad_List, free_List = [], []
    tlist_List, blist_List = [], []
    eplist_List, delist_List = [], []
    mtop_List, mbott_List = [], []
    for i in range(len(tlist)):
        tlist[i] = t_conditional(i, dyad, free, eplist, top, mtop, dyad_offset, 1, 1)
        tlist_List.append(copy.deepcopy(tlist))
        #print 't', i
    for i in range(len(blist)):
        blist[i] = b_conditional(i, dyad, free, delist, bott, mbott, dyad_offset, 1, 1)
        blist_List.append(copy.deepcopy(blist))
        #print 'b', i
    for i in range(len(eplist)):
        eplist[i] = ep_conditional(i, dyad, free, tlist, top, mtop, dyad_offset, 1, 1)
        eplist_List.append(copy.deepcopy(eplist))
        #print 'ep', i
    for i in range(len(delist)):
        delist[i] = de_conditional(i, dyad, free, blist, bott, mbott, dyad_offset, 1, 1)
        delist_List.append(copy.deepcopy(delist))
        #print 'de', i
    for j in range(len(dyad)):
        for i in range(len(dyad[j])):
            dyad[j][i] = D_conditional(j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset, 1, 1)
            dyad_List.append(copy.deepcopy(dyad))
            #print 'dyad', j, i
    for j in range(len(mtop)):
        for i in range(len(mtop[j])):
            if not np.isnan(mtop[j][i]):
                mtop[j][i] = T_conditional(j, i, dyad, free, tlist, eplist, dyad_offset)
                mtop_List.append(copy.deepcopy(mtop))
                #print 'mT', j,i
            if not np.isnan(mbott[j][i]):
                mbott[j][i] = B_conditional(j, i, dyad, free, blist, delist, dyad_offset)
                mbott_List.append(copy.deepcopy(mbott))
                #print 'mB', j,i
    return dyad_List, tlist_List, blist_List, eplist_List, delist_List, mtop_List, mbott_List

ref_length = 225
dyad_offset = 52
cutmaps, mcutmaps, dyadmaps = read_file(['/home/spark159/../../media/spark159/sw/Ascan-5min_S1_L001_R.sort'])
top, bott = [], []
mtop, mbott = [], []
dyad = []
free = []
for key in sorted(cutmaps.keys(), cmp=ID_cmp):
    if not key.startswith('AAA-'):
        break
    top.append(cutmaps[key]['R'])
    bott.append(cutmaps[key]['L'])
    mtop.append(mcutmaps[key]['R'])
    mbott.append(mcutmaps[key]['L'])
    dyad.append(dyadmaps[key])
    free.append([0])
tlist, blist, eplist, delist = [0.5]*ref_length, [0.5]*ref_length, [0.1]*ref_length, [0.1]*ref_length

dyad_List, tlist_List, blist_List, eplist_List, delist_List, mtop_List, mbott_List = Gibbs (1, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset)


"""
alpha = 10
beta = 20
A = 0
B_list = [0.1,0.2]
C_list = [0.3,0.4]
D_list = [20,17]
#print beta_conditional (alpha, beta, A, B_list, C_list, D_list, N=1)
#print gamma_conditional (alpha, beta, A, B_list, C_list, D_list, N=1)
sample1 =  MetroHasting(alpha, beta, A, B_list, C_list, D_list, sigma=0.1, mode=2)
sample2 =  MetroHasting(alpha, beta, A, B_list, C_list, D_list, sigma=1, mode=2)
sample3 =  MetroHasting(alpha, beta, A, B_list, C_list, D_list, sigma=10, mode=2)
fig = plt.figure()
plt.plot(acf(sample1), label='small')
plt.plot(acf(sample2), label='middle')
plt.plot(acf(sample3), label='large')
plt.legend()
plt.show()

X = np.linspace(0.0001,100, num=1000)
#Y = [gamma_dist.pdf(x, a=alpha, scale=1.0/beta) for x in X]
#Y = [prob_func(x) for x in X]
#Z = [env_func(x) for x in X]
#D = [math.exp(-A*x)* product(x) for x in X]
#D = [ gamma_dist.rvs(a=alpha, scale=1.0/beta) for i in range(10000)]
D = [ np.random.beta(alpha,beta) for i in range(10000)]
fig = plt.figure()
sns.kdeplot(MetroHasting(alpha, beta, A, B_list, C_list, D_list, N=10000, sigma=0.1, mode=2), bw=0.01)
sns.kdeplot(D, bw=0.01)
#plt.plot(X,Y)
#plt.plot(X,Z, '--')
#plt.plot(X,D, '--')
plt.xlim([0,1])
#plt.ylim([0,100])
plt.show()

"""
