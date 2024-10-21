import sys
import random
import math
import copy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def rev_comp(seq):
    rev = ""
    for nt in seq[::-1]:
        nt = nt.upper()
        if nt == 'A':
            rev += 'T'
        elif nt == 'T':
            rev += 'A'
        elif nt == 'G':
            rev += 'C'
        elif nt == 'C':
            rev += 'G'
    return rev

def ID_cmp (a, b):
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

def save_data (top,
               bott,
               mtop_list,
               mbott_list,
               dyad_list,
               free_list,
               tlist_list,
               blist_list,
               note = '',
               eplist_list=None,
               delist_list=None):

    def write_file (data_list, name):
        f = open(name, 'w')
        for i in range(len(data_list)):
            print >> f, '>' + str(i)
            s = ""
            data = data_list[i]
            if type(data[0]) == list:
                for j in range(len(data)):
                    for k in range(len(data[j])):
                        s += str(data[j][k]) + ','
            else:
                for j in range(len(data)):
                    s += str(data[j]) + ','                
            print >> f, s[:-1]
        f.close()

    write_file(top, "top" + note + ".txt")
    write_file(bott, "bott" + note + ".txt")
    write_file(mtop_list, "mtop" + note + ".txt")
    write_file(mbott_list, "mbott" + note + ".txt")
    write_file(dyad_list, "dyad" + note + ".txt")
    write_file(free_list, "free" + note + ".txt")
    write_file(tlist_list, "tlist" + note + ".txt")
    write_file(blist_list, "blist" + note + ".txt")

    if eplist_list != None:
        write_file(eplist_list, "eplist" + note + ".txt")
    if delist_list != None:
        write_file(delist_list, "delist" + note + ".txt")

    print >> sys.stderr, "data are saved"
    return

def plot_file (fname, seq_id=0, total=True):    
    sig_list = []
    for line in open(fname):
        if line.startswith('>'):
            count = 0
            continue
        if count == seq_id :
            cols = line.strip().split(',')
            temp = []
            for value in cols:
                try:
                    value = float(value)
                except:
                    value = np.nan
                temp.append(value)
            sig_list.append(temp)
        count +=1
    sig_list = np.asarray(sig_list)
    m,n = sig_list.shape
    #mode_sig = [stats.mode(sig_list[:,i]) for i in range(n)]
    median_sig = [np.median(sig_list[:,i]) for i in range(n)]
    x = range(len(median_sig))
    if total:
        for sig in sig_list:
            plt.plot(x, sig, 'y', linewidth=1, alpha=0.2)
    plt.plot(x, median_sig)
    plt.plot(x, sig_list[0], 'k--')

def plot_data (data_list, seq_id=0, total=True):
    sig_list = []
    for i in range(len(data_list)):
        data = data_list[i]
        if type(data[0]) == list:
            sig_list.append(data[seq_id])
        else:
            sig_list.append(data)
    sig_list = np.asarray(sig_list)
    m,n = sig_list.shape
    #mode_sig = [stats.mode(sig_list[:,i]) for i in range(n)]
    median_sig = [np.median(sig_list[:,i]) for i in range(n)]
    x = range(len(median_sig))
    if total:
        for sig in sig_list:
            plt.plot(x, sig, 'y', linewidth=1, alpha=0.2)
    plt.plot(x, median_sig)
    plt.plot(x, sig_list[0], 'k--')

def dyad_sampling (template, left_bound, right_bound, N):
    dyadmap = [0.0]*len(template)
    st, ed = left_bound, len(template) - 1 - right_bound
    for i in range(N):
        pos = random.randint(st,ed)
        dyadmap[pos] += 1
    return dyadmap

def simulation (id_seq,
                dyad,
                free,
                tlist,
                blist,
                eplist,
                delist,
                dyad_offset):
                
    def get_lambda(template, dyadmap, free_num, tlist, blist, eplist, delist, offset):
        tlamb, blamb = [0.0]*len(template), [0.0]*len(template)
        for i in range(len(template)):
            if i + offset < len(template):
                tlamb[i] += dyadmap[i+offset]*tlist[i]
            if i -  offset >= 0:
                blamb[i] += dyadmap[i-offset]*blist[i]
            total1, total2 = 0.0, 0.0
            for j in range(len(template)):
                if j != i + offset:
                    total1 += dyadmap[j]
                if j != i - offset:
                    total2 += dyadmap[j]
            tlamb[i] += (total1 + free_num)*eplist[i]
            blamb[i] += (total2 + free_num)*delist[i]
        return tlamb, blamb

    def make_cleavage(template, tlamb, blamb):
        tnum_list, bnum_list = [], []
        frags_list = []
        for i in range(len(template)):
            tnum = np.random.poisson(tlamb[i],1)[0]
            bnum = np.random.poisson(blamb[i],1)[0]
            tnum_list.append(tnum)
            bnum_list.append(bnum)
            right_frag = template[i:]
            left_frag = rev_comp(template[:i+1])
            for k in range(tnum):
                frags_list.append(right_frag)
            for k in range(bnum):
                frags_list.append(left_frag)
        return tnum_list, bnum_list, frags_list
    
    def mutations (seq_list, mist_prob, indel_prob):
        def mismatch(seq, prob):
            nts = set(['A','T','C','G'])
            new_seq = ""
            for i in range(len(seq)):
                if random.random() < prob:
                    subs = nts - set(seq[i])
                    new_seq += random.choice(list(subs))
                else:
                    new_seq += seq[i]
            return new_seq
        def indel(seq, prob):
            nts = ['A','T','C','G']
            new_seq = ""
            for i in range(len(seq)):
                if random.random() < prob:
                    if random.random() < 0.5 or i == len(seq):
                        new_seq += random.choice(nts)
                    else:
                        continue
                new_seq += seq[i]
            return new_seq
        new_list = []
        for seq in seq_list:
            new_seq = indel(mismatch(seq, prob=mist_prob), prob=indel_prob)
            new_list.append(new_seq)
        return new_list

    ids = sorted(id_seq.keys(), cmp=ID_cmp)
    top, bott = [], []
    mtop, mbott = [], []
    for i in range(len(ids)):
        id = ids[i]
        template = id_seq[id]
        dyadmap = dyad[i]
        free_num = free[i][0]
        tlamb, blamb = get_lambda(template, dyadmap, free_num, tlist, blist, eplist, delist, offset=dyad_offset)
        tnum_list, bnum_list, frags_list = make_cleavage(template, tlamb, blamb)
        win, loc = id.split('-')
        st, ed = int(loc), int(loc)+len(win)
        top.append(tnum_list[:st] + [np.nan]*(len(template)-st))
        mtop.append([np.nan]*st + tnum_list[st:])
        bott.append([np.nan]*ed +  bnum_list[ed:])
        mbott.append(bnum_list[:ed] + [np.nan]*(len(template)-ed))
        #seqs_list = mutations(frags_list, mist_prob, indel_prob)
        
    print >> sys.stderr, "simulation is done"
    return top, bott, mtop, mbott

def imputation (N,
                id_seq,
                top,
                bott,
                dyad_offset,
                left_bound,
                right_bound,
                noise=False,
                dyad_alphas=None,
                dyad_betas=None,
                t_alphas=None,
                t_betas=None,
                b_alphas=None,
                b_betas=None,
                ep_alphas=None,
                ep_betas=None,
                de_alphas=None,
                de_betas=None):

    # beta-like group coefficients
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

    #gamma-like group coefficients
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

    # beta-like group acceptance ratio funcion
    def accr_beta (x, y, A, B_list, C_list, D_list, alpha, beta, a=0, b=1):
        if x<=a or x>=b:
            return 0
        output = -A*(x-y) + (alpha-1)*(np.log(x-a)-np.log(y-a)) + (beta-1)*(np.log(b-x)-np.log(b-y))
        for i in range(len(B_list)):
            if D_list[i] == 0:
                continue
            output += D_list[i]*(np.log(B_list[i]*x + C_list[i]) - np.log(B_list[i]*y + C_list[i]))
        return np.exp(output)

    # gamma-like group acceptance ratio function
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
        A, Blist, Clist, Dlist = beta_ABCD(i, dyad, free, eplist, top, mtop, dyad_offset)
        newt = MetroHasting (500, b/100.0, b/10.0, accr_beta, A, Blist, Clist, Dlist, alpha, beta, a=a, b=b)
        return newt[-1]

    def b_conditional(i, dyad, free, delist, bott, mbott, dyad_offset, alpha, beta, a=0, b=1):
        A, Blist, Clist, Dlist = beta_ABCD(i, dyad, free, delist, bott, mbott, -dyad_offset)
        newb = MetroHasting (500, b/100.0, b/10.0, accr_beta, A, Blist, Clist, Dlist, alpha, beta, a=a, b=b)
        return newb[-1]

    def ep_conditional(i, dyad, free, tlist, top, mtop, dyad_offset, alpha, beta, a=0, b=0.01):
        A, Blist, Clist, Dlist = beta_ABCD (i, dyad, free, tlist, top, mtop, dyad_offset, mode=2)
        newep = MetroHasting (500, b/100.0, b/10.0, accr_beta, A, Blist, Clist, Dlist, alpha, beta, a=a, b=b)
        return newep[-1]

    def de_conditional(i, dyad, free, blist, bott, mbott, dyad_offset, alpha, beta, a=0, b=0.01):
        A, Blist, Clist, Dlist = beta_ABCD (i, dyad, free, blist, bott, mbott, -dyad_offset, mode=2)
        newde = MetroHasting (500, b/100.0, b/10.0, accr_beta, A, Blist, Clist, Dlist, alpha, beta, a=a, b=b)
        return newde[-1]

    def D_conditional(j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset, alpha, beta):
        A, Blist, Clist, Dlist = gamma_ABCD (j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset)
        if beta + A > 0:
            initial_pt = float(alpha)/(beta + A)
        else:
            initial_pt = 10
        newD = MetroHasting (1000, initial_pt, 100, accr_gamma, A, Blist, Clist, Dlist, alpha, beta)
        return newD[-1]

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
        dyad_list, free_list = [], []
        tlist_list, blist_list = [], []
        mtop_list, mbott_list = [], []
        dyad_list.append(copy.deepcopy(dyad))
        free_list.append(copy.deepcopy(free))
        tlist_list.append(copy.deepcopy(tlist))
        blist_list.append(copy.deepcopy(blist))
        mtop_list.append(copy.deepcopy(mtop))
        mbott_list.append(copy.deepcopy(mbott))
        
        eplist_list, delist_list = [], []
        eplist_list.append(copy.deepcopy(eplist))
        delist_list.append(copy.deepcopy(delist))
        
        st, ed = left_bound, ref_length - right_bound
        for k in range(N):
            for i in range(st-dyad_offset, ed-dyad_offset):
                tlist[i] = t_conditional(i, dyad, free, eplist, top, mtop, dyad_offset, t_alphas[i], t_betas[i])
            tlist_list.append(copy.deepcopy(tlist))

            for i in range(st+dyad_offset, ed+dyad_offset):
                blist[i] = b_conditional(i, dyad, free, delist, bott, mbott, dyad_offset, b_alphas[i], b_betas[i])
            blist_list.append(copy.deepcopy(blist))

            if noise:
                for i in range(len(eplist)):
                    eplist[i] = ep_conditional(i, dyad, free, tlist, top, mtop, dyad_offset, ep_alphas[i], ep_betas[i])
                eplist_list.append(copy.deepcopy(eplist))
                for i in range(len(delist)):
                    delist[i] = de_conditional(i, dyad, free, blist, bott, mbott, dyad_offset, de_alphas[i], de_betas[i])
                delist_list.append(copy.deepcopy(delist))
                
            for j in range(len(dyad)):
                for i in range(st, ed):
                    dyad[j][i] = D_conditional(j, i, dyad, free, tlist, blist, eplist, delist, top, mtop, bott, mbott, dyad_offset, dyad_alphas[j][i], dyad_betas[j][i])
            dyad_list.append(copy.deepcopy(dyad))

            for j in range(len(mtop)):
                for i in range(ref_length):
                    if not np.isnan(mtop[j][i]):
                        mtop[j][i] = T_conditional(j, i, dyad, free, tlist, eplist, dyad_offset)
                    if not np.isnan(mbott[j][i]):
                        mbott[j][i] = B_conditional(j, i, dyad, free, blist, delist, dyad_offset)
            mtop_list.append(copy.deepcopy(mtop))
            mbott_list.append(copy.deepcopy(mbott))

            print >> sys.stderr, "cycle "  + str(k+1)

        return dyad_list, free_list, tlist_list, blist_list, eplist_list, delist_list, mtop_list, mbott_list


    # give initial guess
    pre_mtop, pre_mbott = [], []
    pre_dyad, pre_free = [], []
    ids = sorted(id_seq.keys(), cmp=ID_cmp)
    for i in range(len(ids)):
        id = ids[i]
        win, loc = id.split('-')
        st, ed = int(loc), int(loc) + len(win)
        temp1 = [np.nan]*st + bott[i][st+2*dyad_offset:]
        temp2 = top[i][:max(0, ed-2*dyad_offset)] + [np.nan]*(ref_length-ed)
        pre_mtop.append(temp1 + [0.0]*(ref_length - len(temp1)))
        pre_mbott.append([0.0]*(ref_length - len(temp2)) + temp2)
        full_top = list(np.nan_to_num(top[i]) + np.nan_to_num(pre_mtop[i]))
        full_bott = list(np.nan_to_num(bott[i]) + np.nan_to_num(pre_mbott[i]))
        tdyad = np.asarray([0.0]*dyad_offset + full_top[:ref_length-dyad_offset])
        bdyad = np.asarray(full_bott[dyad_offset:] + [0.0]*dyad_offset)
        pre_dyad.append(list(tdyad + bdyad))
        pre_free.append([0.0])
    pre_tlist, pre_blist = [0.5]*ref_length, [0.5]*ref_length
    pre_eplist, pre_delist = [0.0]*ref_length, [0.0]*ref_length

    # give prior paramters
    if dyad_alphas == None or dyad_betas == None:
        temp_alphas, temp_betas = [], []
    for j in range(len(pre_dyad)):
        var = np.std(pre_dyad[j])**2
        temp1, temp2 = [], []
        for i in range(len(pre_dyad[j])):
            mean = pre_dyad[j][i]
            temp1.append(float(mean**2)/var)
            temp2.append(float(mean)/var)
        temp_alphas.append(temp1)
        temp_betas.append(temp2)
    if dyad_alphas == None:
        dyad_alphas = temp_alphas
    if dyad_betas == None:
        dyad_betas = temp_betas
    if t_alphas == None:
        t_alphas = [3.0]*ref_length
    if t_betas == None:
        t_betas = [3.0]*ref_length
    if b_alphas == None:
        b_alphas = [3.0]*ref_length
    if b_betas == None:
        b_betas = [3.0]*ref_length
    if ep_alphas == None:
        ep_alphas = [1.0]*ref_length
    if ep_betas == None:
        ep_betas = [10.0]*ref_length
    if de_alphas == None:
        de_alphas = [1.0]*ref_length
    if de_betas == None:
        de_betas = [10.0]*ref_length
    
    # imputation by Gibbs sampling    
    return Gibbs (N, pre_dyad, pre_free, pre_tlist, pre_blist, pre_eplist, pre_delist, top, pre_mtop, bott, pre_mbott, dyad_offset)


# read reference
id_seq = {}
for line in open("polyAscanlib.ref"):
    if line.startswith('>AAAA-'):
        break
    if line.startswith('>'):
        id = line[1:].strip()
        continue
    id_seq[id] = line.strip()
    #break
ids = sorted(id_seq.keys(), cmp=ID_cmp)

# basic parameters
ref_length = 225
dyad_offset = 52
left_bound, right_bound = 147/2, 147/2

# true variables
dyad, free = [], []
for i in range(len(id_seq)):
    dyad.append([0.0]*left_bound + [1000]*(ref_length-left_bound-right_bound) + [0.0]*right_bound)
    free.append([0])
tlist, blist = [1.0]*ref_length, [0.5]*ref_length
eplist, delist = [0.0]*ref_length, [0.0]*ref_length

# simulate the data
top, bott, mtop, mbott = simulation (id_seq, dyad, free, tlist, blist, eplist, delist, dyad_offset)

# Gibbs sampling
dyad_list, free_list, tlist_list, blist_list, eplist_list, delist_list, mtop_list, mbott_list = imputation (100, id_seq, top, bott, dyad_offset, left_bound, right_bound)

# save result
save_data (top, bott, mtop_list, mbott_list, dyad_list, free_list, tlist_list, blist_list, note ="")
