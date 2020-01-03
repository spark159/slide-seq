import sys
import random
import math
import copy
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
from scipy.special import comb
from scipy.special import betainc
from scipy.special import beta as betaf

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
                s = s[:-1] +  '\n'
        else:
            for j in range(len(data)):
                s += str(data[j]) + ','                
        print >> f, s[:-1]
    f.close()

def plot_file (fname, seq_id=0, total=True, label=None):    
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
    plt.plot(x, median_sig, label=label)
    plt.plot(x, sig_list[0], 'k--')
    return median_sig

def load_data (fname):    
    sig_list = []
    check = False
    for line in open(fname):
        if line.startswith('>'):
            if check:
                if len(cycle) <= 1:
                    sig_list.append(cycle[0])
                else:
                    sig_list.append(cycle)
            cycle = []
            check = True
            continue
        row = []
        cols = line.strip().split(',')
        for value in cols:
            try:
                value = float(value)
            except:
                value = np.nan
            row.append(value)
        cycle.append(row)
    if check:
        if len(cycle) <= 1:
            sig_list.append(cycle[0])
        else:
            sig_list.append(cycle)
    return sig_list

def plot_data (data_list, seq_id=0, total=True, label=None):
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
    plt.plot(x, median_sig, label=label)
    plt.plot(x, sig_list[0], 'k--')

def rectangular (length, height):
    output = []
    for i in range(length):
        output.append(height)
    return output

def triangular (length, height):
    output = []
    mu = (length-1)/2.0
    for i in range(length):
        if i <= mu:
            value = (height*i)/float(mu)
        else:
            value = -(height*i)/float(mu) + 2*height
        output.append(value)
    return output

def gaussian (length, height):
    sigma = length/6.0
    mu = (length-1)/2.0
    output = []
    for i in range(length):
        value = height*math.exp(-((i-mu)**2)/(2*(sigma**2)))
        output.append(value)
    return output

def random_sig (length, height):
    output = []
    for i in range(length):
        output.append(random.uniform(0, height))
    return output

def simulation (id_seq,
                dyad,
                tlist,
                blist,
                eplist,
                delist,
                dyad_offset):
                
    def get_lambda(template, dyadmap, tlist, blist, eplist, delist, offset):
        tlamb, blamb = [0.0]*len(template), [0.0]*len(template)
        for i in range(len(template)):
            if i + offset < len(template):
                tlamb[i] += dyadmap[i+offset]*tlist[i]
            if i -  offset >= 0:
                blamb[i] += dyadmap[i-offset]*blist[i]
            tlamb[i] += eplist[i]
            blamb[i] += delist[i]
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
        tlamb, blamb = get_lambda(template, dyadmap, tlist, blist, eplist, delist, offset=dyad_offset)
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

def Naive (id_seq,
           top,
           bott,
           ref_length,
           dyad_offset):

    mtop, mbott = [], []
    dyad = []
    ids = sorted(id_seq.keys(), cmp=ID_cmp)
    for i in range(len(ids)):
        id = ids[i]
        win, loc = id.split('-')
        st, ed = int(loc), int(loc) + len(win)
        temp1 = [np.nan]*st + bott[i][st+2*dyad_offset:]
        temp2 = top[i][:max(0, ed-2*dyad_offset)] + [np.nan]*(ref_length-ed)
        mtop.append(temp1 + [0.0]*(ref_length - len(temp1)))
        mbott.append([0.0]*(ref_length - len(temp2)) + temp2)
        full_top = list(np.nan_to_num(top[i]) + np.nan_to_num(mtop[i]))
        full_bott = list(np.nan_to_num(bott[i]) + np.nan_to_num(mbott[i]))
        tdyad = np.asarray([0.0]*dyad_offset + full_top[:ref_length-dyad_offset])
        bdyad = np.asarray(full_bott[dyad_offset:] + [0.0]*dyad_offset)
        dyad.append(list((tdyad + bdyad)*0.5))
        
    print >> sys.stderr, "Naive imputation is done"
    return [dyad], [mtop], [mbott]    

def Linear_reg (id_seq,
                top,
                bott,
                ref_length,
                dyad_offset,
                graph=False):

    def robust_fit(X, Y):
        mX = [[x] for x in X]
        mX = np.asarray(mX)
        mY = np.asarray(copy.deepcopy(Y))
        model = linear_model.HuberRegressor()
        model.fit(mX, mY)
        inmask = np.array([b == False for b in model.outliers_])
        mX = mX[inmask]
        mY = mY[inmask]
        coefs = [model.coef_[0], model.intercept_]
        #RSS = 0.0
        #for i in range(len(mX)):
        #    x, y = mX[i][0], mY[i]
        #    RSS += (y - coefs[0]*x - coefs[1])**2
        rsquare = model.score(mX, mY)
        return coefs, rsquare

    def zero_frac (L):
        num = 0.0
        for e in L:
            if e == 0:
                num +=1
        return num/len(L)

    # collect observable T/B counts pairs
    ids = sorted(id_seq.keys(), cmp=ID_cmp)
    win_pos_TB = {}
    for i in range(len(ids)):
        id = ids[i]
        win, loc = id.split('-')
        size, loc = len(win), int(loc)
        if size <= 0:
            continue
        if win not in win_pos_TB:
            win_pos_TB[win] = {}
        st, ed = loc, loc + size
        Trange = range(max(ed - 2*dyad_offset, 0), min(ref_length - 2*dyad_offset, st))
        for j in Trange:
            if j not in win_pos_TB[win]:
                win_pos_TB[win][j] = {}
                win_pos_TB[win][j]['T'] = []
                win_pos_TB[win][j]['B'] = []
            top_counts = top[i][j]
            bott_counts = bott[i][j + 2*dyad_offset]
            assert not np.isnan(top_counts) and not np.isnan(bott_counts)
            win_pos_TB[win][j]['T'].append(top_counts)
            win_pos_TB[win][j]['B'].append(bott_counts)

    # get linear model of T/B pairs on each position
    win_pos_model = {}
    for win in win_pos_TB:
        for pos in win_pos_TB[win]:
            X, Y = win_pos_TB[win][pos]['T'], win_pos_TB[win][pos]['B']
            if len(X) < 3: # too small data set
                continue
            if zero_frac(X) > 0.7: # X domain is empty
                continue
            if zero_frac(Y) > 0.7: # Y domain is empty
                continue

            try:
                coefs1, rsquare1 = robust_fit(X, Y)  # Top to bottom
            except:
                coefs1, rsquare1 = None, -sys.maxint # fitting error

            try:
                coefs2, rsquare2 = robust_fit(Y, X)  # Bottom to Top
            except:
                coefs2, rsquare2 = None, -sys.maxint # fitting error

            """
            if RSS1 > RSS2:
                if coefs2[0] != 0:
                    coefs1 = [1.0/coefs2[0], -1.0 * coefs2[1]/coefs2[0]]
                else:
                    coefs1 = [0, np.median(Y)]
            else:
                if coefs1[0] != 0:
                    coefs2 = [1.0/coefs1[0], -1.0 * coefs1[1]/coefs1[0]]
                else:
                    coefs2 = [0, np.median(X)]
            """
            
            if rsquare1 >= 0 and rsquare2 >= 0:
                if rsquare1 < rsquare2:
                    coefs1 = [1.0/coefs2[0], -1.0 * coefs2[1]/coefs2[0]]
                if rsquare2 < rsquare1:
                    coefs2 = [1.0/coefs1[0], -1.0 * coefs1[1]/coefs1[0]]
            else:
                if rsquare1 < 0:
                    coefs1 = [0, np.median(Y)]
                if rsquare2 < 0:
                    coefs2 = [0, np.median(X)]
            
            # show fitting graph
            if graph:
                fig = plt.figure()
                plt.plot(X, Y, '.k')
                pY = [np.poly1d(coefs1)(x) for x in X]
                pX = [np.poly1d(coefs2)(y) for y in Y]
                plt.plot(X, pY, 'b', label='Top to Bott fit', alpha=0.5)
                plt.plot(pX, Y, 'r', label='Bott to Top fit', alpha=0.5)
                plt.legend()
                plt.title("Top(" +str(pos) + ")-Bott(" + str(pos+2*dyad_offset) +  ") Pair")
                plt.xlabel("Top counts")
                plt.ylabel("Bott counts")
                plt.gca().set_aspect('equal')
                plt.show()
                plt.close()

            if win not in win_pos_model:
                win_pos_model[win] = {}
            if pos not in win_pos_model[win]:
                win_pos_model[win][pos] = {}
            win_pos_model[win][pos]['T'] = np.poly1d(coefs1)
            pos = pos + 2*dyad_offset
            if pos not in win_pos_model[win]:
                win_pos_model[win][pos] = {}
            win_pos_model[win][pos]['B'] = np.poly1d(coefs2)

    # fill the missing data based on linear model
    mtop, mbott = [], []
    for i in range(len(ids)):
        id = ids[i]
        win, loc = id.split('-')
        size, loc = len(win), int(loc)
        st, ed = loc, loc+size
        mtop.append([np.nan]*st + [0.0]*(ref_length-st))
        mbott.append([0.0]*ed + [np.nan]*(ref_length-ed))
    
    for i in range(len(ids)):
        id = ids[i]
        win, loc = id.split('-')
        size, loc = len(win), int(loc)           
        if size <= 0:
            continue
        st, ed = loc, loc + size
        Brange = range(min(st + 2*dyad_offset, ref_length), ref_length)
        Trange = range(0, max(ed -1 - 2*dyad_offset, -1) + 1)
        for j in Brange:
            try:
                model = win_pos_model[win][j]['B']
                x = bott[i][j]
                y = model(x)
            except:
                continue
            if np.isnan(y) or y < 0:
                continue
            mtop[i][j - 2*dyad_offset] += y
        for j in Trange:
            try:
                model = win_pos_model[win][j]['T']
                x = top[i][j]
                y = model(x)
            except:
                continue
            if np.isnan(y) or y < 0:
                continue
            mbott[i][j + 2*dyad_offset] += y

    dyad = []
    for i in range(len(ids)):
        full_top = list(np.nan_to_num(top[i]) + np.nan_to_num(mtop[i]))
        full_bott = list(np.nan_to_num(bott[i]) + np.nan_to_num(mbott[i]))
        tdyad = np.asarray([0.0]*dyad_offset + full_top[:ref_length-dyad_offset])
        bdyad = np.asarray(full_bott[dyad_offset:] + [0.0]*dyad_offset)
        dyad.append(list((tdyad + bdyad)*0.5))

    print >> sys.stderr, "Linear regression imputation is done"
    return [dyad], [mtop], [mbott]

def Bayesian (N,
              id_seq,
              top,
              bott,
              ref_length,
              dyad_offset,
              left_bound,
              right_bound,
              dyad_alphas=None,
              dyad_betas=None,
              r_alphas=None,
              r_betas=None,
              l_sigmas=None,
              initial="linear",
              note=""):
    
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

    def r_conditional(i, dyad, top, mtop, noise, dyad_offset, alpha, beta):
        assert i - dyad_offset >=0
        assert i + dyad_offset < ref_length
        def prior (x):
            assert x > 0
            output = 0.0
            j = 0
            while j <= beta -1:
                temp = ((-1)**j)*comb(beta-1,j)*(x**(2*alpha+j-1))
                if x <=1:
                    temp *= betaf(((2*alpha+j+1)/float(alpha))+1, beta)
                else:
                    temp *= betainc(((2*alpha+j+1)/float(alpha))+1, beta, (1.0/x)**alpha)*betaf(((2*alpha+j+1)/float(alpha))+1, beta)
                output += temp
                j +=1
            output *= alpha*beta*beta
            return output
        def accr_func (x, y):
            if x <= 0:
                return 0
            if noise[i] < 0:
                for j in range(len(dyad)):
                    if x*dyad[j][i] + noise[i] < 0:
                        return 0
            output = 0.0
            A = 0.0
            for j in range(len(dyad)):
                A += dyad[j][i]
                if not np.isnan(top[j][i-dyad_offset]):
                    cleavage = top[j][i-dyad_offset]
                else:
                    cleavage = mtop[j][i-dyad_offset]
                output += cleavage*(np.log(dyad[j][i]*x + noise[i]) - np.log(dyad[j][i]*y + noise[i]))
            #output += -(A+beta)*(x-y)
            #output += (alpha-1)*(np.log(x) - np.log(y))
            output += -A*(x-y)
            output += np.log(prior(x)) - np.log(prior(y))
            return np.exp(output)
        newr = MetroHasting(1000, 1, 10, accr_func)
        return newr[-1]

    def l_conditional(i, dyad, top, mtop, ratio, dyad_offset, sigma):
        assert i - dyad_offset >= 0
        assert i + dyad_offset < ref_length
        def accr_func (x, y):
            for j in range(len(dyad)):
                if ratio[i]*dyad[j][i] + x < 0:
                    return 0
            output = 0.0
            A = 0.0
            for j in range(len(dyad)):
                A += 1
                if not np.isnan(top[j][i-dyad_offset]):
                    cleavage = top[j][i-dyad_offset]
                else:
                    cleavage = mtop[j][i-dyad_offset]
                output += cleavage*(np.log(x + dyad[j][i]*ratio[i]) - np.log(y + dyad[j][i]*ratio[i]))
            output += -(A)*(x - y)
            output += -0.5*(x**2 - y**2)/float(sigma)
            return np.exp(output)
        newl = MetroHasting(1000, 0, 10, accr_func)
        return newl[-1]

    def D_conditional(j, i, top, mtop, bott, mbott, ratio, noise, dyad_offset, alpha, beta):
        assert i - dyad_offset >= 0
        assert i + dyad_offset < ref_length
        def accr_func (x, y):
            if x <= 0:
                return 0
            if ratio[i]*x + noise[i] < 0:
                return 0
            output = 0.0
            output += -(ratio[i] + beta + 1)*(x-y)
            if not np.isnan(top[j][i-dyad_offset]):
                tcleavage = top[j][i-dyad_offset]
            else:
                tcleavage = mtop[j][i-dyad_offset]
            if not np.isnan(bott[j][i+dyad_offset]):
                bcleavage = bott[j][i+dyad_offset]
            else:
                bcleavage = mbott[j][i+dyad_offset]
            output += (bcleavage + alpha - 1)*(np.log(x) - np.log(y))
            output += tcleavage*(np.log(ratio[i]*x + noise[i]) - np.log(ratio[i]*y + noise[i]))
            return np.exp(output)
        if beta > 0:
            initial_pt = float(alpha)/beta
        else:
            initial_pt = 10
        newD = MetroHasting(1000, initial_pt, 100, accr_func)
        return newD[-1]

    def T_conditional(j, i, dyad, ratio, noise, dyad_offset):
        assert i + 2*dyad_offset < len(dyad[j])
        lamb = dyad[j][i+dyad_offset]*ratio[i+dyad_offset] + noise[i+dyad_offset]
        assert lamb >= 0
        #print noise[i+dyad_offset]
        #print lamb
        #print
        return np.random.poisson(lamb,1)[0]

    def B_conditional(j, i, dyad, dyad_offset):
        assert i - 2*dyad_offset >= 0
        lamb = dyad[j][i-dyad_offset]
        assert lamb >= 0
        return np.random.poisson(lamb,1)[0]

    def Gibbs (N, dyad, ratio, noise, top, mtop, bott, mbott, dyad_offset):
        def save_data (f, data, cycle_num):
            print >> f, '>%s' % (cycle_num)
            if type(data[0]) == list:
                for j in range(len(data)):
                    s = ""
                    for i in range(len(data[j])):
                        s += str(data[j][i]) + ','
                    #print s
                    print >> f, s[:-1]
            else:
                s = ""
                for i in range(len(data)):
                    s += str(data[i]) + ','
                print >> f, s[:-1]
            return None

        t = open("top" + note + ".txt", 'w')
        b = open("bott" + note + ".txt", 'w')
        save_data(t, top, 0)
        save_data(b, bott, 0)
        t.close()
        b.close()

        d = open("dyad" + note + ".txt", 'w')        
        r = open("ratio" + note + ".txt", 'w')
        n = open("noise" + note + ".txt", 'w')
        mt = open("mtop" + note + ".txt", 'w')
        mb = open("mbott" + note + ".txt", 'w')
        save_data(d, dyad, 0)
        save_data(r, ratio, 0)
        save_data(n, noise, 0)
        save_data(mt, mtop, 0)
        save_data(mb, mbott, 0)
                
        st, ed = left_bound, ref_length - right_bound
        print >> sys.stderr, "Gibbs sampling start"
        for k in range(N):            
            for i in range(st, ed):
                ratio[i] = r_conditional(i, dyad, top, mtop, noise, dyad_offset, r_alphas[i], r_betas[i])
                noise[i] = l_conditional(i, dyad, top, mtop, ratio, dyad_offset, l_sigmas[i])
            for j in range(len(dyad)):
                for i in range(st, ed):
                    dyad[j][i] = D_conditional(j, i, top, mtop, bott, mbott, ratio, noise, dyad_offset, dyad_alphas[j][i], dyad_betas[j][i])                
            for j in range(len(mtop)):
                for i in range(st-dyad_offset, ed-dyad_offset):
                    if not np.isnan(mtop[j][i]):
                        mtop[j][i] = T_conditional(j, i, dyad, ratio, noise, dyad_offset)
                for i in range(st+dyad_offset, ed+dyad_offset):
                    if not np.isnan(mbott[j][i]):
                        mbott[j][i] = B_conditional(j, i, dyad, dyad_offset)
    
            save_data(d, dyad, k+1)
            save_data(r, ratio, k+1)
            save_data(n, noise, k+1)
            save_data(mt, mtop, k+1)
            save_data(mb, mbott, k+1)
            print >> sys.stderr, "cycle "  + str(k+1)

        d.close()
        r.close()
        n.close()
        mt.close()
        mb.close()
        return None

    print >> sys.stderr, "\n" + "Bayesian imputation start"

    # give initial guess
    if initial == "linear":
        dyad, mtop, mbott = Linear_reg (id_seq, top, bott, ref_length, dyad_offset)
    else:
        dyad, mtop, mbott = Naive (id_seq, top, bott, ref_length, dyad_offset)
    dyad = dyad[0]
    mtop = mtop[0]
    mbott = mbott[0]
    ratio = [1.0]*ref_length
    noise = [0.0]*ref_length
    
    """
    mtop, mbott = [], []
    dyad = []
    ids = sorted(id_seq.keys(), cmp=ID_cmp)
    for i in range(len(ids)):
        id = ids[i]
        win, loc = id.split('-')
        st, ed = int(loc), int(loc) + len(win)
        temp1 = [np.nan]*st + bott[i][st+2*dyad_offset:]
        temp2 = top[i][:max(0, ed-2*dyad_offset)] + [np.nan]*(ref_length-ed)
        mtop.append(temp1 + [0.0]*(ref_length - len(temp1)))
        mbott.append([0.0]*(ref_length - len(temp2)) + temp2)
        full_top = list(np.nan_to_num(top[i]) + np.nan_to_num(mtop[i]))
        full_bott = list(np.nan_to_num(bott[i]) + np.nan_to_num(mbott[i]))
        tdyad = np.asarray([0.0]*dyad_offset + full_top[:ref_length-dyad_offset])
        bdyad = np.asarray(full_bott[dyad_offset:] + [0.0]*dyad_offset)
        dyad.append(list((tdyad + bdyad)*0.5))
    """

    # give prior paramters
    if dyad_alphas == None or dyad_betas == None:
        temp_alphas, temp_betas = [], []
    for j in range(len(dyad)):
        var = np.std(dyad[j])**2
        temp1, temp2 = [], []
        for i in range(len(dyad[j])):
            mean = dyad[j][i]
            temp1.append(float(mean**2)/var)
            temp2.append(float(mean)/var)
        temp_alphas.append(temp1)
        temp_betas.append(temp2)
    if dyad_alphas == None:
        dyad_alphas = temp_alphas
    if dyad_betas == None:
        dyad_betas = temp_betas
    if r_alphas == None:
        r_alphas = [2.0]*ref_length
    if r_betas == None:
        r_betas = [3.0]*ref_length
    if l_sigmas == None:
        l_sigmas = [0.1]*ref_length
        
    # imputation by Gibbs sampling
    Gibbs (N, dyad, ratio, noise, top, mtop, bott, mbott, dyad_offset)
    dyad_list = load_data("dyad" + note + ".txt")
    ratio_list = load_data("ratio" + note + ".txt")
    noise_list = load_data("noise" + note + ".txt")
    mtop_list = load_data("mtop" + note + ".txt")
    mbott_list = load_data("mbott" + note + ".txt")
    
    print >> sys.stderr, "Bayesian imputation is done"
    return dyad_list, ratio_list, noise_list, mtop_list, mbott_list


"""
# read reference
id_seq = {}
for line in open("polyAscanlib.ref"):
    if line.startswith('>AAAA-'):
        break
    if line.startswith('>'):
        id = line[1:].strip()
        continue
    id_seq[id] = line.strip()
ids = sorted(id_seq.keys(), cmp=ID_cmp)

# basic parameters
ref_length = 225
dyad_offset = 52
left_bound, right_bound = 147/2, 147/2

# true variables
dyad = []
for i in range(len(id_seq)):
    dyad.append([0.0]*left_bound + rectangular(ref_length-left_bound-right_bound, random.randint(800, 1000)) + [0.0]*right_bound)
    #dyad.append([0.0]*left_bound + gaussian(ref_length-left_bound-right_bound, random.randint(800, 1000)) + [0.0]*right_bound)
    #dyad.append([0.0]*left_bound + random_sig(ref_length-left_bound-right_bound, random.randint(800, 1000)) + [0.0]*right_bound)
    
tlist, blist = [], []
eplist, delist = [], []
for i in range(ref_length):
    #tlist.append(random.random())
    #blist.append(random.random())
    tlist.append(1.0)
    blist.append(0.5)
    eplist.append(random.uniform(0,100))
    delist.append(random.uniform(0,100))
    
# simulate the data
top, bott, mtop, mbott = simulation (id_seq, dyad, tlist, blist, eplist, delist, dyad_offset)

# save true values
note = "_rect_test"

write_file([dyad], "dyad_true" + note + ".txt")
write_file([tlist], "tlist_true" + note + ".txt")
write_file([blist], "blist_true" + note + ".txt")
write_file([eplist], "eplist_true" + note + ".txt")
write_file([delist], "delist_true" + note + ".txt")
write_file([top], "top_true" + note + ".txt")
write_file([bott], "bott_true" + note + ".txt")
write_file([mtop], "mtop_true" + note + ".txt")
write_file([mbott], "mbott_true" + note + ".txt")

# Naive imputation
dyad_naive, mtop_naive, mbott_naive = Naive (id_seq, top, bott, dyad_offset)
write_file(dyad_naive, "dyad_naive" + note + ".txt")
write_file(mtop_naive, "mtop_naive" + note + ".txt")
write_file(mbott_naive, "mbott_naive" + note + ".txt")

# Linear regression imputation
dyad_linear, mtop_linear, mbott_linear = Linear_reg (id_seq, top, bott, dyad_offset, graph=False)
write_file(dyad_linear, "dyad_linear" + note + ".txt")
write_file(mtop_linear, "mtop_linear" + note + ".txt")
write_file(mbott_linear, "mbott_linear" + note + ".txt")

# Bayesain imputation
dyad_list, ratio_list, noise_list, mtop_list, mbott_list = Bayesian (1, id_seq, top, bott, dyad_offset, left_bound, right_bound, note='_bayesian' + note)
"""
