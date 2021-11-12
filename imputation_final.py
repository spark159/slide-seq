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

def Jenks_breaks (data_list, number_class):
    data_list.sort()
    mat1 = []
    for i in range(len(data_list) + 1):
        temp = []
        for j in range(number_class + 1):
            temp.append(0)
        mat1.append(temp)
    mat2 = []
    for i in range(len(data_list) + 1):
        temp = []
        for j in range(number_class + 1):
            temp.append(0)
        mat2.append(temp)
    for i in range(1, number_class + 1):
        mat1[1][i] = 1
        mat2[1][i] = 0
        for j in range(2, len(data_list) + 1):
            mat2[j][i] = float('inf')
    v = 0.0
    for l in range(2, len(data_list) + 1):
        s1 = 0.0
        s2 = 0.0
        w = 0.0
        for m in range(1, l + 1):
            i3 = l - m + 1
            val = float(data_list[i3 - 1])
            s2 += val * val
            s1 += val
            w += 1
            v = s2 - (s1 * s1) / w
            i4 = i3 - 1
            if i4 != 0:
                for j in range(2, number_class + 1):
                    if mat2[l][j] >= (v + mat2[i4][j - 1]):
                        mat1[l][j] = i3
                        mat2[l][j] = v + mat2[i4][j - 1]
        mat1[l][1] = 1
        mat2[l][1] = v
    k = len(data_list)
    kclass = []
    for i in range(number_class + 1):
        kclass.append(min(data_list))
    kclass[number_class] = float(data_list[len(data_list) - 1])
    count_num = number_class
    while count_num >= 2:  # print "rank = " + str(mat1[k][count_num])
        idx = int((mat1[k][count_num]) - 2)
        # print "val = " + str(data_list[idx])
        kclass[count_num - 1] = data_list[idx]
        k = int((mat1[k][count_num] - 1))
        count_num -= 1
    return kclass

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

def write_file (data_list, fname):
    f = open(fname, 'w')
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

# make simulated cleavage data
def simulation (wins,
                dyads,
                tlen,
                tlist,
                blist,
                eplist,
                delist,
                dyad_offset):
                
    def get_lambda(dyadmap, tlist, blist, eplist, delist, offset):
        tlamb, blamb = [0.0]*tlen, [0.0]*tlen
        for i in range(tlen):
            if i + offset < tlen:
                tlamb[i] += dyadmap[i+offset]*tlist[i]
            if i -  offset >= 0:
                blamb[i] += dyadmap[i-offset]*blist[i]
            tlamb[i] += eplist[i]
            blamb[i] += delist[i]
        return tlamb, blamb

    def make_cleavage(tlamb, blamb):
        tnum_list, bnum_list = [], []
        for i in range(tlen):
            tnum = np.random.poisson(tlamb[i],1)[0]
            bnum = np.random.poisson(blamb[i],1)[0]
            tnum_list.append(tnum)
            bnum_list.append(bnum)
        return tnum_list, bnum_list

    tops, botts = [], []
    mtops, mbotts = [], []
    for i in range(len(wins)):
        dyadmap = dyads[i]
        tlamb, blamb = get_lambda(dyadmap, tlist, blist, eplist, delist, dyad_offset)
        tnum_list, bnum_list = make_cleavage(tlamb, blamb)
        wst, wed = wins[i]
        tops.append(tnum_list[:wst] + [np.nan]*(tlen-wst))
        mtops.append([np.nan]*wst + tnum_list[wst:])
        botts.append([np.nan]*(wed+1) +  bnum_list[wed+1:])
        mbotts.append(bnum_list[:wed+1] + [np.nan]*(tlen-wed-1))
        
    print >> sys.stderr, "simulation is done"
    return tops, botts, mtops, mbotts

# simple copy and paste imputation
def Naive (wins,
           tops,
           botts,
           dyad_offset,
           silent=False):

    mtops, mbotts = [], []
    dyads = []
    for i in range(len(wins)):
        wst, wed = wins[i]
        top, bott = tops[i], botts[i]
        tlen = len(top)
        mtop, mbott = [], []
        for j in range(tlen):
            if j <= wed:
                if j-2*dyad_offset >=0:
                    mbott.append(top[j-2*dyad_offset])
                else:
                    mbott.append(0)
            else:
                mbott.append(np.nan)
            if j >= wst:
                if j+2*dyad_offset < tlen:
                    mtop.append(bott[j+2*dyad_offset])
                else:
                    mtop.append(0)
            else:
                mtop.append(np.nan)
        mtops.append(mtop)
        mbotts.append(mbott)

        dyad = [0.0]*tlen
        for j in range(dyad_offset, tlen-dyad_offset):
            if j + dyad_offset <= wed:
                dyad[j] += mbott[j+dyad_offset]
            else:
                dyad[j] += bott[j+dyad_offset]
            if j - dyad_offset >= wst:
                dyad[j] += mtop[j-dyad_offset]
            else:
                dyad[j] += top[j-dyad_offset]
        dyads.append(dyad)

    if not silent:
        print >> sys.stderr, "Naive imputation is done"
    return dyads, mtops, mbotts

def get_corr(X, Y):
    assert len(X) == len(Y)
    X_mean = np.mean(X)
    Y_mean = np.mean(Y)
    XYterm = 0.0
    Xterm, Yterm = 0.0, 0.0
    for i in range(len(X)):
        dX = X[i] - X_mean
        dY = Y[i] - Y_mean
        XYterm += dX*dY
        Xterm += dX*dX
        Yterm += dY*dY
    return XYterm / np.sqrt(Xterm * Yterm)
 

# imputation based on linear regression model 
def Linear_reg (wins,
                tops,
                botts,
                tlen,
                dyad_offset,
                graph=False,
                silent=False):

    def get_corr(X, Y):
        assert len(X) == len(Y)
        X_mean = np.mean(X)
        Y_mean = np.mean(Y)
        XYterm = 0.0
        Xterm, Yterm = 0.0, 0.0
        for i in range(len(X)):
            dX = X[i] - X_mean
            dY = Y[i] - Y_mean
            XYterm += dX*dY
            Xterm += dX*dX
            Yterm += dY*dY
        return XYterm / np.sqrt(Xterm * Yterm)

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
        rsquare = model.score(mX, mY)
        return coefs, rsquare

    def zero_frac (L):
        num = 0.0
        for e in L:
            if e == 0:
                num +=1
        return num/len(L)

    # collect observable T/B counts pairs along the template
    pos_TB = {}
    for i in range(len(wins)):
        wst, wed = wins[i]
        Trange = range(max(wed+1 - 2*dyad_offset, 0), min(tlen - 2*dyad_offset, wst))
        for j in Trange:
            if j not in pos_TB:
                pos_TB[j] = {'T':[], 'B':[]}
            top_counts = tops[i][j]
            bott_counts = botts[i][j + 2*dyad_offset]
            assert not np.isnan(top_counts) and not np.isnan(bott_counts)
            pos_TB[j]['T'].append(top_counts)
            pos_TB[j]['B'].append(bott_counts)

    # get linear model of T/B pairs on each position
    pos_model = {}
    for pos in pos_TB:
        X, Y = pos_TB[pos]['T'], pos_TB[pos]['B']

        # too small data set
        if len(X) < 3: 
            continue

        # X or Y domain is mostly empty
        if zero_frac(X) > 0.7 or zero_frac(Y) > 0.7: 
            continue
            #coefs1 = [0, np.median(Y)]
            #coefs2 = [0, np.median(X)]
                    
        # too low correlation for linear fitting 
        if abs(get_corr(X, Y)) < 0.2:
            coefs1 = [0, np.median(Y)]
            coefs2 = [0, np.median(X)]

        # try linear fitting
        else:
            try:
                coefs1, rsquare1 = robust_fit(X, Y)  # Top to bottom
            except:
                coefs1, rsquare1 = None, -sys.maxint # fitting error

            try:
                coefs2, rsquare2 = robust_fit(Y, X)  # Bottom to Top
            except:
                coefs2, rsquare2 = None, -sys.maxint # fitting error

            # one of fittings works
            if rsquare1 >= 0 and rsquare2 >= 0:
                if rsquare1 < rsquare2:
                    coefs1 = [1.0/coefs2[0], -1.0 * coefs2[1]/coefs2[0]]
                if rsquare2 < rsquare1:
                    coefs2 = [1.0/coefs1[0], -1.0 * coefs1[1]/coefs1[0]]
            # fitting is not good, then use median as estimate
            else:
                if rsquare1 < 0:
                    coefs1 = [0, np.median(Y)]
                if rsquare2 < 0:
                    coefs2 = [0, np.median(X)]

            # avoid too risky fitting
            if abs(coefs1[0]) > 100:
                coefs1 = [0, np.median(Y)]
            if abs(coefs2[0]) > 100:
                coefs2 = [0, np.median(X)]


        # show fitting graph
        if graph:
            fig = plt.figure()
            plt.plot(X, Y, '.k')
            pY = [np.poly1d(coefs1)(x) for x in X]
            pX = [np.poly1d(coefs2)(y) for y in Y]
            plt.plot(X, pY, 'b', label='Top to Bott fit', alpha=0.5)
            plt.plot(pX, Y, 'r', label='Bott to Top fit', alpha=0.5)
            plt.text(np.mean(X), np.mean(Y), "Corr:%.2f" %(corr))
            plt.legend()
            plt.title("Top(" +str(pos) + ")-Bott(" + str(pos+2*dyad_offset) +  ") Pair")
            plt.xlabel("Top counts")
            plt.ylabel("Bott counts")
            plt.gca().set_aspect('equal')
            plt.show()
            plt.close()

        if pos not in pos_model:
            pos_model[pos] = {}
        pos_model[pos]['T'] = np.poly1d(coefs1)
        pos = pos + 2*dyad_offset
        if pos not in pos_model:
            pos_model[pos] = {}
        pos_model[pos]['B'] = np.poly1d(coefs2)

    # fill the missing data based on linear model
    mtops, mbotts = [], []
    dyads = []
    for i in range(len(wins)):
        wst, wed = wins[i]
        top, bott = tops[i], botts[i]
        mtop, mbott = [], []
        for j in range(tlen):
            if j < wst:
                mtop.append(np.nan)
            else:
                try:
                    model = pos_model[j+2*dyad_offset]['B']
                    x = bott[j+2*dyad_offset]
                    y = model(x)
                    if np.isnan(y) or y < 0:
                        y = 0
                except:
                    y = 0
                mtop.append(y)
            if j > wed:
                mbott.append(np.nan)
            else:
                try:
                    model = pos_model[j-2*dyad_offset]['T']
                    x = top[j-2*dyad_offset]
                    y = model(x)
                    if np.isnan(y) or y < 0:
                        y = 0
                except:
                    y = 0
                mbott.append(y)
        mtops.append(mtop)
        mbotts.append(mbott)

        dyad = [0.0]*tlen
        for j in range(dyad_offset, tlen-dyad_offset):
            if j + dyad_offset <= wed:
                dyad[j] += mbott[j+dyad_offset]
            else:
                dyad[j] += bott[j+dyad_offset]
            if j - dyad_offset >= wst:
                dyad[j] += mtop[j-dyad_offset]
            else:
                dyad[j] += top[j-dyad_offset]
        dyads.append(dyad)

    if not silent:
        print >> sys.stderr, "Linear regression imputation is done"
    return dyads, mtops, mbotts

# imputation based on Bayesain inference
def Bayesian (N,
              wins,
              tops,
              botts,
              tlen,
              dyad_offset,
              left_bound,
              right_bound,
              dyad_alphas=None,
              dyad_betas=None,
              r_alpha=None,
              r_beta=None,
              l_mu=None,
              l_sigma=None,
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
            #print threshold
            if random.uniform(0,1) <= threshold:
                #success += 1
                samples.append(current)
                prev = copy.deepcopy(current)
            else:
                samples.append(prev)
        #print float(success)/N
        #fig = plt.figure()
        #plt.plot(acf(samples))
        #plt.show()
        #plt.close()
        return samples

    def r_conditional(i, dyads, tops, mtops, botts, mbotts, ratio, noise, dyad_offset, alpha, beta):
        assert i - dyad_offset >=0
        assert i + dyad_offset < tlen
        def accr_func (x, y):
            # impossible case
            if x <= 0 or x >= 1:
                return 0
            # impossible case
            for j in range(len(dyads)):
                if x*dyads[j][i] + noise[i] < 0:
                    return 0
                if (1-x)*dyads[j][i] - noise[i] < 0:
                    return 0
            output = 0.0
            output += (alpha-1)*(np.log(x) - np.log(y))
            output += (beta-1)*(np.log(1-x) - np.log(1-y))
            for j in range(len(dyads)):
                if not np.isnan(tops[j][i-dyad_offset]):
                    tcleavage = tops[j][i-dyad_offset]
                else:
                    tcleavage = mtops[j][i-dyad_offset]
                output += tcleavage*(np.log(dyads[j][i]*x + noise[i]) - np.log(dyads[j][i]*y + noise[i]))
                if not np.isnan(botts[j][i+dyad_offset]):
                    bcleavage = botts[j][i+dyad_offset]
                else:
                    bcleavage = mbotts[j][i+dyad_offset]
                output += bcleavage*(np.log(dyads[j][i]*(1-x) - noise[i]) - np.log(dyads[j][i]*(1-y) - noise[i]))                
            return np.exp(output)
        #newr = MetroHasting(1000, 0.5, 0.1, accr_func)
        #newr = MetroHasting(1000, ratio[i], 0.1, accr_func)
        #print alpha, beta
        std = np.sqrt(float(alpha*beta)/(((alpha+beta)**2)*(alpha+beta+1)))
        #print 'r', std
        newr = MetroHasting(100, ratio[i], std, accr_func)
        return newr[-1]

    def l_conditional(i, dyads, tops, mtops, botts, mbotts, ratio, noise, dyad_offset, mu, sigma):
        assert i - dyad_offset >= 0
        assert i + dyad_offset < tlen
        def accr_func (x, y):
            # impossible case
            for j in range(len(dyads)):
                if ratio[i]*dyads[j][i] + x < 0:
                    return 0
                if (1-ratio[i])*dyads[j][i] - x < 0:
                    return 0
            output = 0.0
            output += -0.5*((x-mu)/sigma)**2 + 0.5*((y-mu)/sigma)**2
            for j in range(len(dyads)):
                if not np.isnan(tops[j][i-dyad_offset]):
                    tcleavage = tops[j][i-dyad_offset]
                else:
                    tcleavage = mtops[j][i-dyad_offset]
                output += tcleavage*(np.log(dyads[j][i]*ratio[i] + x) - np.log(dyads[j][i]*ratio[i] + y))
                if not np.isnan(botts[j][i+dyad_offset]):
                    bcleavage = botts[j][i+dyad_offset]
                else:
                    bcleavage = mbotts[j][i+dyad_offset]
                output += bcleavage*(np.log(dyads[j][i]*(1-ratio[i]) - x) - np.log(dyads[j][i]*(1-ratio[i]) - y))                
            return np.exp(output)
        #newl = MetroHasting(1000, 0, 10, accr_func)
        #newl = MetroHasting(1000, noise[i], 10, accr_func)
        #print 'l', sigma
        newl = MetroHasting(100, noise[i], sigma, accr_func)
        return newl[-1]

    def D_conditional(j, i, dyads, tops, mtops, botts, mbotts, ratio, noise, dyad_offset, alpha, beta):
        assert i - dyad_offset >= 0
        assert i + dyad_offset < tlen
        def accr_func (x, y):
            if x <= 0:
                return 0
            if ratio[i]*x + noise[i] < 0:
                return 0
            if (1 - ratio[i])*x - noise[i] < 0:
                return 0
            output = 0.0
            output += -(1 + beta)*(x-y)
            output += (alpha - 1)*(np.log(x) - np.log(y))
            if not np.isnan(tops[j][i-dyad_offset]):
                tcleavage = tops[j][i-dyad_offset]
            else:
                tcleavage = mtops[j][i-dyad_offset]
            if not np.isnan(botts[j][i+dyad_offset]):
                bcleavage = botts[j][i+dyad_offset]
            else:
                bcleavage = mbotts[j][i+dyad_offset]
            output += tcleavage*(np.log(ratio[i]*x + noise[i]) - np.log(ratio[i]*y + noise[i]))
            output += bcleavage*(np.log((1-ratio[i])*x - noise[i]) - np.log((1-ratio[i])*y - noise[i]))
            return np.exp(output)
        #newD = MetroHasting(1000, float(alpha)/beta, 100, accr_func)
        #newD = MetroHasting(1000, dyads[j][i], 100, accr_func)
        std = np.sqrt(float(alpha)/(beta**2))
        #print 'D', std
        newD = MetroHasting(100, dyads[j][i], std, accr_func)
        return newD[-1]

    def T_conditional(j, i, dyads, ratio, noise, dyad_offset):
        assert i + 2*dyad_offset < len(dyads[j])
        lamb = dyads[j][i+dyad_offset]*ratio[i+dyad_offset] + noise[i+dyad_offset]
        assert lamb >= 0
        return np.random.poisson(lamb,1)[0]

    def B_conditional(j, i, dyads, ratio, noise, dyad_offset):
        assert i - 2*dyad_offset >= 0
        lamb = dyads[j][i-dyad_offset]*(1-ratio[i-dyad_offset]) - noise[i-dyad_offset]
        assert lamb >= 0
        return np.random.poisson(lamb,1)[0]

    def Gibbs (N, dyads, ratio, noise, tops, mtops, botts, mbotts, dyad_offset):
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

        t = open("tops" + note + ".txt", 'w')
        b = open("botts" + note + ".txt", 'w')
        save_data(t, tops, 0)
        save_data(b, botts, 0)
        t.close()
        b.close()

        d = open("dyads" + note + ".txt", 'w')        
        r = open("ratio" + note + ".txt", 'w')
        n = open("noise" + note + ".txt", 'w')
        mt = open("mtops" + note + ".txt", 'w')
        mb = open("mbotts" + note + ".txt", 'w')
        save_data(d, dyads, 0)
        save_data(r, ratio, 0)
        save_data(n, noise, 0)
        save_data(mt, mtops, 0)
        save_data(mb, mbotts, 0)
                
        st, ed = left_bound, tlen - right_bound
        print >> sys.stderr, "Gibbs sampling start"
        for k in range(N):            
            for i in range(st, ed):
                ratio[i] = r_conditional(i, dyads, tops, mtops, botts, mbotts, ratio, noise, dyad_offset, r_alpha[i], r_beta[i])
                noise[i] = l_conditional(i, dyads, tops, mtops, botts, mbotts, ratio, noise, dyad_offset, l_mu[i], l_sigma[i])
            for j in range(len(dyads)):
                for i in range(st, ed):
                    dyads[j][i] = D_conditional(j, i, dyads, tops, mtops, botts, mbotts, ratio, noise, dyad_offset, dyad_alphas[j][i], dyad_betas[j][i])
            for j in range(len(mtops)):
                for i in range(st-dyad_offset, ed-dyad_offset):
                    if not np.isnan(mtops[j][i]):
                        mtops[j][i] = T_conditional(j, i, dyads, ratio, noise, dyad_offset)
                for i in range(st+dyad_offset, ed+dyad_offset):
                    if not np.isnan(mbotts[j][i]):
                        mbotts[j][i] = B_conditional(j, i, dyads, ratio, noise, dyad_offset)
    
            save_data(d, dyads, k+1)
            save_data(r, ratio, k+1)
            save_data(n, noise, k+1)
            save_data(mt, mtops, k+1)
            save_data(mb, mbotts, k+1)
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
        dyads, mtops, mbotts = Linear_reg (wins, tops, botts, tlen, dyad_offset, silent=True)
    else:
        dyads, mtops, mbotts = Naive (wins, tops, botts, dyad_offset, silent=True)
    ratio = [0.5]*tlen
    noise = [0.0]*tlen
    
    # give prior paramters
    # effective nucleosome signal (gamma distribution)
    if dyad_alphas == None or dyad_betas == None:
        dyad_alphas, dyad_betas = [], []
        for j in range(len(dyads)):
            var = np.std(dyads[j])**2
            dyad_alpha, dyad_beta = [], []
            for i in range(len(dyads[j])):
                mean = dyads[j][i]
                dyad_alpha.append(float(mean**2)/var)
                dyad_beta.append(float(mean)/var)
            dyad_alphas.append(dyad_alpha)
            dyad_betas.append(dyad_beta)

    # effective top/bottom signal fraction (beta distribution)
    if r_alpha == None or r_beta == None:
        r_alpha = [3.0]*tlen
        r_beta = [3.0]*tlen
        
        st, ed = left_bound, tlen - right_bound
        for i in range(st, ed):
            fractions = []
            for j in range(len(dyads)):
                if not np.isnan(tops[j][i-dyad_offset]):
                    tcleavage = tops[j][i-dyad_offset]
                else:
                    tcleavage = mtops[j][i-dyad_offset]
                if not np.isnan(botts[j][i+dyad_offset]):
                    bcleavage = botts[j][i+dyad_offset]
                else:
                    bcleavage = mbotts[j][i+dyad_offset]
                if tcleavage == 0 or bcleavage == 0:
                    fraction = 0.5
                else:
                    fraction = float(tcleavage) / (tcleavage + bcleavage)
                fractions.append(fraction)
            mean = np.mean(fractions)
            var = np.std(fractions)**2
            alpha = mean*(mean*(1-mean)/float(var) - 1)
            beta = alpha*(1/float(mean) - 1)
            r_alpha[i] = alpha
            r_beta[i] = beta

        r_alpha = [3.0]*tlen
        r_beta = [3.0]*tlen

        

    # effective noise term (normal distribution)
    if l_mu == None or l_sigma == None:
        def noise_estimation (data_list):
            threshold = Jenks_breaks(data_list, 3)[1]
            noise_list = []
            for data in data_list:
                if data <= threshold:
                    noise_list.append(data)
            return np.mean(noise_list)
        
        st, ed = left_bound, tlen - right_bound
        t_list, b_list = [], []    
        for j in range(len(dyads)):
            for i in range(st, ed):
                if not np.isnan(tops[j][i-dyad_offset]):
                    tcleavage = tops[j][i-dyad_offset]
                else:
                    tcleavage = mtops[j][i-dyad_offset]
                if not np.isnan(botts[j][i+dyad_offset]):
                    bcleavage = botts[j][i+dyad_offset]
                else:
                    bcleavage = mbotts[j][i+dyad_offset]
                t_list.append(tcleavage)
                b_list.append(bcleavage)
                
        tnoise = noise_estimation(t_list)
        bnoise = noise_estimation(b_list)

        fig = plt.figure()
        plt.hist(t_list, bins=100)
        plt.axvline(x=tnoise, color='r')
        plt.show()
        plt.close()

        fig = plt.figure()
        plt.hist(b_list, bins=100)
        plt.axvline(x=bnoise, color='r')
        plt.show()
        plt.close()

        l_mu = [0.0]*tlen
        l_sigma = [0.25*(tnoise + bnoise)]*tlen

        #for i in range(st, ed):
        #    diffs = []
        #    for j in range(len(dyads)):
        #        if not np.isnan(tops[j][i-dyad_offset]):
        #            tcleavage = tops[j][i-dyad_offset]
        #        else:
        #            tcleavage = mtops[j][i-dyad_offset]
        #        if not np.isnan(botts[j][i+dyad_offset]):
        #            bcleavage = botts[j][i+dyad_offset]
        #        else:
        #            bcleavage = mbotts[j][i+dyad_offset]
        #        if tcleavage == 0 or bcleavage == 0:
        #            fraction = 0.5
        #        else:
        #            fraction = float(tcleavage) / (tcleavage + bcleavage)
        #        diff = (1-fraction)*tnoise - fraction*bnoise
        #        diffs.append(diff)
        #    l_mu[i] = np.mean(diffs)
        #    l_sigma[i] = np.std(diffs)
            
        
    # imputation by Gibbs sampling
    Gibbs (N, dyads, ratio, noise, tops, mtops, botts, mbotts, dyad_offset)
    dyads_list = load_data("dyads" + note + ".txt")
    ratio_list = load_data("ratio" + note + ".txt")
    noise_list = load_data("noise" + note + ".txt")
    mtops_list = load_data("mtops" + note + ".txt")
    mbotts_list = load_data("mbotts" + note + ".txt")
    
    print >> sys.stderr, "Bayesian imputation is done"
    return dyads_list, ratio_list, noise_list, mtops_list, mbotts_list

"""
## validate imputation algorithms
# simulated data
# basic parameters
tlen = 225
dyad_offset = 53
left_bound, right_bound = 147/2, 147/2
shape = 'gauss'
note = '_' + shape
libsize = 10

# make random windows
wins = []
for i in range(libsize):
    wst = random.randint(left_bound-20, tlen-right_bound)
    wsize = random.randint(1, 20)
    #wsize = 5
    wed = wst + wsize - 1
    wins.append((wst, wed))

# set true variables
dyads = []
for i in range(libsize):
    dyad = [0.0]*left_bound
    if shape == 'rect':
        dyad += rectangular(tlen-left_bound-right_bound, random.randint(800, 1000))
    elif shape == 'gauss':
        dyad += gaussian(tlen-left_bound-right_bound, random.randint(800, 1000))
    elif shape == 'random':
        dyad += random_sig(tlen-left_bound-right_bound, random.randint(800, 1000))
    dyad += [0.0]*right_bound
    dyads.append(dyad)
    
tlist, blist = [], []
eplist, delist = [], []
for i in range(tlen):
    tlist.append(random.uniform(0.1, 0.9))
    blist.append(random.uniform(0.1, 0.9))
    eplist.append(random.uniform(0,100))
    delist.append(random.uniform(0,100))
    
# simulate the data
tops, botts, mtops, mbotts = simulation (wins, dyads, tlen, tlist, blist, eplist, delist, dyad_offset)

# save true values
write_file([dyads], "dyads_true" + note + ".txt")
write_file([tlist], "tlist_true" + note + ".txt")
write_file([blist], "blist_true" + note + ".txt")
write_file([eplist], "eplist_true" + note + ".txt")
write_file([delist], "delist_true" + note + ".txt")
write_file([tops], "tops_true" + note + ".txt")
write_file([botts], "botts_true" + note + ".txt")
write_file([mtops], "mtops_true" + note + ".txt")
write_file([mbotts], "mbotts_true" + note + ".txt")

# Naive imputation
dyads_naive, mtops_naive, mbotts_naive = Naive (wins, tops, botts, dyad_offset)
write_file([dyads_naive], "dyads_naive" + note + ".txt")
write_file([mtops_naive], "mtops_naive" + note + ".txt")
write_file([mbotts_naive], "mbotts_naive" + note + ".txt")

# Linear regression imputation
dyads_linear, mtops_linear, mbotts_linear = Linear_reg (wins, tops, botts, tlen, dyad_offset, graph=True)
write_file([dyads_linear], "dyads_linear" + note + ".txt")
write_file([mtops_linear], "mtops_linear" + note + ".txt")
write_file([mbotts_linear], "mbotts_linear" + note + ".txt")


# Bayesain imputation
dyad_list, ratio_list, noise_list, mtop_list, mbott_list = Bayesian (1, wins, tops, botts, tlen, dyad_offset, left_bound, right_bound, note='_bayesian' + note)


# plot data
for i in range(libsize):
    fig = plt.figure()
    plt.plot(tops[i], 'r--', label='true top', alpha=0.5)
    plt.plot(mtops[i], 'r--', label='true mtop', alpha=0.5)
    plt.plot(botts[i], 'b--', label='true bott', alpha=0.5)
    plt.plot(mbotts[i], 'b--', label='true mbott', alpha=0.5)
    plot_file("mtops_bayesian" + note + ".txt", seq_id=i, label='bayesian mtop')
    plot_file("mbotts_bayesian" + note + ".txt", seq_id=i, label='bayesian mbott')
    #plt.plot(mtops_linear[i], 'r-', label='linear mtop')
    #plt.plot(mbotts_linear[i], 'b-', label='linear mbott')
    #plt.plot(mtops_naive[i], label='naive mtop')
    #plt.plot(mbotts_naive[i], label='naive mbott')
    wst, wed = wins[i]
    plt.axvspan(wst, wed-1, alpha=0.5, color='red')
    plt.legend()
    plt.show()
    plt.close()

    #fig = plt.figure()
    #plt.plot(dyads[i], linestyle='--', label='true dyad', alpha=0.5)
    #plt.plot(dyads_linear[i], label='linear dyad')
    #plt.plot(dyads_naive[i], label='naive dyad')
    #wst, wed = wins[i]
    #plt.axvspan(wst, wed-1, alpha=0.5, color='red')
    #plt.legend()
    #plt.show()
    #plt.close()
"""
