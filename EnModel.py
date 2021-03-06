import sys
import math
import copy
import numpy as np
from scipy.optimize import curve_fit
from sklearn import linear_model
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
import graph_edit
import random
from SliderClass import Slider

def norm(L):
    total = 0.0
    for value in L:
        total += value
    return [value/total for value in L]

def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq:
        nt = nt.upper()
        if nt == 'A':
            new_seq += 'T'
        elif nt == 'T':
            new_seq += 'A'
        elif nt == 'C':
            new_seq += 'G'
        elif nt == 'G':
            new_seq += 'C'
        else:
            new_seq += nt 
    return new_seq

def is_pal (seq):
    if len(seq) % 2 != 0:
        return False
    for i in range(len(seq)/2):
        if nt[i] != rev_comp(nt[len(seq)-1-i]):
            return False
    return True

def all_path(N, states, arrange=True):
    temp = []
    if N==1:
        temp = list(states)
    else:
        for path in all_path(N-1, states):
            for state in states:
                temp.append(path+state)
    if not arrange:
        return temp
    unique, non_pals = [], []
    for seq in temp:
        if rev_comp(seq) not in unique:
            unique.append(seq)
        else:
            non_pals.append(seq)
    output = unique + non_pals
    assert len(temp) == len(output)
    return output

def get_corr(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = np.average(x)
    avg_y = np.average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    return diffprod / np.sqrt(xdiff2 * ydiff2)

def Amer_len(seq, pos=True):
    num = []
    num_pos = {}
    i = 0
    while i < len(seq):
        if seq[i] in 'AT':
            nt = seq[i]
            count = 1
            j = i + 1
            while j < len(seq):
                if seq[j] != nt:
                    break
                count +=1
                j +=1
            num.append(count)
            if count not in num_pos:
                num_pos[count] = []
            num_pos[count].append(i)
            i = j
        else:
            i +=1
    if pos:
        return num_pos
    if len(num) == 0:
        return 0
    return max(num)

def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
        if nt not in 'ATCG':
            num+=0.5
    return (num/float(len(seq)))*100

#def mask (seq, idxs, type="skips", space=""):
#    new_seq = ""
#    idxs = sorted(idxs)
#    if type == "skips":
#        new_seq += seq[0:idxs[0]] + space
#        for i in range(len(idxs)-1):
#            st, ed = idxs[i]+1, idxs[i+1]
#            if st == ed:
#                continue
#            new_seq += seq[st:ed] + space
#        new_seq += seq[idxs[-1]+1:]
#    elif type == "includes":
#        for i in range(len(idxs)-1):
#            idx = idxs[i]
#            next_idx = idxs[i+1]
#            new_seq += seq[idx]
#            if next_idx - idx > 1:
#                new_seq += space
#        new_seq += seq[next_idx]
#    assert len(new_seq) == len(seq)
#    return new_seq

def mask (seq, idxs, replace='-'):
    nt_list = list(seq)
    for idx in idxs:
        nt_list[idx] = replace
    return "".join(nt_list)
            

class EnergyModel:
    def __init__(self, key_slider, NCPlen=147, templatelen=225, bound=0, mask_idxs=None, shape=False):
        assert NCPlen % 2 != 0
        #self.key_slider = key_slider
        self.NCPlen = NCPlen
        self.templatelen = templatelen
        self.bound = bound
        self.mask_idxs = mask_idxs
        self.nlogprob_list, self.seq_list, self.pos_list = [], [], []
        self.MGW_list, self.HelT_list, self.ProT_list, self.Roll_list = [], [], [], []
        self.key_list = sorted(key_slider.keys())
        for key in self.key_list:
            slider = key_slider[key]
            #print slider.dyadmap
            energy_profile = slider.energy_profile(scale=1000)
            #energy_profile = -np.log(norm(np.asarray(slider.dyadmap)+sys.float_info.min))
            seq = slider.seq
            for i in range(self.NCPlen/2+bound, self.templatelen-self.NCPlen/2-bound):
                nlogprob = energy_profile[i]
                NCPseq = seq[i-self.NCPlen/2:i+self.NCPlen/2+1]

                if self.mask_idxs:
                    NCPseq = mask(NCPseq, self.mask_idxs, replace='-')
                    
                self.nlogprob_list.append(nlogprob)
                self.seq_list.append(NCPseq)
                self.pos_list.append(i)
                if shape:
                    self.MGW_list.append(slider.MGW[i-self.NCPlen/2:i+self.NCPlen/2+1])
                    self.HelT_list.append(slider.HelT[i-self.NCPlen/2:i+self.NCPlen/2+1])
                    self.ProT_list.append(slider.ProT[i-self.NCPlen/2:i+self.NCPlen/2+1])
                    self.Roll_list.append(slider.Roll[i-self.NCPlen/2:i+self.NCPlen/2+1])
            
        self.reg = None
        self.coeff = None
        self.rsquare = None
        self.corr = None
        self.freq = None
        if shape:
            self.shape = True
        else:
            self.shape = False

    def _even_sampling(self, sym=True, shape=False):
        if not shape:
            seq_samples = []
            for i in range(len(self.seq_list)):
                seq = self.seq_list[i]
                seq_samples.append(seq)
                if sym:
                    seq_samples.append(rev_comp(seq))
            return seq_samples
        if shape:
            MGW_samples, HelT_samples, ProT_samples, Roll_samples = [], [], [], []
            for i in range(len(self.seq_list)):
                MGW_samples.append(self.MGW_list[i])
                HelT_samples.append(self.HelT_list[i])
                ProT_samples.append(self.ProT_list[i])
                Roll_samples.append(self.Roll_list[i])
                if sym:
                    MGW_samples.append(self.MGW_list[i][::-1])
                    HelT_samples.append(self.HelT_list[i][::-1])
                    ProT_samples.append(self.ProT_list[i][::-1])
                    Roll_samples.append(self.Roll_list[i][::-1])
            return [MGW_samples, HelT_samples, ProT_samples, Roll_samples]


    def _bias_sampling(self, scale=100, sym=True, shape=False):
        if not shape:
            seq_samples = []
            for i in range(len(self.seq_list)):
                seq = self.seq_list[i]
                nlogprob = self.nlogprob_list[i]
                prob = np.exp(-nlogprob)
                count = int(round(prob*scale))
                for i in range(count):
                    seq_samples.append(seq)
                    if sym:
                        seq_samples.append(rev_comp(seq))
            return seq_samples
        if shape:
            MGW_samples, HelT_samples, ProT_samples, Roll_samples = [], [], [], []
            for i in range(len(self.seq_list)):
                nlogprob = self.nlogprob_list[i]
                prob = np.exp(-nlogprob)
                count = int(round(prob*scale))
                for i in range(count):
                    MGW_samples.append(self.MGW_list[i])
                    HelT_samples.append(self.HelT_list[i])
                    ProT_samples.append(self.ProT_list[i])
                    Roll_samples.append(self.Roll_list[i])
                    if sym:
                        MGW_samples.append(self.MGW_list[i][::-1])
                        HelT_samples.append(self.HelT_list[i][::-1])
                        ProT_samples.append(self.ProT_list[i][::-1])
                        Roll_samples.append(self.Roll_list[i][::-1])

            return [MGW_samples, HelT_samples, ProT_samples, Roll_samples]

             
    def _stat_Markov(self, seq_list, order):
        ntdic = {}
        for nt in all_path(order+1, 'ATCG'):
            ntdic[nt] = 0.0

        sample_num = len(seq_list)
        
        freq = [ copy.deepcopy(ntdic) for i in range(self.NCPlen - order) ]
        for seq in seq_list:
            for i in range(len(seq) - order):
                nt = seq[i:i+1+order]
                try:
                    freq[i][nt] += 1.0 / sample_num
                except:
                    pass
        
        mean, std = [], []
        for ntdic in freq:
            mean.append(np.mean(ntdic.values()))
            std.append(np.std(ntdic.values()))

        stdz_freq = []
        for i in range(len(freq)):
            ntdic = freq[i]
            temp = {}
            for nt, value in ntdic.items():
                if std[i] > 0:
                    temp[nt] = (value - mean[i]) / std[i]
                else:
                    temp[nt] = np.nan
            stdz_freq.append(temp)

        return freq, sample_num, mean, std, stdz_freq
        
    def _stat_Kmer(self, seq_list, knum, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= knum
        
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2

        ntdic = {}
        for nt in all_path(knum, 'ATCG'):
            ntdic[nt] = 0.0
        freq = [ copy.deepcopy(ntdic) for i in range(bnum) ]

        sample_num = len(seq_list)*(seqlen-knum+1)
        
        for seq in seq_list:
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]                
                for i in range(seqlen - knum + 1):
                    nt = bseq[i:i+knum]
                    try:
                        freq[k][nt] += 1.0 / sample_num
                    except:
                        pass

        mean, std = [], []
        for ntdic in freq:
            mean.append(np.mean(ntdic.values()))
            std.append(np.std(ntdic.values()))

        stdz_freq = []
        for i in range(len(freq)):
            ntdic = freq[i]
            temp = {}
            for nt, value in ntdic.items():
                if std[i] > 0:
                    temp[nt] = (value - mean[i]) / std[i]
                else:
                    temp[nt] = np.nan
            stdz_freq.append(temp)

        return freq, sample_num, mean, std, stdz_freq

    def _stat_PolyA(self, seq_list, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2

        freq = [ [] for i in range(bnum) ]
        sample_num = len(seq_list)
        
        for seq in seq_list:
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                num_pos = Amer_len(bseq)
                score = 0.0
                for num, pos in num_pos.items():
                    if num >= 5:
                        score += num*len(pos)
                freq[k].append(score)
                
        mean, std = [], []
        for k in range(len(freq)):
            mean.append(np.mean(freq[k]))
            std.append(np.std(freq[k]))

        return freq, sample_num, mean, std

    def _stat_GC(self, seq_list, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2

        freq = [ [] for i in range(bnum) ]

        sample_num = len(seq_list)
        
        for seq in seq_list:
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                GC = GC_content(bseq)
                freq[k].append(GC)

        mean, std = [], []
        for k in range(len(freq)):
            mean.append(np.mean(freq[k]))
            std.append(np.std(freq[k]))
        return freq, sample_num, mean, std

    def kmer_pair_dist_prob(self, seq_list, knum, max_dist):
        assert self.NCPlen - knum >= max_dist
        pair_count = {}
        pair_dist_prob = {}
        for seq in seq_list:
            pair_dist_count = {}
            for i in range(len(seq)-knum):
                for j in range(i+1, min(len(seq)-knum+1, i+max_dist+1)):
                    kmer1 = seq[i:i+knum]
                    kmer2 = seq[j:j+knum]
                    dist = j - i
                    pair = tuple(sorted([kmer1, kmer2]))
                    if pair not in pair_dist_count:
                        pair_dist_count[pair] = [0.0]*max_dist
                    pair_dist_count[pair][dist-1] += 1
            for pair in pair_dist_count:
                if pair not in pair_dist_prob:
                    pair_dist_prob[pair] = np.asarray([0.0]*max_dist)
                pair_dist_prob[pair] += np.asarray(norm(pair_dist_count[pair]))
                #pair_dist_prob[pair] += np.asarray(pair_dist_count[pair])
                if pair not in pair_count:
                    pair_count[pair] =0
                pair_count[pair] +=1
        for pair in pair_dist_prob:
            pair_dist_prob[pair] = pair_dist_prob[pair]/float(pair_count[pair])
        return pair_dist_prob

    def _stat_shape(self, shape_list):
        shape_list = np.asarray(shape_list)
        #print shape_list.shape
        freq = np.sum(shape_list, axis=0) / len(shape_list)
        return freq

    def report (self,
                MM_orders,
                Kmer_k_b,
                PolyA_b,
                GC_b,
                Harmonic=False,
                PairCorr_k_dist=False,
                scale=100,
                shape=False,
                sym=True,
                graph=False):
        
        freq = {}
        if shape:
            print >> sys.stderr, "sampling data values"        
            even_samples_list = self._even_sampling(sym=sym, shape=True)
            bias_samples_list = self._bias_sampling(scale=scale, sym=sym, shape=True)
            names = ["MGW", "HelT", "ProT", "Roll"] 
            for i in range(len(even_samples_list)):
                even_samples, bias_samples = even_samples_list[i], bias_samples_list[i]
                freq[names[i]] = self._stat_shape(bias_samples) / self._stat_shape(even_samples)
                self.freq = freq
        else:
            # sampling data
            print >> sys.stderr, "sampling data values"        
            even_samples = self._even_sampling(sym=sym)
            bias_samples = self._bias_sampling(scale=scale, sym=sym)
        
        # frequency counts
        print >> sys.stderr, "counting sequence features"
        if MM_orders:
            for order in sorted(MM_orders):
                name = 'MM' + str(order)
                freq1, sample_num, mean, std, stdz_freq = self._stat_Markov(even_samples, order)
                freq2, sample_num, mean, std, stdz_freq = self._stat_Markov(bias_samples, order)
                freq_fold = [{} for i in range(self.NCPlen-order)]
                nts = all_path(order+1, 'ATCG')
                for i in range(self.NCPlen-order):
                    for nt in nts:
                        if freq1[i][nt] > 0:
                            freq_fold[i][nt] = freq2[i][nt] / freq1[i][nt]
                        else:
                            freq_fold[i][nt] = np.NaN
                freq[name] = freq_fold        
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            freq1, sample_num, mean, std, stdz_freq = self._stat_Kmer(even_samples, knum, bnum)
            freq2, sample_num, mean, std, stdz_freq = self._stat_Kmer(bias_samples, knum, bnum)
            nts = all_path(knum, 'ATCG')
            for i in range(bnum):
                name = 'Kmer' + str(i)
                freq_fold = {}
                for nt in nts:
                    if freq1[i][nt] > 0:
                        freq_fold[nt] = freq2[i][nt] / freq1[i][nt]
                    else:
                        freq_fold[nt] = np.NaN
                freq[name] = freq_fold
        if PolyA_b:
            bnum = PolyA_b
            freq1, sample_num, mean1, std = self._stat_PolyA(even_samples, bnum)
            freq2, sample_num, mean2, std = self._stat_PolyA(bias_samples, bnum)
            for i in range(bnum):
                name = 'PolyA' + str(i)
                if mean1[i] > 0:
                    freq[name] = mean2[i]/mean1[i]
                else:
                    freq[name] = np.NaN
        if GC_b:
            bnum = GC_b
            freq1, sample_num, mean1, std = self._stat_GC(even_samples, bnum)
            freq2, sample_num, mean2, std = self._stat_GC(bias_samples, bnum)
            for i in range(bnum):
                name = 'GC' + str(i)
                if mean1[i] > 0:
                    freq[name] = mean2[i]/mean1[i]
                else:
                    freq[name] = np.NaN

        if PairCorr_k_dist:
            knum, max_dist = PairCorr_k_dist
            pair_dist_prob1 = self._kmer_pair_dist_prob(even_samples, knum, max_dist)
            pair_dist_prob2 = self._kmer_pair_dist_prob(bias_samples, knum, max_dist)
            pair_corr = {}
            nts = all_path(knum, 'ATCG')
            for i in range(len(nts)-1):
                for j in range(i, len(nts)):
                    pair = tuple(sorted([nts[i], nts[j]]))
                    if pair not in pair_corr:
                        pair_corr[pair] = [np.nan]*max_dist
                    for k in range(max_dist):
                        try:
                            pair_corr[pair][k] = pair_dist_prob2[pair][k] / pair_dist_prob1[pair][k]
                        except:
                            pass
            freq['PairCorr'] = pair_corr
            
        if Harmonic:
            None

        # To do
        self.freq = freq
        print >> sys.stderr, "Done"
        return None
    

    def _var_Markov (self, seq_list, order):
        nt_pos = {}
        nts = all_path(order+1, 'ATCG')
        for i in range(len(nts)):
            nt = nts[i]
            nt_pos[nt] = i
        if order % 2 == 0:
            palnum = 0
        else:
            palnum = 4**((order+1)/2)
        var_list = []
        for seq in seq_list:
            left = seq
            right = rev_comp(seq)
            row = []
            for i in range((len(seq) - order)/2):
                nt1 = left[i:i+order+1]
                nt2 = right[i:i+order+1]
                temp = [0.0] * len(nts)
                try:
                    temp[nt_pos[nt1]] += 1
                except:
                    pass
                try:
                    temp[nt_pos[nt2]] += 1
                except:
                    pass  
                if sum(temp) > 0:
                    row += temp
                else: # data missing
                    row += [np.nan] * len(nts)
            if (len(seq) - order) % 2 != 0:
                assert order % 2 == 0
                i = len(seq)/2 - (order+1)/2
                nt = left[i:i+order+1]
                temp = [0.0] * ((len(nts)+palnum)/2)
                try:
                    pos = min(nt_pos[nt], nt_pos[rev_comp(nt)])
                    temp[pos] += 1
                except:
                    pass
                if sum(temp) > 0:
                    row += temp
                else: # data missing
                    row += [np.nan] * ((len(nts)+palnum)/2)
            var_list.append(row)
        return var_list

    def _var_Kmer (self, seq_list, knum, bnum):        
        seqlen = self.NCPlen / bnum
        assert seqlen >= knum        
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2
        nt_pos = {}
        nts = all_path(knum, 'ATCG')
        for i in range(len(nts)):
            nt = nts[i]
            nt_pos[nt] = i
        if knum % 2 != 0:
            palnum = 0
        else:
            palnum = 4**(knum/2)
        var_list = []
        for seq in seq_list:
            left = seq
            right = rev_comp(seq)
            row = []
            for k in range(bnum/2):
                temp = [0]*len(nts)
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq1 = left[st:st+seqlen]
                bseq2 = right[st:st+seqlen]
                for i in range(seqlen-knum+1):
                    nt1 = bseq1[i:i+knum]
                    nt2 = bseq2[i:i+knum]
                    try:
                        temp[nt_pos[nt1]] += 1
                    except:
                        pass
                    try:
                        temp[nt_pos[nt2]] += 1
                    except:
                        pass
                if sum(temp) > 0:
                    row += temp
                else: # data missing
                    row += [np.nan] * len(nts)
            if bnum % 2 != 0:
                assert seqlen % 2 != 0
                st = len(seq)/2 - seqlen/2
                bseq = left[st:st+seqlen]
                temp = [0] * ((len(nts)+palnum)/2)
                for i in range(seqlen-knum+1):
                    nt = bseq[i:i+knum]
                    try:
                        pos = min(nt_pos[nt], nt_pos[rev_comp(nt)])
                        temp[pos] += 1
                    except:
                        pass
                if sum(temp) > 0:
                    row += temp
                else: # data missing
                    row += [np.nan] * ((len(nts)+palnum)/2)
            var_list.append(row)
        return var_list

    def _var_PolyA (self, seq_list, bnum):        
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2        
        var_list = []
        for seq in seq_list:
            row = []
            for k in range(bnum):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                num_pos = Amer_len(bseq)
                count = 0
                score = 0.0
                for num, pos in num_pos.items():
                    if num >=5:
                        #count += len(pos)
                        score += (num)*len(pos) 
                #row.append(count)
                row.append(score)
            sym_row = [row[i] + row[::-1][i] for i in range(bnum/2)]
            if bnum % 2 != 0:
                sym_row += [row[bnum/2]]
            var_list.append(sym_row)
        return var_list

    def _var_GC (self, seq_list, bnum):
        seqlen = self.NCPlen / bnum
        assert seqlen >= 1
        extra = self.NCPlen % bnum
        boundoff = extra / 2
        centeroff = extra % 2        
        var_list = []
        for seq in seq_list:
            row = []
            for k in range((bnum+1)/2):
                if k < bnum/2:
                    st = boundoff + k*seqlen
                if k >= bnum/2:
                    st = boundoff + centeroff + k*seqlen
                bseq = seq[st:st+seqlen]
                GC = GC_content(bseq)
                row.append(GC)
            sym_row = [row[i] + row[::-1][i] for i in range(bnum/2)]
            if bnum % 2 != 0:
                sym_row += [row[bnum/2]]
            var_list.append(sym_row)
        return var_list

    def _var_Harmonic (self, pos_list):
        offset = self.templatelen / 2
        var_list = []
        for pos in pos_list:
            row = []
            row.append(0.5*((pos-offset)**2))
            var_list.append(row)
        return var_list

    def _var_dPolyA (self, seq_list, lmin=3, lmax=20):
        var_list = []
        for seq in seq_list:
            row = []
            num_pos = Amer_len(seq)
            for i in range(lmin, lmax+1):
                temp = [0.0]*len(seq)
                try:
                    pos_list = num_pos[i]
                except:
                    row += temp[:self.NCPlen/2 + 1]
                    continue
                for pos in pos_list:
                    for j in range(i):
                        temp[pos+j] += 1
                sym_temp = [temp[i] + temp[::-1][i] for i in range(self.NCPlen/2)]
                sym_temp += [temp[self.NCPlen/2]]
                row += sym_temp
            var_list.append(row)
        return var_list
    
    def train (self,
               MM_orders,
               Kmer_k_b,
               PolyA_b,
               GC_b,
               Harmonic,
               ref_key=None,
               adjust=False,
               dPolyA=False,
               alpha=0.5,
               k_fold=10,
               graph=False):

        self.MM_orders, self.Kmer_k_b, = MM_orders, Kmer_k_b
        self.PolyA_b, self.GC_b, self.Harmonic = PolyA_b, GC_b, Harmonic
        nlogprob_list, seq_list, pos_list = self.nlogprob_list, self.seq_list, self.pos_list
        
        # reading variables
        print >> sys.stderr, "reading data values"
        if MM_orders:
            MM_vars = [ [] for i in range(len(seq_list)) ] 
            for order in sorted(MM_orders):
                temp = self._var_Markov(seq_list, order)
                for i in range(len(temp)):
                    MM_vars[i] += temp[i]
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            Kmer_vars = self._var_Kmer(seq_list, knum, bnum)
        if PolyA_b:
            PolyA_vars = self._var_PolyA(seq_list, PolyA_b)
        if GC_b:
            GC_vars = self._var_GC(seq_list, GC_b)
        if Harmonic:
            Harmonic_vars = self._var_Harmonic(pos_list)
        if dPolyA:
            lmin, lmax = dPolyA
            dPolyA_vars = self._var_dPolyA(seq_list, lmin=lmin, lmax=lmax)

        var_list = [ [] for i in range(len(seq_list))]
        for i in range(len(seq_list)):
            if MM_orders:
                var_list[i] += MM_vars[i]
            if Kmer_k_b:
                var_list[i] += Kmer_vars[i]
            if PolyA_b:
                var_list[i] += PolyA_vars[i]
            if GC_b:
                var_list[i] += GC_vars[i]
            if Harmonic:
                var_list[i] += Harmonic_vars[i]
            if dPolyA:
                var_list[i] += dPolyA_vars[i]

        
        # adjust variables by reference key
        if ref_key:
            new_var_list = []
            new_nlogprob_list = []
            knum = self.key_list.index(ref_key)
            gnum = self.templatelen - self.NCPlen + 1 - 2*self.bound
            st, ed = knum*gnum, (knum+1)*gnum
            ref_var = np.asarray(var_list[st:ed])
            ref_nlogprob = np.asarray(nlogprob_list[st:ed])
            ref_prob = [np.exp(-value) for value in ref_nlogprob]
            #fig = plt.figure()
            #plt.plot(ref_prob)
            
            for i in range(len(self.key_list)):
                if i == knum:
                    continue
                var_group = np.asarray(var_list[i*gnum:(i+1)*gnum])
                nlogprob_group = np.asarray(nlogprob_list[i*gnum:(i+1)*gnum])
                new_var_list += list(var_group - ref_var)
                new_nlogprob_list += list((nlogprob_group - ref_nlogprob)/ref_nlogprob)
                #print self.key_list[i]
                #plt.title(self.key_list[i])
                #prob = [np.exp(-value) for value in nlogprob_group]
                #plt.plot(ref_prob)
                #plt.plot(prob)
                #plt.show()
            var_list = new_var_list
            nlogprob_list = new_nlogprob_list
            
        # adjust variables by reference point
        gnum = self.templatelen - self.NCPlen + 1 - 2*self.bound
        if adjust and gnum > 1:
            feature_list, target_list = [], []
            i = 0
            count = 0
            while i < len(var_list):
                ref_var = copy.deepcopy(var_list[i])
                ref_nlogprob = copy.deepcopy(nlogprob_list[i])
                j = 1
                count += 1
                while j < gnum:
                    row = np.asarray(var_list[i+j]) - np.asarray(ref_var)
                    feature_list.append(row)
                    target_list.append([nlogprob_list[i+j] - ref_nlogprob])
                    j +=1
                    count += 1
                i += j
            assert count == len(var_list)
        else:
            feature_list = var_list
            target_list = [[value] for value in nlogprob_list]

        # drop columns including NaNs (exclude data-missing coefficients)
        column_mask = ~np.any(np.isnan(feature_list), axis=0) 
        feature_list = list(np.asarray(feature_list)[:,column_mask])
        
        # k-fold corss validation
        print >> sys.stderr, "k-fold cross validation"
        feature_list, target_list = shuffle(feature_list, target_list)
        part_size = (len(feature_list) + 1) / k_fold
        corrs = []
        for i in range(k_fold):
            st = i*part_size
            ed = min((i+1)*part_size, len(feature_list))
            train_flist = feature_list[:st] + feature_list[ed:] 
            train_tlist = target_list[:st] + target_list[ed:]
            test_flist, test_tlist = feature_list[st:ed], target_list[st:ed]
            
            print str(i+1) + '-fold'
            reg = linear_model.Ridge(alpha=alpha)
            reg.fit (train_flist, train_tlist)
            #print reg.coef_
            rsquare = reg.score(train_flist, train_tlist)
            print 'r-square: ' + str(rsquare)

            Yexp = [ value[0] for value in test_tlist]
            Ypred = reg.predict(test_flist)
            Ypred = [ value[0] for value in Ypred]
            corr = get_corr(Yexp, Ypred)
            corrs.append(corr)
            print 'correlation: ' + str(corr)

            if graph:
                fig = plt.figure()
                low = min(Yexp + Ypred)
                up = max(Yexp + Ypred)
                mid = np.median(Yexp + Ypred)
                plt.plot(Yexp, Ypred, '.')
                plt.plot([low-mid,up+mid],[low-mid,up+mid], '--')
                plt.xlim([low - mid*0.1, up + mid*0.1])
                plt.ylim([low - mid*0.1, up + mid*0.1])
                plt.show()
                plt.close()

        self.corr = corrs
        print >> sys.stderr, "Mean correlation: " + str(np.mean(self.corr))
        
        # linear fitting all data
        print >> sys.stderr, "linear regression of full data"
        reg = linear_model.Ridge(alpha=alpha)
        reg.fit (feature_list, target_list)
        self.reg = reg
        print reg.coef_
        self.rsquare = reg.score(feature_list, target_list)
        print "r-square: " + str(self.rsquare)

        # putting back data-missing coefficients as NaN
        pt = 0
        full_coef = []
        for boolean in column_mask:
            if boolean:
                full_coef.append(reg.coef_[0][pt])
                pt +=1
            else:
                full_coef.append(np.nan)
        assert pt == len(reg.coef_[0])
        assert len(column_mask) == len(full_coef)

        # retrieve all coefficient 
        pt = 0
        self.coeff = {}
        if MM_orders:
            for order in MM_orders:
                name = 'MM' + str(order)
                self.coeff[name] = [{} for i in range(self.NCPlen-order)]
                nts = all_path(order+1, 'ATCG')
                if order % 2 == 0:
                    palnum = 0
                else:
                    palnum = 4**((order+1)/2)
                for k in range((self.NCPlen - order + 1)/2):
                    count = 0
                    if order % 2 == 0 and k == (self.NCPlen - order + 1)/2 - 1:
                        total = (len(nts)+palnum)/2
                    else:
                        total = len(nts)
                    while count < total:
                        nt = nts[count]
                        self.coeff[name][k][nt] = full_coef[pt]
                        self.coeff[name][self.NCPlen-order-1-k][rev_comp(nt)] = full_coef[pt]
                        count += 1
                        pt += 1
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            if knum % 2 != 0:
                palnum = 0
            else:
                palnum = 4**(knum/2)
            nts = all_path(knum, 'ATCG')
            for i in range((bnum+1)/2):
                name1 = 'Kmer' + str(i)
                name2 = 'Kmer' + str(bnum-1-i)
                self.coeff[name1] = {}
                self.coeff[name2] = {}
                count = 0
                if bnum % 2 != 0 and i == (bnum+1)/2 - 1:
                    total = (len(nts)+palnum)/2
                else:
                    total = len(nts)
                while count < total:
                    nt = nts[count]
                    self.coeff[name1][nt] = full_coef[pt]
                    self.coeff[name2][rev_comp(nt)] = full_coef[pt]
                    count += 1
                    pt += 1
        if PolyA_b:
            bnum = PolyA_b
            for i in range((bnum+1)/2):
                name1 = 'PolyA' + str(i)
                name2 = 'PolyA' + str(bnum-1-i)
                self.coeff[name1] = full_coef[pt]
                self.coeff[name2] = full_coef[pt]
                pt += 1
        if GC_b:
            bnum = GC_b
            for i in range((bnum+1)/2):
                name1 = 'GC' + str(i)
                name2 = 'GC' + str(bnum-1-i)
                self.coeff[name1] = full_coef[pt]
                self.coeff[name2] = full_coef[pt]
                pt += 1
        if Harmonic:
            name = 'Harmonic'
            self.coeff[name] = full_coef[pt]
            pt += 1
        if dPolyA:
            name = 'dPolyA'
            self.coeff[name] = {}
            for i in range(lmin, lmax+1):
                self.coeff[name][i] = {}
                for j in range(self.NCPlen/2 + 1):
                    self.coeff[name][i][j] = full_coef[pt]
                    pt +=1
            
        assert len(full_coef) == pt
        print >> sys.stderr, "Done"
        #print self.coeff
        return None

    def energy_predict(self, seq):
        assert len(seq) == self.NCPlen
        seq_list = [seq]
        pos_list = []

        #for i in range(self.NCPlen/2+bound, len(seq)-self.NCPlen/2-bound):
        #for i in range(self.NCPlen/2+bound, self.templatelen-self.NCPlen/2-bound):
        #    NCPseq = seq[i-self.NCPlen/2:i+self.NCPlen/2+1]
        #    seq_list.append(NCPseq)
        #    pos_list.append(i)
        
        data_dict = self.coeff
        MM_orders = []
        Kmer_k_b = []
        PolyA_b = []
        GC_b = []
        Harmonic = None
        dPolyA = []
        for key in data_dict:
            if key.startswith('MM'):
                MM_orders.append(int(key[2:]))
            if key.startswith('Kmer'):
                knum = len(data_dict[key].keys()[0])
                Kmer_k_b.append(int(key[4:]))
            if key.startswith('PolyA'):
                PolyA_b.append(int(key[5:]))
            if key.startswith('GC'):
                GC_b.append(int(key[2:]))
            if key.startswith('Harmonic'):
                Harmonic = True
            if key.startswith('dPolyA'):
                lmin = min(data_dict['dPolyA'].keys())
                lmax = max(data_dict['dPolyA'].keys())
                dPolyA = [lmin, lmax]
        if MM_orders:
            MM_orders = sorted(MM_orders)
        if Kmer_k_b:
            bnum = max(Kmer_k_b) + 1
            Kmer_k_b = ([knum, bnum])
        if PolyA_b:
            PolyA_b = max(PolyA_b) + 1
        if GC_b:
            GC_b = max(GC_b) + 1

        if MM_orders:
            MM_vars = [ [] for i in range(len(seq_list)) ] 
            for order in sorted(MM_orders):
                temp = self._var_Markov(seq_list, order)
                for i in range(len(temp)):
                    MM_vars[i] += temp[i]
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            Kmer_vars = self._var_Kmer(seq_list, knum, bnum)
        if PolyA_b:
            PolyA_vars = self._var_PolyA(seq_list, PolyA_b)
        if GC_b:
            GC_vars = self._var_GC(seq_list, GC_b)
        if Harmonic:
            Harmonic_vars = self._var_Harmonic(pos_list)
        if dPolyA:
            lmin, lmax = dPolyA
            dPolyA_vars = self._var_dPolyA(seq_list, lmin=lmin, lmax=lmax)

        var_list = [ [] for i in range(len(seq_list))]
        for i in range(len(seq_list)):
            if MM_orders:
                var_list[i] += MM_vars[i]
            if Kmer_k_b:
                var_list[i] += Kmer_vars[i]
            if PolyA_b:
                var_list[i] += PolyA_vars[i]
            if GC_b:
                var_list[i] += GC_vars[i]
            if Harmonic:
                var_list[i] += Harmonic_vars[i]
            if dPolyA:
                var_list[i] += dPolyA_vars[i]

        Ypred = self.reg.predict(var_list)
        Ypred = [ value[0] for value in Ypred]
        #pred_prob = [np.exp(-value) for value in Ypred]
        #pred_prob = norm(pred_prob)
        return Ypred[0]

    def energy_predict_profile(self, seq):
        profile = []
        assert len(seq) >= self.NCPlen
        for i in range(len(seq)):
            if i < self.NCPlen/2:
                profile.append(sys.maxint)
            elif i >= len(seq) - self.NCPlen/2:
                profile.append(sys.maxint)
            else:
                profile.append(self.energy_predict(seq[i-self.NCPlen/2:i+self.NCPlen/2+1]))
        return profile

    def _predict(self, seq, bound=0, scale=1.0):
        assert len(seq) == self.templatelen
        seq_list = []
        pos_list = []
        for i in range(self.NCPlen/2+bound, self.templatelen-self.NCPlen/2-bound):
            NCPseq = seq[i-self.NCPlen/2:i+self.NCPlen/2+1]
            if self.mask_idxs:
                NCPseq = mask(NCPseq, self.mask_idxs, replace='-')
            seq_list.append(NCPseq)
            pos_list.append(i)
        
        data_dict = self.coeff
        MM_orders = []
        Kmer_k_b = []
        PolyA_b = []
        GC_b = []
        Harmonic = None
        dPolyA = []
        for key in data_dict:
            if key.startswith('MM'):
                MM_orders.append(int(key[2:]))
            if key.startswith('Kmer'):
                knum = len(data_dict[key].keys()[0])
                Kmer_k_b.append(int(key[4:]))
            if key.startswith('PolyA'):
                PolyA_b.append(int(key[5:]))
            if key.startswith('GC'):
                GC_b.append(int(key[2:]))
            if key.startswith('Harmonic'):
                Harmonic = (len(seq) == self.templatelen)
            if key.startswith('dPolyA'):
                lmin = min(data_dict['dPolyA'].keys())
                lmax = max(data_dict['dPolyA'].keys())
                dPolyA = [lmin, lmax]
        if MM_orders:
            MM_orders = sorted(MM_orders)
        if Kmer_k_b:
            bnum = max(Kmer_k_b) + 1
            Kmer_k_b = ([knum, bnum])
        if PolyA_b:
            PolyA_b = max(PolyA_b) + 1
        if GC_b:
            GC_b = max(GC_b) + 1

        if MM_orders:
            MM_vars = [ [] for i in range(len(seq_list)) ] 
            for order in sorted(MM_orders):
                temp = self._var_Markov(seq_list, order)
                for i in range(len(temp)):
                    MM_vars[i] += temp[i]
        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            Kmer_vars = self._var_Kmer(seq_list, knum, bnum)
        if PolyA_b:
            PolyA_vars = self._var_PolyA(seq_list, PolyA_b)
        if GC_b:
            GC_vars = self._var_GC(seq_list, GC_b)
        if Harmonic:
            Harmonic_vars = self._var_Harmonic(pos_list)
        if dPolyA:
            lmin, lmax = dPolyA
            dPolyA_vars = self._var_dPolyA(seq_list, lmin=lmin, lmax=lmax)

        var_list = [ [] for i in range(len(seq_list))]
        for i in range(len(seq_list)):
            if MM_orders:
                var_list[i] += MM_vars[i]
            if Kmer_k_b:
                var_list[i] += Kmer_vars[i]
            if PolyA_b:
                var_list[i] += PolyA_vars[i]
            if GC_b:
                var_list[i] += GC_vars[i]
            if Harmonic:
                var_list[i] += Harmonic_vars[i]
            if dPolyA:
                var_list[i] += dPolyA_vars[i]

        column_mask = ~np.any(np.isnan(var_list), axis=0) 
        var_list = list(np.asarray(var_list)[:,column_mask])

        Ypred = self.reg.predict(var_list)
        Ypred = [value[0] for value in Ypred]
        pred_prob = [np.exp(-scale*value) for value in Ypred]
        pred_prob = norm(pred_prob)
        return pred_prob

    def predict(self, key_slider, keys=None, bound=0, scale=1.0):
        pred_key_slider = {}
        padding = [0.0]*(self.NCPlen/2 + bound)
        if keys == None:
            keys = key_slider.keys()
        for key in keys:
            seq = key_slider[key].seq
            prob = padding + self._predict(seq, bound=bound, scale=scale) + padding
            dyadmap = [value*100 for value in prob]
            pred_key_slider[key] = Slider(key, len(seq), len(seq)/2, 52, 52, seq, dyadmap, [], [], None, None, None, None)
        return pred_key_slider
        

    def display(self,
                data_dict,
                vmin=None,
                vmax=None):

        MM_orders = []
        Kmer_k_b = []
        PolyA_b = []
        GC_b = []
        Harmonic = None
        dPolyA = []
        

        for key in data_dict:
            if key.startswith('MM'):
                MM_orders.append(int(key[2:]))
            if key.startswith('Kmer'):
                knum = len(data_dict[key].keys()[0])
                Kmer_k_b.append(int(key[4:]))
            if key.startswith('PolyA'):
                PolyA_b.append(int(key[5:]))
            if key.startswith('GC'):
                GC_b.append(int(key[2:]))
            if key.startswith('Harmonic'):
                Harmonic = True
            if key.startswith('dPolyA'):
                lmin = min(data_dict['dPolyA'].keys())
                lmax = max(data_dict['dPolyA'].keys())
                dPolyA = [lmin, lmax]

        if MM_orders:
            MM_orders = sorted(MM_orders)
        if Kmer_k_b:
            bnum = max(Kmer_k_b) + 1
            Kmer_k_b = ([knum, bnum])
        if PolyA_b:
            PolyA_b = max(PolyA_b) + 1
        if GC_b:
            GC_b = max(GC_b) + 1

        if MM_orders:
            for order in sorted(MM_orders):
                name = 'MM' + str(order)
                img = []
                nts = all_path(order+1, 'ATCG')
                labels = []
                for nt in nts:
                    labels.append(nt)
                    row = []
                    for i in range(self.NCPlen - order):
                        row.append(data_dict[name][i][nt])
                    img.append(row)
                fig = plt.figure()
                plt.imshow(img, interpolation='none', aspect='auto', vmin=vmin, vmax=vmax)
                plt.colorbar()
                plt.show()
                plt.close()

                fig = plt.figure()
                for i in range(len(img)):
                    row = img[i]
                    label = labels[i]
                    x_axis = [ i - len(row)*0.5 + 0.5  for i in range(len(row))]
                    plt.plot(x_axis, row, label=label)
                plt.legend()
                plt.show()
                plt.close()
                
        if 1 in MM_orders:
            name = 'MM1'
            AT = ['AA', 'AT', 'TA', 'TT']
            GC = ['GG', 'GC', 'CG', 'CC']
            Y1, Y2 = np.zeros(self.NCPlen-1), np.zeros(self.NCPlen-1)
            for i in range(self.NCPlen - 1):
                for key in AT:
                    Y1[i] += data_dict[name][i][key] / len(AT)
                for key in GC:
                    Y2[i] += data_dict[name][i][key] / len(GC)
            fig = plt.figure()
            x_axis = [ i - len(Y1)*0.5 + 0.5  for i in range(len(Y1))]
            plt.plot(x_axis, Y1, label='AA/AT/TA/TT')
            plt.plot(x_axis, Y2, label='GG/GC/CG/CC')
            plt.legend()
            plt.show()
            plt.close()

        if Kmer_k_b:
            knum, bnum = Kmer_k_b
            nts = all_path(knum, 'ATCG')                

            #for nt in nts:
            #    row = []
            #    for i in range(bnum):
            #        name = 'Kmer' + str(i)
            #        row.append(data_dict[name][nt])
            #    img.append(row)
            #print img
            #fig = plt.figure()
            #plt.imshow(img, interpolation='none', aspect='auto')
            #plt.colorbar()
            #plt.show()
            #plt.close()

            nt_value_list = []
            for i in range(bnum):
                nt_value = []
                for nt in nts:
                    name = 'Kmer' + str(i)
                    nt_value.append((nt,data_dict[name][nt]))
                nt_value_list.append(nt_value)

            for i in range(bnum):
                nt_value = nt_value_list[i]
                X, Y = [], []
                for j in range(len(nt_value)):
                    X.append(nt_value[j][0])
                    Y.append(nt_value[j][1])
                fig = plt.figure()
                plt.plot(range(len(Y)), Y, '.')
                plt.show()
                plt.close()

        if PolyA_b:
            bnum = PolyA_b
            Y = []
            for i in range(bnum):
                name = 'PolyA' + str(i)
                Y.append(data_dict[name])
            fig = plt.figure()
            plt.plot(Y)
            plt.xlabel('Bin')
            plt.show()
            plt.close()

        if GC_b:
            bnum = GC_b
            Y = []
            for i in range(bnum):
                name = 'GC' + str(i)
                Y.append(data_dict[name])
            fig = plt.figure()
            plt.plot(Y)
            plt.xlabel('Bin')
            plt.show()
            plt.close()

        if Harmonic:
            k = data_dict['Harmonic']
            X = np.asarray(range(self.templatelen))
            Y = [0.5*k*(x-self.templatelen/2)**2 for x in X]
            fig = plt.figure()
            plt.plot(X,Y)
            plt.show()
            plt.close()

        if dPolyA:
            name = 'dPolyA'
            lmin, lmax = dPolyA
            img = np.zeros((lmax-lmin+1, self.NCPlen))
            for i in range(lmin, lmax+1):
                for j in range(self.NCPlen/2 + 1):
                    img[i-lmin][j] = data_dict[name][i][j]
                    if j < self.NCPlen/2:
                        img[i-lmin][self.NCPlen-1-j] = data_dict[name][i][j]
            fig = plt.figure()
            plt.imshow(img, interpolation='none', aspect='auto', vmin=vmin, vmax=vmax)
            plt.colorbar()
            plt.show()
            plt.close()

            fig = plt.figure()
            for i in range(len(img)):
                row = img[i]
                plt.plot(row, label=str(i+lmin))
            plt.legend()
            plt.show()
            plt.close()
            
            self.img = img
            #fig = plt.figure()
            #plt.plot(row, label=str(i+lmin))
            #plt.legend()
            #plt.show()
            #plt.close()

        if self.shape:
            names = ["MGW", "HelT", "ProT", "Roll"]
            fig = plt.figure()
            for name in names:
                freq = data_dict[name]
                plt.plot(freq, label=name)
            plt.legend()
            plt.show()
            plt.close()
