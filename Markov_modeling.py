import sys
import copy
import math
import numpy as np
import matplotlib.pyplot as plt
import SliderClass
from SliderClass import Slider
import analysis
import load
import HMM
import pickle
import seaborn as sns
import random
import graph_final
import matplotlib.ticker as tck

def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq:
        nt = nt.upper()
        if nt == 'A':
            new_seq += 'T'
        if nt == 'T':
            new_seq += 'A'
        if nt == 'C':
            new_seq += 'G'
        if nt == 'G':
            new_seq += 'C'
    return new_seq

def all_path (N, states='ATCG'):
    if N==0:
        return ['']
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def read_NCPprob (sort_fname):
    id_type_count = {}
    for line in open(sort_fname):
        if line.strip():
            read_id, type, mapped_id, cut_loc, read_seq = line.strip().split()
            if mapped_id == '*':
                continue
            #id = int(mapped_id)
            id = mapped_id
            if type != 'valid':
                type = 'invalid'
            if id not in id_type_count:
                id_type_count[id] = {'valid':0, 'invalid':0}
            id_type_count[id][type] +=1

    id_NCPprob = {}
    for id in id_type_count:
        success = id_type_count[id]['valid']
        fail = id_type_count[id]['invalid']
        id_NCPprob[id] = float(success)/(success + fail)
    return id_NCPprob

# files and parameters
sort_fname1 = "Plslib-HS_S1_L001_R.sort"
sort_fname2 = "Plslib-HS-30min_S2_L001_R.sort"
ref_fname = "plusonelib.ref"
ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
NCPlen = 145
scale = 200
sym = True
order = 0 # MarkovModel order
klen = 5 # kmer length

# load data
#pickle_fname1 = "plusonelib_new:corrected_0.pickle"
#pickle_fname2 = "plusonelib_new:corrected_30.pickle"

pickle_fname1 = "plusonelib_new_0.pickle"
pickle_fname2 = "plusonelib_new_30.pickle"

key_slider1 = pickle.load(open(pickle_fname1, 'rb'))
key_slider2 = pickle.load(open(pickle_fname2, 'rb'))

#key_slider1 = load.load_files([sort_fname1], ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None, load_ref=ref_fname)
#key_slider2 = load.load_files([sort_fname2], ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None, load_ref=ref_fname

keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))


# get NCP recontitution probability for each sequence in library
key_NCPprob1 = read_NCPprob(sort_fname1)
key_NCPprob2 = read_NCPprob(sort_fname2)
key_NCPprob = {}
for key in keys:
    NCPprob1 = key_NCPprob1[key]
    NCPprob2 = key_NCPprob2[key]
    NCPprob = 0.5*(NCPprob1+NCPprob2)
    key_NCPprob[key] = NCPprob
    #key_NCPprob[key] = 1.0

print "NCP prob is calculated"


# get sampling based on positioning probability
seq_sample0 = []
seq_sample1, seq_sample2 = [], []
for key in keys:
    seq = key_slider1[key].seq
    psig1 = analysis.norm(key_slider1[key].dyadmap)
    psig2 = analysis.norm(key_slider2[key].dyadmap)

    for i in range(NCPlen/2, ref_length-NCPlen/2):
        subseq = seq[i-NCPlen/2:i+NCPlen/2+1]
        if sym:
            rev_subseq = rev_comp(seq[i-NCPlen/2:i+NCPlen/2+1])
        count1 = key_NCPprob[key]*psig1[i]*scale
        count1 = int(round(count1))
        for k in range(count1):
            seq_sample1.append(subseq)
            if sym:
                seq_sample1.append(rev_subseq)

        count2 = key_NCPprob[key]*psig2[i]*scale
        count2 = int(round(count2))
        for k in range(count2):
            seq_sample2.append(subseq)
            if sym:
                seq_sample2.append(rev_subseq)

        seq_sample0.append(subseq)
        if sym:
            seq_sample0.append(rev_subseq)

print "sampling is done"


# get Markov transition probability
def get_markov_prob (seq_list, order):
    # prob[order][position][previous_bases + current_base]
    prob = []
    for k in range(order+1):
        prob_pos = []
        for i in range(NCPlen):
            prob_i = {}
            if i >= k:
                for prev in all_path(k):
                    prob_i[prev] = {nt:0 for nt in 'ATCG'}
            prob_pos.append(prob_i)
        prob.append(prob_pos)
    
    #prob = [[{} for j in range(NCPlen)] for i in range(order+1)]
    for seq in seq_list:
        for i in range(len(seq)):
            for j in range(order + 1):
                if i < j:
                    continue
                prev = seq[i-j:i]
                cur = seq[i]
                assert cur in "ACGT"
                prob_i = prob[j][i]
                if prev not in prob_i:
                    prob_i[prev] = {}
                prob_prev = prob_i[prev]
                if cur not in prob_prev:
                    prob_prev[cur] = 1
                else:
                    prob_prev[cur] += 1

    # Normalize the data
    for prob_order in prob:
        for prob_pos in prob_order:
            for prev, prob_prev in prob_pos.items():
                total_count = sum(prob_prev.values())
                total_count = float(total_count)
                for cur, count in prob_prev.items():
                    prob_prev[cur] = count / total_count

    return prob

# get normalized Markov transition probability
def get_markov_nprob (prob, bg_prob):
    nprob = []
    for k in range(len(prob)):
        nprob_pos = []
        for i in range(len(prob[k])):
            nprob_i = {}
            if i >= k:
                for prev in all_path(k):
                    nprob_i[prev] = {nt:0 for nt in 'ATCG'}
            nprob_pos.append(nprob_i)
        nprob.append(nprob_pos)

    #nprob = [[{} for j in range(NCPlen)] for i in range(order+1)]
    # get fold change compared to background
    for order in range(len(prob)):
        for pos in range(len(prob[order])):
            if pos < order:
                continue
            trans_prob = prob[order][pos]
            trans_prob0 = bg_prob[order][pos]
            for prev in trans_prob:
                for cur in trans_prob[prev]:
                    prob_value = trans_prob[prev][cur]
                    prob_value0 = trans_prob0[prev][cur]
                    trans_nprob = nprob[order][pos]
                    trans_nprob[prev][cur] = prob_value/prob_value0

    # Normalize the data
    for nprob_order in nprob:
        for nprob_pos in nprob_order:
            for prev, nprob_prev in nprob_pos.items():
                total_count = sum(nprob_prev.values())
                total_count = float(total_count)
                for cur, count in nprob_prev.items():
                    nprob_prev[cur] = count / total_count

    return nprob

# get kmer probability
def get_kmer_prob (seq_list, klen):
    #prob[prev][cur]
    prob = {}
    for prev in all_path(klen-1):
        for cur in 'ATCG':
            if prev not in prob:
                prob[prev] = {}
            if cur not in prob[prev]:
                prob[prev][cur] = 0
                
    for seq in seq_list:
        for i in range(len(seq)):
            if i < klen-1:
                continue
            prev = seq[i-klen+1:i]
            cur = seq[i]
            prob[prev][cur] +=1

    for prev in prob:
        total = float(sum(prob[prev].values()))
        for cur in prob[prev]:
            prob[prev][cur] /= total

    return prob

# get kmer normalized probability
def get_kmer_nprob (prob, bg_prob):
    nprob = {}
    # get fold change compared to background
    for prev in prob:
        for cur in prob[prev]:
            if prev not in nprob:
                nprob[prev] = {}
            nprob[prev][cur] = float(prob[prev][cur]/bg_prob[prev][cur])

    # normalize the data
    for prev in nprob:
        total = float(sum(nprob[prev].values()))
        for cur in nprob[prev]:
            nprob[prev][cur] /= total
    return nprob

markov_prob0 = get_markov_prob(seq_sample0, order)
markov_prob1 = get_markov_prob(seq_sample1, order)
markov_prob2 = get_markov_prob(seq_sample2, order)

markov_nprob1 = get_markov_nprob(markov_prob1, markov_prob0)
markov_nprob2 = get_markov_nprob(markov_prob2, markov_prob0)

kmer_prob0 = get_kmer_prob(seq_sample0, klen)
kmer_prob1 = get_kmer_prob(seq_sample1, klen)
kmer_prob2 = get_kmer_prob(seq_sample2, klen)

kmer_nprob1 = get_kmer_nprob(kmer_prob1, kmer_prob0)
kmer_nprob2 = get_kmer_nprob(kmer_prob2, kmer_prob0)

print "training is done"
    

# Markov Model prediction
def MM_predict (seq, prob, order):
    logscore = 0.0
    for i in range(len(seq)):
        if i < order:
            continue
        prev = seq[i-order:i]
        cur = seq[i]
        prob_value = prob[order][i][prev][cur]
        logscore += np.log(prob_value)
    return logscore

# interpolated between Markov and Kmer
def Kmer_predict(seq, kprob, klen):
    logscore = 0.0
    for i in range(len(seq)):
        if i < klen-1:
            continue
        prev = seq[i-klen+1:i]
        cur = seq[i]
        prob_value = kprob[prev][cur]
        logscore += np.log(prob_value)
    return logscore


# Interpolated Markov Model prediction
def IMM_predict (seq, prob, order1, order2, lamb):
    logscore = 0.0
    for i in range(len(seq)):
        if i < min(order1, order2):
            continue
        if i < max(order1, order2):
            order = min(order1, order2)
            prev = seq[i-order:i]
            cur = seq[i]
            prob_value = prob[order][i][prev][cur]
            logscore += np.log(prob_value)
            continue
        prev1 = seq[i-order1:i]
        prev2 = seq[i-order2:i]
        cur = seq[i]
        prob_value1 = prob[order1][i][prev1][cur]
        prob_value2 = prob[order2][i][prev2][cur]
        prob_value = lamb*prob_value1 + (1-lamb)*prob_value2
        logscore += np.log(prob_value)
    return logscore

# interpolated between Markov and Kmer
def MM_Kmer_predict(seq, mprob, kprob, order, klen, lamb):
    logscore = 0.0
    for i in range(len(seq)):
        if i < order and i < klen-1:
            continue
        if i >= order and i < klen-1:
            prev = seq[i-order:i]
            cur = seq[i]
            prob_value = mprob[order][i][prev][cur]
            logscore += np.log(prob_value)
            continue
        if i < order and i >= klen-1:
            prev = seq[i-klen+1:i]
            cur = seq[i]
            prob_value = kprob[prev][cur]
            logscore += np.log(prob_value)
            continue
        prev1 = seq[i-order:i]
        cur1 = seq[i]
        prob_value1 = mprob[order][i][prev1][cur1]
        prev2 = seq[i-klen+1:i]
        cur2 = seq[i]
        prob_value2 = kprob[prev2][cur2]
        prob_value = lamb*prob_value1 + (1-lamb)*prob_value2
        logscore += np.log(prob_value)
    return logscore


# predict scores along sequence
def predict_profile (seq, predict_func, *args, **kargs):
    prob_profile = []
    for i in range(len(seq)):
        if i < NCPlen/2:
            prob_profile.append(0.0)
            continue
        elif i >= len(seq)-NCPlen/2:
            prob_profile.append(0.0)
            continue
        subseq = seq[i-NCPlen/2:i+NCPlen/2+1]
        logscore = predict_func(subseq, *args, **kargs)
        #print logscore
        prob_profile.append(np.exp(logscore))
    total = float(sum(prob_profile))
    prob_profile = [value/total for value in prob_profile]
    return prob_profile

# Experiment VS prediction
X1, Y1 = [], []
X2, Y2 = [], []
pred_key_slider1 = {}
pred_key_slider2 = {}
for key in keys:
    slider1 = key_slider1[key]
    slider2 = key_slider2[key]
    psig1, psig2 = analysis.norm(slider1.dyadmap), analysis.norm(slider2.dyadmap)
    seq = slider1.seq
    assert len(seq) == 225
    #pred_psig1 = predict_profile(seq, MM_predict, markov_nprob1, order)
    #pred_psig2 = predict_profile(seq, MM_predict, markov_nprob2, order)

    #pred_psig1 = predict_profile(seq, IMM_predict, markov_nprob1, order, order-1, 0.5)
    #pred_psig2 = predict_profile(seq, IMM_predict, markov_nprob2, order, order-1, 0.5)

    #pred_psig1 = predict_profile(seq, Kmer_predict, kmer_nprob1, klen)
    #pred_psig2 = predict_profile(seq, Kmer_predict, kmer_nprob2, klen)

    pred_psig1 = predict_profile(seq, MM_Kmer_predict, markov_nprob1, kmer_nprob1, order, klen, 0.9)
    pred_psig2 = predict_profile(seq, MM_Kmer_predict, markov_nprob2, kmer_nprob1, order, klen, 0.9)

    #fig = plt.figure()
    #plt.plot(psig2, 'k-')
    #plt.plot(pred_psig2, 'r-')
    #plt.title(key)
    #plt.show()
    #plt.close()

    X1 += psig1[NCPlen/2:ref_length-NCPlen/2]
    Y1 += pred_psig1[NCPlen/2:ref_length-NCPlen/2]
    X2 += psig2[NCPlen/2:ref_length-NCPlen/2]
    Y2 += pred_psig2[NCPlen/2:ref_length-NCPlen/2]

    slider1.dyadmap = psig1
    slider2.dyadmap = psig2

    newslider1 = copy.deepcopy(slider1)
    newslider1.dyadmap = pred_psig1
    pred_key_slider1[key] = newslider1

    newslider2 = copy.deepcopy(slider2)
    newslider2.dyadmap = pred_psig2
    pred_key_slider2[key] = newslider2

    
# before
print analysis.get_corr(X1, Y1)
fig = plt.figure(figsize=(2, 1.8))
plt.plot(X1, Y1, 'k.', markersize=2, alpha=0.1)
sns.kdeplot(X1, Y1, cmap='Reds', shade=False, shade_lowest=False, zorder=100, linewidth=0.1, alpha=0.5)
plt.xlabel("Experiment", fontsize=6)
plt.ylabel("Prediction", fontsize=6)
plt.title("HS data", fontsize=8)
plt.gca().tick_params(axis='both', which='major', labelsize=5)
plt.gca().tick_params(axis='both', which='minor', labelsize=5)
plt.xlim([0, 0.03])
plt.ylim([0, 0.03])
plt.savefig("ExVsPre1.png", dpi=500, bbox_inches='tight')
plt.close()

# after
print analysis.get_corr(X2, Y2)
fig = plt.figure(figsize=(2, 1.8))
plt.plot(X2, Y2, 'k.', markersize=2, alpha=0.1)
sns.kdeplot(X2, Y2, cmap='Reds', shade=False, shade_lowest=False, zorder=100, linewidth=0.1, alpha=0.5)
plt.xlabel("Experiment", fontsize=6)
plt.ylabel("Prediction", fontsize=6)
plt.title("Chd1 data", fontsize=8)
plt.gca().tick_params(axis='both', which='major', labelsize=5)
plt.gca().tick_params(axis='both', which='minor', labelsize=5)
plt.xlim([0, 0.03])
plt.ylim([0, 0.03])
plt.savefig("ExVsPre2.png", dpi=500, bbox_inches='tight')
plt.close()



tlen = 225
random.seed(0)
sample_list = random.sample(key_slider1.keys(), 5)

for key in sample_list:
    fig = plt.figure(figsize=(1.85, 1.2))
    plt.plot(key_slider1[key].dyadmap, 'k--', label='Exp', alpha=0.6, lw=0.8)
    plt.plot(pred_key_slider1[key].dyadmap, 'r-', label='Model', alpha=0.6, lw=0.8)
    plt.xticks([i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)], fontsize=5)
    plt.yticks(fontsize=5)
    plt.title('%s (HS)' % key, fontsize=8)
    plt.gca().xaxis.set_minor_locator(tck.AutoMinorLocator(2))
    leg = plt.legend(prop={'size': 5})
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
        lh.set_linewidth(1.0)
    plt.savefig('%s_1.svg' % (key), format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

for key in sample_list:
    fig = plt.figure(figsize=(1.85, 1.2))
    plt.plot(key_slider2[key].dyadmap, 'k--', label='Exp', alpha=0.6, lw=0.8)
    plt.plot(pred_key_slider2[key].dyadmap, 'r-', label='Model', alpha=0.6, lw=0.8)
    plt.xticks([i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)], fontsize=5)
    plt.yticks(fontsize=5)
    plt.title('%s (Chd1)' % key, fontsize=8)
    plt.gca().xaxis.set_minor_locator(tck.AutoMinorLocator(2))
    leg = plt.legend(prop={'size': 5})
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
        lh.set_linewidth(1.0)
    plt.savefig('%s_2.svg' % (key), format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

sys.exit(1)

graph_final.plot_map(key_slider1, Slider.KDE, ids=sample_list, thickness=[8, 2, 2], cmap='YlGnBu', save=True, note='_before', xticks=[[i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)]], figscale=150, fontsize=5)
graph_final.plot_map(pred_key_slider1, Slider.KDE, ids=sample_list, thickness=[8, 2, 2], cmap='YlGnBu', save=True, note='_before_pred', xticks=[[i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)]], figscale=150, fontsize=5, norm=False)
graph_final.plot_map(key_slider2, Slider.KDE, ids=sample_list, thickness=[8, 2, 2], cmap='YlGnBu', save=True, note='_after', xticks=[[i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)]], figscale=150, fontsize=5)
graph_final.plot_map(pred_key_slider2, Slider.KDE, ids=sample_list, thickness=[8, 2, 2], cmap='YlGnBu', save=True, note='_after_pred', xticks=[[i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)]], figscale=150, fontsize=5, norm=False)

