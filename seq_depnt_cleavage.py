import math
import random
import pickle
import EnModel
import LinModel
import SliderClass
import numpy as np
import matplotlib.pyplot as plt

def key_cmp (key1, key2):
    if float(key1) < float(key2):
        return -1
    else:
        return 1

def tuple_cmp(a, b):
    if a[0] < b[0]:
        return -1
    else:
        return 1
    
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

def norm(signal_list):
    total = sum(signal_list)
    return [float(value)/total for value in signal_list]

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

# Nucleosomal DNA length
NCPlen = 147

# k-mer length and kmers
k = 3

# offset length
offset = 52

# read slider data
key_slider1 = pickle.load(open("/home/../../media/spark159/sw/dataforslide/slider1.p", "rb"))
key_slider2 = pickle.load(open("/home/../../media/spark159/sw/dataforslide/slider2.p", "rb"))

keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
keys = sorted(keys, cmp=key_cmp)

# read cleavage count for each kmers
kmer_back_top, kmer_back_bott = {}, {}
kmer_count1_top, kmer_count1_bott = {}, {}
kmer_count2_top, kmer_count2_bott = {}, {}
kmer_ncount1_top, kmer_ncount1_bott = {}, {}
kmer_ncount2_top, kmer_ncount2_bott = {}, {}
for key in keys:
    slider1 = key_slider1[key]
    slider2 = key_slider2[key]
    seq = slider1.seq
    top_cutmap1, top_cutmap2 = slider1.right_cutmap, slider2.right_cutmap
    bott_cutmap1, bott_cutmap2 = slider1.left_cutmap, slider2.left_cutmap

    for i in range(NCPlen/2, len(top_cutmap1)-NCPlen/2):
        count1, count2 = top_cutmap1[i], top_cutmap2[i]
        kmer = seq[i-1-k/2:i-1+k/2+1]
        #kmer = seq[i-1:i-1+2]
        if kmer not in kmer_back_top:
            kmer_back_top[kmer] = 0
        kmer_back_top[kmer] += 1
        if kmer not in kmer_count1_top:
            kmer_count1_top[kmer] = 0
            kmer_ncount1_top[kmer] = 0
        kmer_count1_top[kmer] += float(count1)
        kmer_ncount1_top[kmer] += float(count1) / sum(top_cutmap1)
        if kmer not in kmer_count2_top:
            kmer_count2_top[kmer] = 0
            kmer_ncount2_top[kmer] = 0
        kmer_count2_top[kmer] += float(count2)
        kmer_ncount2_top[kmer] += float(count2) / sum(top_cutmap2)
    
    for i in range(NCPlen/2, len(bott_cutmap1)-NCPlen/2):
        count1, count2 = bott_cutmap1[i], bott_cutmap2[i]
        kmer = rev_comp(seq[i+1-k/2:i+1+k/2+1])
        #kmer = rev_comp(seq[i:i+2])
        if kmer not in kmer_back_bott:
            kmer_back_bott[kmer] = 0
        kmer_back_bott[kmer] += 1
        if kmer not in kmer_count1_bott:
            kmer_count1_bott[kmer] = 0
            kmer_ncount1_bott[kmer] = 0
        kmer_count1_bott[kmer] += float(count1)
        kmer_ncount1_bott[kmer] += float(count1) / sum(bott_cutmap1)
        if kmer not in kmer_count2_bott:
            kmer_count2_bott[kmer] = 0
            kmer_ncount2_bott[kmer] = 0
        kmer_count2_bott[kmer] += float(count2)
        kmer_ncount2_bott[kmer] += float(count2) / sum(bott_cutmap2)

kmers = list(set(all_path(k, 'ATCG')) - set(kmer_back_top.keys()))
for kmer in kmers:
    kmer_back_top[kmer] = 0
    kmer_count1_top[kmer] = 0
    kmer_ncount1_top[kmer] = 0
    kmer_count2_top[kmer] = 0
    kmer_ncount2_top[kmer] = 0

kmers = list(set(all_path(k, 'ATCG')) - set(kmer_back_bott.keys()))
for kmer in kmers:
    kmer_back_bott[kmer] = 0
    kmer_count1_bott[kmer] = 0
    kmer_ncount1_bott[kmer] = 0
    kmer_count2_bott[kmer] = 0
    kmer_ncount2_bott[kmer] = 0
    
kmer_count_list = [kmer_count1_top, kmer_count1_bott, kmer_back_top, kmer_back_bott]

total_kmer = []
for kmer in all_path(k, 'ATCG'):
    total = 0
    for kmer_count in kmer_count_list:
        total += kmer_count[kmer]
    total_kmer.append([total, kmer])
total_kmer = sorted(total_kmer, cmp=tuple_cmp, reverse=True)

X = []
Y1, Y2 = [], []
Y3, Y4 = [], []
for count, kmer in total_kmer:
    X.append(kmer)
    Y1.append(kmer_count1_top[kmer])
    Y2.append(kmer_count1_bott[kmer])
    Y3.append(kmer_back_top[kmer])
    Y4.append(kmer_back_bott[kmer])
    
fig = plt.figure()
plt.plot(range(len(X)), Y1, label="HS top")
plt.plot(range(len(X)), Y2, label="HS bott")
plt.plot(range(len(X)), Y3, label="Back top")
plt.plot(range(len(X)), Y4, label="Back bott")
plt.xticks(range(len(X)), X, rotation=90)
plt.legend()
#plt.show()
plt.close()

kmer_count_list = [kmer_count2_top, kmer_count2_bott, kmer_back_top, kmer_back_bott]
kmers = set(kmer_back_top.keys()) | set(kmer_back_bott.keys())

total_kmer = []
for kmer in all_path(k, 'ATCG'):
    total = 0
    for kmer_count in kmer_count_list:
        total += kmer_count[kmer]
    total_kmer.append([total, kmer])
total_kmer = sorted(total_kmer, cmp=tuple_cmp, reverse=True)

X = []
Y1, Y2 = [], []
for count, kmer in total_kmer:
    X.append(kmer)
    Y1.append(kmer_count2_top[kmer])
    Y2.append(kmer_count2_bott[kmer])

fig = plt.figure()
plt.plot(range(len(X)), Y1, label="Chd1 top")
plt.plot(range(len(X)), Y2, label="Chd1 bott")
plt.plot(range(len(X)), Y3, label="Back top")
plt.plot(range(len(X)), Y4, label="Back bott")
plt.xticks(range(len(X)), X, rotation=90)
plt.legend()
#plt.show()
plt.close()

# estimate bott/top cleavage ratio
g1 = sum(kmer_count1_bott.values()) / sum(kmer_count1_top.values())
g2 = sum(kmer_count2_bott.values()) / sum(kmer_count2_top.values())
top_total = sum(kmer_count1_top.values()) + sum(kmer_count2_top.values())
bott_total = sum(kmer_count1_bott.values()) + sum(kmer_count2_bott.values())
g = float(bott_total) / top_total


# estimate kmer dependent cleavage efficiency
def fold (a, b):
    if b > 0:
        return float(a)/float(b)
    assert a == 0 and b == 0 
    return 1.0

kmers = all_path(k, 'ATCG')
kmer_top_fold1, kmer_bott_fold1 = {}, {}
kmer_top_fold2, kmer_bott_fold2 = {}, {}
kmer_fold1, kmer_fold2 = {}, {}
kmer_fold = {}
for kmer in kmers:
    top_fold1 = fold(kmer_ncount1_top[kmer], kmer_back_top[kmer])
    bott_fold1 = fold(kmer_ncount1_bott[kmer], kmer_back_bott[kmer])
    top_fold2 = fold(kmer_ncount2_top[kmer], kmer_back_top[kmer])
    bott_fold2 = fold(kmer_ncount2_bott[kmer], kmer_back_bott[kmer])
    fold1 = top_fold1 + bott_fold1
    fold2 = top_fold2 + bott_fold2
    if kmer not in kmer_top_fold1:
        kmer_top_fold1[kmer] = 0
    kmer_top_fold1[kmer] += top_fold1
    if kmer not in kmer_bott_fold1:
        kmer_bott_fold1[kmer] = 0
    kmer_bott_fold1[kmer] += bott_fold1
    if kmer not in kmer_top_fold2:
        kmer_top_fold2[kmer] = 0
    kmer_top_fold2[kmer] += top_fold2
    if kmer not in kmer_bott_fold2:
        kmer_bott_fold2[kmer] = 0
    kmer_bott_fold2[kmer] += bott_fold2
    if kmer not in kmer_fold1:
        kmer_fold1[kmer] = 0
    kmer_fold1[kmer] += fold1
    if kmer not in kmer_fold2:
        kmer_fold2[kmer] = 0
    kmer_fold2[kmer] += fold2
    if kmer not in kmer_fold:
        kmer_fold[kmer] = 0
    kmer_fold[kmer] += fold1 + fold2

kmer_top_freq1, kmer_top_freq2 = {}, {}
kmer_bott_freq1, kmer_bott_freq2 = {}, {}
kmer_freq1, kmer_freq2 = {}, {}
kmer_freq = {}
for kmer in kmers:
    kmer_top_freq1[kmer] = float(kmer_top_fold1[kmer]) / sum(kmer_top_fold1.values())
    kmer_bott_freq1[kmer] = float(kmer_bott_fold1[kmer]) / sum(kmer_bott_fold1.values())
    kmer_top_freq2[kmer] = float(kmer_top_fold2[kmer]) / sum(kmer_top_fold2.values())
    kmer_bott_freq2[kmer] = float(kmer_bott_fold2[kmer]) / sum(kmer_bott_fold2.values())
    kmer_freq1[kmer] = float(kmer_fold1[kmer]) / sum(kmer_fold1.values())
    kmer_freq2[kmer] = float(kmer_fold2[kmer]) /sum(kmer_fold2.values())
    kmer_freq[kmer] = float(kmer_fold[kmer]) / sum (kmer_fold.values())

topY1, bottY1 = [], []
topY2, bottY2 = [], []
Y1, Y2 = [], []
Y = []
for kmer in kmers:
    topY1.append(kmer_top_freq1[kmer])
    bottY1.append(kmer_bott_freq1[kmer])
    topY2.append(kmer_top_freq2[kmer])
    bottY2.append(kmer_bott_freq2[kmer])
    Y1.append(kmer_freq1[kmer])
    Y2.append(kmer_freq2[kmer])
    Y.append(kmer_freq[kmer])

fig = plt.figure()
plt.plot(range(len(kmers)), topY1, label='HS top')
plt.plot(range(len(kmers)), bottY1, label='HS bott')
plt.plot(range(len(kmers)), topY2, label='Chd1 top')
plt.plot(range(len(kmers)), bottY2, label='Chd1 bott')
plt.xticks(range(len(kmers)), kmers, rotation=90)
plt.legend()
#plt.show()
plt.close()

fig = plt.figure()
plt.plot(range(len(kmers)), Y1, label='HS')
plt.plot(range(len(kmers)), Y2, label='Chd1')
plt.plot(range(len(kmers)), Y, label='all')
plt.xticks(range(len(kmers)), kmers, rotation=90)
plt.legend()
#plt.show()
plt.close()


# build linear model from collected kmers
seq_list = []
score_list = []
count_list = []
for kmer, freq in kmer_freq.items():
    seq_list.append(kmer)
    score_list.append(freq)
    count_list.append(freq)

m = LinModel.SeqLinearModel(seq_list, score_list, count_list)
m.train(MM_orders=[1], Kmer_k_b=None, PolyA_b=False, GC_b=False, Harmonic=False, sym=False)


# estimate the corrected dyad signal
for key in keys:
    slider1 = key_slider1[key]
    slider2 = key_slider2[key]
    seq = slider1.seq
    top_cutmap1, top_cutmap2 = slider1.right_cutmap, slider2.right_cutmap
    bott_cutmap1, bott_cutmap2 = slider1.left_cutmap, slider2.left_cutmap
    new_dyadmap1, new_dyadmap2 = [], []
    for i in range(len(top_cutmap1)):
        if i < offset+k:
            signal1, signal2 = 0.0, 0.0
        elif i >= len(top_cutmap1) - (offset+k):
            signal1, signal2 = 0.0, 0.0
        else:
            top_kmer = seq[i-offset-1-k/2:i-offset-1+k/2+1]
            bott_kmer = rev_comp(seq[i+offset+1-k/2:i+offset+1+k/2+1])
            #top_kmer = seq[i-1:i-1+2]
            #bott_kmer = rev_comp(seq[i:i+2])
            t = kmer_freq[top_kmer]
            b = kmer_freq[bott_kmer]
            #t = m.score_predict(top_kmer)
            #b = m.score_predict(bott_kmer)
            signal1 = float(top_cutmap1[i-offset])/t + float(bott_cutmap1[i+offset])/(g*b)
            signal2 = float(top_cutmap2[i-offset])/t + float(bott_cutmap2[i+offset])/(g*b)  
        new_dyadmap1.append(signal1)
        new_dyadmap2.append(signal2)

    new_dyadmap1 = norm(new_dyadmap1)
    new_dyadmap2 = norm(new_dyadmap2)
    dyadmap1 = norm(slider1.dyadmap)
    dyadmap2 = norm(slider2.dyadmap)

    #key_slider1[key].dyadmap = new_dyadmap1
    #key_slider2[key].dyadmap = new_dyadmap2
    key_slider1[key].dyadmap = dyadmap1
    key_slider2[key].dyadmap = dyadmap2


    fig = plt.figure()
    plt.plot(dyadmap1, label='HS old')
    plt.plot(new_dyadmap1, label='HS new')
    plt.legend()
    #plt.show()
    plt.close()

    fig = plt.figure()
    plt.plot(dyadmap2, label='Chd1 old')
    plt.plot(new_dyadmap2, label='Chd1 new')
    plt.legend()
    #plt.show()
    plt.close()
    
    
# check the correction by checking nucleosome positioning sequences
seq_list = []
score_list1, count_list1 = [], []
score_list2, count_list2 = [], []

for key in keys:
    dyadmap1 = key_slider1[key].dyadmap
    dyadmap2 = key_slider2[key].dyadmap
    seq = key_slider1[key].seq
    for i in range(NCPlen/2, len(seq)-NCPlen/2):
        NCPseq = seq[i - NCPlen/2 : i + NCPlen/2 + 1]
        score1, score2 = dyadmap1[i], dyadmap2[i]
        count1, count2 = dyadmap1[i], dyadmap2[i]
        seq_list.append(NCPseq)
        score_list1.append(score1)
        score_list2.append(score2)
        count_list1.append(count1)
        count_list2.append(count2)

m1 = LinModel.SeqLinearModel(seq_list, score_list1, count_list1)
m1.report(MM_orders=[1], Kmer_k_b=False, PolyA_b=False, GC_b=False, Harmonic=False, sym=True)
m2 = LinModel.SeqLinearModel(seq_list, score_list2, count_list2)
m2.report(MM_orders=[1], Kmer_k_b=False, PolyA_b=False, GC_b=False, Harmonic=False, sym=True)

"""
group_freq1 = m1.spectrum(MM_orders=[1], Kmer_k_b=False, PolyA_b=False, GC_b=False, Harmonic=False, sym=True)
for i in range(len(group_freq1)):
    freq = group_freq1[i]
    m1.display(freq)
    
group_freq2 = m2.spectrum(MM_orders=[1], Kmer_k_b=False, PolyA_b=False, GC_b=False, Harmonic=False, sym=True)
for i in range(len(group_freq2)):
    freq = group_freq2[i]
    m2.display(freq)
"""
