import numpy as np
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import random
import pickle

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

def tuple_cmp(a,b):
    if a[1] < b[1]:
        return -1
    else:
        return 0

output1 = open("m1shape.pkl", "rb")
output2 = open("m2shape.pkl", "rb")
m1 = pickle.load(output1)
m2 = pickle.load(output2)

"""
# MM comparison
order=1
NCPlen= 147
for order in range(1,2):
    name = 'MM' + str(order)
    img1, img2, img3 = [], [], []
    nts = all_path(order+1, 'ATCG')
    labels = []
    for nt in nts:
        labels.append(nt)
        row1, row2, row3 = [], [], []
        for i in range(NCPlen - order):
            row1.append(m1[name][i][nt])
            row2.append(m2[name][i][nt])
            fold = (m2[name][i][nt] - m1[name][i][nt])/m1[name][i][nt]
            row3.append(fold)
        img1.append(row1)
        img2.append(row2)
        img3.append(row3)
    fig = plt.figure()
    plt.imshow(img3, interpolation='none', aspect='auto', vmin=-0.1, vmax=0.1)
    plt.yticks(range(len(nts)), labels)
    plt.colorbar()
    plt.show()
    plt.close()

    fig = plt.figure()
    for i in range(len(img3)):
        plt.plot(img3[i], label=labels[i])
    plt.legend()
    plt.show()
    plt.close()

# K-mer comparison
Kmer_k_b = [5,1]
knum, bnum = Kmer_k_b
nts = all_path(knum, 'ATCG')                
nt_value_list = []
for i in range(bnum):
    nt_value = []
    for nt in nts:
        name = 'Kmer' + str(i)
        fold = (m2[name][nt] - m1[name][nt]) / m1[name][nt]
        nt_value.append((nt,fold))
    nt_value = sorted(nt_value, cmp=tuple_cmp)
    nt_value_list.append(nt_value)

for i in range(bnum):
    nt_value = nt_value_list[i]
    X, Y = [], []
    label_loc = {}
    for j in range(len(nt_value)):
        nt = nt_value[j][0]
        if nt == 'AAAAA':
            label_loc[nt] = j
        if nt == 'GGGGG':
            label_loc[nt] = j
        value = nt_value[j][1]
        X.append(nt)
        Y.append(value)
    fig = plt.figure()
    plt.plot(Y, '.')
    locs, labels = [], []
    for label, loc in label_loc.items():
        locs.append(loc)
        labels.append(label)
    plt.xticks(locs, labels)
    plt.axhline(0, color='k', linestyle='--')
    plt.show()
    plt.close()
"""

## shape compare
names = ["MGW", "HelT", "ProT", "Roll"]
fig = plt.figure()
for name in names:
    freq = (np.asarray(m2[name]) - np.asarray(m1[name])) #/ abs(np.asarray(m1[name]))
    plt.plot(freq, label=name)
plt.legend()
plt.show()
plt.close()

names = ["MGW", "HelT", "ProT", "Roll"]
for name in names:
    fig = plt.figure()
    freq = (np.asarray(m2[name]) - np.asarray(m1[name])) #/ abs(np.asarray(m1[name]))
    plt.plot(freq, label=name)
    plt.legend()
    plt.show()
    plt.close()
