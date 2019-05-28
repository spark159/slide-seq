import math
import random
import pickle
import EnModel
import SliderClass
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

NCPlen = 147

def key_cmp (key1, key2):
    if float(key1) < float(key2):
        return -1
    else:
        return 1

def tuple_cmp (a, b):
    if a[0] < b[0]:
        return -1
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

# read slider data
key_slider1 = pickle.load(open("slider1.p", "rb"))
key_slider2 = pickle.load(open("slider2.p", "rb"))

keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
keys = sorted(keys, cmp=key_cmp)


# JJ model vs Linear model

m1 = EnModel.EnergyModel(key_slider1, NCPlen=NCPlen)
m1.train(MM_orders=[1], Kmer_k_b=False, PolyA_b=False,  GC_b=False, Harmonic=False)
m2 = EnModel.EnergyModel(key_slider2, NCPlen=NCPlen)
m2.train(MM_orders=[1], Kmer_k_b=False, PolyA_b=False,  GC_b=False, Harmonic=False)

coeff1 = m1.coeff['MM1']
coeff2 = m2.coeff['MM1']

# check dinucleotides
all_din = all_path(2, 'ATCG')

din_count1 = {}
for din in all_din:
    if din not in din_count1:
        din_count1[din] = 0
    total = 0
    for i in range(len(coeff1)):
        try:
            count = coeff1[i][din]
        except:
            continue
        din_count1[din] += count
        total += 1
    din_count1[din] = float(din_count1[din]) / total

din_count2 = {}
for din in all_din:
    if din not in din_count2:
        din_count2[din] = 0
    total = 0
    for i in range(len(coeff2)):
        try:
            count = coeff2[i][din]
        except:
            continue
        din_count2[din] += count
        total += 1
    din_count2[din] = float(din_count2[din]) / total

dinCount1, dinCount2 = [], []

for din in all_din:
    if din not in din_count1:
        count1 = 0
    else:
        count1 = din_count1[din]
    if din not in din_count2:
        count2 = 0
    else:
        count2 = din_count2[din]
    dinCount1.append([count1, din])
    dinCount2.append([count2, din])

dinCount1 = sorted(dinCount1, cmp=tuple_cmp)
dinCount2 = sorted(dinCount2, cmp=tuple_cmp)

X1, Y1 = [], []
X2, Y2 = [], []
for i in range(len(dinCount1)):
    din = dinCount2[i][1]
    X1.append(din)
    X2.append(din)
    Y1.append(din_count1[din])
    Y2.append(din_count2[din])

fig = plt.figure()
plt.plot(range(len(Y1)), Y1, 'x-', label='Heat Shift')
plt.plot(range(len(Y2)), Y2, 'x-', label='Chd1 Sliding')
plt.xticks(range(len(X1)), X1)
plt.title("Average along position")
plt.legend()
plt.show()
plt.close()


"""
fig = plt.figure()
X1, Y1 = [], []
#X2, Y2 = [], []
for i in range(len(dinCount1)):
    X1.append(dinCount1[i][1])
    Y1.append(dinCount1[i][0])
    #X2.append(dinCount2[0])
    #Y2.append(dinCount2[1])
plt.plot(range(len(Y1)), Y1, 'x-')
plt.xticks(range(len(Y1)), X1)
plt.show()
plt.close()

fig = plt.figure()
X, Y = [], []
for i in range(len(dinCount2)):
    X.append(dinCount2[0])
    Y.append(dinCount2[1])
plt.plot(range(len(Y)), Y, 'x-')
plt.xticks(range(len(Y)), X)
plt.show()
plt.close()
"""
