import math
import random
import pickle
import LinModel
import SliderClass
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

def key_cmp (key1, key2):
    if float(key1) < float(key2):
        return -1
    else:
        return 1

def norm(L):
    total = 0.0
    for value in L:
        total += value
    return [value/total for value in L]

def mask (seq, idxs, type="skips", space=""):
    new_seq = ""
    idxs = sorted(idxs)
    if type == "skips":
        new_seq += seq[0:idxs[0]] + space
        for i in range(len(idxs)-1):
            st, ed = idxs[i]+1, idxs[i+1]
            if st == ed:
                continue
            new_seq += seq[st:ed] + space
        new_seq += seq[idxs[-1]+1:]
    elif type == "includes":
        for i in range(len(idxs)-1):
            idx = idxs[i]
            next_idx = idxs[i+1]
            new_seq += seq[idx]
            if next_idx - idx > 1:
                new_seq += space
        new_seq += seq[next_idx] 
    return new_seq

# parameters
NCPlen = 147
scale = 100
mask_type, mask_idxs = None, None
#mask_type, mask_idxs = "includes", range(NCPlen/2-70-3, NCPlen/2-70+3) + range(NCPlen/2+70-3, NCPlen/2+70+3)
#space = 'G'
#mask_type, mask_idxs = "skips", range(NCPlen/2-70-3, NCPlen/2-70+3) + range(NCPlen/2+70-3, NCPlen/2+70+3)
#space = 'G'


# read slider data
key_slider1 = pickle.load(open("slider1.p", "rb"))
key_slider2 = pickle.load(open("slider2.p", "rb"))

keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
keys = sorted(keys, cmp=key_cmp)

seq_list = []
score_list1, count_list1 = [], []
score_list2, count_list2 = [], []

for key in keys:
    slider1 = key_slider1[key]
    slider2 = key_slider2[key]
    energy_profile1 = slider1.energy_profile()
    energy_profile2 = slider2.energy_profile()
    prob_profile1 = norm(slider1.dyadmap)
    prob_profile2 = norm(slider2.dyadmap)
    seq = slider1.seq
    for i in range(NCPlen/2, len(seq)-NCPlen/2):
        NCPseq = seq[i-NCPlen/2:i+NCPlen/2+1]
        if mask_type:
            NCPseq = mask(NCPseq, mask_idxs, type=mask_type, space='G')
        #score1, score2 = energy_profile1[i], energy_profile2[i]
        score1, score2 = prob_profile1[i], prob_profile2[i]
        count1, count2 = prob_profile1[i], prob_profile2[i]
        seq_list.append(NCPseq)
        score_list1.append(score1)
        score_list2.append(score2)
        count_list1.append(count1)
        count_list2.append(count2)

# subsample of data by score


        

m1 = LinModel.SeqLinearModel(seq_list, score_list1, count_list1)
m1.train([1], None, None, None, None)
m2 = LinModel.SeqLinearModel(seq_list, score_list2, count_list2)
m2.train([1], None, None, None, None)

predict_list1, predict_list2 = [], []
for i in range(len(seq_list)):
    NCPseq = seq_list[i]
    predict1 = m1.predict(NCPseq)[0]
    predict2 = m2.predict(NCPseq)[0]
    predict_list1.append(predict1)
    predict_list2.append(predict2)

fig = plt.figure()
plt.title("Linear model prediction")
plt.plot(score_list1, predict_list1, '.', alpha=0.2, label="Heat Shift")
plt.plot(score_list2, predict_list2, '.', alpha=0.2, label="Chd1 Sliding")
plt.plot(np.linspace(0,0.2,100), np.linspace(0,0.2,100), 'k--')
plt.xlabel("Measured Probability")
plt.ylabel("Predicted Probability")
leg = plt.legend()

for lh in leg.legendHandles:
    lh._legmarker.set_alpha(1)                    

plt.savefig("Linmodel.png", bbox_inches='tight')

#plt.show()
plt.close()
