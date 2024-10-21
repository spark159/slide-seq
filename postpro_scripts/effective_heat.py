import math
import random
import pickle
import EnModel
import SliderClass
import numpy as np
import matplotlib.pyplot as plt

def norm (data_list):
    total = sum(data_list)
    assert total > 0
    return [float(data)/total for data in data_list]

def rescale (prob_list, scale):
    new_prob_list = []
    for i in range(len(prob_list)):
        prob = prob_list[i]
        new_prob = prob**scale
        new_prob_list.append(new_prob)
    return norm(new_prob_list)

# read slider data
key_slider1 = pickle.load(open("slider1.p", "rb"))
key_slider2 = pickle.load(open("slider2.p", "rb"))

keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))

X, Y = [], []
for key in keys:
    slider1 = key_slider1[key]
    slider2 = key_slider2[key]
    prob1 = norm(slider1.dyadmap)
    prob2 = norm(slider2.dyadmap)
    X.append(slider1.entropy())
    Y.append(slider2.entropy())
    prob3 = rescale (prob1, 1.011562225135828)
    fig = plt.figure()
    plt.plot(prob1, label='Heat Shift', alpha=0.8, linewidth=1)
    plt.plot(prob2, label='Chd1 Sliding', alpha=0.8, linewidth=1)
    #plt.plot(prob3, label='Chd1 Sliding (rescale)', alpha=0.8, linewidth=1, linestyle='--')
    plt.legend()
    plt.show()
    plt.close()

fig = plt.figure()
plt.plot(X,Y,'.')
plt.xlabel("Entropy (Heat Shift)")
plt.ylabel("Entropy (Chd1 Sliding)")
plt.xlim([4,5])
plt.ylim([4,5])
plt.plot([4,5],[4,5],'--')
#plt.savefig("entropy.png")
plt.show()
plt.close()
