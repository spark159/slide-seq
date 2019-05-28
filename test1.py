import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

def read(fname, line_num=0, neg=False):    
    sig_list = []
    #line_num = 125
    for line in open(fname):
        if line.startswith('>'):
            count = 0
            continue
        if count == line_num :
            cols = line.strip().split(',')
            cols = [float(value) for value in cols]
            sig_list.append(cols)
        count +=1

    sig_list = np.asarray(sig_list)
    m,n = sig_list.shape
    #mode_sig = [stats.mode(sig_list[:,i]) for i in range(n)]
    median_sig = [np.median(sig_list[:,i]) for i in range(n)]

    if not neg:
        x = range(len(median_sig))
    else:
        x = [i-len(median_sig) for i in range(len(median_sig))]
    
    #fig = plt.figure()
    for sig in sig_list:
        plt.plot(x, sig, 'y', linewidth=1, alpha=0.2)
    #plt.plot(x, median_sig, 'k')
    plt.plot(x, median_sig)
    plt.plot(x, sig_list[0], 'k--')
    #plt.show()

"""
sig_list = []
line_num = 0
for line in open("eplist.txt"):
    if line.startswith('>'):
        count = 0
        continue
    if count == line_num :
        cols = line.strip().split(',')
        cols = [float(value) for value in cols]
        sig_list.append(cols)
    count +=1

sig_list = np.asarray(sig_list)
m,n = sig_list.shape
#mode_sig = [stats.mode(sig_list[:,i]) for i in range(n)]
median_sig = [np.median(sig_list[:,i]) for i in range(n)]

fig = plt.figure()
for sig in sig_list:
    plt.plot(sig, 'y', linewidth=1, alpha=0.2)
plt.plot(median_sig, 'k')
plt.plot(sig_list[0], 'r--')
plt.show()
"""
