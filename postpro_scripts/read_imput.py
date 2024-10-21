import math
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

note1 = "_naive"
#note1 = "_linear"
#note1 = "_bayesian"

note2 = "_rect"
#note2 = "_gauss"
#note2 = "_random"

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

def read_mean (fname, iter_max=sys.maxint):
    mean_list = []
    for line in open(fname):
        if line.startswith('>'):
            if int(line[1:]) > iter_max:
                break
            if int(line[1:]) <= 0:
                data = []
                continue
            try:
                mean_data += np.asarray(data)
                mean_list.append(copy.deepcopy(mean_data))
            except:
                mean_data = copy.deepcopy(np.asarray(data))
                mean_list.append(copy.deepcopy(mean_data))
            data = []
            continue
        cols = line.strip().split(',')
        temp = []
        for value in cols:
            try:
                value = float(value)
            except:
                value = np.nan
            temp.append(value)
        data.append(temp)
    try:
        mean_data += np.asarray(data)
        mean_list.append(copy.deepcopy(mean_data))
    except:
        mean_data = copy.deepcopy(np.asarray(data))
        mean_list.append(copy.deepcopy(mean_data))
    mean_list = [mean_list[i]/float(i+1) for i in range(len(mean_list))]
    return mean_list

def plot_strand (note1, note2):
    for i in range(0, 1000, 3):
        plot_file("mtop" + note1 + note2 + ".txt", i, total=False, label="Top predict")
        plot_file("top_true" + note2 + ".txt", i)
        plot_file("mbott" + note1 + note2 + ".txt", i, total=False, label="Bott predict")
        plot_file("bott_true" + note2 + ".txt", i)
        plot_file("mtop_true" + note2 + ".txt", i)
        plot_file("mbott_true" + note2 + ".txt", i)
        plt.legend()
        plt.show()

def get_residuals (note1, note2, iter_max=sys.maxint):
    residuals = []
    true_top = read_mean("mtop_true" + note2 + ".txt", iter_max=iter_max)[-1]
    true_bott = read_mean("mbott_true" + note2 + ".txt", iter_max=iter_max)[-1]
    pre_top = read_mean("mtop" + note1 + note2 + ".txt", iter_max=iter_max)[-1]
    pre_bott = read_mean("mbott" + note1 + note2 + ".txt", iter_max=iter_max)[-1]
    true_top = true_top.ravel()
    true_bott = true_bott.ravel()
    pre_top = pre_top.ravel()
    pre_bott = pre_bott.ravel()
    pre_top = pre_top[~np.isnan(pre_top)]
    pre_bott = pre_bott[~np.isnan(pre_bott)]
    true_top = true_top[~np.isnan(true_top)]
    true_bott = true_bott[~np.isnan(true_bott)]
    residuals += list((pre_top - true_top)**2)
    residuals += list((pre_bott - true_bott)**2)
    return residuals

def R_trace (note1, note2):
    true_top = read_mean("mtop_true" + note2 + ".txt")[-1]
    true_bott = read_mean("mbott_true" + note2 + ".txt")[-1]
    true_top = true_top.ravel()
    true_bott = true_bott.ravel()
    true_top = true_top[~np.isnan(true_top)]
    true_bott = true_bott[~np.isnan(true_bott)]
    pre_top_trace = read_mean("mtop" + note1 + note2 + ".txt")
    pre_bott_trace = read_mean("mbott" + note1 + note2 + ".txt")
    total_len = min(len(pre_top_trace), len(pre_bott_trace))
    trace = []
    for i in range(total_len):
        pre_top = pre_top_trace[i]
        pre_bott = pre_bott_trace[i]
        pre_top = pre_top.ravel()
        pre_bott = pre_bott.ravel()
        pre_top = pre_top[~np.isnan(pre_top)]
        pre_bott = pre_bott[~np.isnan(pre_bott)]
        R = sum((pre_top - true_top)**2)
        R += sum((pre_bott - true_bott)**2)
        trace.append(R)
    #fig = plt.figure()
    plt.plot(range(len(trace)), trace)
    return trace


"""
a = get_residuals("_linear", "_random")
b = get_residuals("_bayesian", "_random")
fig = plt.figure()
plt.plot(a, b, '.', markersize=3)
plt.xlabel("linear R")
plt.ylabel("Bayes R")
plt.xlim([-5000, 300000])
plt.ylim([-5000, 300000])
plt.title("Random Signal")
plt.savefig("linearVSbayes_random.png", bbox_inches='tight')


a = np.asarray(get_residuals("_naive", "_random"))
b = np.asarray(get_residuals("_linear", "_random"))
c = np.asarray(get_residuals("_bayesian", "_random"))

fig = plt.figure()
plt.title("Random Signal")
plt.plot([0,1,2], [np.mean(np.log(a+1)),np.mean(np.log(b+1)),np.mean(np.log(c+1))], 'rx', markersize=8, zorder=10)
sns.stripplot(data=[np.log(a+1),np.log(b+1),np.log(c+1)], alpha=0.5, size=1, jitter=True)
plt.xticks([0,1,2], ["Naive", "Linear", "Bayes"])
plt.ylabel("log(residual+1)")
plt.savefig("residuals_random.png", bbox_inches='tight')
plt.show()
"""
