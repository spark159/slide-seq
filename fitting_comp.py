import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func1(x, a, b, c):
    return a*np.exp(-b*x)+c

def func2(x, a, b, c, d, e):
    return a*np.exp(-b*x)+c*np.exp(-d*x)+e

def read_data (fname):
    sample_times = {}
    sample_values = {}
    sample_errors = {}
    for line in open(fname):
        if line.startswith('@'):
            name = line[1:].strip()
            continue
        cols = line.strip().split()
        try:
            cols = [float(col) for col in cols]
        except:
            continue
        time = int(cols[0])
        value = float(cols[1])
        error = float(cols[2])
        if name not in sample_times:
            sample_times[name] = []
        sample_times[name].append(time)
        if name not in sample_values:
            sample_values[name] = []
        sample_values[name].append(value)
        if name not in sample_errors:
            sample_errors[name] = []
        sample_errors[name].append(error)
    return sample_times, sample_values, sample_errors

#sample_times, sample_values, sample_errors = read_data("/home/spark159/../../media/spark159/sw/dataforslide/0N80_mismatch_insertion_comp.csv")
sample_times, sample_values, sample_errors = read_data("competition_2nd.csv")


#names = ['601 0N80', 'A1L 0N80', 'A1R 0N80', 'A1LM 0N80', 'A1RM 0N80', 'A2L 0N80', 'A2R 0N80', 'A2LM 0N80', 'A2RM 0N80', 'A1IL_typeI', 'A1IL_typeII', 'A1IR_typeI', 'A1IR_typeII']
#names = ['601 0N80', 'A1L 0N80', 'A1R 0N80', 'A1LM 0N80', 'A1RM 0N80', 'A2L 0N80', 'A2R 0N80', 'A2LM 0N80', 'A2RM 0N80']

names1 = ['601 0N80', 'A1L 0N80', 'A1R 0N80', 'A1LM 0N80', 'A1RM 0N80']
names2 = ['601 0N80', 'A2L 0N80', 'A2R 0N80', 'A2LM 0N80', 'A2RM 0N80']
names3 = ['601 0N80', 'A1IL_typeI', 'A1IL_typeII', 'A1IR_typeI', 'A1IR_typeII']
names4 = ['601 0N80', 'AP control', 'Top1AP', 'Bott1AP', 'Both1AP']
names5 = ['601 0N80', 'AP control', 'Top2AP', 'Bott2AP', 'Both2AP']
names6 = ['601 0N80', 'A1R 0N80', 'A2R 0N80', 'A1RM 0N80', 'A2RM 0N80', 'A1IR_typeI', 'A1IR_typeII']
names7 = ['601 0N80', 'A1R 0N80', 'A1RM 0N80', 'A1IR_typeI', 'A1IR_typeII', 'Top2AP', 'Bott2AP', 'Both2AP']
names_list = [names1, names2, names3, names4, names5]
#names_list = [names7]

colors = ['k', 'r', 'g', 'b', 'orange', 'm', 'y', 'cyan', 'steelblue']
color_list = np.linspace(0.01, 1, num=len(sample_times))
cmap = mpl.cm.get_cmap("Set1")

sample_ks = {}

for k in range(len(names_list)):
    names = names_list[k]
    fig = plt.figure()
    for i in range(len(names)):
        name = names[i]
        if name not in sample_ks:
            sample_ks[name] = []
        X = sample_times[name]
        Y = sample_values[name]
        Z = sample_errors[name]
        #popt, pcov = curve_fit(func1, X, Y)#, bounds=([-np.inf, -np.inf, 0], [0, 0, np.inf]))
        popt, pcov = curve_fit(func2, X, Y, bounds=([-np.inf, 0, -np.inf, 0, 0], [0, np.inf, 0, np.inf, np.inf]))
        #sample_ks[name].append(popt[1])
        print name
        print popt
        print
        xx = np.linspace(min(X), max(X), num=100000)
        #yy = func1(xx, *popt)
        #print yy
        #print 
        yy = func2(xx, *popt)
        plt.plot(X, Y, '.', markersize=10, color=colors[i], label=name)
        plt.errorbar(X, Y, yerr=Z, fmt='.', color=colors[i])
        #plt.plot(xx, yy, color=cmap(color_list[i]), alpha=0.6, label = name)

    #plt.xscale("log")
    plt.title("Nucleosome competition assay")
    plt.ylim([0, 1])
    plt.xlabel("Time (min)")
    plt.ylabel("Two-bound fraction")
    leg = plt.legend()
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    plt.savefig("0N80_comp_" + str(k+1) + ".png", bbox_inches='tight')
    plt.show()
    plt.close()

"""
k_means, k_stds = [], []
t_means, t_stds = [], []
for name in names:
    ks = sample_ks[name]
    k_mean = np.mean(ks)
    k_means.append(k_mean)
    k_std = np.std(ks)
    k_stds.append(k_std)
    ts = [1.0/k for k in ks]
    t_mean = np.mean(ts)
    t_means.append(t_mean)
    t_std = np.std(ts)
    t_stds.append(t_std)
    
fig = plt.figure()
k_means[-1] = 0
k_stds[-1] = 0
plt.bar(range(len(k_means)), k_means, width=0.5, yerr=k_stds, color=colors)
plt.title("0N80 Sliding")
plt.yscale("log")
plt.xticks(range(len(k_means)), names, rotation=30)
plt.ylabel("$k_{on} (min^{-1})$")
plt.ylim([0.1, 10])
plt.savefig("0N80_sliding_kons.png", bbox_inches='tight')
#plt.show()
plt.close()

fig = plt.figure()
t_means[-1] = 0
t_stds[-1] = 0
plt.bar(range(len(t_means)), t_means, width=0.5, yerr=t_stds, color=colors)
plt.title("0N80 Sliding")
plt.yscale("log")
plt.xticks(range(len(t_means)), names, rotation=30)
plt.ylabel("Sliding time (min)")
plt.ylim([0.1, 10])
plt.savefig("0N80_sliding_times.png", bbox_inches='tight')
#plt.show()
plt.close()

        

"""
