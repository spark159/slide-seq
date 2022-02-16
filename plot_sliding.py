import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func1(x, a, b, c):
    return a*np.exp(-b*x)+c

def func2(x, a, b, c, d, e):
    return a*np.exp(-b*x)+c*np.exp(-d*x)+e

def read_data (fname_list):
    sample_rep_values = {}
    sample_times = {}
    for fname in fname_list:
        for line in open(fname):
            if line.startswith('@'):
                name = line[1:].strip()
                continue
            cols = line.strip().split()
            try:
                cols = [float(col) for col in cols]
            except:
                continue
            time = cols[0]
            values = cols[1:]
            if name not in sample_rep_values:
                sample_rep_values[name] = [[] for i in range(len(values))]
            if name not in sample_times:
                sample_times[name] = []
            for k in range(len(values)):
                value = values[k]
                sample_rep_values[name][k].append(value)
            sample_times[name].append(time)
    return sample_rep_values, sample_times

path = ""
fname_list = ["0N80_NCP_sliding_1AP.csv", 
              "0N80_NCP_sliding_A1.csv", 
              "0N80_NCP_sliding_A2.csv", 
              "0N80_NCP_sliding_Ainsert.csv"]

sample_rep_values, sample_times = read_data([path + fname for fname in fname_list])
maxtime = np.max(sample_times.values())

#names = ['601 0N80', 'Top1AP 0N80', 'Bott1AP 0N80']
#names = ['601 0N80', 'Top1AP 0N80', 'Bott1AP 0N80', 'Both1AP 0N80']
#names = ['601 0N80', 'Top2AP 0N80', 'Bott2AP 0N80', 'Both2AP 0N80']
#names = ['601 0N80']
#names = ['601 0N80', 'A1L 0N80', 'A1R 0N80', 'A1LM 0N80', 'A1RM 0N80']
#names = ['601 0N80', 'A2L 0N80', 'A2R 0N80', 'A2LM 0N80', 'A2RM 0N80']
#names = ['601 0N80', 'A1LI typeI 0N80', 'A1LI typeII 0N80', 'A1RI typeI 0N80', 'A1RI typeII 0N80']
#names = ['601 0N80', 'Top 1dU 0N80', 'Bott 1dU 0N80', 'Top 2dU 0N80', 'Bott 2dU 0N80', 'Top 1AP 0N80', 'Bott 1AP 0N80', 'Top 2AP 0N80', 'Bott 2AP 0N80']
#names = ['601 0N80', 'Top 1dU 0N80', 'Bott 1dU 0N80', 'Top 2dU 0N80', 'Bott 2dU 0N80', 'Top 1AP 0N80', 'Bott 1AP 0N80', 'Top 2AP 0N80', 'Bott 2AP 0N80']
sample_names = ['601 0N80', 'A1R 0N80', 'A1RM 0N80', 'A2R 0N80', 'A2RM 0N80', 'A1RI typeI 0N80', 'A1RI typeII 0N80']
figure_names = ['WT 601', '1bp dA:dT', '2bp dA:dT', '1bp mismatch', '2bp mismatch', 'A-insertion(guide)', 'A-insertion(tracking)']
#colors = ['k', 'r', 'g', 'b', 'orange', 'm', 'y', 'cyan', 'steelblue']

colors = ['k', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']


sample_ks = {}
fig = plt.figure()
for i in range(len(sample_names)):
    name = sample_names[i]
    if name not in sample_ks:
        sample_ks[name] = []
    rep_values = sample_rep_values[name]
    color = colors[i]
    X = sample_times[name]
    for j in range(len(rep_values)):
        Y = rep_values[j]
        popt, pcov = curve_fit(func1, X, Y, bounds=([-np.inf, 0, 0], [0, np.inf, np.inf]))
        #popt, pcov = curve_fit(func1, X, Y)
        sample_ks[name].append(popt[1])
        #popt, pcov = curve_fit(func2, X, Y, bounds=([-np.inf, 0, -np.inf, 0, 0], [0, np.inf, 0, np.inf, np.inf]))
        print name, j
        print popt
        print
        #xx = np.linspace(min(X), maxtime, num=100000)
        xx = np.linspace(min(X), max(X), num=100000)
        yy = func1(xx, *popt)
        #yy = func2(xx, *popt)
        plt.plot(X, Y, '.', color=color)
        if j == 0:
            plt.plot(xx, yy, color, alpha=0.6, label=figure_names[i])
        else:
            plt.plot(xx, yy, color, alpha=0.6)

plt.xscale("log")
plt.title("0N80 Sliding")
plt.ylim([0, 1])
plt.xlabel("Time (min)")
plt.ylabel("Slided fraction")
leg = plt.legend()
for lh in leg.legendHandles:
    lh.set_alpha(1)
#plt.savefig("0N80_Sliding.png", bbox_inches='tight')
plt.savefig("0N80_Sliding_log.png", bbox_inches='tight')
plt.show()
plt.close()

k_means, k_stds = [], []
t_means, t_stds = [], []
for name in sample_names:
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
    
#fig = plt.figure()
#k_means[-1] = 0
#k_stds[-1] = 0 
##k_means[-2], k_means[-1] = 0, 0
##k_stds[-2], k_stds[-1] = 0, 0 
##plt.bar(range(len(k_means)), k_means, width=0.5, yerr=k_stds, color=colors)
#plt.barh(range(len(k_means)), k_means, xerr=k_stds, color=colors, height=0.5)
#plt.title("0N80 Sliding")
##plt.yscale("log")
#plt.xscale("log")
##plt.xticks(range(len(k_means)), names, rotation=30)
#plt.yticks(range(len(k_means)), figure_names)
##plt.ylabel("$k_{on} (min^{-1})$")
#plt.xlabel("$k_{on} (min^{-1})$")
#plt.gca().invert_yaxis()
#plt.xlim([None, 10])
#plt.tight_layout()
##plt.ylim([0.1, 10])
#plt.savefig("0N80_sliding_kons.png", bbox_inches='tight')
##plt.show()
#plt.close()

fig = plt.figure()
k_means[-1] = 0
k_stds[-1] = 0 
plt.bar(range(len(k_means)), k_means, yerr=k_stds, color=colors, width=0.5)
plt.title("0N80 Sliding")
plt.yscale("log")
plt.xticks(range(len(k_means)), figure_names, rotation=30, ha="right", rotation_mode="anchor")
plt.ylabel("$k_{on} (min^{-1})$")
plt.tight_layout()
plt.savefig("0N80_sliding_kons.png", bbox_inches='tight')
#plt.show()
plt.close()


fig = plt.figure()
t_means[-1] = 0
t_stds[-1] = 0
plt.bar(range(len(t_means)), t_means, width=0.5, yerr=t_stds, color=colors)
plt.title("0N80 Sliding")
plt.yscale("log")
plt.xticks(range(len(t_means)), figure_names, rotation=30, ha="right", rotation_mode="anchor")
plt.ylabel("Sliding time (min)")
plt.ylim([0.1, 10])
plt.tight_layout()
plt.savefig("0N80_sliding_times.png", bbox_inches='tight')
#plt.show()
plt.close()

        

"""
def read_file (fname, skip_num=2):
    X, Y = [], []
    Z = []
    count = -1
    for line in open(fname):
        count +=1
        if count < skip_num:
            continue
        cols = line.strip().split(',')
        X.append(float(cols[0]))
        Y.append(float(cols[1]))
        Z.append(float(cols[2]))
    return X, Y, Z

fnames = ["data1.csv", "data2.csv", "data3.csv"]
colors = ['r', 'g', 'b']
names = ["601", "A3LM", "A3RM"]
3
sliding_time = []
fig = plt.figure()
for i in range(len(fnames)):
    fname = fnames[i]
    color = colors[i]
    name = names[i]
    X, Y, Z = read_file(fname)
    popt, pcov = curve_fit(func1, X, Y)
    sliding_time.append(float(1)/popt[1])
    xx = np.linspace(0, int(max(X)), num=1000)
    yy = func1(xx, *popt)
    plt.plot(X, Y, '.'+color)
    plt.errorbar(X, Y, yerr=Z, fmt='.', color=color)
    plt.plot(xx, yy, color, label=name)
plt.legend()
plt.title("80N0 Sliding")
plt.xlabel("Time (min)")
plt.ylabel("Slided fraction")
plt.show()
plt.close()
    

fig = plt.figure()
plt.bar(range(3), sliding_time)
plt.xticks(range(3), names)
plt.ylabel("Sliding time (min)")
plt.show()
plt.close()
"""
