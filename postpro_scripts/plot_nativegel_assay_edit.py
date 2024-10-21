import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# sliding 
def func1(x, a, b, c):
    return a*np.exp(-b*x)+c

# competition
def func2(x, a, b, c, d, e):
    return a*np.exp(-b*x)+c*np.exp(-d*x)+e

# binding
def func3(x, a, b, c):
    return a*x/(x+b) + c

def read_data (fname_list):
    sample_rep_values = {}
    sample_times = {}
    for fname in fname_list:
        for line in open(fname):
            if line.startswith('@'):
                name = line[1:].strip()
                if name in sample_times:
                    repeat = True
                else:
                    repeat = False
                continue
            if repeat:
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

def curve_fitting (sample_times, sample_rep_values, exp, graph=False, sample_names=None):
    sample_ks = {}

    if sample_names == None:
        sample_names = sample_times.keys()
    
    if graph:
        fig = plt.figure(figsize=(5,4))

    for i in range(len(sample_names)):
        name = sample_names[i]
        color = colors[i]

        if name not in sample_ks:
            sample_ks[name] = []

        rep_values = sample_rep_values[name]

        X = sample_times[name]
        for j in range(len(rep_values)):
            Y = rep_values[j]
            if exp == 'sliding':
                #popt, pcov = curve_fit(func1, X, Y, bounds=([-np.inf, 0, 0], [0, np.inf, np.inf]), p0=[-0.82260531, 7.14972408, 0.94119315])
                popt, pcov = curve_fit(func1, X, Y, bounds=([-10, 0, 0], [0, 10, 10]))
                sample_ks[name].append(popt[1])
            elif exp == 'competition':
                popt, pcov = curve_fit(func2, X, Y, bounds=([-np.inf,0,-np.inf,0, 0], [np.inf,np.inf,np.inf,np.inf, np.inf]), p0=[ 8.49028742e-01, 4.02463916e-03, -8.03394545e-02, 1.09613822e+01, 6.66899789e-11])
                idx = np.argmax([abs(popt[0]), abs(popt[2])])
                sample_ks[name].append(popt[idx+1])
                #k = (abs(popt[0])*popt[1] + abs(popt[2])*popt[3])/float(abs(popt[0]) + abs(popt[2]))
                #sample_ks[name].append(k)
            elif exp == 'binding':
                popt, pcov = curve_fit(func3, X, Y, bounds=([0,0,0], [np.inf,np.inf,np.inf]))
                sample_ks[name].append(popt[1])

            if graph:
                if exp == 'sliding':
                    xx = np.linspace(min(X), max(X), num=100000)
                    yy = func1(xx, *popt)
                elif exp == 'competition':
                    xx = np.linspace(min(X), max(X), num=100000)
                    yy = func2(xx, *popt)
                elif exp == 'binding':
                    xx = np.linspace(min(X), max(X), num=100000)
                    yy = func3(xx, *popt)
                    
                plt.plot(X, Y, '.', color=color)
                if j == 0:
                    plt.plot(xx, yy, color, alpha=0.6, label=legend_names[i])
                else:
                    plt.plot(xx, yy, color, alpha=0.6)

    if graph:
        if exp == "sliding":
            plt.xlabel("Time (min)")
            plt.ylabel("Fraction shifted", fontsize=6)
            #plt.title("Nucleoosme Sliding assay")
        elif exp == 'competition':
            plt.xlabel("Time (min)")
            plt.ylabel("Fraction two-bound", fontsize=6)
            #plt.title("Nucleosome Competition assay")
        elif exp == 'binding':
            plt.xlabel("Chd1 (nM)")
            plt.ylabel("Bound fraction", fontsize=6)
            
        plt.xscale("log")
        leg = plt.legend()
        plt.ylim([0, 1])
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        plt.savefig("Curve_fitting_" + exp + ".svg", format='svg', bbox_inches='tight')
        #plt.show()
        plt.close()
    
    return sample_ks

def combine_replicates (sample_rep_values):
    sample_means, sample_stds = {}, {}
    for sample in sample_rep_values:
        rep_values = np.asarray(sample_rep_values[sample])
        means = np.mean(rep_values, axis=0)
        stds = np.std(rep_values, axis=0)
        sample_means[sample] = means
        sample_stds[sample] = stds
    return sample_means, sample_stds

def plot_data (sample_times, sample_means, sample_stds, sample_names, exp):
    fig = plt.figure(figsize=(2.4,1.7))
    for i in range(len(sample_names)):
        name = sample_names[i]
        X = sample_times[name]
        Y = sample_means[name]
        Z = sample_stds[name]

        if exp == 'sliding':
            popt, pcov = curve_fit(func1, X, Y, bounds=([-np.inf, 0, 0], [0, np.inf, np.inf]), p0=[-0.82260531, 7.14972408, 0.94119315])
            #if name == '601 0N80':
            #    print popt
            xx = np.linspace(min(X), max(X), num=100000)
            yy = func1(xx, *popt)
        elif exp == 'competition':
            popt, pcov = curve_fit(func2, X, Y, bounds=([-np.inf,0,-np.inf,0, 0], [np.inf,np.inf,np.inf,np.inf, np.inf]), p0=[ 8.49028742e-01, 4.02463916e-03, -8.03394545e-02, 1.09613822e+01, 6.66899789e-11])
            xx = np.linspace(min(X), max(X), num=100000)
            yy = func2(xx, *popt)
        elif exp == 'binding':
            popt, pcov = curve_fit(func3, X, Y, bounds=([0,0,0], [np.inf,np.inf,np.inf]))
            xx = np.linspace(min(X), max(X), num=100000)
            yy = func3(xx, *popt)

        plt.plot(X, Y, '.', markersize=3, color=colors[i], label=legend_names[i])
        plt.errorbar(X, Y, yerr=Z, fmt='.', markersize=3, lw=1, color=colors[i])
        plt.plot(xx, yy, '-', color=colors[i], alpha=0.6)

    if exp == "sliding":
        plt.xlabel("Time (min)", fontsize=6)
        plt.ylabel("Fraction shifted", fontsize=6)
        #plt.title("Nucleosome Sliding assay")
    elif exp == 'competition':
        plt.xlabel("Time (min)", fontsize=6)
        plt.ylabel("Fraction two-bound", fontsize=6)
        #plt.title("Nucleosome Competition assay")
    elif exp == 'binding':
        plt.xlabel("Chd1 (nM)", fontsize=6)
        plt.ylabel("Fraction bound", fontsize=6)

    plt.xscale("log")
    plt.ylim([0, 1])
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    leg = plt.legend(frameon=False, loc='upper left', fontsize=5)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        lh._legmarker.set_markersize(5)
    plt.savefig("mean_data_graph_" + exp + ".svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

def plot_barplot (sample_ks, exp, clip_off=False):
    k_means = [np.mean(sample_ks[name]) for name in sample_names]
    k_stds = [np.std(sample_ks[name]) for name in sample_names]
    fig = plt.figure(figsize=(2.6, 1.6))
    if clip_off:
        k_means[-1] = 0
        k_stds[-1] = 0 
    plt.bar(range(len(k_means)), k_means, yerr=k_stds, color=colors, width=0.5)
    if exp == 'sliding':
        #plt.title("Apparent Nucleosome Sliding Rate")
        plt.ylabel("Sliding rate (min$^{-1}$)", fontsize=6)
    elif exp == 'competition':
        #plt.title("Apparent Chd1 Dissociation Rate")
        plt.ylabel("Dissociation rate (min$^{-1}$)", fontsize=6)
    elif exp == 'binding':
        plt.ylabel("Effective $k_{d}$", fontsize=6)

    plt.yscale("log")
    plt.xticks(range(len(k_means)), legend_names, rotation=30, ha="right", rotation_mode="anchor",
               fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig("rate_" + exp + ".svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    
# load sliding data
path = ""
#sliding_fname_list = ["0N80_NCP_sliding_A1.csv", 
#                      "0N80_NCP_sliding_A2.csv", 
#                      "0N80_NCP_sliding_Ainsert.csv",
#                      "0N80_NCP_sliding_1AP.csv"]

sliding_fname_list = ["New_80N0_native_sliding.csv"]

sliding_rep_values, sliding_times = read_data([path + fname for fname in sliding_fname_list])

# load competition data
#path = ""
#compet_fname_list = ['0N80_competition_assay.csv']
#compet_rep_values, compet_times = read_data([path + fname for fname in compet_fname_list])

# load binding data
path = ""
binding_fname_list = ["New_80N0_Chd1_binding_assay_with_AMPPNP.csv"]
binding_rep_values, binding_times = read_data([path + fname for fname in binding_fname_list])


# set parameters
#sample_names = ['601 0N80', 'A1R 0N80', 'A1RM 0N80', 'A2R 0N80', 'A2RM 0N80', 'A1RI typeI 0N80', 'A1RI typeII 0N80']
#legend_names = ['WT 601', '1bp dA:dT', '2bp dA:dT', '1bp mismatch', '2bp mismatch', 'A-insertion(guide)', 'T-insertion(tracking)']
#sample_names = ['601 0N80', 'Top1AP 0N80', 'Bott1AP 0N80', 'Both1AP 0N80']
#legend_names = ['WT 601', '1AP(guide)', '1AP(tracking)', '1AP(both)']

sample_names = ['601 80N0', '2bp mm SHL2', '2bp mm SHL2.7', '1AP SHL2', '1AP SHL2.7']
legend_names = ['WT 601', '2bp mismatch SHL2 ', '2bp mismatch SHL2.7', '1AP SHL2', '1AP SHL2.7']

colors = ['k', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

# curve fitting
sliding_ks = curve_fitting(sliding_times, sliding_rep_values, exp='sliding', graph=True, sample_names=sample_names)
#compet_ks = curve_fitting(compet_times, compet_rep_values, exp='competition', graph=True, sample_names=sample_names)
binding_ks = curve_fitting(binding_times, binding_rep_values, exp='binding', graph=True, sample_names=sample_names)

# collapse the replicates
sliding_means, sliding_stds = combine_replicates(sliding_rep_values)
#compet_means, compet_stds = combine_replicates(compet_rep_values)
binding_means, binding_stds = combine_replicates(binding_rep_values)

# plot mean data
plot_data(sliding_times, sliding_means, sliding_stds, sample_names, exp='sliding')
#plot_data(compet_times, compet_means, compet_stds, sample_names, exp='competition')
plot_data(binding_times, binding_means, binding_stds, sample_names, exp='binding')

# plot rate bargraph
plot_barplot(sliding_ks, exp='sliding')
#plot_barplot(compet_ks, exp='competition')
plot_barplot(binding_ks, exp='binding')

# scater plot (sliding rate VS dissociation rate)
k_means1 = [np.mean(sliding_ks[name]) for name in sample_names]
k_stds1 = [np.std(sliding_ks[name]) for name in sample_names]
#k_means2 = [np.mean(compet_ks[name]) for name in sample_names]
#k_stds2 = [np.std(compet_ks[name]) for name in sample_names]

k_means2 = [np.mean(binding_ks[name]) for name in sample_names]
k_stds2 = [np.std(binding_ks[name]) for name in sample_names]


#k_means1[-1] = 0
#k_stds1[-1] = 0



fig = plt.figure(figsize=(2.4,1.7))
plt.scatter(k_means1, k_means2, s=3, c='k', edgecolor='k')
plt.errorbar(k_means1, k_means2, xerr=k_stds1, yerr=k_stds2,  mec='k', ms=2, lw=1, fmt='.')
for i in range(len(sample_names)):
    plt.annotate(' ' + legend_names[i], (k_means1[i], k_means2[i]), ha='left', va='bottom', fontsize=6)
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Sliding rate (min$^{-1}$)", fontsize=6)
#plt.ylabel("Dissociation rate (min$^{-1}$)", fontsize=6)
plt.ylabel("Effective $k_{d}$", fontsize=6)
#plt.title("Nucleosome Sliding VS Chd1 dissociation rate")
plt.xticks(fontsize=5)
plt.yticks(fontsize=5)
plt.savefig("sliding_VS_dissociation.svg", format='svg', bbox_inches='tight')
#plt.show()
plt.close()
