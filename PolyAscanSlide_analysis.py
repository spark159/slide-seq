import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from SliderClass import Slider
import load
import graph
import graph_edit
import analysis
import sample
import pickle
from sklearn import linear_model

ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
NCP_len = 147
mtype_choice = ['M']
size_st, size_ed = 3, 25

# load 601 control data
fnames1 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_before_.combined.sort"]
fnames2 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_after_.combined.sort"]
Control1 = load.load_files(fnames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)['601']
Control2 = load.load_files(fnames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)['601'] 

# pickle initialization
try:
    f1 = open("temp1.pickle", "rb")
    f2 = open("temp2.pickle", "rb")

    key_slider1 = pickle.load(f1)
    key_slider2 = pickle.load(f2)

except:
    fnames1 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/Ascan0_S1_L001_R.oldsort"]
    fnames2 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/Ascan-5min_S1_L001_R.oldsort"]
    #fnames1 = ["/home/spark159/../../media/spark159/sw/all_slide_seq_data/Ascan0_S1_L001_R.sort"]
    #fnames2 = ["/home/spark159/../../media/spark159/sw/all_slide_seq_data/Ascan-5min_S1_L001_R.sort"]
    
    key_slider1 = load.load_files(fnames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill="linear", mtype_choice=mtype_choice)
    key_slider2 = load.load_files(fnames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill="linear", mtype_choice=mtype_choice)
    
    f1 = open("temp1.pickle", "wb")
    pickle.dump(key_slider1, f1)
    f1.close()

    f2 = open("temp2.pickle", "wb")
    pickle.dump(key_slider2, f2)
    f2.close()

# define keys
keys = []
for key in list(set(key_slider1.keys()) & set(key_slider2.keys())):
    loc, mtype, nts = key.split('-')
    if mtype not in mtype_choice:
        continue
    if len(nts) < size_st:
        continue
    if len(nts) > size_ed:
        continue
    keys.append(key)


# dyad map plot (601, polyA, logratio)
key_logratio1, key_logratio2 = {}, {}
Size_logratio1, Size_logratio2 = {}, {}
control_map1 = analysis.norm(Control1.dyadmap)
control_map2 = analysis.norm(Control2.dyadmap)
for i in range(len(keys)):
    key = keys[i]
    dyadmap1 = analysis.norm(key_slider1[key].dyadmap)
    dyadmap2 = analysis.norm(key_slider2[key].dyadmap)
    logratio1 = [np.log2(float(dyadmap1[i]+1)/(control_map1[i]+1)) for i in range(ref_length)]
    logratio2 = [np.log2(float(dyadmap2[i]+1)/(control_map2[i]+1)) for i in range(ref_length)]
    loc, mtype, nts = key.split('-')
    Size, loc = len(nts), int(loc)
    ref_pt = loc + Size/2
    if Size not in Size_logratio1:
        Size_logratio1[Size] = []
    if Size not in Size_logratio2:
        Size_logratio2[Size] = []
    Size_logratio1[Size].append(logratio1[ref_pt-40:ref_pt+41])
    Size_logratio2[Size].append(logratio2[ref_pt-40:ref_pt+41])
    assert key not in key_logratio1
    key_logratio1[key] = Slider(key,
                               ref_length,
                               dyad_axis,
                               dyad_offset,
                               dyad_offset,
                               "",
                               logratio1,
                               [],
                               [],
                               None,
                               None,
                               None,
                               None)
    assert key not in key_logratio2
    key_logratio2[key] = Slider(key,
                               ref_length,
                               dyad_axis,
                               dyad_offset,
                               dyad_offset,
                               "",
                               logratio2,
                               [],
                               [],
                               None,
                               None,
                               None,
                               None)


Size_list = sorted(Size_logratio1.keys())
sample_mode = "polyA:" + "-".join([str(Size) for Size in Size_list])
sample_list = sample.sampling(key_logratio1, sample_mode)
#graph.plot_map(key_slider1, sample_list, norm_choice=True, draw_key=True, draw_vert=False, note="_polyA_before")
#graph.plot_map(key_slider2, sample_list, norm_choice=True, draw_key=True, draw_vert=False, note="_polyA_after")
#graph_edit.plot_map(key_logratio1, sample_list, norm_choice=False, obs_func = Slider.get_dyadmap, draw = "key", note="_polyA_logratio_before")
#graph_edit.plot_map(key_logratio2, sample_list, norm_choice=False, obs_func = Slider.get_dyadmap, draw = "key", note="_polyA_logratio_after")


#sys.exit(0)

# log ratio signal around Poly-A
Size_logratios = [Size_logratio1, Size_logratio2]
names = ['before', 'after']
for k in range(len(Size_logratios)):
    Size_logratio = Size_logratios[k]
    for i in range(size_st, size_ed+1):
        fig = plt.figure()
        X = range(-40, 41)
        for trace in Size_logratio[i]:
            plt.plot(X, trace, alpha=1)
        plt.title("PolyA " + str(i) + "bp")
        plt.xlabel("Position w.r.t polyA (bp)")
        plt.ylabel("log2 (PolyA/Control)")
        #plt.savefig("polyA_align_" + names[k] + '_' + str(i) + ".png", bbox_inches='tight')
        #plt.show()
        plt.close()

for u in range(len(Size_logratios)):
    Size_logratio = Size_logratios[u]
    fig1 = plt.figure(1)
    fig2 = plt.figure(2)
    fig3 = plt.figure(3)
    for i in range(size_st, size_ed+1):
        X = range(-40, 41)
        max_trace = [0.0]*81
        min_trace = [0.0]*81
        mean_trace = np.asarray([0.0]*81)
        for trace in Size_logratio[i]:
            for k in range(len(trace)):
                if max_trace[k] < trace[k]:
                    max_trace[k] = trace[k]
                if min_trace[k] > trace[k]:
                    min_trace[k] = trace[k]
            mean_trace += np.asarray(trace)
        mean_trace = mean_trace/ float(len(Size_logratio[i]))
        plt.figure(1)
        plt.plot(X, max_trace, alpha=1, label=str(i) + 'bp')
        plt.figure(2)
        plt.plot(X, min_trace, alpha=1, label=str(i) + 'bp')
        plt.figure(3)
        plt.plot(X, mean_trace, alpha=1, label=str(i) + 'bp')
    plt.figure(1)
    plt.title("Max log2 (PolyA/Control)")
    plt.figure(2)
    plt.title("Min log2 (PolyA/Control)")
    plt.figure(3)
    plt.title("Mean log2 (PolyA/Control)")
    for i in range(3):
        plt.figure(i+1)
        plt.xlabel("Position w.r.t polyA (bp)")
        plt.ylabel("log2 (PolyA/Control)")
        plt.legend()
        #plt.savefig("polyA_sum_" + names[u] + str(i) + ".png", bbox_inches='tight')
    #plt.show()
    plt.close('all')


# left/right scattering plot
left, right = [], []
for key in keys:
    loc, mtype, nts = key.split('-')
    size, loc = len(nts), int(loc)
    if size % 2 != 0:
        loc = loc + size/2
    else:
        loc = loc + size/2 - 0.5
    if loc < dyad_axis:
        left.append(key)
    elif loc > dyad_axis:
        right.append(key)
sample_list = [left, right]
graph.plot_corr2(key_slider1, Slider.Amer_len, Slider.median_pos, xlabel='Poly-A length', ylabel='Mean position', sample_labels = ['PolyA in Left','PolyA in Right'], sample_list=sample_list, note="before")
graph.plot_corr2(key_slider2, Slider.Amer_len, Slider.median_pos, xlabel='Poly-A length', ylabel='Mean position', sample_labels = ['PolyA in Left','PolyA in Right'], sample_list=sample_list, note="after")


# calculate mean/KL-div sensitivity map
Size_loc_mlist1, Size_loc_mlist2 = {}, {}
Size_loc_mlist3, Size_loc_mlist4 = {}, {}
Size_loc_KL1, Size_loc_KL2 = {}, {}

for key in keys:
    loc, mtype, nts = key.split('-')
    size, loc = len(nts), int(loc)
            
    slider1 = key_slider1[key]
    slider2 = key_slider2[key]

    mean_pos1 = slider1.median_pos()
    mean_pos2 = slider2.median_pos()
    mean_pos3 = mean_pos1 - Control1.median_pos()
    mean_pos4 = mean_pos2 - Control2.median_pos()

    KL1 = analysis.KL_div(analysis.norm(slider1.get_dyadmap()), analysis.norm(Control1.get_dyadmap()))
    KL2 = analysis.KL_div(analysis.norm(slider2.get_dyadmap()), analysis.norm(Control2.get_dyadmap()))

    st = loc
    ed = st + size
    assert st >= 0
    assert ed <= ref_length 

    if size not in Size_loc_mlist1:
        Size_loc_mlist1[size] = {}
    if size not in Size_loc_mlist2:
        Size_loc_mlist2[size] = {}
    if size not in Size_loc_mlist3:
        Size_loc_mlist3[size] = {}
    if size not in Size_loc_mlist4:
        Size_loc_mlist4[size] = {}

    if size not in Size_loc_KL1:
        Size_loc_KL1[size] = {}
    if size not in Size_loc_KL2:
        Size_loc_KL2[size] = {}

    for i in range(st, ed):
        if i not in Size_loc_mlist1[size]:
            Size_loc_mlist1[size][i] = []
        if i not in Size_loc_mlist2[size]:
            Size_loc_mlist2[size][i] = []
        if i not in Size_loc_mlist3[size]:
            Size_loc_mlist3[size][i] = []
        if i not in Size_loc_mlist4[size]:
            Size_loc_mlist4[size][i] = []

        Size_loc_mlist1[size][i].append(mean_pos1)
        Size_loc_mlist2[size][i].append(mean_pos2)
        Size_loc_mlist3[size][i].append(mean_pos3)
        Size_loc_mlist4[size][i].append(mean_pos4)

        if i not in Size_loc_KL1[size]:
            Size_loc_KL1[size][i] = []
        Size_loc_KL1[size][i].append(KL1)

        if i not in Size_loc_KL2[size]:
            Size_loc_KL2[size][i] = []
        Size_loc_KL2[size][i].append(KL2)


def plot_MKscatter (Size_loc_data, ylabel="", labels=[None], note="", ylim=None, hline=True):
    fig = plt.figure()
    color_list = np.linspace(0.01, 1, num=len(Size_loc_data))
    cmap = cm.get_cmap("jet")
    for i in range(len(Size_loc_data)):
        size = Size_loc_data.keys()[i]
        loc_data = Size_loc_data[size]
        X, Y = [], []
        for loc in loc_data:
            X += [loc for k in range(len(loc_data[loc]))]
            Y += loc_data[loc]
        plt.plot(X, Y, '.', markersize=3, alpha=0.5, label=labels[i], color=cmap(color_list[i]))
    plt.axvline(x=ref_length/2, linestyle='--', color='k')
    if hline:
        plt.axhline(y=0, linestyle='--', color='k')
    Lticks = [-60, -50, -40, -30, -20, -10]
    Rticks = [10, 20, 30, 40, 50, 60]
    for tick in Lticks:
        plt.axvline(x=tick+ref_length/2, linestyle='--', color='r', alpha=0.2)
    for tick in Rticks:
        plt.axvline(x=tick+ref_length/2, linestyle='--', color='b', alpha=0.2)
    plt.xticks([ref_length/2 + value for value in range(-60, 61, 10)], [str(value) for value in range(-60, 61, 10)])
    plt.xlabel("DNA template (bp)")
    plt.ylabel(ylabel)
    if ylim != None:
        plt.ylim(ylim)
    leg = plt.legend(loc='best', numpoints=1, prop={'size': 10})
    for lh in leg.legendHandles:
        lh._legmarker.set_markersize(10)
        lh._legmarker.set_alpha(1)
    plt.savefig('MKscatter' + note + '.png', bbox_inches='tight')
    plt.close()

def plot_MKmean (Size_loc_data, ylabel="", labels=[None], note="", ylim=None, hline=True):
    fig = plt.figure()
    color_list = np.linspace(0.01, 1, num=len(Size_loc_data))
    cmap = cm.get_cmap("jet")
    for i in range(len(Size_loc_data)):
        size = Size_loc_data.keys()[i]
        loc_data = Size_loc_data[size]
        X, Y = [], []
        Z = []
        for loc in sorted(loc_data.keys()):
            X.append(loc)
            Y.append(np.mean(loc_data[loc]))
            #Y.append(np.median(loc_data[loc]))
        plt.plot(X, Y, '.', markersize=5, label=labels[i], color=cmap(color_list[i]))
        #plt.plot(X, Y, label=labels[i])
    plt.axvline(x=ref_length/2, linestyle='--', color='k')
    if hline:
        plt.axhline(y=0, linestyle='--', color='k')
    #Lticks = [-60, -40, -20, -10]
    #Rticks = [60, 40, 20, 10]
    Lticks = [-60, -50, -40, -30, -20, -10]
    Rticks = [10, 20, 30, 40, 50, 60]
    for tick in Lticks:
        plt.axvline(x=tick+ref_length/2, linestyle='--', color='r', alpha=0.2)
    for tick in Rticks:
        plt.axvline(x=tick+ref_length/2, linestyle='--', color='b', alpha=0.2)
    plt.xticks([ref_length/2 + value for value in range(-60, 61, 10)], [str(value) for value in range(-60, 61, 10)])
    plt.xlabel("DNA template (bp)")
    plt.ylabel(ylabel)
    if ylim != None:
        plt.ylim(ylim)
    leg = plt.legend(loc='best', numpoints=1, prop={'size': 10})
    for lh in leg.legendHandles:
        lh._legmarker.set_markersize(10)
        lh._legmarker.set_alpha(1)
    plt.savefig('MKmean' + note + '.png', bbox_inches='tight')
    plt.close()

plot_MKscatter(Size_loc_mlist3, labels=[str(size) + " bp" for size in Size_loc_mlist3], ylabel="Mean position (bp)", note='mean_before')
plot_MKmean(Size_loc_mlist3, labels=[str(size) + " bp" for size in Size_loc_mlist3], ylabel="Mean position (bp)", note='mean_before')
plot_MKscatter(Size_loc_mlist4, labels=[str(size) + " bp" for size in Size_loc_mlist4], ylabel="Mean position (bp)", note='mean_after')
plot_MKmean(Size_loc_mlist4, labels=[str(size) + " bp" for size in Size_loc_mlist4], ylabel="Mean position (bp)", note='mean_after')

plot_MKscatter(Size_loc_KL1, labels=[str(size) + " bp" for size in Size_loc_KL1], ylabel="KL-divergence (Bits)", note='KL_before', hline=False)
plot_MKmean(Size_loc_KL1, labels=[str(size) + " bp" for size in Size_loc_KL1], ylabel="KL-divergence (Bits)", note='KL_before', hline=False)
plot_MKscatter(Size_loc_KL2, labels=[str(size) + " bp" for size in Size_loc_KL2], ylabel="KL-divergence (Bits)", note='KL_after', hline=False)
plot_MKmean(Size_loc_KL2, labels=[str(size) + " bp" for size in Size_loc_KL2], ylabel="KL-divergence (Bits)", note='KL_after', hline=False)


# calculate energy/force sensitivity map
try:
    f1 = open("temp_energy_before.pickle", "rb")
    f2 = open("temp_force_before.pickle", "rb")
    f3 = open("temp_energy_after.pickle", "rb")
    f4 = open("temp_force_after.pickle", "rb")
    
    Size_loc_energy1 = pickle.load(f1)
    Size_loc_force1 = pickle.load(f2)
    Size_loc_energy2 = pickle.load(f3)
    Size_loc_force2 = pickle.load(f4)

except:
    Size_loc_energy1, Size_loc_force1 = {}, {}
    Size_loc_energy2, Size_loc_force2 = {}, {}

    for key in keys:
        loc, mtype, nts = key.split('-')
        size, loc = len(nts), int(loc)

        slider1 = key_slider1[key]
        slider2 = key_slider2[key]

        energy1 = slider1.energy_profile() - Control1.energy_profile()
        force1 = slider1.force_profile() - Control1.force_profile()
        energy2 = slider2.energy_profile() - Control2.energy_profile()
        force2 = slider2.force_profile() - Control2.force_profile()

        if size not in Size_loc_energy1:
            Size_loc_energy1[size] = {}
        if size not in Size_loc_force1:
            Size_loc_force1[size] = {}        
        if size not in Size_loc_energy2:
            Size_loc_energy2[size] = {}
        if size not in Size_loc_force2:
            Size_loc_force2[size] = {}
        

        for i in range(NCP_len/2, ref_length-NCP_len/2):
            st = max(loc - i + NCP_len/2, 0)
            ed = min(loc - i + NCP_len/2 + size, NCP_len)
            if ed - st <= 0:
                continue
            for k in range(st, ed):
                if k not in Size_loc_energy1[size]:
                    Size_loc_energy1[size][k] = []
                Size_loc_energy1[size][k].append(energy1[i])
                if k not in Size_loc_force1[size]:
                    Size_loc_force1[size][k] = []
                Size_loc_force1[size][k].append(force1[i])

                if k not in Size_loc_energy2[size]:
                    Size_loc_energy2[size][k] = []
                Size_loc_energy2[size][k].append(energy2[i])
                if k not in Size_loc_force2[size]:
                    Size_loc_force2[size][k] = []
                Size_loc_force2[size][k].append(force2[i])
                
            # symmetrize
            st = max(-(loc-i + size -1) + NCP_len/2, 0)
            ed = min(-(loc-i-1) + NCP_len/2, NCP_len)
            if ed - st <= 0:
                continue
            for k in range(st, ed):
                if k not in Size_loc_energy1[size]:
                    Size_loc_energy1[size][k] = []
                Size_loc_energy1[size][k].append(energy1[i])
                if k not in Size_loc_force1[size]:
                    Size_loc_force1[size][k] = []
                Size_loc_force1[size][k].append(-force1[i])

                if k not in Size_loc_energy2[size]:
                    Size_loc_energy2[size][k] = []
                Size_loc_energy2[size][k].append(energy2[i])
                if k not in Size_loc_force2[size]:
                    Size_loc_force2[size][k] = []
                Size_loc_force2[size][k].append(-force2[i])


    f1 = open("temp_energy_before.pickle", "wb")
    pickle.dump(Size_loc_energy1, f1)
    f1.close()

    f2 = open("temp_force_before.pickle", "wb")
    pickle.dump(Size_loc_force1, f2)
    f2.close()

    f3 = open("temp_energy_after.pickle", "wb")
    pickle.dump(Size_loc_energy2, f3)
    f3.close()

    f4 = open("temp_force_after.pickle", "wb")
    pickle.dump(Size_loc_force2, f4)
    f4.close()

            
def plot_EFscatter (Size_loc_data, ylabel="", labels=[None], note="", ylim=None):
    fig = plt.figure()
    alpha_list = np.linspace(0.02, 1, num=len(Size_loc_data))
    color_list = np.linspace(0.01, 1, num=len(Size_loc_data))
    cmap = cm.get_cmap("jet")
    for i in range(len(Size_loc_data)):
        size = sorted(Size_loc_data.keys())[i]
        loc_data = Size_loc_data[size]
        X, Y = [], []
        for loc in loc_data:
            data = loc_data[loc]
            X += [loc - NCP_len/2 for k in range(len(data))]
            Y += data
        plt.plot(X, Y, '.', markersize=1, label=labels[i], color=cmap(color_list[i]), alpha=alpha_list[i])    
    plt.axvline(x=0, linestyle='--', color='k')
    plt.axhline(y=0, linestyle='--', color='k')
    Lticks = [-60, -40, -20, -10]
    Rticks = [60, 40, 20, 10]
    for tick in Lticks:
        plt.axvline(x=tick, linestyle='--', color='r', alpha=0.2)
    for tick in Rticks:
        plt.axvline(x=tick, linestyle='--', color='b', alpha=0.2)
    #plt.title(titles[i])
    plt.xlabel("Position w.r.t Dyad (bp)")
    plt.ylabel(ylabel)
    plt.legend()
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig('EFscatter' + note + '.png', bbox_inches='tight')
    plt.close()

def plot_EFmean (Size_loc_data, ylabel="", labels=[None], note="", ylim=None):
    fig = plt.figure()
    color_list = np.linspace(0.01, 1, num=len(Size_loc_data))
    cmap = cm.get_cmap("jet")
    for i in range(len(Size_loc_data)):
        size = sorted(Size_loc_data.keys())[i]
        loc_data = Size_loc_data[size]
        X, Y = [], []
        for loc in loc_data:
            data = loc_data[loc]
            X.append(loc - NCP_len/2)
            Y.append(np.mean(data))
        #plt.plot(X, Y, '.', markersize=5, label=labels[i], color=cmap(color_list[i]))
        plt.plot(X, Y,label=labels[i], color=cmap(color_list[i]))
    plt.axvline(x=0, linestyle='--', color='k')
    plt.axhline(y=0, linestyle='--', color='k')
    Lticks = [-60, -40, -20, -10]
    Rticks = [60, 40, 20, 10]
    for tick in Lticks:
        plt.axvline(x=tick, linestyle='--', color='r', alpha=0.2)
    for tick in Rticks:
        plt.axvline(x=tick, linestyle='--', color='b', alpha=0.2)
    #plt.title(titles[i])
    plt.xlabel("Position w.r.t Dyad (bp)")
    plt.ylabel(ylabel)
    plt.legend()
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig('EFmean' + note + '.png', bbox_inches='tight')
    plt.close()
    
plot_EFmean(Size_loc_energy1, labels=[str(size) + " bp" for size in sorted(Size_loc_energy1.keys())], ylabel="Energy (A.U.)", note="energy_before")
plot_EFscatter(Size_loc_energy1, labels=[str(size) + " bp" for size in sorted(Size_loc_energy1.keys())], ylabel="Energy (A.U.)", note="energy_before")
plot_EFmean(Size_loc_energy2, labels=[str(size) + " bp" for size in sorted(Size_loc_energy2.keys())], ylabel="Energy (A.U.)", note="energy_after")
plot_EFscatter(Size_loc_energy2, labels=[str(size) + " bp" for size in sorted(Size_loc_energy2.keys())], ylabel="Energy (A.U.)", note="energy_after")


plot_EFmean(Size_loc_force1, labels=[str(size) + " bp" for size in sorted(Size_loc_force1.keys())], ylabel="Force (A.U.)", note="force_before", ylim=[-0.1, 0.1])
plot_EFscatter(Size_loc_force1, labels=[str(size) + " bp" for size in sorted(Size_loc_force1.keys())], ylabel="Force (A.U.)", note="force_before")
plot_EFmean(Size_loc_force2, labels=[str(size) + " bp" for size in sorted(Size_loc_force2.keys())], ylabel="Force (A.U.)", note="force_after", ylim=[-0.1, 0.1])
plot_EFscatter(Size_loc_force2, labels=[str(size) + " bp" for size in sorted(Size_loc_force2.keys())], ylabel="Force (A.U.)", note="force_after")
