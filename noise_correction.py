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
from scipy import signal
import matplotlib.backends.backend_pdf
import analysis
import random

def key_cmp (key1, key2):
    loc1, mtype1, nts1 = key1.split('-')
    loc2, mtype2, nts2 = key2.split('-')
    if len(nts1) < len(nts2):
        return -1
    elif len(nts1) > len(nts2):
        return 1
    else:
        if int(loc1) < int(loc2):
            return -1
        elif int(loc1) > int(loc2):
            return 1
        else:
            return 0

ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
NCP_len = 147

# load data
name_key_slider = {}

# 601 control data
fnames1 = ["/home/spark159/../../media/spark159/sw/polyAlibFinal/601_before_.combined.sort"]
fnames2 = ["/home/spark159/../../media/spark159/sw/polyAlibFinal/601_after_.combined.sort"]
Control1 = load.load_files(fnames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)
Control2 = load.load_files(fnames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)
name_key_slider['Control1'] = Control1
name_key_slider['Control2'] = Control2

# PolyA library
for condition in []:
    if condition == 'old':
        path = "/home/spark159/../../media/spark159/sw/AscanlibFinal/"
    elif condition == 'new':
        #path = "/home/spark159/../../media/spark159/sw/all_slide_seq_data/"
        path = "/home/spark159/../../media/spark159/sw/AscanlibFinal/"
    for time in [0, 5]:
        fname = "%slib_%s_%s_%srep" % ('polyA', condition, time, 1)
        if condition =='old' and time == 0:
            sort_fname = "Ascan0_S1_L001_R.oldsort"
        elif condition == 'old' and time == 5:
            sort_fname = "Ascan-5min_S1_L001_R.oldsort"
        elif condition == 'new' and time == 0:
            sort_fname = "Ascan0_S1_L001_R.combined.sort"
        elif condition == 'new' and time == 5:
            sort_fname = "Ascan-5min_S1_L001_R.combined.sort"
        try:
            with open(fname + ".pickle", "rb") as f:
                key_slider = pickle.load(f)
        except:
            key_slider = load.load_files([path +  sort_fname], ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear')
            with open(fname + ".pickle", "wb") as f:
                pickle.dump(key_slider, f)
        assert fname not in name_key_slider
        name_key_slider[fname] = key_slider

# Mismatch/Indel library
path = "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/"
for mtype in []:
    if mtype == 'M':
        library_type = 'mm'
    elif mtype in ['I', 'D']:
        library_type = 'ID'
    for condition in ['bubble']:
        for time in [0, 5]:
            for rep in [1, 2]:
                fname = "%slib_%s_%s_%srep" % (library_type, condition, time, rep)
                try:
                    with open(fname + ".pickle", "rb") as f:
                        key_slider = pickle.load(f)
                except:
                    key_slider = load.load_files([path + fname + "_.combined.sort"], ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', mtype_choice=[mtype])
                    with open(fname + ".pickle", "wb") as f:
                        pickle.dump(key_slider, f)
                assert fname not in name_key_slider
                name_key_slider[fname] = key_slider

# Plusone library
path = "/home/spark159/../../media/spark159/sw/all_slide_seq_data/"
for condition in ['new:corrected']:
    for time in [0, 30]:
        fname = "%slib_%s_%s" % ('plusone', condition, time)
        if condition =='old' and time == 0:
            sort_fname = "plusone-0_S1_L001_R.sort"
        elif condition == 'old' and time == 30:
            sort_fname = "plusone-30_S2_L001_R.sort"
        elif condition == 'new' and time == 0:
            sort_fname = "Plslib-HS_S1_L001_R.sort"
        elif condition == 'new' and time == 30:
            sort_fname = "Plslib-HS-30min_S2_L001_R.sort"
        try:
            with open(fname + ".pickle", "rb") as f:
                key_slider = pickle.load(f)
        except:
            key_slider = load.load_files([path +  sort_fname], ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None, load_ref="/home/spark159/../../media/spark159/sw/plusonelibFinal/plusonelib.ref")
            with open(fname + ".pickle", "wb") as f:
                pickle.dump(key_slider, f)
        assert fname not in name_key_slider
        name_key_slider[fname] = key_slider

# key selection
name_keys = {}
for name in name_key_slider:
    assert name not in name_keys
    key_slider = name_key_slider[name]
    keys = []
    if name.startswith("mm"):
        for key in key_slider:
            loc, mtype, nts = key.split('-')
            if len(nts) < 1:
                continue
            if len(nts) > 5:
                continue
            keys.append(key)
    elif name.startswith('polyAlib_old'):
        for key in key_slider:
            loc, mtype, nts = key.split('-')
            if len(nts) < 1:
                continue
            if len(nts) > 15:
                continue
            keys.append(key)
    else:
        keys = key_slider.keys()
    if not name.startswith("plusone"):
        keys = sorted(keys, cmp=key_cmp)
    name_keys[name] = keys

# noise subtraction
name_strand_threshold = {"Control2":{"top":0.05, "bott":0.02}, "mmlib_control_5_1rep":{"top":0.02, "bott":0.06}, "mmlib_control_5_2rep":{"top":0.04, "bott":0.085}, "polyAlib_old_5_1rep":{"top":0.02, "bott":0.06}}
for name in name_key_slider:
    if name not in name_strand_threshold:
        continue
    key_slider = name_key_slider[name]
    top_frac = name_strand_threshold[name]['top']
    bott_frac = name_strand_threshold[name]['bott']
    for key in key_slider:
        top_cutmap = analysis.sub_background(key_slider[key].right_cutmap, frac=top_frac)
        bott_cutmap = analysis.sub_background(key_slider[key].left_cutmap, frac=bott_frac)
        dyad_map = []
        for i in range(ref_length):
            temp = 0
            try:
                temp += top_cutmap[i-dyad_offset]
            except:
                pass
            try:
                temp += bott_cutmap[i+dyad_offset]
            except:
                pass
            dyad_map.append(temp)
        key_slider[key].right_cutmap = top_cutmap
        key_slider[key].left_cutmap = bott_cutmap
        key_slider[key].dyadmap = dyad_map


# get average signal
#graph.plot_signal(Control1, note='601_before')
#graph.plot_signal(Control2, note='601_after')
#graph.plot_signal(key_slider1, note='mmlib_control_0_1rep')
#graph.plot_signal(key_slider2, note='mmlib_control_0_2rep')
#graph.plot_signal(key_slider3, note='mmlib_control_5_1rep')
#graph.plot_signal(key_slider4, note='mmlib_control_5_2rep')
#graph.plot_signal(key_slider5, note='mmlib_bubble_0_1rep')
#graph.plot_signal(key_slider6, note='mmlib_bubble_0_2rep')
#graph.plot_signal(key_slider7, note='mmlib_bubble_5_1rep')
#graph.plot_signal(key_slider8, note='mmlib_bubble_5_2rep')

#sample_mode = "polyA:" + "-".join([str(size) for size in [3]])
#sample_list = sample.sampling(key_slider1, sample_mode)
#graph.plot_map(key_slider1, sample_list, norm_choice=True, draw_key=True, note='_mmlib_control_0_1rep')
#graph.plot_map(key_slider2, sample_list, norm_choice=True, draw_key=True, note='_mmlib_control_0_2rep')
#graph.plot_map(key_slider3, sample_list, norm_choice=True, draw_key=True, note='_mmlib_control_5_1rep')
#graph.plot_map(key_slider4, sample_list, norm_choice=True, draw_key=True, note='_mmlib_control_5_2rep')
#graph.plot_map(key_slider5, sample_list, norm_choice=True, draw_key=True, note='_mmlib_bubble_0_1rep')
#graph.plot_map(key_slider6, sample_list, norm_choice=True, draw_key=True, note='_mmlib_bubble_0_2rep')
#graph.plot_map(key_slider7, sample_list, norm_choice=True, draw_key=True, note='_mmlib_bubble_5_1rep')
#graph.plot_map(key_slider8, sample_list, norm_choice=True, draw_key=True, note='_mmlib_bubble_5_2rep')

sample_mode = "r:30"
sample_list = sample.sampling(name_key_slider['plusonelib_new:corrected_0'], sample_mode)
sub_sample_list = []
for index in random.sample(range(30),5):
    print index
    sub_sample_list = [[sample_list[0][index]]]

for name, key_slider in name_key_slider.items():
    if name.startswith('Control'):
        continue
    #sample_mode = "polyA:" + "-".join([str(size) for size in range(3, 16)])
    graph_edit.plot_map(key_slider, sample_list, True, Slider.peak_signal, slicing=0, note='_' + name)
    #sample_list = [[choice] for choice in random.sample(sample_list[0], 10)]
    graph_edit.plot_signal (key_slider, sub_sample_list, True, Slider.get_dyadmap, mean_choice=False, slicing = 0, note = "_" + name)

"""
# draw data on pdf
m, n = 5, 2 # number of rows and colums
for name, key_slider in name_key_slider.items():
    print name
    print
    pdf = matplotlib.backends.backend_pdf.PdfPages(name + ".pdf")
    keys = name_keys[name]
    page_nums = int(math.ceil(len(keys)/float(m*n))) # number of pages
    #page_nums = 10
    for i in range(page_nums):
        fig = plt.figure(figsize=(15,20))
        #fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
        j = 0
        while j < min(m*n, len(keys)-m*n*i):
            key = keys[m*n*i + j]
            top_cutmap, bott_cutmap = key_slider[key].right_cutmap, key_slider[key].left_cutmap
            dyad_map = [value*0.5 for value in key_slider[key].dyadmap]
            plt.subplot(m, n, j+1)
            plt.plot(top_cutmap, '-', label='Top', alpha=0.5)
            plt.plot(bott_cutmap, '-', label='Bottom', alpha=0.5)
            plt.plot(dyad_map, label='Dyad')
            #plt.hist(dyad_map, 100)
            loc, mtype, nts = key.split('-')
            st = int(loc)
            ed = st+len(nts)
            plt.axvspan(st, ed-1, alpha=0.5, color='red')
            plt.title(key)
            plt.legend(loc='upper right', frameon=False)
            plt.xlim([0,225])
            j +=1
        pdf.savefig(fig)
        plt.close()
    pdf.close()
"""
