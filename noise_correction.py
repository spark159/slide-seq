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

ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
NCP_len = 147
mtype_choice = ["M"]
size_st, size_ed = 1, 5

# pickle initialization
try:
    f1 = open("temp1.pickle", "rb")
    f2 = open("temp2.pickle", "rb")
    f3 = open("temp3.pickle", "rb")
    f4 = open("temp4.pickle", "rb")
    f5 = open("temp5.pickle", "rb")
    f6 = open("temp6.pickle", "rb")
    f7 = open("temp7.pickle", "rb")
    f8 = open("temp8.pickle", "rb")


    #f1 = open("/home/spark159/scripts/slide-seq/pickle_files/mmlib_control_before.pickle", "rb")
    #f2 = open("/home/spark159/scripts/slide-seq/pickle_files/mmlib_bubble_before.pickle", "rb")
    #f3 = open("/home/spark159/scripts/slide-seq/pickle_files/mmlib_control_after.pickle", "rb")
    #f4 = open("/home/spark159/scripts/slide-seq/pickle_files/mmlib_bubble_after.pickle", "rb")

    #f1 = open("/home/spark159/scripts/slide-seq/pickle_files/Insertion_control_before.pickle", "rb")
    #f2 = open("/home/spark159/scripts/slide-seq/pickle_files/Insertion_bubble_before.pickle", "rb")
    #f3 = open("/home/spark159/scripts/slide-seq/pickle_files/Insertion_control_after.pickle", "rb")
    #f4 = open("/home/spark159/scripts/slide-seq/pickle_files/Insertion_bubble_after.pickle", "rb")

    key_slider1 = pickle.load(f1)
    key_slider2 = pickle.load(f2)
    key_slider3 = pickle.load(f3)
    key_slider4 = pickle.load(f4)
    key_slider5 = pickle.load(f5)
    key_slider6 = pickle.load(f6)
    key_slider7 = pickle.load(f7)
    key_slider8 = pickle.load(f8)


except:
    fnames1 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_control_0_1rep_.combined.sort"]
    fnames2 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_control_0_2rep_.combined.sort"]

    fnames3 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_control_5_1rep_.combined.sort"]
    fnames4 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_control_5_2rep_.combined.sort"]

    fnames5 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_bubble_0_1rep_.combined.sort"]
    fnames6 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_bubble_0_2rep_.combined.sort"]

    fnames7 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_bubble_5_1rep_.combined.sort"]
    fnames8 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_bubble_5_2rep_.combined.sort"]


    #fnames1 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_control_0_1rep_.combined.sort", "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_control_0_2rep_.combined.sort"]
    #fnames2 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_bubble_0_1rep_.combined.sort", "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_bubble_0_2rep_.combined.sort"]
    #fnames3 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_control_5_1rep_.combined.sort", "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_control_5_2rep_.combined.sort"]
    #fnames4 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_bubble_5_1rep_.combined.sort", "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/mmlib_bubble_5_2rep_.combined.sort"]

    #fnames1 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/IDlib_control_0_1rep_.combined.sort", "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/IDlib_control_0_2rep_.combined.sort"]
    #fnames2 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/IDlib_bubble_0_1rep_.combined.sort", "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/IDlib_bubble_0_2rep_.combined.sort"]
    #fnames3 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/IDlib_control_5_1rep_.combined.sort", "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/IDlib_control_5_2rep_.combined.sort"]
    #fnames4 = ["/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/IDlib_bubble_5_1rep_.combined.sort", "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/IDlib_bubble_5_2rep_.combined.sort"]

    key_slider1 = load.load_files(fnames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', mtype_choice=mtype_choice)
    key_slider2 = load.load_files(fnames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', mtype_choice=mtype_choice)
    key_slider3 = load.load_files(fnames3, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', mtype_choice=mtype_choice)
    key_slider4 = load.load_files(fnames4, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', mtype_choice=mtype_choice)
    key_slider5 = load.load_files(fnames5, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', mtype_choice=mtype_choice)
    key_slider6 = load.load_files(fnames6, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', mtype_choice=mtype_choice)
    key_slider7 = load.load_files(fnames7, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', mtype_choice=mtype_choice)
    key_slider8 = load.load_files(fnames8, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', mtype_choice=mtype_choice)    


    f1 = open("temp1.pickle", "wb")
    pickle.dump(key_slider1, f1)
    f1.close()
    f2 = open("temp2.pickle", "wb")
    pickle.dump(key_slider2, f2)
    f2.close()
    f3 = open("temp3.pickle", "wb")
    pickle.dump(key_slider3, f3)
    f3.close()
    f4 = open("temp4.pickle", "wb")
    pickle.dump(key_slider4, f4)
    f4.close()
    f5 = open("temp5.pickle", "wb")
    pickle.dump(key_slider5, f5)
    f5.close()
    f6 = open("temp6.pickle", "wb")
    pickle.dump(key_slider6, f6)
    f6.close()
    f7 = open("temp7.pickle", "wb")
    pickle.dump(key_slider7, f7)
    f7.close()
    f8 = open("temp8.pickle", "wb")
    pickle.dump(key_slider8, f8)
    f8.close()


# define keys
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
keys = []
for key in list(set(key_slider1.keys()) & set(key_slider2.keys()) & set(key_slider3.keys()) & set(key_slider4.keys())):
    loc, mtype, nts = key.split('-')
    if mtype not in mtype_choice:
        continue
    if len(nts) < size_st:
        continue
    if len(nts) > size_ed:
        continue
    keys.append(key)
keys = sorted(keys, cmp=key_cmp)


# noise subtraction
key_slider_list = [key_slider3, key_slider4]
top_frac_list = [0.02, 0.06]
bott_frac_list = [0.04, 0.085]
for key_slider, top_frac, bott_frac in zip(key_slider_list, top_frac_list, bott_frac_list):
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
#graph.plot_signal(key_slider1, note='mmlib_control_0_1rep')
#graph.plot_signal(key_slider2, note='mmlib_control_0_2rep')
#graph.plot_signal(key_slider3, note='mmlib_control_5_1rep')
#graph.plot_signal(key_slider4, note='mmlib_control_5_2rep')
#graph.plot_signal(key_slider5, note='mmlib_bubble_0_1rep')
#graph.plot_signal(key_slider6, note='mmlib_bubble_0_2rep')
#graph.plot_signal(key_slider7, note='mmlib_bubble_5_1rep')
#graph.plot_signal(key_slider8, note='mmlib_bubble_5_2rep')

sample_mode = "polyA:" + "-".join([str(size) for size in [3]])
sample_list = sample.sampling(key_slider1, sample_mode)
#graph.plot_map(key_slider1, sample_list, norm_choice=True, draw_key=True, note='_mmlib_control_0_1rep')
#graph.plot_map(key_slider2, sample_list, norm_choice=True, draw_key=True, note='_mmlib_control_0_2rep')
#graph.plot_map(key_slider3, sample_list, norm_choice=True, draw_key=True, note='_mmlib_control_5_1rep')
#graph.plot_map(key_slider4, sample_list, norm_choice=True, draw_key=True, note='_mmlib_control_5_2rep')
#graph.plot_map(key_slider5, sample_list, norm_choice=True, draw_key=True, note='_mmlib_bubble_0_1rep')
#graph.plot_map(key_slider6, sample_list, norm_choice=True, draw_key=True, note='_mmlib_bubble_0_2rep')
#graph.plot_map(key_slider7, sample_list, norm_choice=True, draw_key=True, note='_mmlib_bubble_5_1rep')
#graph.plot_map(key_slider8, sample_list, norm_choice=True, draw_key=True, note='_mmlib_bubble_5_2rep')

graph_edit.plot_map(key_slider1, sample_list, True, Slider.peak_signal, draw = "key", slicing=0, note='_mmlib_control_0_1rep_KDE')



"""
# draw data on pdf
fname_list = ["mmlib_control_0_1rep", "mmlib_control_0_2rep", "mmlib_control_5_1rep", "mmlib_control_5_2rep", "mmlib_bubble_0_1rep", "mmlib_bubble_0_2rep", "mmlib_bubble_5_1rep", "mmlib_bubble_5_2rep"]
key_slider_list = [key_slider1, key_slider2, key_slider3, key_slider4, key_slider5, key_slider6, key_slider7, key_slider8]
m, n = 5, 2 # number of rows and colums
page_nums = int(math.ceil(len(keys)/float(m*n))) # number of pages
#page_nums = 3

for key_slider, fname in zip(key_slider_list, fname_list):
    pdf = matplotlib.backends.backend_pdf.PdfPages(fname + ".pdf")
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
