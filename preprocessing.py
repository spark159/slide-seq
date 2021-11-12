import sys
import code
from argparse import ArgumentParser, FileType
import numpy as np
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from SliderClass_final import Slider
import math
import sample
import analysis
import load
import random
import readLoopseq
from sklearn import linear_model
import pickle
import os, sys, subprocess, re
import matplotlib.backends.backend_pdf
import graph_final as graph
import analysis_final as analysis

def size_loc_cmp (a, b):
    loc1, mtype1, nts1 = a.split('-')
    loc2, mtype2, nts2 = b.split('-')
    loc1, loc2 = int(loc1), int(loc2)
    size1, size2 = len(nts1), len(nts2)
    if size1 < size2:
        return -1
    elif size1 == size2:
        if loc1 < loc2:
            return -1
        else:
            return 1
    else:
        return 1


# basic pre-processing of slide-seq data
names = []

# 601
inpath = "/home/spark159/../../media/spark159/sw/polyAlibFinal/"
outpath = "/home/spark159/../../media/spark159/sw/601analysis_20210714/"
ref_fname = "601"

for time in [5,0]:
    if time == 0:
        read_pair = ["601_before_1.fastq", "601_before_2.fastq"]
    elif time == 5:
        read_pair = ["601_after_1.fastq", "601_after_2.fastq"]
    name = "%s_%smin" % ("601", time)
    print name
    names.append(name)
    
    # sorting reads
    #sort_cmd = ['python', 'sort_final.py', inpath+read_pair[0], inpath+read_pair[1], outpath+ref_fname, '-o', outpath+name]
    #subprocess.call(sort_cmd)

    # make data file
    data_cmd = ['python', 'make_data_final.py', outpath+name+'.sort', '-o', outpath+name,  '--oe', 'data', '-x', inpath+ref_fname+'.ref']
    subprocess.call(data_cmd)
    data_cmd = ['python', 'make_data_final.py', outpath+name+'.sort', '-o', outpath+name,  '--oe', 'pickle', '-x', inpath+ref_fname+'.ref']
    subprocess.call(data_cmd)



# window 1nt padding off / multiHit off
# PolyAlib
inpath = "/home/spark159/../../media/spark159/sw/polyAlibFinal/"
outpath = "/home/spark159/../../media/spark159/sw/polyAlibanalysis_20210713/"
ref_fname = "polyAscanlib_reindexed"

for time in [5,0]:
    if time == 0:
        read_pair = ["Ascan0_S1_L001_R1_001.fastq", "Ascan0_S1_L001_R2_001.fastq"]
        read_fname = "Ascan0_S1_L001_R.combined.fastq.gz"
    elif time == 5:
        read_pair = ["Ascan-5min_S1_L001_R1_001.fastq", "Ascan-5min_S1_L001_R2_001.fastq"]
        read_fname = "Ascan-5min_S1_L001_R.combined.fastq"
    name = "%s_%smin" % ("polyAscanlib", time)
    print name
    names.append(name)
    
    # sorting reads
    #sort_cmd = ['python', 'sort_final.py', inpath+read_fname, outpath+ref_fname, '-o', outpath+name]
    #subprocess.call(sort_cmd)

    # make data file
    data_cmd = ['python', 'make_data_final.py', outpath+name+'.sort', '-o', outpath+name,  '--fill', 'linear', '--oe', 'data', '-x', outpath+ref_fname+'.ref']
    subprocess.call(data_cmd)
    data_cmd = ['python', 'make_data_final.py', outpath+name+'.sort', '-o', outpath+name,  '--fill', 'linear', '--oe', 'pickle', '-x', outpath+ref_fname+'.ref']
    subprocess.call(data_cmd)



# Mismatch and Indel library
inpath = "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/"
outpath = "/home/spark159/../../media/spark159/sw/mmlibIDlibanalysis_20210713/"

for mtype in ['M', 'ID']:
    if mtype == 'M':
        library_type = 'mm'
        ref_fname = "polyAMismatch"
        cutpad = 3
    elif mtype == 'ID':
        library_type = 'ID'
        ref_fname = "singleInDel"
        cutpad = 3
    for condition in ['bubble', 'control']:
        for time in [0, 5]:
            for rep in [1, 2]:
                name = "%slib_%s_%smin_%srep" % (library_type, condition, time, rep)
                print name
                names.append(name)
                fastq_name = "%slib_%s_%s_%srep" % (library_type, condition, time, rep)
                read_pair = [fastq_name + "_1_trim.fastq", fastq_name + "_2_trim.fastq"]
                
                # sorting reads
                #sort_cmd = ['python', 'sort_final.py', inpath+read_pair[0], inpath+read_pair[1], inpath+ref_fname, '-o', outpath+name, '--direct']
                #subprocess.call(sort_cmd)

                # make data file
                data_cmd = ['python', 'make_data_final.py', outpath+name+'.sort', '-o', outpath+name,  '--fill', 'linear', '--oe', 'data', '--pad', str(cutpad), '-x', inpath+ref_fname+'.ref']
                subprocess.call(data_cmd)
                data_cmd = ['python', 'make_data_final.py', outpath+name+'.sort', '-o', outpath+name,  '--fill', 'linear', '--oe', 'pickle', '--pad', str(cutpad), '-x', inpath+ref_fname+'.ref']
                subprocess.call(data_cmd)


    
"""
# Basic QC of slide-seq data
names = ["IDlib_control_0min_1rep", "IDlib_control_5min_1rep", "IDlib_control_5min_2rep", "IDlib_bubble_0min_1rep", "IDlib_bubble_5min_1rep"]
for name in names:
    if name.startswith("polyAscanlib"):
        path = "/home/spark159/../../media/spark159/sw/polyAlibanlaysis_20210701/"
    elif name.startswith("mm") or name.startswith("ID"):
        path = "/home/spark159/../../media/spark159/sw/mmlibIDlibanalysis_20210701/"

    # make heatmaps
    key_slider = pickle.load(open(path+name+'.pickle'))

    size_keys = {}
    for key in key_slider.keys():
        loc, mtype, nts = key.split('-')
        if mtype != 'I':
            continue
        loc = int(loc)
        size = len(nts)
        if size not in size_keys:
            size_keys[size] = []
        size_keys[size].append(key)

    sample_list = []
    for size in sorted(size_keys.keys()):
        sample_list.append(sorted(size_keys[size], cmp=size_loc_cmp))

    graph.plot_map(key_slider, sample_list, norm_choice=True, draw_key=True, note='_' + name)
    #graph_edit.plot_map(key_slider, sample_list, True, Slider.peak_signal, draw = "key", note='_' + name)
    
    # make PDF files
    m, n = 5, 2
    pdf = matplotlib.backends.backend_pdf.PdfPages(path+name+".pdf")
    keys = sorted(key_slider.keys(), cmp=size_loc_cmp)
    page_nums = int(math.ceil(len(keys)/float(m*n))) # number of pages
    #page_nums = 1
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
            #plt.plot(dyad_map, 'k-', label=key)
            loc, mtype, nts = key.split('-')
            st = int(loc)
            ed = st+len(nts)
            plt.axvspan(st, ed-1, alpha=0.5, color='red')
            if name.startswith("polyA"):
                plt.title('-'.join([loc, mtype, str(len(nts))]))
            else:
                plt.title(key)
            leg = plt.legend(loc='upper right', frameon=False)
            plt.xlim([0,225])
            j +=1
        pdf.savefig(fig)
        plt.close()
    pdf.close()
"""    

    
"""
# make data file

outfname = "/home/spark159/../../media/spark159/sw/mmlibIDlibanalysis_20210701/IDlib_bubble_0_1rep"
outfname = "/home/spark159/../../media/spark159/sw/polyAlibanlaysis_20210701/polyAscanlib_5_1rep"
data_cmd = ['python', 'make_data_final.py', outfname+'.sort', '-o', 'test',  '--fill', 'linear', '--oe', 'pickle', '--pad', str(0)]
subprocess.call(data_cmd)

# make heatmaps
key_slider = pickle.load(open('test.pickle'))

size_keys = {}
for key in key_slider.keys():
    loc, mtype, nts = key.split('-')
    loc = int(loc)
    size = len(nts)
    if size not in size_keys:
        size_keys[size] = []
    size_keys[size].append(key)

sample_list = []
for size in sorted(size_keys.keys()):
    sample_list.append(sorted(size_keys[size], cmp=size_loc_cmp))

#graph.plot_map(key_slider, sample_list, norm_choice=True, draw_key=True)

# make PDF files
#m, n = 10, 4 # number of rows and colums
m, n = 5, 2
pdf = matplotlib.backends.backend_pdf.PdfPages('test' + ".pdf")
keys = sorted(key_slider.keys(), cmp=size_loc_cmp)
#keys = sorted(size_keys[15], cmp=size_loc_cmp)
page_nums = int(math.ceil(len(keys)/float(m*n))) # number of pages
#page_nums = 1
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
        #plt.plot(dyad_map, 'k-', label=key)
        loc, mtype, nts = key.split('-')
        st = int(loc)
        ed = st+len(nts)
        plt.axvspan(st, ed-1, alpha=0.5, color='red')
        #plt.title(key)
        plt.title('-'.join([loc, mtype, str(len(nts))]))
        leg = plt.legend(loc='upper right', frameon=False)
        plt.xlim([0,225])
        j +=1
    pdf.savefig(fig)
    plt.close()
pdf.close()



# preprocessing pipelines
path = "/home/spark159/../../media/spark159/sw/polyAlibanlaysis_20210701/"

lib_input = {"polyAscanlib":
             {"ref":path + "polyAscanlib_reindexed",
              "before":(path+"Ascan0_S1_L001_R1_001.fastq", path+"Ascan0_S1_L001_R2_001.fastq"),
              "after":(path+"Ascan-5min_S1_L001_R1_001.fastq", path+"Ascan-5min_S1_L001_R2_001.fastq")}}
                  
for libname in ['polyAscanlib']:
    ref_fname = lib_input[libname]["ref"]
    for time in [0]:
        for rep in [1]:
            outfname = path + "%s_%s_%srep" % (libname, time, rep)
            print outfname
            # sorting read files
            #read_pair = lib_input[libname][time]
            #sort_cmd = ['python', 'sort_final.py', read_pair[0], read_pair[1], ref_fname, '-o', outfname]
            #subprocess.call(sort_cmd)

            # basic QC for sorting
            #qc_cmd = ['python', 'libqc.py', outfname+'.sort', '-x', ref_fname+'.ref', '-c', 'valid', '-o', outfname]
            #subprocess.call(qc_cmd)

            # make data file
            data_cmd = ['python', 'make_data_final.py', outfname+'.sort', '-x', ref_fname+'.ref', '--fill', 'linear', '--oe', 'pickle']
            subprocess.call(data_cmd)

            # make heatmaps
            key_slider = pickle.load(open(outfname+'.pickle'))

            size_keys = {}
            for key in key_slider.keys():
                loc, mtype, nts = key.split('-')
                loc = int(loc)
                size = len(nts)
                if size not in size_keys:
                    size_keys[size] = []
                size_keys[size].append(key)

            sample_list = []
            for size in sorted(size_keys.keys()):
                sample_list.append(sorted(size_keys[size], cmp=size_loc_cmp))

            graph.plot_map(key_slider, sample_list, norm_choice=True, draw_key=True)

            
            # make PDF files
            #m, n = 10, 4 # number of rows and colums
            m, n = 5, 2
            pdf = matplotlib.backends.backend_pdf.PdfPages(outfname + ".pdf")
            #keys = sorted(key_slider.keys(), cmp=size_loc_cmp)
            keys = sorted(size_keys[15], cmp=size_loc_cmp)
            page_nums = int(math.ceil(len(keys)/float(m*n))) # number of pages
            #page_nums = 1
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
                    #plt.plot(dyad_map, 'k-', label=key)
                    loc, mtype, nts = key.split('-')
                    st = int(loc)
                    ed = st+len(nts)
                    plt.axvspan(st, ed-1, alpha=0.5, color='red')
                    plt.title(key)
                    leg = plt.legend(loc='upper right', frameon=False)
                    plt.xlim([0,225])
                    j +=1
                pdf.savefig(fig)
                plt.close()
            pdf.close()
            

    
    

#filenames = ['polyAMismatch_simulated.sort']
#filenames = ["polyAscanlib_reindexed_simulated.sort"]
#filenames = ["singleInDel_simulated.sort"]
#key_slider = load.load_files(filenames, ref_length=225, dyad_axis=225/2, dyad_offset=53, filter_num = 0,  fill='linear', mtype_choice='I')


#fname = "polyAMismatch_simulated.pickle"
#fname = "polyAscanlib_reindexed_simulated.pickle"
fname = "singleInDel_simulated.pickle"
key_slider = pickle.load(open(fname))

key_list = []
for key in key_slider:
    if key != 'BACKBONE':
        key_list.append(key)

size_keys = {}
for key in key_list:
    loc, mtype, nts = key.split('-')
    loc = int(loc)
    size = len(nts)
    if size not in size_keys:
        size_keys[size] = []
    size_keys[size].append(key)

sample_list = []
for size in sorted(size_keys.keys()):
    sample_list.append(sorted(size_keys[size], cmp=size_loc_cmp))
        
#graph_edit.plot_map(key_slider, sample_list = [sample_list], norm_choice=True, obs_func = Slider.get_cutmap, draw = 'key')

graph.plot_map(key_slider, sample_list, norm_choice=True, draw_key=True)



#filenames = ['polyAscanlib_reindexed_simulated']
filenames = ["polyAMismatch_simulated"]
#filenames = ["singleInDel_simulated"]
for fname in filenames:
    # make data file
    data_cmd = ['python', 'make_data_final.py', fname+'.sort', '-o', fname,  '--fill', 'Naive', '--oe', 'pickle']
    subprocess.call(data_cmd)

    id_slider = pickle.load(open(fname+'.pickle'))

    target_ids = []
    for id in id_slider:
        loc, mtype, nts = id.split('-')

        #if mtype == 'I':
        #    target_ids.append(id)

        if len(nts) == 3:
            target_ids.append(id)

    size_ids = {}
    for id in target_ids:
        loc, mtype, nts = id.split('-')
        loc = int(loc)
        size = len(nts)
        if size not in size_ids:
            size_ids[size] = []
        size_ids[size].append(id)

    #sample_list = []
    #for size in sorted(size_ids.keys()):
    #    sample_list.append(sorted(size_ids[size], cmp=analysis.cmp_wid_st))

            
    graph.plot_map(id_slider, [Slider.get_top_cutmap, Slider.get_bottom_cutmap], ids=sorted(target_ids, cmp=analysis.wid_cmp_len), mark='wid', cmap=['Reds', 'Blues'], thickness=[8, 0, 2], save=True, note='_cut')
    graph.plot_map(id_slider, Slider.get_dyadmap, ids=sorted(target_ids, cmp=analysis.wid_cmp_len), mark='wid', cmap=['jet'], thickness=[8, 0, 2], save=True, note='_dyad')
    graph.plot_sig(id_slider, [Slider.get_top_cutmap, Slider.get_bottom_cutmap, Slider.get_dyadmap], ids=sorted(target_ids, cmp=analysis.wid_cmp_len), mark='wid', label = ['top', 'bottom', 'dyad'], alpha=[0.5, 0.5, 1], save=True, note='')

    #graph.plot_map(id_slider, [Slider.get_top_cutmap, Slider.get_bottom_cutmap], ids=sample_list, mark='wid', cmap=['Reds', 'Blues'], thickness=[8, 0, 2], save=True, note='_cut')
    #graph.plot_map(id_slider, Slider.get_dyadmap, ids=sample_list, mark='wid', cmap=['jet'], thickness=[8, 0, 2], save=True, note='_dyad')
    #graph.plot_sig(id_slider, [Slider.get_top_cutmap, Slider.get_bottom_cutmap, Slider.get_dyadmap], ids=sample_list, mark='wid', label = ['top', 'bottom', 'dyad'], alpha=[0.5, 0.5, 1], save=True, note='')
"""
