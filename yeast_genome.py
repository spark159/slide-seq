# Yeast genome data, Saccharomyces cerevisiae, UCSC, sacCer1 and sacCer2
# Nucleosome data, A map of nucleosome positions in yeast at base-pair resolution, Brogaard et al, Nature 2012.
# TSS data, Bidirectional promoters generate pervasive transcription in yeast, Xu et al, Nature 2009.x

import sys, os
import random, math
from argparse import ArgumentParser, FileType
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import load
from SliderClass import Slider
import HMM
import analysis
import graph
from termcolor import colored
import pickle
from pyliftover import LiftOver

#lo = LiftOver('/home/spark159/../../media/spark159/sw/plusonelibFinal/sacCer1ToSacCer2.over.chain')
lo = LiftOver("sacCer1", "sacCer2")

def change (chr, pos):
    try:
        new_chr, new_pos, _, _ = lo.convert_coordinate(chr, pos)[0]
        return new_chr, new_pos
    except:
        return None

roman_num = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16}
num_roman = {value:key for key, value in roman_num.items()}

def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq:
        nt = nt.upper()
        if nt == 'A':
            new_seq += 'T'
        if nt == 'T':
            new_seq += 'A'
        if nt == 'C':
            new_seq += 'G'
        if nt == 'G':
            new_seq += 'C'
    return new_seq

def tuple_cmp (a, b):
    if a[0] < b[0]:
        return -1
    elif a[0] > b[0]:
        return 1
    else:
        return 0

def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

def Amer_len(seq, pos=False):
    num = []
    num_pos = {}
    i = 0
    while i < len(seq):
        if seq[i] in 'AT':
            nt = seq[i]
            count = 1
            j = i + 1
            while j < len(seq):
                if seq[j] != nt:
                    break
                count +=1
                j +=1
            num.append(count)
            if count not in num_pos:
                num_pos[count] = []
            num_pos[count].append(i)
            i = j
        else:
            i +=1
    if pos:
        return num_pos
    if len(num) == 0:
        return 0
    return max(num)

# read sequence from genome
def get_seq(ref_fname, chr, st_ID, win):
    seq = ""
    pt = -1
    k = 0
    left = []
    Find = False
    stID = [[st, st_ID[st]] for st in sorted(st_ID.keys())]
    pos, ID = stID[k]
    seq_list = []
    ID_seq = {}
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith(">"):
            if Find:
                break
            if line[1:] == chr:
                Find = True
            continue
        if Find:
            if len(left) == 0 and pt + len(line) < pos:
                pt += len(line)
                continue
            for i in range(len(left)):
                leftID, seq = left.pop(0)
                ed = min(len(line), win-len(seq))
                seq += line[:ed]
                if len(seq) == win:
                    ID_seq[leftID] = seq
                else:
                    left.append([leftID, seq])
            while pt + len(line) >= pos and k < len(stID):
                loc = pos - pt - 1
                seq = line[loc:min(loc+win,len(line))]
                if len(seq) == win:
                    ID_seq[ID] = seq
                else:
                    left.append([ID, seq])
                k += 1
                try:
                    pos, ID = stID[k]
                except:
                    None
            if len(left) == 0 and len(ID_seq) == len(stID):
                break
            pt += len(line)
    while len(left) > 0:
        leftID, seq = left.pop(0)
        ID_seq[leftID] = seq
    assert len(ID_seq) == len(stID)
    return ID_seq

# read out library information
def read_library(fname):
    id_info = {}
    id = 0
    First = True
    for line in open(fname):
        if First:
            key_list = line.strip().split()
            First = False
            continue
        cols = line.strip().split()
        info = {}
        for i in range(len(key_list)):
            key = key_list[i]
            if key == 'pos':
                value = int(cols[i])
            elif key in ['score', 'snr']:
                value = float(cols[i])
            else:
                value = cols[i]
            info[key] = value
        id_info[id] = info
        id +=1
    return id_info
id_info = read_library("/home/spark159/../../media/spark159/sw/plusonelibFinal/plusonelib_table.txt")

# read out yeast genome reference (sacCer2)
def read_genome(fname):
    chr_dic = {}
    chr_name, sequence = "", ""
    for line in open(fname):
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
    return chr_dic
#chr_seq = read_genome("/home/spark159/../../media/spark159/sw/plusonelibFinal/sacCer2.fa")

# read Widom yeast nucleosome map (mapped on sacCer2)
def read_Nucmap(fname):
    chr_pos_score_snr = {}
    for line in open(fname):
        cols = line.strip().split()
        chr, pos, score, snr = cols
        pos, score, snr = int(pos) - 1, float(score), float(snr)
        if chr not in chr_pos_score_snr:
            chr_pos_score_snr[chr] = {}
        if pos not in chr_pos_score_snr[chr]:
            chr_pos_score_snr[chr][pos] = {}
        chr_pos_score_snr[chr][pos] = (score, snr)
    return chr_pos_score_snr
chr_pos_score_snr = read_Nucmap("/home/spark159/../../media/spark159/sw/plusonelibFinal/widom_yeast_map.txt")

# read Steinmetz TSS table (mapped on sacCer1)
def read_TSS(fname):
    gname_TSS = {}
    First = True
    for line in open(fname):
        if First:
            key_list = line.strip().split()[1:]
            First = False
            continue
        cols = line.strip().split()
        ID, chrnum, strand, txStart, txEnd, type, gname = cols[:7]
        if type != "ORF-T":
            continue
        chrnum = "chr" + chrnum
        if strand == '+':
            txStart, txEnd = int(txStart)-1, int(txEnd)-1
        else:
            txStart, txEnd = int(txEnd)-1, int(txStart)-1
        chr, txStart = change(chrnum, txStart)
        chr, txEnd = change(chrnum, txEnd)
        TSS = {"chrom":chr, "strand":strand, "txStart":txStart, "txEnd":txEnd}
        while gname in gname_TSS:
            gname += "-alt"
        gname_TSS[gname] = TSS
    return gname_TSS
gname_TSS = read_TSS("/home/spark159/../../media/spark159/sw/plusonelibFinal/TSS_anot.csv")

# read SGD gene hgTable
def read_hgTable(fname):
    gname_info = {}
    First = True
    for line in open(fname):
        if First:
            key_list = line.strip().split()[2:6]
            First = False
            continue
        cols = line.strip().split()
        bin, name, cols = cols[0], cols[1], cols[2:6]
        info = {}
        for i in range(len(key_list)):
            key = key_list[i]
            try:
                info[key] = int(cols[i])
            except:
                info[key] = cols[i]
        gname_info[name] = info
    return gname_info
gname_info = read_hgTable("/home/spark159/../../media/spark159/sw/plusonelibFinal/hgTable_SGD.txt")

# read RNA-seq data
def read_FPKM (fname):
    gname_FPKM = {}
    First = True
    for line in open(fname):
        cols = line.strip().split()
        if First:
            First = False
            continue
        _, gname, _, FPKM = cols[:4]
        FPKM = float(FPKM)
        if gname == '-':
            continue
        assert gname not in gname_FPKM
        gname_FPKM[gname] = FPKM
    return gname_FPKM
gname_FPKM = read_FPKM("/home/spark159/../../media/spark159/sw/plusonelibFinal/FPKM_NR_CR.txt")

# sanity check: all library data in widom data
for id in id_info:
    info = id_info[id]
    chr, pos, score = info['chr'], info['pos'], info['score']
    try:
        w_score, w_snr = chr_pos_score_snr[chr][pos]
        assert score == w_score
    except:
        raise ValueError

# find the +1 nucleosomes
plusone_chr_pos_info = {}
for gname in gname_TSS:
    TSS_info = gname_TSS[gname]
    chr, strand, txStart, txEnd = TSS_info['chrom'], TSS_info['strand'], TSS_info['txStart'], TSS_info['txEnd']
    try:
        pos_list = sorted(chr_pos_score_snr[chr].keys())
    except:
        continue
    if strand == '+':
        idx = np.searchsorted(pos_list, txStart, side='left')
    elif strand == '-':
        idx = np.searchsorted(pos_list, txStart, side='right') - 1
    if idx < 0 or idx >= len(pos_list):
        continue
    if strand == '+' and not (pos_list[idx] >= txStart and pos_list[idx] <= txEnd):
        continue
    if strand == '-' and not (pos_list[idx] <= txStart and pos_list[idx] >= txEnd):
        continue
    if strand == '+':
        assert txStart <= pos_list[idx]
        if idx > 0:
            assert txStart > pos_list[idx-1]
    else:
        assert txStart >= pos_list[idx]
        if idx < len(pos_list) - 1:
            assert txStart < pos_list[idx+1]     
    pos = pos_list[idx]
    info = {'gname':gname}
    if chr not in plusone_chr_pos_info:
        plusone_chr_pos_info[chr] = {}
    plusone_chr_pos_info[chr][pos] = info

# sanity check: all library data are +1 nucleosomes
out_count = 0
for id in id_info:
    info = id_info[id]
    chr, pos, score = info['chr'], info['pos'], info['score']
    try:
        gname = plusone_chr_pos_info[chr][pos]['gname']
        assert id_info[id]['strand'] == gname_TSS[gname]['strand'] 
        id_info[id]['gname'] = gname
    except:
        #print id
        id_info[id]['gname'] = np.nan
        out_count +=1

print out_count
        
# assign FPKM value to library
out_count = 0
for id in id_info:
    gname = id_info[id]['gname']
    try:
        FPKM = gname_FPKM[gname]
        id_info[id]['FPKM'] = FPKM
    except:
        id_info[id]['FPKM'] = np.nan
        out_count +=1

print out_count    

X, Y = [], []
for id in id_info:
    score = id_info[id]['score']
    FPKM = id_info[id]['FPKM']
    X.append(score)
    Y.append(FPKM)

with open("plusonelib_info.pickle", "wb") as f:
    pickle.dump(id_info, f)


"""
fig = plt.figure()
plt.plot(X, Y, '.')
plt.show()
plt.close()
"""
