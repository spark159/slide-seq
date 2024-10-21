#!/usr/bin/env python

# Nucleosome_positioning by Sangwoo Park and Daehwan Kim, September 2016
# Calculate probability of nucleosome positioning on given DNA sequence
# Yeast genome data Based on Pubmed, Saccharomyces cerevisiae (UCSC -SAC2)
#    http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/
# Positioning data Based on Kristin, et al, Nature 2012
# Algorithm Based on Segal, et al, Nature 2006

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

def Nuc_positioning(nucleosome_dna_len,
        markov_order,
        background,
        verbose):

    def color_A (seq):
        text = ''
        for nt in seq:
            if nt == 'A':
                text += colored(nt, 'red')
            else:
                text += nt
        return text


        # get reverse complementary sequence
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

    # read out yeast genome reference
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

    # read out library reference
    def read_library(fname):
        id_seq = {}
        for line in open(fname):
            if line.startswith('>'):
                id = line[1:].strip()
            else:
                seq = line.strip()
                id_seq[id] = seq
        return id_seq

    # read out nucleosome positionin data
    def read_NCP_list(NCP_fname):
        NCP_list = []
        for line in open(NCP_fname):
            line = line.strip()
            try:
                chr, pos, score, noise= line.split()
                pos, score, noise = int(pos) - 1, float(score), float(noise)
                if noise <= 0.0:
                    continue
                NCP_list.append([chr, pos, "+",  score / noise])
                NCP_list.append([chr, pos, "-",  score / noise])
            except ValueError:
                continue
        return NCP_list

    # read out nucleosome positioning sequences
    def read_NCP_sequences(genome, NCP_list):
        def read_NCP_sequence(genome, NCP):
            chr, pos, strand, score = NCP
            st = pos - nucleosome_dna_len / 2; en = pos + nucleosome_dna_len / 2
            if st < 0 or en >= len(genome[chr]):
                return None
            seq = genome[chr][st:en+1]
            if strand == "+":
                return seq
            elif strand == "-":
                return rev_comp(seq)
        NCP_seq_list = []
        for NCP in NCP_list:
            seq = read_NCP_sequence(genome, NCP)
            if seq == None:
                continue
            NCP_seq_list.append(seq)
            #NCP_seq_list.append(rev_comp(seq))
        return NCP_seq_list

    # convert slide-seq data into nucleosome positioning data
    def Slider_to_NCP_list(key_slider, scale=1000, pick_frac=0.5, duplicate=True):
        NCP_list = []
        for key in key_slider:
            dyadmap = key_slider[key].dyadmap
            total = sum(dyadmap)
            countmap = [ int(round(float(value)*scale/total)) for value in dyadmap]
            count_pos = [(countmap[i], i) for i in range(nucleosome_dna_len/2,len(countmap)-nucleosome_dna_len/2)]
            count_pos = sorted(count_pos, cmp=tuple_cmp, reverse=True)[:int(round(len(count_pos)*pick_frac))]
            for count, pos in count_pos:
                if duplicate:
                    for i in range(count):
                        NCP_list.append([key, pos, "+", count])
                        NCP_list.append([key, pos, "-", count])
                else:
                    NCP_list.append([key, pos, "+", count])
                    NCP_list.append([key, pos, "-", count])
        return NCP_list


    # k-mer analysis
    def kmer_freq (NCP_seq_list, kmer_len, states='ATCG'):
        def all_path(N, states):
            if N==1:
                return list(states)
            output=[]
            for path in all_path(N-1, states):
                for state in states:
                    output.append(path+state)
            return output
        assert NCP_seq_list[0] >= kmer_len
        kmers = all_path(kmer_len, states)
        kmer_freq = {}
        for kmer in kmers:
            kmer_freq[kmer] = 0.0            
        for seq in NCP_seq_list:
            for i in range(len(seq)-kmer_len+1):
                kmer = seq[i:i+kmer_len]
                kmer_freq[kmer] += 1.0
        #print kmer_freq
        #print
        total = sum(kmer_freq.values())
        for kmer in kmer_freq:
            kmer_freq[kmer] = kmer_freq[kmer] / total
        return kmer_freq

    
    # dinucleotide step analysis (A vs. G)
    def AG_freq (NCP_seq_list):
        Afreq=np.zeros(nucleosome_dna_len - 1); Gfreq=np.zeros(nucleosome_dna_len - 1)
        for seq in NCP_seq_list:
            for i in range(len(seq)-1):
                dint = seq[i:i+2]
                if dint in ['AA','AT','TA','TT']:
                    Afreq[i] += 1.0
                elif dint in ['CC','CG','GC','GG']:
                    Gfreq[i] += 1.0
        return Afreq / len(NCP_seq_list), Gfreq / len(NCP_seq_list)

    # dinucleotide step analysis (full)
    def din_freq (NCP_seq_list):
        freq = {}
        for seq in NCP_seq_list:
            for i in range(len(seq)-1):
                dint = seq[i:i+2]
                if dint not in freq:
                    freq[dint] = np.zeros(nucleosome_dna_len)
                freq[dint][i] += 1.0 / len(NCP_seq_list)
        return freq

    # get transition matrix of NCP positining
    def get_NCP_m (NCP_seq_list):    
        freq=din_freq(NCP_seq_list, nucleosome_dna_len);
        result=[]; base=["A", "T", "C", "G"]
        for i in range(nucleosome_dna_len - 1):
            matrix=np.zeros([4,4])
            for j in range(4):
                norm=0.0; row=np.zeros(4)
                for k in range(4):
                    key=base[j]+base[k]
                    row[k]= freq[key][i]; norm += freq[key][i]
                matrix[j]=row/float(norm)
            result.append(matrix)
        return result

    # get NCP positioning probabiliy profile for a long sequence
    def NCPprob_profile (seq):
        F=get_forward(seq, len(seq)); R=get_reverse(seq, len(seq))
        profile=[]
        for i in range(1, len(seq)+1):
            if len(seq[i-1:]) < nucleosome_dna_len:
                prob=0.0;
            else:
                prob=(F[i-1] * get_NCP_prob (seq[i-1:i+nucleosome_dna_len-1]) * R[len(seq)-(i+nucleosome_dna_len-1)]) / R[-1]
            profile.append(prob)
        return profile

    # get NCP occupancy profile of given sequence
    def NCPoccupancy (profile):
        result=[]
        for i in range(len(profile)):
            prob=0.0
            for j in range(nucleosome_dna_len):
                pos=i-j
                if pos <0:
                    break
                else:
                    prob += profile[pos]
            result.append(prob)
        return result

    # read out yeast genome data
    #ygenome = read_genome("data/scerevisiae.fa")    
    #NCP_list = read_NCP_list("data/nature11142-s2.txt")
    #NCP_seq_list = read_NCP_sequences(ygenome, NCP_list)
    
    # read pluse-one library
    ref_length = 225
    dyad_axis = ref_length/2
    dyad_offset = 52
    filenames1 = ["../../Illumina/plusoneHS/data/Plslib-HS_S1_L001_R.sort"]
    filenames2 = ["../../Illumina/plusoneHS/data/Plslib-HS-30min_S2_L001_R.sort"]
    key_slider1 = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)
    key_slider2 = load.load_files(filenames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)
    id_seq = read_library("../../Illumina/plusoneHS/plusonelib.ref")

    # background
    #NCP_listb = []
    #for key in key_slider1:
    #    for u in range(nucleosome_dna_len/2, ref_length-nucleosome_dna_len/2):
    #        NCP_listb.append([key, u, "+", 1])
    #        NCP_listb.append([key, u, "-", 1])
    #NCP_seq_listb = read_NCP_sequences(id_seq, NCP_listb)
    #Afreqb, Gfreqb = AG_freq(NCP_seq_listb)
    #Dinfreqb = din_freq(NCP_seq_listb)
    #Afreqb = np.asarray(Afreqb)
    #Gfreqb = np.asarray(Gfreqb)
    #Kmerfreqb = kmer_freq (NCP_seq_listb, kmer_len=5, states='ATCG')
    #print Kmerfreqb

    #pickle.dump(Afreqb, open("Afreqb.p", "wb"))
    #pickle.dump(Gfreqb, open("Gfreqb.p", "wb"))
    #pickle.dump(Kmerfreqb, open("Kmerfreqb.p", "wb"))
    #pickle.dump(Dinfreqb, open("Dinfreqb.p", "wb"))

    Afreqb = pickle.load(open("Afreqb.p", "rb"))
    Gfreqb = pickle.load(open("Gfreqb.p", "rb"))
    Kmerfreqb = pickle.load(open("Kmerfreqb.p", "rb"))
    Dinfreqb = pickle.load(open("Dinfreqb.p", "rb"))
    

    x_axis = [ i - len(Afreqb)*0.5 + 0.5  for i in range(len(Afreqb))]
    fig = plt.figure()
    plt.plot(x_axis, Afreqb, label='AA/AT/TA/TT')
    plt.plot(x_axis, Gfreqb, label='CC/CG/GC/GG')
    #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
    plt.legend()
    #plt.show()
    plt.ylim([0,0.7])
    plt.xlim([-80,80])
    plt.ylabel("Frequency")
    plt.xlabel("Location from dyad (bp)")
    plt.savefig("Background_ATGC.png")
    plt.close()

    """
    # nucleosome positioning signal
    for pick_frac in [1]:
        pick_name = str(pick_frac * 100) + ' %'
        
        #NCP_list1 = Slider_to_NCP_list(key_slider1, pick_frac=pick_frac)
        #NCP_list2 = Slider_to_NCP_list(key_slider2, pick_frac=pick_frac)

        #NCP_seq_list1 = read_NCP_sequences(id_seq, NCP_list1)
        #NCP_seq_list2 = read_NCP_sequences(id_seq, NCP_list2)
        #Afreq1, Gfreq1 = AG_freq(NCP_seq_list1)
        #Afreq2, Gfreq2 = AG_freq(NCP_seq_list2)
        #Dinfreq1 = din_freq(NCP_seq_list1)
        #Dinfreq2 = din_freq(NCP_seq_list2)
        
        #Kmerfreq1 = kmer_freq (NCP_seq_list1, kmer_len=5, states='ATCG')
        #Kmerfreq2 = kmer_freq (NCP_seq_list2, kmer_len=5, states='ATCG')

        #pickle.dump(Afreq1, open("Afreq1.p", "wb"))
        #pickle.dump(Gfreq1, open("Gfreq1.p", "wb"))
        #pickle.dump(Kmerfreq1, open("Kmerfreq1.p", "wb"))
        #pickle.dump(Dinfreq1, open("Dinfreq1.p", "wb"))

        Afreq1 = pickle.load(open("Afreq1.p", "rb"))
        Gfreq1 = pickle.load(open("Gfreq1.p", "rb"))
        Kmerfreq1 = pickle.load(open("Kmerfreq1.p", "rb"))
        Dinfreq1 = pickle.load(open("Dinfreq1.p", "rb"))

        #pickle.dump(Afreq2, open("Afreq2.p", "wb"))
        #pickle.dump(Gfreq2, open("Gfreq2.p", "wb"))
        #pickle.dump(Kmerfreq2, open("Kmerfreq2.p", "wb"))
        #pickle.dump(Dinfreq2, open("Dinfreq2.p", "wb"))

        Afreq2 = pickle.load(open("Afreq2.p", "rb"))
        Gfreq2 = pickle.load(open("Gfreq2.p", "rb"))
        Kmerfreq2 = pickle.load(open("Kmerfreq2.p", "rb"))
        Dinfreq2 = pickle.load(open("Dinfreq2.p", "rb"))
        
        kmers = Kmerfreq1.keys()
        kmers_fold1 = {}
        kmers_fold2 = {}
        for kmer in kmers:
            kmers_fold1[kmer] = Kmerfreq1[kmer] / Kmerfreqb[kmer]
            kmers_fold2[kmer] = Kmerfreq2[kmer] / Kmerfreqb[kmer]

        temp1 = [(fold, kmer) for kmer, fold in kmers_fold1.items()]
        temp1 = sorted(temp1, reverse=True)
        fold1, fold2 = [], []
        kmers = []
        mark = np.zeros(2)
        clist= []
        for u in range(len(temp1)):
            fold, kmer = temp1[u]
            kmers.append(kmer)
            if kmer == 'AAAAA':
                mark[0] = u*20
            if kmer == 'GGGGG':
                mark[1] = u*20
            fold1.append(kmers_fold1[kmer])
            #fold2.append(kmers_fold2[kmer])
            #clist.append(GC_content(kmer))
            clist.append(GC_content(kmer))

        fig = plt.figure()
        #barWidth = 0.25
        r1 = [ i*20 for i in range(len(kmers))]
        #r2 = [x + barWidth for x in r1]
        #plt.bar(r1, freq1_fold, color='#7f6d5f', width=barWidth, edgecolor='white', label='Before')
        #plt.bar(r2, freq2_fold, color='#557f2d', width=barWidth, edgecolor='white', label='After')
        plt.scatter(r1, fold1 ,c=clist, cmap="jet", s=3)
        #plt.scatter(r1, fold2, label='After',c=clist, cmap="Greens", s=3)
        for k in range(len(mark)):
            plt.axvline(x=mark[k], color='r')
        #plt.xlabel('group', fontweight='bold')
        #plt.xticks([r + barWidth/2 for r in range(len(freq1_fold))], kmers)
        #plt.xticks(r1, kmers, rotation=90, fontsize=3)
        plt.xticks(mark, ['AAAAA', 'GGGGG'], rotation=45, fontsize=8)
        plt.ylabel('Fold change')
        plt.xlabel("5-mers")
        #plt.legend()
        plt.savefig("K-merfreq_Before_" + pick_name + ".png")
        #plt.show()
        plt.close()

        temp2 = [(fold, kmer) for kmer, fold in kmers_fold2.items()]
        temp2 = sorted(temp2, reverse=True)
        fold1, fold2 = [], []
        kmers = []
        mark = np.zeros(2)
        clist= []
        for u in range(len(temp2)):
            fold, kmer = temp2[u]
            kmers.append(kmer)
            if kmer == 'AAAAA':
                mark[0] = u*20
            if kmer == 'GGGGG':
                mark[1] = u*20
            #fold1.append(kmers_fold1[kmer])
            fold2.append(kmers_fold2[kmer])
            #clist.append(GC_content(kmer))
            clist.append(GC_content(kmer))


        fig = plt.figure()
        #barWidth = 0.25
        r1 = [ i*20 for i in range(len(kmers))]
        #r2 = [x + barWidth for x in r1]
        #plt.bar(r1, freq1_fold, color='#7f6d5f', width=barWidth, edgecolor='white', label='Before')
        #plt.bar(r2, freq2_fold, color='#557f2d', width=barWidth, edgecolor='white', label='After')
        #plt.scatter(r1, fold1, label='Before',c=clist, cmap="Blues", s=3)
        plt.scatter(r1, fold2, c=clist, cmap="jet", s=3)
        for k in range(len(mark)):
            plt.axvline(x=mark[k], color='r')
        #plt.xlabel('group', fontweight='bold')
        #plt.xticks([r + barWidth/2 for r in range(len(freq1_fold))], kmers)
        #plt.xticks(r1, kmers, rotation=90, fontsize=3)
        plt.xticks(mark, ['AAAAA', 'GGGGG'], rotation=45, fontsize=8)
        plt.ylabel('Fold change')
        plt.xlabel("5-mers")
        #plt.legend()
        plt.savefig("K-merfreq_After_" + pick_name + ".png")
        #plt.show()
        plt.close()

        
        fig = plt.figure()
        for kmer in kmers:
            x = kmers_fold1[kmer]
            y = kmers_fold2[kmer]
            plt.plot(x, y, '.')
        small = min(kmers_fold1.values() + kmers_fold2.values())
        large = max(kmers_fold1.values() + kmers_fold2.values())
        plt.plot([0,2], [0,2], '--')
        plt.xlim([small-0.1, large+0.1])
        plt.ylim([small-0.1, large+0.1])
        plt.xlabel("Before")
        plt.ylabel("After")
        plt.savefig("K-merfreq_BvsA_" + pick_name + ".png")
        #plt.show()    
        plt.close()

        
        kmerchange = []
        for kmer in kmers:
            change = 100 * (kmers_fold2[kmer] - kmers_fold1[kmer]) / kmers_fold1[kmer]
            kmerchange.append((change, kmer))
        kmerchange = sorted(kmerchange, cmp=tuple_cmp)

        changes = []
        mark = np.zeros(2)
        clist = []
        GC = []
        for i in range(len(kmerchange)):
            change, kmer = kmerchange[i]
            changes.append(change)
            clist.append(GC_content(kmer))
            if kmer == 'AAAAA':
                mark[0] = i*20
            if kmer == 'GGGGG':
                mark[1] = i*20
            GC.append(GC_content(kmer))
            print "%s\t%0.2f" % (color_A(kmer), change)

        fig = plt.figure()
        r1 = [ i*20 for i in range(len(changes))]
        plt.scatter(r1, changes, c=clist, cmap="jet", s=3)
        for k in range(len(mark)):
            plt.axvline(x=mark[k], color='r')
        plt.xticks(mark, ['AAAAA', 'GGGGG'], rotation=45, fontsize=8)
        plt.ylabel('Fold change (%)')
        plt.xlabel("5-mers")
        #plt.legend()
        plt.savefig("K-merfreq_change_" + pick_name + ".png")
        #plt.show()
        plt.close()

        fig = plt.figure()
        #r1 = [ i*20 for i in range(len(changes))]
        plt.scatter(GC, changes, s=3)
        #plt.axvline(x=x[0], color='r')
        #plt.xticks(x, ['AAAAA'], rotation=45, fontsize=5)
        plt.ylabel('Fold change (%)')
        plt.xlabel("GC content (%)")
        #plt.legend()
        plt.savefig("K-merfreq_GC_" + pick_name + ".png")
        #plt.show()
        plt.close()

        
        hmap1, hmap2 = [], []
        hmap3 = []
        din_list = ['AA', 'AT', 'TA', 'TT', 'GT', 'TG', 'AC', 'CA', 'AG', 'GA', 'TC', 'CT', 'GG', 'GC', 'CG', 'CC']
        for din in din_list:
            hmap1.append(np.asarray(Dinfreq1[din]) / np.asarray(Dinfreqb[din]))
            hmap2.append(np.asarray(Dinfreq2[din]) / np.asarray(Dinfreqb[din]))
            hmap3.append(100 * (Dinfreq2[din]-Dinfreq1[din]) / Dinfreq1[din])

        #x_axis = [ i - len(hmap1[0])*0.5 + 0.5  for i in range(len(hmap1[0]))]
        tick_loc = []
        tick_name = []
        for i in range(0, len(hmap1[0])/2, 10):
            if i == 0:
                tick_loc.append(i + len(hmap1[0])/2)
                tick_name.append(str(i))
            else:
                tick_loc.append(i + len(hmap1[0])/2)
                tick_loc.append(-i + len(hmap1[0])/2)
                tick_name.append(str(i))
                tick_name.append(str(-i))
        
        fig = plt.figure()
        plt.imshow(hmap1, aspect='auto', interpolation='none', cmap='Purples', vmin=0.8, vmax=1.2)
        plt.yticks(range(len(din_list)), din_list, fontsize=10)
        plt.xticks(tick_loc, tick_name)
        plt.xlabel("Location from dyad (bp)")
        plt.colorbar()
        plt.savefig("Dinfreq_before_" + pick_name + ".png")
        #plt.show()
        plt.close()

        fig = plt.figure()
        plt.imshow(hmap2, aspect='auto', interpolation='none', cmap='Purples', vmin=0.8, vmax=1.2)
        plt.yticks(range(len(din_list)), din_list, fontsize=10)
        plt.xticks(tick_loc, tick_name)
        #plt.xticks(range(0, len(hmap1[0])+3, 10), [ i - 150/2 for i in range(0, len(hmap1[0])+3, 10)])
        plt.xlabel("Location from dyad (bp)")
        plt.colorbar()
        plt.savefig("Dinfreq_after_" + pick_name + ".png")
        plt.colorbar()
        plt.close()
        
        fig = plt.figure()
        plt.imshow(hmap3, aspect='auto', interpolation='none', cmap='jet')
        plt.yticks(range(len(din_list)), din_list, fontsize=10)
        plt.xticks(tick_loc, tick_name)
        #plt.xticks(range(0, len(hmap1[0])+3, 10), [ i - 150/2 for i in range(0, len(hmap1[0])+3, 10)])
        plt.xlabel("Location from dyad (bp)")
        plt.colorbar()
        plt.savefig("Dinfreq_change_" + pick_name + ".png")
        #plt.show()
        plt.close()

        fig = plt.figure()
        plt.plot(x_axis, Afreq1, label='AA/AT/TA/TT')
        plt.plot(x_axis, Gfreq1, label='CC/CG/GC/GG')
        #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
        plt.legend()
        #plt.show()
        plt.ylim([0,0.7])
        plt.xlim([-80,80])
        plt.ylabel("Frequency")
        plt.xlabel("Location from dyad (bp)")
        plt.savefig("Before_ATGC_" + pick_name + ".png")
        plt.close()

        fig = plt.figure()
        Afreq1 = np.asarray(Afreq1)
        Gfreq1 = np.asarray(Gfreq1)    
        plt.plot(x_axis, Afreq1/Afreqb, label='AA/AT/TA/TT')
        plt.plot(x_axis, Gfreq1/Gfreqb, label='CC/CG/GC/GG')
        #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
        plt.legend()
        #plt.show()
        #plt.ylim([0,0.7])
        plt.xlim([-80,80])
        plt.ylabel("Fold change")
        plt.xlabel("Location from dyad (bp)")
        plt.savefig("Before_ATGC_fold_" + pick_name + ".png")
        plt.close()


        fig = plt.figure()
        plt.plot(x_axis, Afreq2, label='AA/AT/TA/TT')
        plt.plot(x_axis, Gfreq2, label='CC/CG/GC/GG')
        plt.legend()
        #plt.show()
        plt.ylim([0,0.7])
        plt.xlim([-80,80])
        plt.ylabel("Frequency")
        plt.xlabel("Location from dyad (bp)")
        plt.savefig("After_ATGC_" + pick_name + ".png")
        plt.close()

        fig = plt.figure()
        Afreq2 = np.asarray(Afreq2)
        Gfreq2 = np.asarray(Gfreq2)    
        plt.plot(x_axis, Afreq2/Afreqb, label='AA/AT/TA/TT')
        plt.plot(x_axis, Gfreq2/Gfreqb, label='CC/CG/GC/GG')
        #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
        plt.legend()
        #plt.show()
        #plt.ylim([0,0.7])
        plt.xlim([-80,80])
        plt.ylabel("Fold change")
        plt.xlabel("Location from dyad (bp)")
        plt.savefig("After_ATGC_fold_" + pick_name + ".png")
        plt.close()

        
        fig = plt.figure()
        plt.plot(x_axis, 100*(Afreq2-Afreq1)/Afreq1, label='AA/AT/TA/TT')
        plt.plot(x_axis, 100*(Gfreq2-Gfreq1)/Gfreq1, label='CC/CG/GC/GG')
        #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
        plt.legend()
        #plt.show()
        #plt.ylim([0,0.7])
        plt.xlim([-80,80])
        plt.ylabel("Fold change (%)")
        plt.xlabel("Location from dyad (bp)")
        plt.savefig("change_ATGC_fold_" + pick_name + ".png")
        plt.close()

        
    
    
    # Before vs After differntial positioning signal
    keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
    before_NCP_list = []
    after_NCP_list = []
    for key in keys:
        slider1 = key_slider1[key]
        slider2 = key_slider2[key]
        dyadmap1 = slider1.dyadmap
        total1 = sum(dyadmap1)
        dyadmap1 = [value / total1 for value in dyadmap1]
        dyadmap2 = slider2.dyadmap
        total2 = sum(dyadmap2)
        dyadmap2 = [value / total2 for value in dyadmap2]
        before_ov_after = []
        after_ov_before = [] 
        for i in range(len(dyadmap1)):
            sig1 = dyadmap1[i]
            sig2 = dyadmap2[i]
            #if sig1 <= 0:
            #    sig1 = 1.0
            #if sig2 <= 0:
            #    sig2 = 1.0
            #ratio1 = sig1/sig2
            ratio1 = sig1 - sig2
            before_ov_after.append(ratio1)
            #ratio2 = sig2/sig1
            ratio2 = sig2 - sig1
            after_ov_before.append(ratio2)
        #print before_ov_after
        #print
        peaksig1 = analysis.find_peaks(before_ov_after, num=50)
        #print peaksig1
        peaksig2 = analysis.find_peaks(after_ov_before, num=50)
        scale = 10000
        for i in range(len(peaksig1)):
            pos1, sig1 = peaksig1[i]
            num = int(round(sig1*scale))
            if num < 0 :
                continue
            for k in range(num):
            #for k in range(1):
                before_NCP_list.append([key, pos1, "+", sig1])
                before_NCP_list.append([key, pos1, "-", sig1])
        for i in range(len(peaksig2)):
            pos2, sig2 = peaksig2[i]
            num = int(round(sig2*scale))
            if num < 0:
                continue
            for k in range(num):
            #for k in range(1):
                after_NCP_list.append([key, pos2, "+", sig2])
                after_NCP_list.append([key, pos2, "-", sig2])

    NCP_seq_list1 = read_NCP_sequences(id_seq, before_NCP_list)
    NCP_seq_list2 = read_NCP_sequences(id_seq, after_NCP_list)
    Afreq1, Gfreq1 = AG_freq(NCP_seq_list1)
    Afreq2, Gfreq2 = AG_freq(NCP_seq_list2)

    fig = plt.figure()
    plt.plot(x_axis, Afreq1, label='AA/AT/TA/TT')
    plt.plot(x_axis, Gfreq1, label='CC/CG/GC/GG')
    #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
    plt.legend()
    #plt.show()
    plt.ylim([0,0.7])
    plt.xlim([-80,80])
    plt.ylabel("Frequency")
    plt.xlabel("Location from dyad (bp)")
    plt.savefig("Before_ATGC_" + "enriched" + ".png")
    plt.close()

    fig = plt.figure()
    Afreq1 = np.asarray(Afreq1)
    Gfreq1 = np.asarray(Gfreq1)    
    plt.plot(x_axis, Afreq1/Afreqb, label='AA/AT/TA/TT')
    plt.plot(x_axis, Gfreq1/Gfreqb, label='CC/CG/GC/GG')
    #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
    plt.legend()
    #plt.show()
    #plt.ylim([0,0.7])
    plt.xlim([-80,80])
    plt.ylabel("Fold change")
    plt.xlabel("Location from dyad (bp)")
    plt.savefig("Before_ATGC_fold_" + "enriched" + ".png")
    plt.close()


    fig = plt.figure()
    plt.plot(x_axis, Afreq2, label='AA/AT/TA/TT')
    plt.plot(x_axis, Gfreq2, label='CC/CG/GC/GG')
    plt.legend()
    #plt.show()
    plt.ylim([0,0.7])
    plt.xlim([-80,80])
    plt.ylabel("Frequency")
    plt.xlabel("Location from dyad (bp)")
    plt.savefig("After_ATGC_" + "enriched" + ".png")
    plt.close()

    fig = plt.figure()
    Afreq2 = np.asarray(Afreq2)
    Gfreq2 = np.asarray(Gfreq2)    
    plt.plot(x_axis, Afreq2/Afreqb, label='AA/AT/TA/TT')
    plt.plot(x_axis, Gfreq2/Gfreqb, label='CC/CG/GC/GG')
    #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
    plt.legend()
    #plt.show()
    #plt.ylim([0,0.7])
    plt.xlim([-80,80])
    plt.ylabel("Fold change")
    plt.xlabel("Location from dyad (bp)")
    plt.savefig("After_ATGC_fold_" + "enriched" + ".png")
    plt.close()
    
    """
   
    
    # sliding motif analysis
    keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
    
    #Right = [[] for i in range(3)]
    #Left = [[] for i in range(3)]
    #for key in keys:
    #    slider1 = key_slider1[key]
    #    slider2 = key_slider2[key]
    #    eqm_flux1 = slider1.eqm_flux()
    #    eqm_flux2 = slider2.eqm_flux()
    #    #diff = [ (eqm_flux2[i] - eqm_flux1[i]) / abs(eqm_flux1[i]) for i in range(len(eqm_flux1)) ]
    #    diff = [ (eqm_flux2[i] - eqm_flux1[i]) for i in range(len(eqm_flux1)) ]
    #    for u in range(3):
    #        if u == 0:
    #            sig_list = eqm_flux1
    #            scale = 100
    #        elif u == 1:
    #            sig_list = eqm_flux2
    #            scale = 100
    #        else:
    #            sig_list = diff
    #            #scale = 10
    #            scale = 100
    #        sig_list = sig_list[10:len(sig_list)-10]
    #        mag = [ abs(value) for value in sig_list]
    #        peaksig = analysis.find_peaks(mag, num=5)
    #        #print peaksig
    #        offset = nucleosome_dna_len/2 + 10
    #        for i in range(len(peaksig)):
    #            pos = peaksig[i][0]
    #            sig = sig_list[pos]
    #            dup_num = int(round(abs(sig)*scale))
    #            if sig > 0:
    #                for k in range(dup_num):
    #                    Right[u].append([key, pos+offset, "+", sig])
    #                    Left[u].append([key, pos+offset, "-", sig])
    #            elif sig < 0:
    #                for k in range(dup_num):
    #                    Left[u].append([key, pos+offset, "+", sig])
    #                    Right[u].append([key, pos+offset, "-", sig])

    #NCP_seq_list1 = read_NCP_sequences(id_seq, Right[0])
    #NCP_seq_list2 = read_NCP_sequences(id_seq, Right[1])
    #NCP_seq_list3 = read_NCP_sequences(id_seq, Right[2])

    #Afreq1, Gfreq1 = AG_freq(NCP_seq_list1)
    #Afreq2, Gfreq2 = AG_freq(NCP_seq_list2)
    #Afreq3, Gfreq3 = AG_freq(NCP_seq_list3)

    #Dinfreq1 = din_freq(NCP_seq_list1)
    #Dinfreq2 = din_freq(NCP_seq_list2)
    #Dinfreq3 = din_freq(NCP_seq_list3)

    #Kmerfreq1 = kmer_freq (NCP_seq_list1, kmer_len=5, states='ATCG')
    #Kmerfreq2 = kmer_freq (NCP_seq_list2, kmer_len=5, states='ATCG')
    #Kmerfreq3 = kmer_freq (NCP_seq_list3, kmer_len=5, states='ATCG')

    #pickle.dump(Afreq1, open("Afreq1.p", "wb"))
    #pickle.dump(Gfreq1, open("Gfreq1.p", "wb"))
    #pickle.dump(Kmerfreq1, open("Kmerfreq1.p", "wb"))
    #pickle.dump(Dinfreq1, open("Dinfreq1.p", "wb"))

    Afreq1 = pickle.load(open("Afreq1.p", "rb"))
    Gfreq1 = pickle.load(open("Gfreq1.p", "rb"))
    Kmerfreq1 = pickle.load(open("Kmerfreq1.p", "rb"))
    Dinfreq1 = pickle.load(open("Dinfreq1.p", "rb"))

    #pickle.dump(Afreq2, open("Afreq2.p", "wb"))
    #pickle.dump(Gfreq2, open("Gfreq2.p", "wb"))
    #pickle.dump(Kmerfreq2, open("Kmerfreq2.p", "wb"))
    #pickle.dump(Dinfreq2, open("Dinfreq2.p", "wb"))

    Afreq2 = pickle.load(open("Afreq2.p", "rb"))
    Gfreq2 = pickle.load(open("Gfreq2.p", "rb"))
    Kmerfreq2 = pickle.load(open("Kmerfreq2.p", "rb"))
    Dinfreq2 = pickle.load(open("Dinfreq2.p", "rb"))

    #pickle.dump(Afreq3, open("Afreq3.p", "wb"))
    #pickle.dump(Gfreq3, open("Gfreq3.p", "wb"))
    #pickle.dump(Kmerfreq3, open("Kmerfreq3.p", "wb"))
    #pickle.dump(Dinfreq3, open("Dinfreq3.p", "wb"))

    Afreq3 = pickle.load(open("Afreq3.p", "rb"))
    Gfreq3 = pickle.load(open("Gfreq3.p", "rb"))
    Kmerfreq3 = pickle.load(open("Kmerfreq3.p", "rb"))
    Dinfreq3 = pickle.load(open("Dinfreq3.p", "rb"))

    kmers = Kmerfreq1.keys()
    kmers_fold1 = {}
    kmers_fold2 = {}
    kmers_fold3 = {}
    for kmer in kmers:
        kmers_fold1[kmer] = Kmerfreq1[kmer] / Kmerfreqb[kmer]
        kmers_fold2[kmer] = Kmerfreq2[kmer] / Kmerfreqb[kmer]
        kmers_fold3[kmer] = Kmerfreq3[kmer] / Kmerfreqb[kmer]

    temp1 = [(fold, kmer) for kmer, fold in kmers_fold1.items()]
    temp1 = sorted(temp1, reverse=True)
    fold1 = []
    kmers = []
    mark = np.zeros(2)
    clist= []
    for u in range(len(temp1)):
        fold, kmer = temp1[u]
        kmers.append(kmer)
        if kmer == 'AAAAA':
            mark[0] = u*20
        if kmer == 'GGGGG':
            mark[1] = u*20
        fold1.append(kmers_fold1[kmer])
        #fold2.append(kmers_fold2[kmer])
        #clist.append(GC_content(kmer))
        clist.append(GC_content(kmer))

    fig = plt.figure()
    #barWidth = 0.25
    r1 = [ i*20 for i in range(len(kmers))]
    #r2 = [x + barWidth for x in r1]
    #plt.bar(r1, freq1_fold, color='#7f6d5f', width=barWidth, edgecolor='white', label='Before')
    #plt.bar(r2, freq2_fold, color='#557f2d', width=barWidth, edgecolor='white', label='After')
    plt.scatter(r1, fold1 ,c=clist, cmap="jet", s=3)
    #plt.scatter(r1, fold2, label='After',c=clist, cmap="Greens", s=3)
    for k in range(len(mark)):
        plt.axvline(x=mark[k], color='r')
    #plt.xlabel('group', fontweight='bold')
    #plt.xticks([r + barWidth/2 for r in range(len(freq1_fold))], kmers)
    #plt.xticks(r1, kmers, rotation=90, fontsize=3)
    plt.xticks(mark, ['AAAAA', 'GGGGG'], rotation=45, fontsize=8)
    plt.ylabel('Fold change')
    plt.xlabel("5-mers")
    #plt.legend()
    plt.savefig("K-merfreq_Right_before.png")
    #plt.show()
    plt.close()

    temp2 = [(fold, kmer) for kmer, fold in kmers_fold2.items()]
    temp2 = sorted(temp2, reverse=True)
    fold2 = []
    kmers = []
    mark = np.zeros(2)
    clist= []
    for u in range(len(temp2)):
        fold, kmer = temp2[u]
        kmers.append(kmer)
        if kmer == 'AAAAA':
            mark[0] = u*20
        if kmer == 'GGGGG':
            mark[1] = u*20
        #fold1.append(kmers_fold1[kmer])
        fold2.append(kmers_fold2[kmer])
        #clist.append(GC_content(kmer))
        clist.append(GC_content(kmer))


    fig = plt.figure()
    #barWidth = 0.25
    r1 = [ i*20 for i in range(len(kmers))]
    #r2 = [x + barWidth for x in r1]
    #plt.bar(r1, freq1_fold, color='#7f6d5f', width=barWidth, edgecolor='white', label='Before')
    #plt.bar(r2, freq2_fold, color='#557f2d', width=barWidth, edgecolor='white', label='After')
    #plt.scatter(r1, fold1, label='Before',c=clist, cmap="Blues", s=3)
    plt.scatter(r1, fold2, c=clist, cmap="jet", s=3)
    for k in range(len(mark)):
        plt.axvline(x=mark[k], color='r')
    #plt.xlabel('group', fontweight='bold')
    #plt.xticks([r + barWidth/2 for r in range(len(freq1_fold))], kmers)
    #plt.xticks(r1, kmers, rotation=90, fontsize=3)
    plt.xticks(mark, ['AAAAA', 'GGGGG'], rotation=45, fontsize=8)
    plt.ylabel('Fold change')
    plt.xlabel("5-mers")
    #plt.legend()
    plt.savefig("K-merfreq_Right_after.png")
    #plt.show()
    plt.close()

    
    temp3 = [(fold, kmer) for kmer, fold in kmers_fold3.items()]
    temp3 = sorted(temp3, reverse=True)
    #fold1, fold2 = [], []
    fold3 = []
    kmers = []
    mark = np.zeros(2)
    clist= []
    for u in range(len(temp3)):
        fold, kmer = temp3[u]
        kmers.append(kmer)
        if kmer == 'AAAAA':
            mark[0] = u*20
        if kmer == 'GGGGG':
            mark[1] = u*20
        #fold1.append(kmers_fold1[kmer])
        fold3.append(kmers_fold3[kmer])
        #clist.append(GC_content(kmer))
        clist.append(GC_content(kmer))



    fig = plt.figure()
    #barWidth = 0.25
    r1 = [ i*20 for i in range(len(kmers))]
    #r2 = [x + barWidth for x in r1]
    #plt.bar(r1, freq1_fold, color='#7f6d5f', width=barWidth, edgecolor='white', label='Before')
    #plt.bar(r2, freq2_fold, color='#557f2d', width=barWidth, edgecolor='white', label='After')
    #plt.scatter(r1, fold1, label='Before',c=clist, cmap="Blues", s=3)
    plt.scatter(r1, fold3, c=clist, cmap="jet", s=3)
    for k in range(len(mark)):
        plt.axvline(x=mark[k], color='r')
    #plt.xlabel('group', fontweight='bold')
    #plt.xticks([r + barWidth/2 for r in range(len(freq1_fold))], kmers)
    #plt.xticks(r1, kmers, rotation=90, fontsize=3)
    plt.xticks(mark, ['AAAAA', 'GGGGG'], rotation=45, fontsize=8)
    plt.ylabel('Fold change')
    plt.xlabel("5-mers")
    #plt.legend()
    plt.savefig("K-merfreq_Right_diff.png")
    #plt.show()
    plt.close()


    fig = plt.figure()
    for kmer in kmers:
        x = kmers_fold1[kmer]
        y = kmers_fold2[kmer]
        plt.plot(x, y, '.')
    small = min(kmers_fold1.values() + kmers_fold2.values())
    large = max(kmers_fold1.values() + kmers_fold2.values())
    plt.plot([0,2], [0,2], '--')
    plt.xlim([small-0.1, large+0.1])
    plt.ylim([small-0.1, large+0.1])
    plt.xlabel("Before")
    plt.ylabel("After")
    plt.savefig("K-merfreq_Right_BvsA.png")
    #plt.show()    
    plt.close()


    kmerchange = []
    for kmer in kmers:
        change = 100 * (kmers_fold2[kmer] - kmers_fold1[kmer]) / kmers_fold1[kmer]
        kmerchange.append((change, kmer))
    kmerchange = sorted(kmerchange, cmp=tuple_cmp)

    changes = []
    mark = np.zeros(2)
    clist = []
    GC = []
    for i in range(len(kmerchange)):
        change, kmer = kmerchange[i]
        changes.append(change)
        clist.append(GC_content(kmer))
        if kmer == 'AAAAA':
            mark[0] = i*20
        if kmer == 'GGGGG':
            mark[1] = i*20
        GC.append(GC_content(kmer))
        print "%s\t%0.2f" % (color_A(kmer), change)

    fig = plt.figure()
    r1 = [ i*20 for i in range(len(changes))]
    plt.scatter(r1, changes, c=clist, cmap="jet", s=3)
    for k in range(len(mark)):
        plt.axvline(x=mark[k], color='r')
    plt.xticks(mark, ['AAAAA', 'GGGGG'], rotation=45, fontsize=8)
    plt.ylabel('Fold change (%)')
    plt.xlabel("5-mers")
    #plt.legend()
    plt.savefig("K-merfreq_Right_change.png")
    #plt.show()
    plt.close()

    hmap1, hmap2, hmap3 = [], [], []
    hmap4 = []
    din_list = ['AA', 'AT', 'TA', 'TT', 'GT', 'TG', 'AC', 'CA', 'AG', 'GA', 'TC', 'CT', 'GG', 'GC', 'CG', 'CC']
    for din in din_list:
        hmap1.append(Dinfreq1[din] / Dinfreqb[din])
        hmap2.append(Dinfreq2[din] / Dinfreqb[din])
        hmap3.append(Dinfreq3[din] / Dinfreqb[din])
        hmap4.append(100 * (Dinfreq2[din]-Dinfreq1[din]) / Dinfreq1[din])

    #x_axis = [ i - len(hmap1[0])*0.5 + 0.5  for i in range(len(hmap1[0]))]
    tick_loc = []
    tick_name = []
    for i in range(0, len(hmap1[0])/2, 10):
        if i == 0:
            tick_loc.append(i + len(hmap1[0])/2)
            tick_name.append(str(i))
        else:
            tick_loc.append(i + len(hmap1[0])/2)
            tick_loc.append(-i + len(hmap1[0])/2)
            tick_name.append(str(i))
            tick_name.append(str(-i))

    fig = plt.figure()
    plt.imshow(hmap1, aspect='auto', interpolation='none', cmap='Purples', vmin=0.8, vmax=1.2)
    plt.yticks(range(len(din_list)), din_list, fontsize=10)
    plt.xticks(tick_loc, tick_name)
    plt.xlabel("Location from dyad (bp)")
    plt.colorbar()
    plt.savefig("Dinfreq_Right_before.png")
    #plt.show()
    plt.close()

    fig = plt.figure()
    plt.imshow(hmap2, aspect='auto', interpolation='none', cmap='Purples', vmin=0.8, vmax=1.2)
    plt.yticks(range(len(din_list)), din_list, fontsize=10)
    plt.xticks(tick_loc, tick_name)
    #plt.xticks(range(0, len(hmap1[0])+3, 10), [ i - 150/2 for i in range(0, len(hmap1[0])+3, 10)])
    plt.xlabel("Location from dyad (bp)")
    plt.colorbar()
    plt.savefig("Dinfreq_Right_after.png")
    plt.colorbar()
    plt.close()

    fig = plt.figure()
    plt.imshow(hmap3, aspect='auto', interpolation='none', cmap='Purples', vmin=0.8, vmax=1.2)
    plt.yticks(range(len(din_list)), din_list, fontsize=10)
    plt.xticks(tick_loc, tick_name)
    #plt.xticks(range(0, len(hmap1[0])+3, 10), [ i - 150/2 for i in range(0, len(hmap1[0])+3, 10)])
    plt.xlabel("Location from dyad (bp)")
    plt.colorbar()
    plt.savefig("Dinfreq_Right_diff.png")
    plt.colorbar()
    plt.close()

    fig = plt.figure()
    plt.imshow(hmap4, aspect='auto', interpolation='none', cmap='jet')
    plt.yticks(range(len(din_list)), din_list, fontsize=10)
    plt.xticks(tick_loc, tick_name)
    #plt.xticks(range(0, len(hmap1[0])+3, 10), [ i - 150/2 for i in range(0, len(hmap1[0])+3, 10)])
    plt.xlabel("Location from dyad (bp)")
    plt.colorbar()
    plt.savefig("Dinfreq_Right_change.png")
    #plt.show()
    plt.close()        

    x_axis = [ i - len(Afreq1)*0.5 + 0.5  for i in range(len(Afreq1))]

    fig = plt.figure()
    plt.plot(x_axis, Afreq1/Afreqb, label='AA/AT/TA/TT')
    plt.plot(x_axis, Gfreq1/Gfreqb, label='CC/CG/GC/GG')
    #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
    plt.legend()
    for n in [-40,-20,0,20,40]:
        plt.axvline(x=n, linestyle='--')
    plt.axhline(y=1, linestyle='--')
    #plt.show()
    #plt.ylim([0,0.7])
    plt.xlim([-80,80])
    plt.ylabel("Frequency")
    plt.xlabel("Location from dyad (bp)")
    plt.savefig("AGfreq_Right_before.png")
    plt.close()

    fig = plt.figure()
    plt.plot(x_axis, Afreq2/Afreqb, label='AA/AT/TA/TT')
    plt.plot(x_axis, Gfreq2/Gfreqb, label='CC/CG/GC/GG')
    plt.legend()
    for n in [-40,-20,0,20,40]:
        plt.axvline(x=n, linestyle='--')
    plt.axhline(y=1, linestyle='--')
    #plt.show()
    #plt.ylim([0,0.7])
    plt.xlim([-80,80])
    plt.ylabel("Frequency")
    plt.xlabel("Location from dyad (bp)")
    plt.savefig("AGfreq_Right_after.png")
    plt.close()

    fig = plt.figure()
    plt.plot(x_axis, Afreq3/Afreqb, label='AA/AT/TA/TT')
    plt.plot(x_axis, Gfreq3/Gfreqb, label='CC/CG/GC/GG')
    #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
    plt.legend()
    for n in [-40,-20,0,20,40]:
        plt.axvline(x=n, linestyle='--')
    plt.axhline(y=1, linestyle='--')
    #plt.show()
    #plt.ylim([0,0.7])
    plt.xlim([-80,80])
    plt.ylabel("Fold change")
    plt.xlabel("Location from dyad (bp)")
    plt.savefig("AGfreq_Right_diff.png")
    plt.close()

    fig = plt.figure()
    plt.plot(x_axis, Afreq2/Afreq1, label='AA/AT/TA/TT')
    plt.plot(x_axis, Gfreq2/Gfreq1, label='CC/CG/GC/GG')
    #xticks = [ i - nucleosome_dna_len/2 for i in range(0,nucleosome_dna_len,5)]
    plt.legend()
    for n in [-40,-20,0,20,40]:
        plt.axvline(x=n, linestyle='--')
    plt.axhline(y=1, linestyle='--')
    #plt.show()
    #plt.ylim([0,0.7])
    plt.xlim([-80,80])
    plt.ylabel("Fold change")
    plt.xlabel("Location from dyad (bp)")
    plt.savefig("AGfreq_Right_change.png")
    plt.close()


    #graph.plot_weblogo(NCP_seq_list1, note='Right_' + name)
    #graph.plot_weblogo(NCP_seq_list2, note='Left_' + name)

    """
    # Randomly select 99% of NCPs as a training set and the rest, 1%, as a validation set
    NCP_random_list = NCP_list[:] # deep copy
    random.shuffle(NCP_random_list)
    num99 = len(NCP_random_list) * 99 / 100
    NCP_tlist, NCP_vlist = NCP_random_list[:num99], NCP_random_list[num99:]
    NCP_seq_tlist, NCP_seq_vlist = read_NCP_sequences(ygenome, NCP_tlist), read_NCP_sequences(ygenome, NCP_vlist)

    mm = HMM.MarkovModel(markov_order)
    mm.train(ygenome, NCP_seq_tlist)
    num_test, num_correct = 0, 0
    for NCP in NCP_vlist:
        chr, pos = NCP[:2]
        max_i, max_score = -1, -sys.float_info.max
        for i in range(max(nucleosome_dna_len / 2, pos - 50), pos + 50):
            seq = read_NCP_sequence(ygenome, [chr, i, 0.0])
            cur_score = mm.predict(seq, background)
            
            # comp_seq = get_comp(seq)
            # cur_score = max(mm.predict(seq), mm.predict(comp_seq))
            if max_score < cur_score:
                max_i = i
                max_score = cur_score

        num_test += 1
        if pos == max_i:
            num_correct += 1

        # DK - for debugging purposes
        if pos != max_i and False:
            print NCP
            print "predicted: %d, score: %f" % (max_i, max_score)
            print "true: %d, score: %f" % (pos, mm.predict(read_NCP_sequence(ygenome, [chr, pos, 0.0])))
            for i in range(max(nucleosome_dna_len / 2, pos - 50), pos + 50):
                seq = read_NCP_sequence(ygenome, [chr, i, 0.0])
                cur_score = mm.predict(seq)
                print "\t%f at %d" % (cur_score, i)

            chr_seq = ygenome[chr]
            chr_len = len(chr_seq)
            left = max(0, pos - 5000)
            right = min(chr_len, pos + 5000)
            seq = chr_seq[left:right]
            F = mm.get_forward(seq)
            R = mm.get_reverse(seq)

            max_i, max_score = -1, -sys.float_info.max
            for i in range(max(nucleosome_dna_len / 2, pos - 50), pos + 50):
                cur_score = mm.logprob_i(seq, i - left, F, R)
                print "genome wide: %f at %d" % (cur_score, i)
                if max_score < cur_score:
                    max_i = i
                    max_score = cur_score

            print "Global: %d, score: %f (%.2f%%)" % (max_i, max_score, math.exp(max_score) * 100)

        
            sys.exit(1)
            
    print "%d-order Markov Model: %.2f%% (%d/%d)" % (markov_order, float(num_correct)/num_test*100, num_correct, num_test)
    # mm.help()
    #sys.exit(1)
    """

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Nucleosome positioning analysis')
    parser.add_argument('--nucleosome-dna-len',
                        dest='nucleosome_dna_len',
                        type=int,
                        default=147,
                        help='Flanking sequence length (both sides including the center, default: 147')
    parser.add_argument('--markov-order',
                        dest='markov_order',
                        type=int,
                        default=1,
                        help='Markov Model order (default: 1)')    
    parser.add_argument('--seed',
                        dest='seed',
                        type=int,
                        default=1,
                        help='Random seeding value (default: 1)')
    parser.add_argument('--no-background',
                        dest='background',
                        action='store_false',
                        help='No background noise for calculating score')
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    if args.nucleosome_dna_len % 2 == 0:
        print >> sys.stderr, "Error: please use an odd number for --nucleosome-dna-len, perhaps %d instead of %d" % (args.nucleosome_dna_len + 1, args.nucleosome_dna_len)
        sys.exit(1)

    random.seed(args.seed)

    Nuc_positioning(args.nucleosome_dna_len,
        args.markov_order,
        args.background,
        args.verbose)
