import math
import random
import pickle
import load
import EnModel
import LinModel
import SliderClass
import numpy as np
import matplotlib.pyplot as plt

def key_cmp (key1, key2):
    if float(key1) < float(key2):
        return -1
    else:
        return 1

def tuple_cmp(a, b):
    if a[0] < b[0]:
        return -1
    else:
        return 1
    
def rev_comp(seq):
    rev = ""
    for nt in seq[::-1]:
        nt = nt.upper()
        if nt == 'A':
            rev += 'T'
        elif nt == 'T':
            rev += 'A'
        elif nt == 'G':
            rev += 'C'
        elif nt == 'C':
            rev += 'G'
    return rev

def normalize (data):
    if type(data) == list:
        total = sum(data)
        return [float(value)/total for value in data]
    elif type(data) == dict:
        total = sum(data.values())
        return {key:float(value)/total for key, value in data.items()}

def combine (dict_list):
    keys = set([])
    for i in range(len(dict_list)):
        if i == 0:
            keys |= set(dict_list[i].keys())
            continue
        keys &= set(dict_list[i].keys())
    keys = list(keys)

    new_dict = {}
    for key in keys:
        value = 0.0
        for i in range(len(dict_list)):            
            try:
                value += dict_list[i][key]
            except:
                pass
        new_dict[key] = value
    return new_dict

def all_path(N, states, arrange=True):
    temp = []
    if N==1:
        temp = list(states)
    else:
        for path in all_path(N-1, states):
            for state in states:
                temp.append(path+state)
    if not arrange:
        return temp
    unique, non_pals = [], []
    for seq in temp:
        if rev_comp(seq) not in unique:
            unique.append(seq)
        else:
            non_pals.append(seq)
    output = unique + non_pals
    assert len(temp) == len(output)
    return output


### parameters
klen = 3 # k-mer length and kmers
NCPlen = 145
ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 53
padding = klen/2
dyad_start = dyad_offset + padding
dyad_end = ref_length - dyad_offset - padding

### load library
path = "./"
ref_name = 'plusonelib.ref'
names = ['SD', 'HS', 'HS+Chd1']
data_names = ['Plusone_SD.pickle', 'Plusone_HS.pickle', 'Plusone_HS+Chd1.pickle']

name_key_slider = {}
for name, data_name in zip(names, data_names):
    key_slider = pickle.load(open(data_name, "rb"))
    name_key_slider[name] = key_slider

    
### set common keys
keys = set([])
for i in range(len(names)):
    name = names[i]
    key_slider = name_key_slider[name]
    if i == 0:
        keys |= set(key_slider.keys())
        continue
    keys &= set(key_slider.keys())
keys = sorted(list(keys), cmp=key_cmp)


### estimate kmer dependent cleavage efficiency
# read cleavage count for each kmers
def get_kmer_count (key_slider, keys=None):
    if keys == None:
        keys = key_slider.keys()

    kmer_back_top, kmer_back_bott = {}, {}
    kmer_count_top, kmer_count_bott = {}, {}
    for kmer in all_path(klen, 'ATCG'):
        kmer_back_top[kmer] = 0
        kmer_back_bott[kmer] = 0
        kmer_count_top[kmer] = 0
        kmer_count_bott[kmer] = 0
    
    for key in keys:
        slider = key_slider[key]
        seq = slider.seq
        top_cutmap, bott_cutmap = slider.right_cutmap, slider.left_cutmap

        for i in range(dyad_start, dyad_end):
            top_pos = i - dyad_offset
            count = top_cutmap[top_pos]
            kmer = seq[top_pos-klen/2:top_pos+klen/2+1]

            kmer_back_top[kmer] += 1
            kmer_count_top[kmer] += float(count)
            
            bott_pos = i + dyad_offset
            count = bott_cutmap[bott_pos]
            kmer = rev_comp(seq[bott_pos-klen/2:bott_pos+klen/2+1])

            kmer_back_bott[kmer] += 1
            kmer_count_bott[kmer] += float(count)

    return kmer_count_top, kmer_count_bott, kmer_back_top, kmer_back_bott

# get fold change w.r.t. background
def get_kmer_fold (kmer_count, kmer_back):
    def fold (a, b):
        if b > 0:
            return float(a)/float(b)
        assert a == 0 and b == 0 
        return 1.0

    kmer_fold = {}
    for kmer in all_path(klen, 'ATCG'):
        kmer_fold[kmer] = fold(kmer_count[kmer], kmer_back[kmer])
    return kmer_fold

name_kmer_fold = {}
name_kmer_freq = {}
name_kmer_freq_top = {}
name_kmer_freq_bott = {}
for name in names:
    key_slider = name_key_slider[name]
    kmer_count_top, kmer_count_bott, kmer_back_top, kmer_back_bott = get_kmer_count(key_slider,
                                                                                    keys = keys)

    kmer_fold_top = get_kmer_fold(normalize(kmer_count_top), kmer_back_top)
    kmer_fold_bott = get_kmer_fold(normalize(kmer_count_bott), kmer_back_bott)
    kmer_fold = combine([kmer_fold_top, kmer_fold_bott])
    name_kmer_fold[name] = kmer_fold

    kmer_freq = normalize(kmer_fold)
    kmer_freq_top = normalize(kmer_fold_top)
    kmer_freq_bott = normalize(kmer_fold_bott)

    name_kmer_freq[name] = kmer_freq
    name_kmer_freq_top[name] = kmer_freq_top
    name_kmer_freq_bott[name] = kmer_freq_bott
    
total_kmer_fold = combine(name_kmer_fold.values())
total_kmer_freq = normalize(kmer_fold)


### plot kmer-frequency
kmers = sorted(all_path(klen, 'ATCG'))

fig = plt.figure(figsize=(12, 5))
Y = [total_kmer_freq[kmer] for kmer in kmers]
plt.plot(range(len(kmers)), Y, '.-')
plt.xticks(range(len(kmers)), kmers, fontsize=6, rotation=45)
plt.title("kmer cleavage frequency")
plt.savefig("kmer_freq_total.png", bbox_inches='tight')
#plt.show()
plt.close()

fig = plt.figure(figsize=(12, 5))
for name in names:
    kmer_freq = name_kmer_freq[name]
    Y = [kmer_freq[kmer] for kmer in kmers]
    plt.plot(range(len(kmers)), Y, '.-', label=name)
plt.xticks(range(len(kmers)), kmers, fontsize=6, rotation=45)
plt.legend()
plt.savefig("kmer_freq_by_exp.png", bbox_inches='tight')
#plt.show()
plt.close()

for name in names:
    fig = plt.figure(figsize=(12, 5))
    kmer_freq_top = name_kmer_freq_top[name]
    kmer_freq_bott = name_kmer_freq_bott[name]
    Y1 = [kmer_freq_top[kmer] for kmer in kmers]
    Y2 = [kmer_freq_bott[kmer] for kmer in kmers]
    plt.plot(range(len(kmers)), Y1, '.-', label='top')
    plt.plot(range(len(kmers)), Y2, '.-', label='bottom')
    plt.xticks(range(len(kmers)), kmers, fontsize=6, rotation=45)
    plt.legend()
    plt.title("Kmer frequency (%s)" % (name))
    plt.savefig("kmer_freq_by_strand_%s.png" % (name), bbox_inches='tight')
    #plt.show()
    plt.close()


### estimate the corrected dyad signal
def correct_signal (key_slider,
                    kmer_freq,
                    keys=None,
                    pickle_name=''):

    if keys == None:
        keys = key_slider.keys()
    
    for key in keys:
        slider = key_slider[key]
        seq = slider.seq
        top_cutmap, bott_cutmap = slider.right_cutmap, slider.left_cutmap
        new_map = []

        for i in range(ref_length):
            if i < dyad_start:
                signal = 0.0
            elif i >= dyad_end:
                signal = 0.0
            else:
                top_pos = i - dyad_offset
                bott_pos = i + dyad_offset
                top_kmer = seq[top_pos-klen/2:top_pos+klen/2+1]
                bott_kmer = rev_comp(seq[bott_pos-klen/2:bott_pos+klen/2+1])
                t = total_kmer_freq[top_kmer]
                b = total_kmer_freq[bott_kmer]
                signal = float(top_cutmap[top_pos])/t + float(bott_cutmap[bott_pos])/b
            new_map.append(signal)

        key_slider[key].dyadmap = normalize(new_map)

    pickle.dump(key_slider, open(pickle_name + ".pickle", "wb"))
    return key_slider

for name in names:
    key_slider = name_key_slider[name]
    correct_signal (key_slider, kmer_freq, keys=keys, pickle_name='Plusone_' + name + '_corrected')
