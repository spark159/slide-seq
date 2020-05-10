import HMM
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
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import seaborn as sns
import sklearn
from sklearn.decomposition import PCA
from pydiffmap import diffusion_map as dm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import ImageGrid
import EnModel
import copy

def read_ref (ref_fname):
    id_seq = {}
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith('>'):
            id = line[1:]
            continue
        if line:
            assert id not in id_seq
            id_seq[id] = line
    return id_seq

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

def tuple_cmp (a, b, priority=[0,1]):
    pos1, pos2 = priority
    if a[pos1] < b[pos1]:
        return -1
    elif a[pos1] > b[pos1]:
        return 1
    else:
        if a[pos2] < b[pos2]:
            return -1
        elif a[pos2] > b[pos2]:
            return 1
        else:
            return 0

def dict_sort (dict, sort_by='value'):
    def tuple_cmp (a, b):
        if a[0] < b[0]:
            return -1
        elif a[0] > b[1]:
            return 1
        else:
            if a[1] < b[1]:
                return -1
            elif a[1] > b[1]:
                return 1
            else:
                return 0
    if sort_by == 'value':
        temp = [(value, key) for key, value in dict.items()]
        temp = sorted(temp, cmp=tuple_cmp)
    return [temp[i][1] for i in range(len(temp))]

ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
NCP_len = 147

# load model
with open("plusonelib_new_0_model.pickle", "rb") as f:
    m1 = pickle.load(f)
with open("plusonelib_new_30_model.pickle", "rb") as f:
    m2 = pickle.load(f)

# load data
name_key_slider = {}

# 601 control data
fnames1 = ["/home/spark159/../../media/spark159/sw/polyAlibFinal/601_before_.combined.sort"]
fnames2 = ["/home/spark159/../../media/spark159/sw/polyAlibFinal/601_after_.combined.sort"]
seq = "ATCCGACTGGCACCGGCAAGGTCGCTGTTCGCCACATGCGCAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCAGGGCGTCCTCGTATAGGGTCCATCACATAAGGGATGAACT"
Control1 = load.load_files(fnames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)
Control2 = load.load_files(fnames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)
Control1['601'].seq = seq
Control2['601'].seq = seq
name_key_slider['Control1'] = Control1
name_key_slider['Control2'] = Control2

Control1_predict = copy.deepcopy(Control1)
pred_energy = m1.predict(seq)
predict = [np.exp(-3.0*value) for value in pred_energy]
predict = [value*100 for value in analysis.norm(predict)]
predict = [0.0]*(NCP_len/2) + predict + [0.0]*(NCP_len/2)
Control1_predict['601'].dyadmap = predict

fig = plt.figure()
plt.plot(analysis.norm(Control1['601'].dyadmap), label='exp')
plt.plot(analysis.norm(Control1_predict['601'].dyadmap), label='predict')
plt.legend()
plt.show()
plt.close()

with open("Control1-predict.pickle", "wb") as f:
    pickle.dump(Control1_predict, f)
    
Control2_predict = copy.deepcopy(Control2)
pred_energy = m2.predict(seq)
predict = [np.exp(-3.0*value) for value in pred_energy]
predict = [value*100 for value in analysis.norm(predict)]
predict = [0.0]*(NCP_len/2) + predict + [0.0]*(NCP_len/2)
Control2_predict['601'].dyadmap = predict

fig = plt.figure()
plt.plot(analysis.norm(Control2['601'].dyadmap), label='exp')
plt.plot(analysis.norm(Control2_predict['601'].dyadmap), label='predict')
plt.legend()
plt.show()
plt.close()

with open("Control2-predict.pickle", "wb") as f:
    pickle.dump(Control2_predict, f)


# PolyA library
polyA_key_seq = read_ref("/home/spark159/../../media/spark159/sw/polyAlibFinal/polyAscanlib.ref")
for condition in ['old']:
    if condition == 'old':
        path = "/home/spark159/../../media/spark159/sw/polyAlibFinal/"
    elif condition == 'new':
        #path = "/home/spark159/../../media/spark159/sw/all_slide_seq_data/"
        path = "/home/spark159/../../media/spark159/sw/polyAlibFinal/"
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
            key_slider = load.load_files([path +  sort_fname], ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', load_ref="/home/spark159/../../media/spark159/sw/polyAlibFinal/polyAscanlib.ref")
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
    for condition in ['control', 'bubble']:
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
plusone_key_seq = read_ref("/home/spark159/../../media/spark159/sw/plusonelibFinal/plusonelib.ref")
for condition in ['new']:
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

# train the model with yeast library data
def sample_seq (key_slider, key_seq, scale=100):
    seq_list = []
    for key in key_slider:
        dyadmap = analysis.norm(key_slider[key].dyadmap)
        for i in range(NCP_len/2, ref_length-NCP_len/2):
            count = int(dyadmap[i]*scale)
            seq = key_seq[key][i - NCP_len/2 : i + NCP_len/2 + 1]
            seq_list += [seq]*count
    return seq_list

def logprob_to_prob (logprob_profile, scale=1.0):
    prob_profile = [ np.exp(scale*value) for value in logprob_profile ]
    total = sum(prob_profile)
    prob_profile = [ float(value)/total for value in prob_profile ]
    return prob_profile

for name in name_key_slider:
    if not name.startswith("polyAlib"):
        continue
    library_type, condition, time, rep = name.split('_')
    if time == '0':
        m = m1
    elif time == '5':
        m = m2

    print name
    key_slider = name_key_slider[name]
    keys = name_keys[name]

    #fname = name + "_" + "predict"
    fname = "_".join([library_type, condition + "-" + "predict", time, rep])
    try:
        with open(fname + ".pickle", "rb") as f:
            key_pslider = pickle.load(f)
    except:
        key_pslider = {}
        for key in keys:
            slider = key_slider[key]
            seq = slider.seq
            pred_energy = m.predict(seq)
            predict = [np.exp(-3.0*value) for value in pred_energy]
            predict = [value*100 for value in analysis.norm(predict)]
            predict = [0.0]*(NCP_len/2) + predict + [0.0]*(NCP_len/2)
            assert len(predict) == ref_length
            key_pslider[key] = Slider(key,
                                      ref_length,
                                      dyad_axis,
                                      dyad_offset,
                                      dyad_offset,
                                      seq,
                                      predict,
                                      [],
                                      [],
                                      None,
                                      None,
                                      None,
                                      None)
        
        with open(fname + ".pickle", "wb") as f:
            pickle.dump(key_pslider, f)

    sample_mode = "polyA:" + "-".join([str(size) for size in [5, 12]])
    sample_list = sample.sampling(key_pslider, sample_mode)
    graph_edit.plot_map(key_pslider, sample_list, True, Slider.peak_signal, draw = "key", slicing=0, note='_' + fname)


"""
try:
    with open("temp.pickle", "rb") as f:
        pred_key_slider = pickle.load(f)
except:
    model1 = EnModel.EnergyModel(name_key_slider['plusonelib_new_30'], NCPlen=NCP_len)
    model1.train(MM_orders=[1], Kmer_k_b=[5,1], PolyA_b=False,  GC_b=False, Harmonic=True, k_fold=3)
    pred_key_slider = model1.predict(name_key_slider['polyAlib_old_5_1rep'])
    with open("temp.pickle", "wb") as f:
        pickle.dump(pred_key_slider, f)


key_slider = name_key_slider['polyAlib_old_5_1rep']
for key in key_slider.keys():
    dyadmap = analysis.norm(key_slider[key].dyadmap)
    pred_dyadmap = pred_key_slider[key].dyadmap
    fig = plt.figure()
    plt.plot(dyadmap, label='experiment')
    plt.plot(pred_dyadmap, label='prediction')
    loc, mtype, nts = key.split('-')
    st = int(loc)
    ed = st+len(nts)
    plt.axvspan(st, ed-1, alpha=0.1, color='red')
    plt.legend()
    plt.show()
"""    



#predict = logprobz_to_prob(model1.energy_predict_profile(seq), scale=-1.0)

#plusone_seq_list1 = sample_seq (name_key_slider['plusonelib_new_0'], plusone_key_seq)
#model1 = HMM.MarkovModel(2)
#model1.train(plusone_key_seq, plusone_seq_list1)
#predict = logprob_to_prob(model1.single_logprob_profile(seq))

#plusone_seq_list2 = sample_seq (name_key_slider['plusonelib_new_30'], plusone_key_seq)
#model2 = HMM.MarkovModel(1)
#model2.train(plusone_key_seq, plusone_seq_list2)

"""
fig = plt.figure()
plt.plot(predict)
plt.show()
plt.close()

key_pslider = {}
for key in polyA_key_seq:
    if key == "BACKBONE":
        continue
    seq = polyA_key_seq[key]
    predict = logprob_to_prob(model1.energy_predict_profile(seq), scale=-1.0)
    key_pslider[key] = Slider(key,
                              ref_length,
                              dyad_axis,
                              dyad_offset,
                              dyad_offset,
                              "",
                              predict,
                              [],
                              [],
                              None,
                              None,
                              None,
                              None)

sample_mode = "polyA:" + "-".join([str(size) for size in [3]])
sample_list = sample.sampling(key_pslider, sample_mode)
graph_edit.plot_map(key_pslider, sample_list, True, Slider.get_dyadmap, draw = "key", slicing=0, note='_predict')
"""
