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
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
import matplotlib as mpl
import LinModel
import EnModel
from scipy import signal
from scipy import fft
from scipy.spatial.distance import squareform
from sklearn import preprocessing
import random
import graph_final
import HMM
import copy

def get_AG_sig (freq):
    freq_MM0 = freq['MM0']
    length = len(freq_MM0)
    nts = LinModel.all_path(1, 'ATCG')
    A_sig, G_sig = np.zeros(length), np.zeros(length)
    for nt in nts:
        row = [ freq_MM0[i][nt] for i in range(length)]
        if nt in ['A', 'T']:
            A_sig += np.asarray(row)
        if nt in ['G', 'C']:
            G_sig += np.asarray(row)
    for i in range(length):
        A_sig[i] = float(A_sig[i]) / sum(freq_MM0[i].values())
        G_sig[i] = float(G_sig[i]) / sum(freq_MM0[i].values())
    return A_sig, G_sig

def get_ATGC_sig (freq):
    freq_MM1 = freq['MM1']
    length = len(freq_MM1)
    nts = LinModel.all_path(2, 'ATCG')
    AT_sig, GC_sig = np.zeros(length), np.zeros(length)
    for nt in nts:
        row = [ freq_MM1[i][nt] for i in range(length)]
        if nt in ['AA', 'AT', 'TA', 'TT']:
            AT_sig += np.asarray(row)
        if nt in ['GG', 'GC', 'CG', 'CC']:
            GC_sig += np.asarray(row)
    for i in range(length):
        AT_sig[i] = float(AT_sig[i]) / sum(freq_MM1[i].values())
        GC_sig[i] = float(GC_sig[i]) / sum(freq_MM1[i].values())
    return AT_sig, GC_sig

def stat_Markov(seq_list, NCPlen, order, sym=False):
    ntdic = {}
    for nt in all_path(order+1, 'ATCG'):
        ntdic[nt] = 0.0

    sample_num = len(seq_list)

    if sym:
        sample_num = 2*sample_num

    freq = [ copy.deepcopy(ntdic) for i in range(NCPlen - order) ]
    for seq in seq_list:
        if sym:
            rev_seq = rev_comp(seq)
        for i in range(len(seq) - order):
            nt = seq[i:i+1+order]
            freq[i][nt] += 1.0 / sample_num
            if sym:
                nt = rev_seq[i:i+1+order]
                freq[i][nt] += 1.0 / sample_num

    mean, std = [], []
    for ntdic in freq:
        mean.append(np.mean(ntdic.values()))
        std.append(np.std(ntdic.values()))

    stdz_freq = []
    for i in range(len(freq)):
        ntdic = freq[i]
        temp = {}
        for nt, value in ntdic.items():
            temp[nt] = (value - mean[i]) / std[i]
        stdz_freq.append(temp)

    return freq, sample_num, mean, std, stdz_freq


def GC_content (seq):
    count = 0
    for nt in seq:
        if nt in ['G', 'C']:
            count += 1
    return 100*float(count)/len(seq)

def poly_score (seq, nts='AT', pos=False):
    num = []
    num_pos = {}
    i = 0
    while i < len(seq):
        if seq[i] in nts:
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
    score = max(num)
    if score < 3:
        return 0
    return score

def get_dincount(seq, din=None):
    if din:
        count = 0
        for i in range(len(seq)-1):
            if seq[i:i+2].upper() == din.upper():
                count +=1
        return count
    din_count = {}
    seq = seq.upper()
    for i in range(len(seq)-1):
        din = seq[i:i+2]
        if 'N' in din:
            continue
        if din not in din_count:
            din_count[din] = 0
        din_count[din] += 1
    return din_count


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

def get_descendants (node_children, root):
    output = set([])
    for node in node_children[root]:
        output.add(node)
        output |= get_descendants(node_children, node)
    return output

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
        path = "/home/spark159/../../media/spark159/sw/polyAlibFinal/"
    elif condition == 'new':
        #path = "/home/spark159/../../media/spark159/sw/all_slide_seq_data/"
        path = "/home/spark159/../../media/spark159/sw/polyAlibFinal/"
    for time in [5]:
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

        # PolyA libary reindexing
        polyAlib_id_newid = pickle.load(open("polyAlib_id_newid.p", "rb"))
        polyAlib_newid_id = pickle.load(open("polyAlib_newid_id.p", "rb"))

        newkey_slider = {}
        for key, slider in key_slider.items():
            newkey = polyAlib_id_newid[key]
            newkey_slider[newkey] = slider
                
        assert fname not in name_key_slider
        #name_key_slider[fname] = key_slider
        name_key_slider[fname] = newkey_slider



# Mismatch/Indel library
path = "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/"
for mtype in []:
    if mtype == 'M':
        library_type = 'mm'
    elif mtype in ['I', 'D']:
        library_type = 'ID'
    for condition in ['bubble']:
        for time in [0, 5]:
            for rep in [1]:
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

# load prediction data
#fname_list = ['Control1-predict', 'Control2-predict',  "polyAlib_old-predict_0_1rep", "polyAlib_old-predict_5_1rep"]
#for fname in fname_list:
#    with open(fname+".pickle", "rb") as f:
#        name_key_slider[fname] = pickle.load(f)

# Plusone library
path = "/home/spark159/../../media/spark159/sw/all_slide_seq_data/"
for condition in ['new:corrected']:
#for condition in []:
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
    elif name.startswith('polyAlib'):
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
#name_strand_threshold = {}
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


# reproduciblity check
if False:
    names = []
    data_mtype = "M"
    for name in name_key_slider:
        if name.startswith("Control"):
            continue
        library_type, condition, time, rep = name.split('_')
        if data_mtype == "M":
            if library_type == "mmlib":
                names.append(name)
        elif data_mtype == "I":
            if library_type == "IDlib":
                names.append(name)
    names = sorted(names)

    common_keys = set([])
    for name in names:
        common_keys |= set(name_key_slider[name].keys())
    common_keys = sorted(list(common_keys), cmp=key_cmp)

    pair_matrix = {}
    for i in range(len(names)-1):
        for j in range(i+1, len(names)):
            name1, name2 = names[i], names[j]
            matrix = np.zeros((len(common_keys), len(common_keys)))
            matrix[:] = np.nan
            for k in range(len(common_keys)-1):
                for l in range(k, len(common_keys)):
                    key1, key2 = common_keys[k], common_keys[l]
                    dyadmap1 = analysis.norm(name_key_slider[name1][key1].dyadmap)
                    dyadmap2 = analysis.norm(name_key_slider[name2][key2].dyadmap)
                    JS = analysis.JS_dist(dyadmap1, dyadmap2)
                    matrix[k][l] = np.exp(-10*(JS))
                    matrix[l][k] = np.exp(-10*(JS))
            assert (name1, name2) not in pair_matrix
            pair_matrix[(name1, name2)] = matrix

    pair_value = {}
    for i in range(len(names)-1):
        for j in range(i+1, len(names)):
            name1, name2 = names[i], names[j]
            value = np.nanmean([pair_matrix[(name1,name2)][k] for k in range(len(common_keys))])
            pair_value[(name1, name2)] = value

    fig = plt.figure()
    grid = ImageGrid(fig, 111, nrows_ncols=(len(names), len(names)), axes_pad=0.1, cbar_mode='single')

    for i in range(len(names)):
        for j in range(len(names)):
            idx = len(names)*i + j
            name1, name2 = names[i], names[j]
            try:
                matrix = pair_matrix[(name1, name2)]
                img = grid[idx].imshow(matrix, cmap="YlOrRd", interpolation='kaiser')
            except:
                matrix = np.zeros((len(common_keys), len(common_keys)))
                if i == j:
                    matrix[:] = np.nan
                    grid[idx].imshow(matrix)
                    library_type, condition, time, rep = name1.split('_')
                    s = library_type + " " + condition + "\n" + time + "min " + rep
                    grid[idx].text(len(common_keys)/2, len(common_keys)/2, s, ha="center", va="center", fontsize=9, weight='bold')
                else:
                    value = pair_value[(name2, name1)]
                    matrix[:] = value
                    grid[idx].imshow(matrix, cmap="GnBu", vmin=min(pair_value.values()), vmax=max(pair_value.values()))
                    if value < max(pair_value.values())/2.0:
                        color = "black"
                    else:
                        color = "white"
                    grid[idx].text(len(common_keys)/2, len(common_keys)/2, str(round(value,2)), ha="center", va="center", fontsize=10, color=color, weight='bold')

    cbar = grid.cbar_axes[0].colorbar(img)
    cbar.ax.set_ylabel('Similarity value', rotation=-90, va="bottom")
    plt.tight_layout()
    if data_mtype == 'M':
        title = "Mismatch library"
    elif data_mtype == 'I':
        title = "Insertion library"
    plt.title(title)
    plt.show()
    plt.close()

# clustering analysis
if False:
    for name in name_key_slider:
        if name.startswith('Control'):
            continue
        print name
        print
        library_type, condition, time, rep = name.split('_')

        prediction = ('predict' in condition.split('-'))

        key_slider = name_key_slider[name]
        keys = name_keys[name]

        # get features
        X = []
        key_KL = {}
        key_num, num_key = {}, []
        for i in range(len(keys)):
            key = keys[i]
            key_num[key] = i
            num_key.append(key)
            dyadmap = analysis.norm(key_slider[key].dyadmap)
            X.append(analysis.norm(dyadmap[NCP_len/2:ref_length-NCP_len/2]))
            if time == '0':
                if prediction:
                    control_map = analysis.norm(name_key_slider['Control1-predict']['601'].dyadmap)
                else:
                     control_map = analysis.norm(Control1['601'].dyadmap)
            elif time == '5':
                if prediction:
                    control_map = analysis.norm(name_key_slider['Control2-predict']['601'].dyadmap)
                else:
                    control_map = analysis.norm(Control2['601'].dyadmap)
            #diffmap = [dyadmap[i] - control_map[i] for i in range(len(dyadmap))]
            #X.append(diffmap[NCP_len/2:ref_length-NCP_len/2])
            key_KL[key] = analysis.KL_div(dyadmap, control_map)

        # Dimensionality reduction by PCA
        # no standardization
        #pca = PCA(n_components=10).fit(X)

        # with standardization
        #scaler = preprocessing.StandardScaler().fit(X)
        scaler = preprocessing.StandardScaler(with_std=False).fit(X)
        X_scaled = scaler.transform(X)
        pca = PCA(n_components=10).fit(X_scaled)
        
        variance_ratio_list = 100* pca.explained_variance_ratio_

        num = 1
        while sum(variance_ratio_list[:num]) <= 90: # can explain 90% of variance
            num +=1

        print num
        #num=10
        #num=3        

        fig = plt.figure()
        plt.plot(range(1, len(variance_ratio_list)+1), variance_ratio_list, '.-')
        plt.axvline(x=num, color = 'red', linestyle='--', alpha=0.5)
        plt.xlabel("PCA component")
        plt.ylabel("Variance (%)")
        plt.xticks(range(1, len(variance_ratio_list)+1))
        #plt.show()
        plt.close()

        Xr = [vect[:num] for vect in pca.transform(X)]

        #sys.exit(1)

        fig = plt.figure()
        for i in range(len(pca.components_)):
            plt.subplot(pca.n_components_, 1, i+1) 
            plt.plot(pca.components_[i])
        #plt.show()
        plt.close()

        # distance matrix
        metric = 'JS'
        fname = name + "_" + metric + "_" + 'matrix'
        try:
            with open(fname + ".pickle", "rb") as f:
                key1_key2_dist = pickle.load(f)
        except:
            key1_key2_dist = {}
            for i in range(len(keys)):
                key1 = keys[i]
                for j in range(i, len(keys)):
                    key2 = keys[j]
                    dyadmap1 = analysis.norm(key_slider[key1].dyadmap)
                    dyadmap2 = analysis.norm(key_slider[key2].dyadmap)
                    dist = analysis.JS_dist(dyadmap1, dyadmap2)
                    if key1 not in key1_key2_dist:
                        key1_key2_dist[key1] = {}
                    key1_key2_dist[key1][key2] = dist
                    if key2 not in key1_key2_dist:
                        key1_key2_dist[key2] = {}
                    key1_key2_dist[key2][key1] = dist

            with open(fname + ".pickle", "wb") as f:
                pickle.dump(key1_key2_dist, f)

        ## Pick clustering algorithm and number of clusters
        algorithm = 'hierarchy'
        cnum = 7
        color_list = ['r', 'g', 'b', 'orange', 'm', 'c', 'darksalmon', 'navy', 'seagreen', 'skyblue']

        if algorithm == 'Kmeans':
            #Kmeans clustering
            key_xr = {}
            for i in range(len(Xr)):
                key = num_key[i]
                key_xr[key] = Xr[i]
            key_cID, cID_keys = analysis.Kmeans(key_xr, cluster_num=cnum, sample_list=None, type_targets=[None, []])
            #color_list = ['r', 'g', 'b', 'y', 'm', 'c', 'darksalmon', 'seagreen', 'skyblue', 'navy']

        
        elif algorithm == 'hierarchy':
            # hierarichal clustering by PCA components
            #Z = linkage(Xr, 'ward')
            #Z = linkage(Xr, 'average', optimal_ordering=True)

            # hierarchial clustering with JS distance metric of raw data
            dist_matrix = np.zeros((len(keys), len(keys))) # similarity matrix
            for i in range(len(keys)-1):
                for j in range(i+1, len(keys)):
                    key1, key2 = keys[i], keys[j]
                    JS = key1_key2_dist[key1][key2]
                    dist_matrix[i][j] = JS
                    dist_matrix[j][i] = JS

            y = squareform(dist_matrix)
            Z = linkage(y, 'average', optimal_ordering=True)
            #Z = linkage(y, 'single', optimal_ordering=True)

            
            #Z = linkage(X, 'weighted', metric=analysis.JS_dist)
            #Z = linkage(X, 'average', metric=analysis.JS_dist)

            node_children = {i:{} for i in range(len(num_key))}
            node_keys = {i:{num_key[i]} for i in range(len(num_key))}
            for i in range(len(Z)):
                node1, node2 = int(Z[i][0]), int(Z[i][1])
                new_node = max(node_keys.keys()) + 1
                node_children[new_node] = set([node1, node2])
                node_keys[new_node] = node_keys[node1] | node_keys[node2]

            #num_cID = [cID-1 for cID in fcluster(Z, 8, 'maxclust')]
            #num_cID = [cID-1 for cID in fcluster(Z, cnum, 'maxclust')]
            num_cID = [cID-1 for cID in fcluster(Z, t=0.4, criterion='distance')]
            cID_keys, key_cID = {}, {}
            for i in range(len(num_cID)):
                cID = num_cID[i]
                key = num_key[i]
                if cID not in cID_keys:
                    cID_keys[cID] = set([])
                cID_keys[cID].add(key)
                key_cID[key] = cID

            cID_nodes = {}
            for cID, keys1 in cID_keys.items():
                for node, keys2 in node_keys.items():
                    if keys1 == keys2:
                        if cID not in cID_nodes:
                            cID_nodes[cID] = set([])
                        cID_nodes[cID].add(node)
                        cID_nodes[cID] |= get_descendants(node_children, node)
                        break

            # KL-div color map
            #color_list = np.linspace(0, 1, num=len(keys))
            cmap = cm.get_cmap('RdPu')
            key_color = {}
            key_list = dict_sort(key_KL)
            for i in range(len(key_list)):
                key = key_list[i]
                #color = cmap(color_list[i])
                scaled_KL = 0.01 + (key_KL[key]-min(key_KL.values()))*float(1.0 - 0.01)/(max(key_KL.values())-min(key_KL.values()))
                color = cmap(scaled_KL)
                key_color[key] = color

            # cluster color map
            #color_list = np.linspace(0.01, 1, num=len(cID_nodes))
            #cmap = cm.get_cmap("jet")
            color_list = ['r', 'g', 'b', 'y', 'm', 'c', 'darksalmon', 'seagreen', 'skyblue', 'navy']*10
            node_color = ['black'] * (2 * len(Xr) - 1)
            for i in range(len(cID_nodes)):
                cID = cID_nodes.keys()[i]
                #color = cmap(color_list[i])
                color = color_list[i]
                for node in cID_nodes[cID]:
                    node_color[node] = color
                    
            #node_color = ['black'] * (2 * len(Xr) - 1)
            
            # plot dendogram
            fig = plt.figure()
            dn = dendrogram(Z, link_color_func=lambda k: node_color[k])

            ax = plt.gca()
            xlbls = ax.get_xmajorticklabels()
            for lbl in xlbls:
                    lbl.set_color(key_color[num_key[int(lbl.get_text())]])

            #plt.savefig("dendrogram.png")
            plt.show()
            plt.close()

            # plot signal for each cluster
            fig = plt.figure()
            for i in range(len(cID_keys)):
                key_list = cID_keys[i]
                for key in key_list:
                    plt.subplot(len(cID_keys),1,i+1)
                    plt.plot(analysis.norm(key_slider[key].dyadmap), alpha=0.2)
                    loc, mtype, nts = key.split('-')
                    st = int(loc)
                    ed = st+len(nts)
                    plt.axvspan(st, ed-1, alpha=0.1, color='red')
            #plt.show()
            plt.close()

        elif algorithm == 'Spectral':
            #Spectral clustering
            A = np.zeros((len(keys), len(keys))) # similarity matrix
            for i in range(len(keys)-1):
                for j in range(i+1, len(keys)):
                    key1, key2 = keys[i], keys[j]
                    dist = key1_key2_dist[key1][key2]
                    #dyadmap1 = analysis.norm(key_slider[key1].dyadmap)
                    #dyadmap2 = analysis.norm(key_slider[key2].dyadmap)
                    #vect1, vect2 = Xr[i], Xr[j]
                    #dist = np.sqrt(np.sum((vect1-vect2)**2))
                    #score = np.exp(-dist**2)
                    #JS = analysis.JS_dist(dyadmap1, dyadmap2)
                    score = np.exp(-10*dist)
                    #score = np.exp(-100*dist**2)
                    A[i][j] = score
                    A[j][i] = score

            idx_cID = sklearn.cluster.spectral_clustering(affinity=A, n_clusters=cnum)
            key_cID, cID_keys = {}, {}
            for i in range(len(idx_cID)):
                key = keys[i]
                cID = idx_cID[i]
                key_cID[key] = cID
                if cID not in cID_keys:
                    cID_keys[cID] = []
                cID_keys[cID].append(key)

            D = np.diag(A.sum(axis=0))
            L = D - A # graph laplacian
            w, v = eigsh(L, 5, which='SM')

            # plot of eigenvector
            color_list = ['r', 'g', 'b', 'orange', 'm', 'c', 'darksalmon', 'navy', 'seagreen', 'skyblue']*10
            eigen_key = sorted([(v[i,1], keys[i]) for i in range(len(keys))], cmp=tuple_cmp)
            fig = plt.figure()
            for i in range(len(eigen_key)):
                eigen, key = eigen_key[i]
                cID = key_cID[key]
                plt.plot(i, eigen, '.', color=color_list[cID])
            plt.ylabel("Eigenvector")
            #plt.show()
            plt.close()

        # plot of reordered similarity matrix
        img = np.zeros((len(keys), len(keys)))
        img[:] = np.nan
        for i in range(len(keys)-1):
            for j in range(i+1, len(keys)):
                #img[i][j] = A[i][j]
                key1, key2 = keys[i], keys[j]
                dist = key1_key2_dist[key1][key2]
                #vect1, vect2 = Xr[i], Xr[j]
                #dist = np.sqrt(np.sum((vect1-vect2)**2))
                #score = np.exp(-dist**2)
                score = np.exp(-10*dist)
                img[i][j] = score

        new_keys = []
        for cID in sorted(cID_keys.keys()):
            #new_keys += cID_keys[cID]
            temp = list(cID_keys[cID])
            random.shuffle(temp)
            new_keys += temp

        #new_keys = [ key for value, key in eigen_key ]

        for i in range(len(new_keys)-1): 
           for j in range(i+1, len(new_keys)):
                key1, key2 = new_keys[i], new_keys[j]
                #idx1, idx2 = key_num[key1], key_num[key2]
                #vect1, vect2 = Xr[i], Xr[j]
                #dist = np.sqrt(np.sum((vect1-vect2)**2))
                dist = key1_key2_dist[key1][key2]
                #score = np.exp(-dist**2)
                score = np.exp(-10*dist)
                img[j][i] = score

        fig = plt.figure()
        #plt.imshow(img, cmap="YlOrRd")
        plt.imshow(img, cmap="jet")
        plt.colorbar()
        #plt.show()
        plt.close()
            
        # PCA plot
        fig = plt.figure()
        for i in range(len(Xr)):
            key = num_key[i]
            plt.plot(Xr[i][0], Xr[i][1], '.', color=color_list[key_cID[key]])
        plt.title("PCA plot")
        #plt.show()
        plt.close()

        # tSNE plot
        perp=23
        tsne = sklearn.manifold.TSNE(n_components=2, perplexity=perp, init='pca', random_state=0)
        trans_data = tsne.fit_transform(X).T

        fig = plt.figure()
        for i in range(len(trans_data[0])):
            key = num_key[i]
            plt.plot(trans_data[0][i], trans_data[1][i], '.', color=color_list[key_cID[key]])
        plt.title("tSNE plot")
        #plt.show()
        plt.close()

        # mean position for each cluster
        cID_meanpos = {}
        for cID in cID_keys:
            temp = []
            for key in cID_keys[cID]:
                meanpos = key_slider[key].median_pos()
                temp.append(meanpos)
            cID_meanpos[cID] = np.mean(temp)

        # read counts for each cluster
        cID_readnums = {}
        for cID in cID_keys:
            for key in cID_keys[cID]:
                readnum = sum(key_slider[key].dyadmap)
                if cID not in cID_readnums:
                    cID_readnums[cID] = []
                cID_readnums[cID].append(readnum)

        fig = plt.figure()
        cIDs = dict_sort(cID_meanpos)
        for i in range(len(cID_readnums)):
                cID = cIDs[i]
                readnums = cID_readnums[cID]
                plt.plot([i]*len(readnums), readnums, '.')
        #plt.show()
        plt.close()

        # plot clustering map
        fig = plt.figure()
        for key in keys:
            loc, mtype, nts = key.split('-')
            st = int(loc)
            ed = st+len(nts)
            length = len(nts)
            pos = (st + ed)/2
            cID = key_cID[key]
            #alpha = 0.01 + (key_KL[key]-min(key_KL.values()))*float(1.0 - 0.01)/(max(key_KL.values())-min(key_KL.values()))
            #plt.plot(range(st, ed), [len(nts)]*len(nts), '.', color=color_list[cID], alpha=alpha)
            plt.plot([st], [len(nts)], 'o', color=color_list[cID])
        plt.show()
        plt.close()

        # plot heatmap for each cluster
        #graph_edit.plot_map(key_slider, [cID_keys[cID] for cID in cIDs], norm_choice=True, obs_func = Slider.get_dyadmap, draw = 'key', slicing=0, note='Clustering')

        #graph_edit.plot_map(key_slider, [cID_keys[cID] for cID in cIDs], norm_choice=True, obs_func = Slider.peak_signal, draw='key', slicing=0, note='Clustering')

        
        # plot KL divergence map
        fig = plt.figure()
        for key in keys:
            loc, mtype, nts = key.split('-')
            st = int(loc)
            ed = st+len(nts)
            length = len(nts)
            pos = (st + ed)/2
            cID = key_cID[key]
            scaled_KL = 0.01 + (key_KL[key]-min(key_KL.values()))*float(1.0 - 0.01)/(max(key_KL.values())-min(key_KL.values()))
            #plt.plot(range(st, ed), [len(nts)]*len(nts), '.', color=color_list[cID], alpha=alpha)
            plt.plot([st], [len(nts)], 'o', color= cm.get_cmap('RdPu')(scaled_KL))
            #plt.plot(range(st, ed), [len(nts)]*len(nts), 'o', color= cm.get_cmap('RdPu')(alpha), alpha=0.5)
        #plt.show()
        plt.close()
        

        # make Widom 601 perturbation map
        size_loc_KLs = {}
        for key in keys:
            loc, mtype, nts = key.split('-')
            size = len(nts)
            if size not in size_loc_KLs:
                size_loc_KLs[size] = {}
            st = int(loc) # - (ref_length - NCP_len)/2
            ed = st+len(nts)
            for i in range(st, ed):
                if i not in size_loc_KLs[size]:
                    size_loc_KLs[size][i] = []
                size_loc_KLs[size][i].append(key_KL[key])

        min_size, max_size = min(size_loc_KLs.keys()), max(size_loc_KLs.keys())

        size_loc_meanKL = {}
        meanKLs = []
        for size in size_loc_KLs:
            if size not in size_loc_meanKL:
                size_loc_meanKL[size] = {}
            for loc in size_loc_KLs[size]:
                meanKL = np.mean(size_loc_KLs[size][loc])
                size_loc_meanKL[size][loc] = meanKL
                meanKLs.append(meanKL)

        min_KL, max_KL = min(meanKLs), max(meanKLs)

        #group_num = 5
        group_num = 7
        #group_num = 10
        page_num = int(math.ceil(float(len(cID_keys)) / group_num))
        for u in range(page_num):
            fig = plt.figure()
            cIDs = dict_sort(cID_meanpos)[group_num*u:min(group_num*(u+1), len(cID_keys))]
            for i in range(len(cIDs)):
                X, Y = [], []
                C = []
                plt.subplot(group_num, 1, i+1)
                ax1 = plt.gca()
                #ax1.text(0, 0, "Cluster" + str(i+1), ha="center", va="center", fontsize=9, weight='bold')
                ax2 = ax1.twinx()
                ax2.plot(control_map, 'k--', linewidth=2)
                sig_list = []
                for key in cID_keys[cIDs[i]]:
                    sig_list.append(analysis.norm(key_slider[key].dyadmap))
                    loc, mtype, nts = key.split('-')
                    st = int(loc)
                    ed = st+len(nts)
                    for k in range(st, ed):
                        X.append(k)
                        Y.append(len(nts))
                        C.append(size_loc_meanKL[len(nts)][k])
                mean_sig = np.mean(sig_list, axis=0)
                error_sig = np.std(sig_list, axis=0)
                ax2.plot(mean_sig, color=color_list[group_num*u+i], linewidth=2)
                ax2.fill_between(range(len(mean_sig)), mean_sig-error_sig, mean_sig+error_sig, color=color_list[group_num*u+i], alpha=0.5)
                ax2.set_xlim([(ref_length - NCP_len)/2, ref_length - (ref_length - NCP_len)/2])
                ax = ax1.scatter(X, Y, c=C, vmin=min_KL, vmax=max_KL, cmap='RdPu', alpha=0.25)
                ax1.set_xticks(range((ref_length - NCP_len)/2, ref_length - (ref_length - NCP_len)/2 + 1, 7))
                ax1.set_ylim([min_size-0.5, max_size+0.5])
                if i < len(cIDs) - 1:
                    ax1.set_xticklabels([])
                else:
                    ax1.set_xticklabels([str(u) for u in range(1, NCP_len+2, 7)])
                ax1.set_yticks(range(min_size, max_size+1, 1))
                ax1.set_ylabel('Size (bp)')
                cbar = plt.colorbar(ax)
                cbar.ax.set_ylabel('KL-div', rotation=-90, va="bottom")
            plt.suptitle("Widom 601 DNA perturbation map")
            plt.show()
            plt.close()

        #sys.exit(1)

# some basic statistics
if False:
    for name in name_key_slider:
        if name.startswith("Control"):
            continue
        library_type, condition, time, rep = name.split('_')
        key_slider = name_key_slider[name]
        keys = name_keys[name]
        size_counts = {}
        for key in keys:
            loc, mtype, nts = key.split('-')
            size, loc = len(nts), int(loc)
            ref_pt = loc + size/2
            if size not in size_counts:
                size_counts[size] = [np.nan]*ref_length
            size_counts[size][ref_pt] = sum(key_slider[key].dyadmap)

        color_list = np.linspace(0, 1, num=len(size_counts))
        cmap = cm.get_cmap("jet")
        fig = plt.figure()
        for i in range(len(size_counts)):
            size = sorted(size_counts.keys())[i]
            counts = size_counts[size]
            plt.plot(range(1, ref_length+1), counts, '.', color=cmap(color_list[i]), label=str(size))
        plt.legend()
        plt.show()
        plt.close()
        
    
    # Align by perturbated sites
    for name in name_key_slider:
        if name.startswith("Control"):
            continue
        library_type, condition, time, rep = name.split('_')
        if time == '0':
            control_map = analysis.norm(Control1['601'].dyadmap)
        elif time == '5':
            control_map = analysis.norm(Control2['601'].dyadmap)

        key_slider = name_key_slider[name]
        keys = name_keys[name]

        key_logratio = {}
        size_logratios = {}
        for key in keys:
            dyadmap = analysis.norm(key_slider[key].dyadmap)
            loc, mtype, nts = key.split('-')
            size, loc = len(nts), int(loc)
            ref_pt = loc + size/2
            logratio = [np.log2(float(dyadmap[i]+1)/(control_map[i]+1)) for i in range(ref_length)]

            key_logratio[key] = Slider(key,
                                       ref_length,
                                       dyad_axis,
                                       dyad_offset,
                                       dyad_offset,
                                       "",
                                       logratio,
                                       [],
                                       [],
                                       None,
                                       None,
                                       None,
                                       None)

            if size not in size_logratios:
                size_logratios[size] = []
            size_logratios[size].append(logratio[ref_pt-40:ref_pt+41])

        for size in size_logratios:
            fig = plt.figure()
            for logratio in size_logratios[size]:
                plt.plot(range(-40, 41), logratio, alpha=1)
            plt.title("Bubble " + str(size) + "bp")
            plt.xlabel("Position w.r.t bubble (bp)")
            #plt.show()
            plt.close()

        #sample_mode = "polyA:" + "-".join([str(size) for size in sorted(size_logratios.keys())])
        #sample_list = sample.sampling(key_logratio, sample_mode)
        #graph_edit.plot_map(key_logratio, sample_list, True, Slider.get_dyadmap, draw=None, slicing=0, note='_signal_' + name)

# energy/force profile analysis
if False:
    name_profile = {}
    for profile in ["energy"]:
        for control_type in ["601"]:
            for name in name_key_slider:
                if name.startswith("Control"):
                    continue
                library_type, condition, time, rep = name.split('_')
                if condition == "control":
                    continue

                prediction = ('predict' in condition.split('-'))
                if prediction:
                    control_type = '601p'

                if time == '0':
                    WT_map = analysis.norm(Control1['601'].dyadmap)[NCP_len/2:ref_length-NCP_len/2]
                    if control_type == '601':
                        control = Control1
                    elif control_type == 'polyA':
                        control = name_key_slider["_".join([library_type, 'control', time, '1rep'])]
                    elif control_type == '601p':
                        control = name_key_slider["Control1-predict"]
                        WT_map = analysis.norm(control['601'].dyadmap)[NCP_len/2:ref_length-NCP_len/2]
                elif time == '5':
                    WT_map = analysis.norm(Control2['601'].dyadmap)[NCP_len/2:ref_length-NCP_len/2]
                    if control_type == '601':
                        control = Control2
                    elif control_type == 'polyA':
                        control = name_key_slider["_".join([library_type, 'control', time, '1rep'])]
                    elif control_type == '601p':
                        control = name_key_slider["Control2-predict"]
                        WT_map = analysis.norm(control['601'].dyadmap)[NCP_len/2:ref_length-NCP_len/2]

                key_slider = name_key_slider[name]
                keys = name_keys[name]

                # read energy file
                fname = name + "_" + profile + "_wrt_" + control_type
                print fname
                try:
                    with open(fname + ".pickle", "rb") as f:
                        size_dyad_shl_values = pickle.load(f)
                except:
                    size_dyad_shl_values = {}
                    for key in keys:
                        loc, mtype, nts = key.split('-')
                        size, loc = len(nts), int(loc)

                        if profile == "energy":
                            test_profile = key_slider[key].energy_profile()
                            if control_type.startswith('601'):
                                control_profile = control['601'].energy_profile()
                            elif control_type == 'polyA':
                                control_profile = control[key].energy_profile()
                            value_profile = test_profile - control_profile
                        elif profile == "force":
                            test_profile = key_slider[key].force_profile()
                            if control_type.startswith('601'):
                                control_profile = control['601'].force_profile()
                            elif control_type == 'polyA':
                                control_profile = control[key].force_profile()
                            value_profile = test_profile - control_profile

                        if size not in size_dyad_shl_values:
                            size_dyad_shl_values[size] = {}

                        for i in range(NCP_len/2, ref_length-NCP_len/2):
                            st = loc - i
                            ed = loc - i + size
                            if ed - st <= 0:
                                continue
                            if i not in size_dyad_shl_values[size]:
                                size_dyad_shl_values[size][i] = {}
                            for k in range(st, ed):
                                if k not in size_dyad_shl_values[size][i]:
                                    size_dyad_shl_values[size][i][k] = []
                                size_dyad_shl_values[size][i][k].append(value_profile[i])

                    with open(fname + ".pickle", "wb") as f:
                        pickle.dump(size_dyad_shl_values, f)
                name_profile[fname] = size_dyad_shl_values

                # check the signal variation for each positions
                size_pos_changes = {}
                all_pos_changes = {}
                if profile == None:
                    for size in size_dyad_shl_values:
                        if size not in size_pos_changes:
                            size_pos_changes[size] = {}
                        for i in range(NCP_len/2, ref_length-NCP_len/2):
                            pos = i - NCP_len/2
                            changes = []
                            for values in size_dyad_shl_values[size][i].values():
                                changes += values
                            size_pos_changes[size][pos] = changes
                            if pos not in all_pos_changes:
                                all_pos_changes[pos] = []
                            all_pos_changes[pos] += changes

                elif profile in ['energy', 'force']:
                    for key in keys:
                        loc, mtype, nts = key.split('-')
                        size, loc = len(nts), int(loc)
                        dyadmap = analysis.norm(key_slider[key].dyadmap)[NCP_len/2:ref_length-NCP_len/2]
                        if size not in size_pos_changes:
                            size_pos_changes[size] = {}
                        for pos in range(len(WT_map)):
                            if pos not in size_pos_changes[size]:
                                size_pos_changes[size][pos] = []
                            change = float(dyadmap[pos] - WT_map[pos])
                            size_pos_changes[size][pos].append(change)
                            if pos not in all_pos_changes:
                                all_pos_changes[pos] = []
                            all_pos_changes[pos].append(change)

                # plot data variation by size
                for size in sorted(size_pos_changes.keys()):
                    fig = plt.figure()
                    ax1 = plt.gca()
                    ax1.plot(WT_map, 'k--')
                    ax2 = ax1.twinx()
                    pos_changes = size_pos_changes[size]
                    for pos in sorted(pos_changes.keys()):
                        changes = pos_changes[pos]
                        ax2.plot([pos]*len(changes), changes, '.', alpha=0.5)
                    ax1.set_ylabel("WT probability")
                    ax2.set_ylabel("Probability change")
                    ax1.set_xlabel("Position (bp)")
                    plt.title("size " + str(size))
                    plt.tight_layout()
                    #plt.show()
                    plt.close()

                # plot total data variation for each position
                #color_list = np.linspace(0, 1, num=len(all_pos_changes))
                #cmap = cm.get_cmap("viridis")
                pos_var = {}
                for pos in sorted(all_pos_changes.keys()):
                    var = np.var(all_pos_changes[pos])
                    pos_var[pos] = var

                def rescale (value_list, old_st, old_ed, new_st, new_ed):
                    output = []
                    for value in value_list:
                        new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
                        output.append(new_value)
                    return output

                var_list = [pos_var[pos] for pos in sorted(pos_var.keys())]
                color_list = rescale(var_list, min(var_list), max(var_list), 0.3, 1)
                cmap = cm.get_cmap('Wistia')

                fig = plt.figure()
                ax1 = plt.gca()
                ax2 = ax1.twinx()
                ax2.plot(WT_map, 'k--')
                for i in range(len(all_pos_changes)):
                    pos = sorted(all_pos_changes.keys())[i]
                    changes = all_pos_changes[pos]
                    ax1.scatter([pos]*len(changes), changes, s=2, c=changes, cmap='Spectral', vmin=-0.1, vmax=0.1)
                    #ax2.plot([pos]*len(changes), changes, '.', color=cmap(color_list[i]), alpha=0.1)
                    #ax2.plot([pos]*len(changes), changes, '.', alpha=0.5, color=cmap(color_list[i]))
                ax2.set_ylabel("WT probability")
                ax1.set_ylabel("Probability change")
                ax2.set_xlabel("Position (bp)")
                plt.title("Signal variation by position")
                #norm = mpl.colors.Normalize(vmin=0.3,vmax=1)
                #sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
                #sm.set_array([])
                #cbar = plt.colorbar(sm)
                #cbar.ax.set_ylabel('Variance', rotation=-90, va="bottom")
                plt.tight_layout()
                #plt.show()
                plt.close()

                # plot total data variance for each position
                fig = plt.figure()
                ax1 = plt.gca()
                ax1.plot(WT_map, 'k--')
                ax2 = ax1.twinx()
                X, Y = [], []
                for pos in sorted(all_pos_changes.keys()):
                    variance = np.var(all_pos_changes[pos])
                    X.append(pos)
                    Y.append(variance)
                ax2.plot(X, Y)
                ax1.set_ylabel("WT probability")
                ax2.set_ylabel("Variance")
                ax1.set_xlabel("Position (bp)")
                plt.title("Signal variance by position")
                plt.tight_layout()
                #plt.show()
                plt.close()    

                # sorting the position by data variance
                size_ordered_pos = {}
                for size in size_pos_changes:
                    pos_changes = size_pos_changes[size]
                    variance_pos = []
                    for pos in size_pos_changes[size]:
                        variance = np.var(size_pos_changes[size][pos])
                        variance_pos.append((variance, pos))
                    variance_pos = sorted(variance_pos, cmp=tuple_cmp, reverse=True)
                    size_ordered_pos[size] = [pos for var, pos in variance_pos]

                variance_pos = [(np.var(changes), pos) for pos, changes in all_pos_changes.items()]
                all_ordered_pos = [pos for var, pos in sorted(variance_pos, cmp=tuple_cmp, reverse=True)]

                # plot nucleosome positioning energy/force perturbation map
                color_list = np.linspace(0, 1, num=len(size_dyad_shl_values))
                cmap = cm.get_cmap("jet")
                marker_list = ["o", "x", "d", ">", "s"]
                for pos in sorted(all_ordered_pos[:3]):
                    #pos = 101 - NCP_len/2
                    fig = plt.figure()
                    for i in range(len(size_dyad_shl_values)):
                        size = sorted(size_dyad_shl_values.keys())[i]
                        shl_values = size_dyad_shl_values[size][pos + NCP_len/2]
                        X, Y = [], []
                        for shl in sorted(shl_values.keys()):
                            X.append(shl)
                            Y.append(np.mean(shl_values[shl]))
                        p = plt.plot(X, Y, '.', markersize=4, color=cmap(color_list[i]), label="size " + str(size))
                        #p = plt.plot(X, Y, '.', markersize=4, label="size " + str(size))
                        plt.plot(X, Y, '-', alpha=0.5, color=p[0].get_color())
                    plt.axvline(x=0, linestyle='--', color='k')
                    plt.axhline(y=0, linestyle='--', color='k')
                    for tick in [-60, -40, -20, -10, 0, 10, 20, 40, 60]:
                        plt.axvline(x=tick, linestyle = '--', linewidth=1, color='k', alpha=0.25)
                    plt.title("Position " + str(pos + NCP_len/2))
                    plt.xlim([-80, 80])
                    plt.xlabel("SHL coordinate")
                    plt.ylabel("$\Delta$ E")
                    leg = plt.legend()
                    for lh in leg.legendHandles:
                        lh._legmarker.set_markersize(5)
                    plt.tight_layout()
                    plt.show()
                    plt.close()

                # Fourier analysis (need to be corrected)
                """
                from astroML.time_series import lomb_scargle
                #from astropy.timeseries import LombScargle

                pos_size_pgram = {}
                T = np.linspace(2, 41, 80)
                F = [2*np.pi/t for t in T]
                for pos in sorted(all_ordered_pos[:3]):
                    if pos not in pos_size_pgram:
                        pos_size_pgram[pos] = {}
                    for i in range(len(size_dyad_shl_values)):
                        size = sorted(size_dyad_shl_values.keys())[i]
                        shl_values = size_dyad_shl_values[size][pos + NCP_len/2]
                        X, Y = [], []
                        for shl in sorted(shl_values.keys()):
                            X.append(shl)
                            Y.append(np.mean(shl_values[shl]))

                        #model = LombScargle(fit_offset=True).fit(T, Y)
                        #pgram = model.score(periods)
                        #pgram = LombScargle(X, Y).power(F)
                        #pgram = lomb_scargle(X, Y, [10**-10]*len(X), F)
                        pgram = signal.lombscargle(X, Y, F, normalize=True, precenter=True)
                        pos_size_pgram[pos][size] = pgram
                        #fig = plt.figure()
                        #ax1 = plt.gca()
                        #ax2 = ax1.twinx()
                        #ax1.plot(X, Y, color='red')
                        #plt.plot(T, pgram)
                        #plt.title(str(size))
                        #plt.show()
                        #plt.close()


                pos_img = {}
                for pos in sorted(all_ordered_pos[:3]):
                    img = []
                    fig = plt.figure()
                    for i in range(len(pos_size_pgram[pos])):
                        #fig = plt.figure()
                        size = sorted(pos_size_pgram[pos].keys())[i]
                        pgram = pos_size_pgram[pos][size]
                        plt.plot(T, pgram, color=cmap(color_list[i]), alpha=0.5, label="size " + str(size))
                        img.append(pgram)
                    pos_img[pos] = img
                    plt.xlabel("Period (bp)")
                    plt.ylabel("Power (a.u.)")
                    plt.title("Position " + str(pos + NCP_len/2))
                    plt.legend()
                    #plt.show()
                    plt.close()

                for i in range(len(pos_img.keys())):
                    fig = plt.figure()
                    pos = sorted(pos_img.keys())[i]
                    img = pos_img[pos]
                    plt.imshow(img, vmin=np.min(img), vmax=np.max(img), interpolation='kaiser', cmap='jet', aspect='auto')
                    plt.yticks(range(len(pos_size_pgram[pos])), sorted(pos_size_pgram[pos].keys()))
                    plt.xticks(rescale(range(5, 41, 5), 2, 41, 0, len(T)), range(5, 41, 5))
                    for tick in rescale(range(5, 41, 5), 2, 41, 0, len(T)):
                        plt.axvline(x=tick, linestyle = '--', linewidth=1, color='white', alpha=0.5)
                    plt.xlabel("Period (bp)")
                    plt.ylabel("Size (nt)")
                    cbar = plt.colorbar()
                    cbar.ax.set_ylabel('Power', rotation=-90, va="bottom")
                    plt.title("Position " + str(pos + NCP_len/2))

                    plt.show()
                    plt.close()
                """

#sys.exit(1)
                
# yeast library analysis
key_slider1 = name_key_slider["plusonelib_new:corrected_0"]
key_slider2 = name_key_slider["plusonelib_new:corrected_30"]

with open("plusonelib_info.pickle", "rb") as f:
    id_info = pickle.load(f)


# check all averaged profile
if True:
    keys = sorted(list(set(key_slider1.keys()) & set(key_slider2.keys())))
    dyadmap_list1 =[]
    dyadmap_list2 = []
    energy_list1 = []
    energy_list2 = []
    for key in keys:
        dyadmap1 = analysis.norm(key_slider1[key].dyadmap)
        dyadmap2 = analysis.norm(key_slider2[key].dyadmap)
        dyadmap_list1.append(dyadmap1)
        dyadmap_list2.append(dyadmap2)
        energy1 = key_slider1[key].energy_profile(scale=1000)
        energy2 = key_slider2[key].energy_profile(scale=1000)
        energy_list1.append(energy1)
        energy_list2.append(energy2)

    mean_dyadmap1 = np.asarray(dyadmap_list1).mean(axis=0)
    mean_dyadmap2 = np.asarray(dyadmap_list2).mean(axis=0)

    mean_energy1 = np.asarray(energy_list1).mean(axis=0)
    mean_energy2 = np.asarray(energy_list2).mean(axis=0)

    fig = plt.figure()
    plt.plot(mean_dyadmap1, label='HS')
    plt.plot(mean_dyadmap2, label='Chd1')
    plt.legend()
    plt.show()

    sys.exit(1)

    
# check Entropy change
if False:
    keys = sorted(list(set(key_slider1.keys()) & set(key_slider2.keys())))
    X, Y = [], []
    dS = []
    Z = []
    F = []
    for key in keys:
        entropy1 = key_slider1[key].entropy(scale=1000)
        entropy2 = key_slider2[key].entropy(scale=1000)
        X.append(entropy1)
        Y.append(entropy2)
        dS.append(entropy2-entropy1)
        dyadmap1 = analysis.norm(key_slider1[key].dyadmap)
        dyadmap2 = analysis.norm(key_slider2[key].dyadmap)
        strand = id_info[int(key)]['strand']
        if strand == '-':
            key_slider1[key].dyadmap = dyadmap1[::-1]
            key_slider2[key].dyadmap = dyadmap2[::-1]
        KL = analysis.KL_div(dyadmap2, dyadmap1)
        Z.append(KL)
        FPKM = np.log2(id_info[int(key)]['FPKM'])
        F.append(FPKM)

    fig = plt.figure(figsize=(2,1.4))
    plt.scatter(X, Y, c=Z, s=2, cmap='rainbow')
    #plt.plot(X, Y, '.')
    plt.plot([min(X+Y)*(1-0.2), max(X+Y)*(1+0.2)], [min(X+Y)*(1-0.2), max(X+Y)*(1+0.2)], 'k--', linewidth=1)
    #plt.xscale("log", basex=2)
    #plt.yscale("log", basey=2)
    plt.gca().tick_params(axis='both', which='major', labelsize=5)
    plt.gca().tick_params(axis='both', which='minor', labelsize=5)
    plt.xlim([4.2, 4.8])
    plt.ylim([4.2, 4.8])
    plt.xlabel("Before", fontsize=6)
    plt.ylabel("After", fontsize=6)
    plt.title("Entropy change", fontsize=6)
    cbar = plt.colorbar(ticks=[min(Z), max(Z)], shrink=0.6, aspect=20)
    cbar.ax.set_yticklabels([round(min(Z), 2), round(max(Z), 2)], fontsize=5)
    cbar.ax.set_ylabel('KL-divergence', rotation=-90, va="bottom", fontsize=5, labelpad=-5)
    plt.savefig("Entropy_scatter.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    fig = plt.figure()
    plt.plot(dS, Z, '.', markersize=4)
    #plt.plot([min(dS+Z), max(dS+Z)], [min(dS+Z), max(dS+Z)], '--')
    plt.axvline(x=0, linestyle='--', color='k')
    plt.axhline(y=0, linestyle='--', color='k')
    plt.xlabel("Entropy change")
    plt.ylabel("KL-divergence")
    plt.xscale("log", basex=2)
    plt.yscale("log", basey=2)
    #plt.show()
    plt.close()

    fig = plt.figure()
    plt.scatter(Y, F, c=Z, s=4, cmap='jet')
    #plt.plot(X, Y, '.')
    #plt.plot([min(X+Y), max(X+Y)], [min(X+Y), max(X+Y)], '--')
    plt.xlabel("Before")
    plt.ylabel("After")
    plt.title("Entropy change")
    cbar = plt.colorbar()
    cbar.set_label("KL-divergence")
    #plt.show()
    plt.close()


    
    sys.exit(1)
    
# check by strand
if False:
    X1, Y1, Z1 = [], [], []
    X2, Y2, Z2 = [], [], []
    for id in id_info:
        mean1 = key_slider1[str(id)].mean_pos()
        mean2 = key_slider2[str(id)].mean_pos()
        entropy1 = key_slider1[str(id)].entropy(scale=1000)
        entropy2 = key_slider2[str(id)].entropy(scale=1000)
        entropy_change = entropy2 - entropy1
        mean_change = mean2 - mean1
        FPKM = np.log2(id_info[id]['FPKM'])
        if np.isnan(FPKM):
            continue
        strand = id_info[id]['strand']
        if strand == '+':
            X1.append(mean1)
            Y1.append(mean2)
            Z1.append(FPKM)
        elif strand == '-':
            X2.append(mean1)
            Y2.append(mean2)
            Z2.append(FPKM)

    fig = plt.figure()
    plt.plot(X1, Y1, '.', color='tab:red', alpha=0.8, label='+ strand')
    plt.plot(X2, Y2, '.', color='tab:blue', alpha=0.8, label='- strand')
    ax = plt.gca()
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    plt.plot(lims, lims, 'k--', alpha=0.5)
    plt.xlim(lims)
    plt.ylim(lims)
    plt.xlabel("Before")
    plt.ylabel("After")
    plt.title('Mean position of nucleosome')
    plt.legend()
    plt.show()
    plt.close()

    fig = plt.figure()
    positions1 = [pos - 0.2 for pos in range(2)]
    positions2 = [pos + 0.2 for pos in range(2)]
    box1 = plt.boxplot([X1, Y1], positions=positions1, widths=0.25, notch=True, patch_artist=True)
    for patch in box1['boxes']:
        patch.set_facecolor('tab:red')
    box2 = plt.boxplot([X2, Y2], positions=positions2, widths=0.25, notch=True, patch_artist=True)
    for patch in box2['boxes']:
        patch.set_facecolor('tab:blue')
    plt.xticks(range(2), ['Before', 'After'])
    plt.xlim([-0.5, 2-0.5])
    plt.legend([box1["boxes"][0], box2["boxes"][0]], ['+ strand', '- strand'], loc='best')
    plt.ylabel("Mean position")
    plt.title('Mean position of nucleosome')
    plt.show()
    plt.close()

    #sys.exit(1)

# group by KL divergence value
if False:

    keys = sorted(list(set(key_slider1.keys()) & set(key_slider2.keys())))
    key_KL, key_input = {}, {}
    for key in keys:
        dyadmap1 = analysis.norm(key_slider1[key].dyadmap)
        dyadmap2 = analysis.norm(key_slider2[key].dyadmap)
        strand = id_info[int(key)]['strand']
        if strand == '-':
            key_slider1[key].dyadmap = dyadmap1[::-1]
            key_slider2[key].dyadmap = dyadmap2[::-1]
        KL = analysis.KL_div(dyadmap2, dyadmap1)
        key_KL[key] = KL
        key_input[key] = [KL]

    KL_key = sorted([(KL, key) for key, KL in key_KL.items()], cmp=tuple_cmp)
    key_order = {}
    for i in range(len(KL_key)):
        KL, key = KL_key[i]
        key_order[key] = i

    # group by k-means clustering
    #gnum = 3
    #key_cdx, cdx_keys = analysis.Kmeans(key_input, cluster_num=gnum)
    #temp = sorted([(key_KL[cdx_keys[i][0]], i) for i in range(len(cdx_keys))], cmp=tuple_cmp)
    #temp = [cdx for KL, cdx in temp]

    #cdx_group = {}
    #for i in range(gnum):
    #    cdx = temp[i]
    #    cdx_group[cdx] = i

    #key_group, group_keys = {}, [[] for i in range(gnum)]
    #for key in key_cdx:
    #    key_group[key] = cdx_group[key_cdx[key]]

    #for cdx in cdx_keys:
    #    group = cdx_group[cdx]
    #    temp = sorted([(key_KL[key], key) for key in cdx_keys[cdx]], cmp=tuple_cmp)
    #    group_keys[group] += [key for KL, key in temp]

    # group by equal partition
    gnum = 4
    gsize = int(math.ceil(float(len(keys))/gnum))
    group_keys = []
    for i in range(gnum):
        st = i*gsize
        ed = min((i+1)*gsize, len(keys))
        temp = []
        for KL, key in KL_key[st:ed]:
            temp.append(key)
        group_keys.append(temp)
    assert len(keys) == sum([len(key_list) for key_list in group_keys])

    # group by equal value range
    #gnum = 3
    #min_KL, max_KL = min(key_KL.values()), max(key_KL.values())
    #step_size = float(max_KL - min_KL)/gnum
    #bins = []
    #for i in range(gnum):
    #    st = min_KL + i*step_size
    #    if i != gnum-1:
    #        ed = min_KL + (i+1)*step_size
    #    else:
    #        ed = max_KL + 0.1
    #    bins.append((st,ed))
    #group_keys = [[] for i in range(gnum)]
    #for KL, key in KL_key:
    #    for i in range(gnum):
    #        st, ed = bins[i]
    #        if KL >= st and KL < ed:
    #            group_keys[i].append(key)
    #            break

    color_list = ['r', 'g', 'b', 'orange']

    fig = plt.figure(figsize=(2.2,1.2))
    for i in range(len(group_keys)):
        X = [key_order[key] + 1 for key in group_keys[i]]
        Y = [key_KL[key] for key in group_keys[i]]
        plt.plot(X, Y, '.', markersize=2, label='Cluster ' + str(i+1), color=color_list[i])
    plt.xlabel("Sequences", fontsize=6)
    plt.ylabel("KL-divergence", fontsize=6)
    plt.title("Before VS After", fontsize=6)
    plt.gca().tick_params(axis='both', which='major', labelsize=5)
    plt.gca().tick_params(axis='both', which='minor', labelsize=5)
    leg = plt.legend(frameon=False, fontsize=6, markerscale=2)
    #plt.savefig("group_byKL.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    sample_num = 10
    sample_list = []
    for key_list in group_keys:
        step_size = len(key_list) / sample_num
        for i in range(sample_num):
            idx = step_size*i
            key = key_list[idx]
            sample_list.append(key)
    sample_list = [sample_list]

    #sample_list = [[key for KL, key in KL_key]]

    #graph_edit.plot_map(key_slider1, sample_list, True, Slider.KDE, draw=None, slicing=0, note='_before')
    #graph_edit.plot_map(key_slider2, sample_list, True, Slider.KDE, draw=None, slicing=0, note='_after')
    #tlen =225
    #graph_final.plot_map(key_slider1, Slider.KDE, ids=sample_list[0], thickness=[8, 2, 2], cmap='YlGnBu', save=True, note='_before', xticks=[[i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)]], figscale=150, fontsize=5)

    #graph_final.plot_map(key_slider2, Slider.KDE, ids=sample_list[0], thickness=[8, 2, 2], cmap='YlGnBu', save=True, note='_after', xticks=[[i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)]], figscale=150, fontsize=5)



    entropies_list1, entropies_list2 = [], []
    means_list1, means_list2 = [], []
    FPKMs_list1, FPKMs_list2 = [], []
    for i in range(gnum):
        entropies1, entropies2 = [], []
        means1, means2 = [], []
        FPKMs = []
        for key in group_keys[i]:
            entropy1 = key_slider1[key].entropy(scale=100)
            entropy2 = key_slider2[key].entropy(scale=100)
            entropies1.append(entropy1)
            entropies2.append(entropy2)
            mean1 = key_slider1[key].mean_pos()
            mean2 = key_slider2[key].mean_pos()
            means1.append(mean1)
            means2.append(mean2)
            FPKM = np.log2(id_info[int(key)]['FPKM'])
            if np.isnan(FPKM):
                continue
            FPKMs.append(FPKM)
        entropies_list1.append(entropies1)
        entropies_list2.append(entropies2)
        means_list1.append(means1)
        means_list2.append(means2)
        FPKMs_list1.append(FPKMs)
        FPKMs_list2.append(FPKMs)

    fig = plt.figure()
    for i in range(gnum):
        X, Y = entropies_list1[i], entropies_list2[i]
        plt.plot(X, Y, '.', alpha=0.75, label="Cluster " + str(i+1), color=color_list[i])
    ax = plt.gca()
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    plt.plot(lims, lims, 'k--', alpha=0.5)
    plt.legend()
    plt.xlabel("Before")
    plt.ylabel("After")
    plt.xlim(lims)
    plt.ylim(lims)
    plt.title("Entropy change")
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(2, 1.4))
    positions1 = [pos - 0.2 for pos in range(gnum)]
    positions2 = [pos + 0.2 for pos in range(gnum)]
    box1 = plt.boxplot(entropies_list1, positions=positions1, widths=0.25, notch=True, patch_artist=True, showfliers=False)
    for patch in box1['boxes']:
        patch.set_facecolor('tab:red')
    box2 = plt.boxplot(entropies_list2, positions=positions2, widths=0.25, notch=True, patch_artist=True, showfliers=False)
    for patch in box2['boxes']:
        patch.set_facecolor('tab:blue')
    plt.xticks(range(gnum), ["Cluster " + str(i+1) for i in range(gnum)], fontsize=6, rotation=40)
    plt.xlim([-0.5, gnum-0.5])
    plt.legend([box1["boxes"][0], box2["boxes"][0]], ['Before', 'After'], loc='lower left', fontsize=6)
    plt.ylabel("Entropy", fontsize=6)
    plt.title('Entropy by cluster', fontsize=6)
    plt.gca().tick_params(axis='y', which='major', labelsize=5)
    plt.gca().tick_params(axis='y', which='minor', labelsize=5)
    #plt.savefig("Entropy_by_group_boxr.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(2, 1.4))
    positions1 = [pos - 0.2 for pos in range(gnum)]
    positions2 = [pos + 0.2 for pos in range(gnum)]
    box1 = plt.boxplot(means_list1, positions=positions1, widths=0.25, notch=True, patch_artist=True, showfliers=False)
    for patch in box1['boxes']:
        patch.set_facecolor('tab:red')
    box2 = plt.boxplot(means_list2, positions=positions2, widths=0.25, notch=True, patch_artist=True, showfliers=False)
    for patch in box2['boxes']:
        patch.set_facecolor('tab:blue')
    plt.xticks(range(gnum), ["Cluster " + str(i+1) for i in range(gnum)], fontsize=6, rotation=40)
    plt.xlim([-0.5, gnum-0.5])
    plt.legend([box1["boxes"][0], box2["boxes"][0]], ['Before', 'After'], loc='lower left', fontsize=6)
    plt.ylabel("Mean position", fontsize=6)
    plt.title('Mean position by cluster', fontsize=6)
    plt.gca().tick_params(axis='y', which='major', labelsize=5)
    plt.gca().tick_params(axis='y', which='minor', labelsize=5)
    #plt.savefig("MPosition_by_group_boxr.svg", format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    fig = plt.figure()
    positions1 = [pos - 0.2 for pos in range(gnum)]
    positions2 = [pos + 0.2 for pos in range(gnum)]
    box1 = plt.boxplot(FPKMs_list1, positions=positions1, widths=0.25, notch=True, patch_artist=True)
    for patch in box1['boxes']:
        patch.set_facecolor('tab:red')
    box2 = plt.boxplot(FPKMs_list2, positions=positions2, widths=0.25, notch=True, patch_artist=True)
    for patch in box2['boxes']:
        patch.set_facecolor('tab:blue')
    plt.xticks(range(gnum), ["Cluster " + str(i+1) for i in range(gnum)])
    plt.xlim([-0.5, gnum-0.5])
    plt.legend([box1["boxes"][0], box2["boxes"][0]], ['Before', 'After'], loc='best')
    plt.ylabel("FPKMs")
    plt.title('Gene expression by cluster')
    #plt.show()
    plt.close()


    #sys.exit(1)

    # Mononucleotide periodicity analysis
    M1 = LinModel.LinSliderModel(key_slider1)
    M1.report (MM_orders=[0], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None)
    M2 = LinModel.LinSliderModel(key_slider2)
    M2.report (MM_orders=[0], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None)

    AT_sig1, GC_sig1 = get_AG_sig(M1.freq)
    AT_sig2, GC_sig2 = get_AG_sig(M2.freq)


    #fig = plt.figure(figsize=(2.6, 2))
    fig = plt.figure(figsize=(3, 2.2))
    plt.plot(AT_sig1, color='hotpink', lw=1.2, label='A/T (HS)')
    plt.plot(GC_sig1, color='lightskyblue', lw=1.2, label='G/C (HS)')
    plt.plot(AT_sig2, '--', color='red', lw=1.2, label='A/T (Chd1)')
    plt.plot(GC_sig2, '--', color='blue', lw=1.2, label='G/C (Chd1)')
    line_list = [147/2 - i*20 for i in range(4)] + [147/2 + i*20 for i in range(1, 4)]
    for line in line_list:
        plt.axvline(x=line, color='k', linestyle='--', alpha=0.25, lw=1)
    plt.xticks([147/2 + 10*i for i in range(-7, 8, 2)], [str(10*i) for i in range(-7, 8, 2)], fontsize=5)
    plt.xlabel("Super Helical Location", fontsize=6)
    plt.ylabel("Relative frequency", fontsize=6)
    plt.gca().tick_params(axis='both', which='major', labelsize=5)
    plt.gca().tick_params(axis='both', which='minor', labelsize=5)
    leg = plt.legend(frameon=True, fontsize=5)

    #plt.ylim([0.235, 0.275])
    #plt.ylim([0.22, 0.28])
    #plt.savefig('ATGCperiod_' + "slide" + '.png')
    plt.title("Mononucleotide periodicity", fontsize=8)
    plt.savefig('AGperiod_' + "slide" + '.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    sys.exit(1)

    # homopolymer analysis
    nts_list = ['AT', 'GC']
    M1 = LinModel.LinSliderModel(key_slider1)
    M1.report (MM_orders=[], Kmer_k_b=[], PolyA_b=None, GC_b=None, Harmonic=None, Polylen=nts_list)
    M2 = LinModel.LinSliderModel(key_slider2)
    M2.report (MM_orders=[], Kmer_k_b=[], PolyA_b=None, GC_b=None, Harmonic=None, Polylen=nts_list)

    for nts in nts_list:
        mname = 'Polylen' + '/'.join(nts)
        profile1 = M1.freq[mname]
        profile2 = M2.freq[mname]

        fig = plt.figure(figsize=(3, 2.2))
        plt.plot(profile1, '--', color='tab:orange', lw=1.2, label='HS')
        plt.plot(profile2, 'k-', lw=1.2, label='Chd1')
        line_list = [147/2 - i*20 for i in range(4)] + [147/2 + i*20 for i in range(1, 4)]
        for line in line_list:
            plt.axvline(x=line, color='k', linestyle='--', alpha=0.25, lw=1)
        plt.xticks([147/2 + 10*i for i in range(-7, 8, 2)], [str(10*i) for i in range(-7, 8, 2)], fontsize=5)
        plt.xlabel("Super Helical Location", fontsize=6)
        plt.ylabel("Fold change", fontsize=6)
        plt.gca().tick_params(axis='both', which='major', labelsize=5)
        plt.gca().tick_params(axis='both', which='minor', labelsize=5)
        leg = plt.legend(frameon=True, fontsize=5)

        #plt.ylim([0.235, 0.275])
        #plt.ylim([0.22, 0.28])
        #plt.savefig('ATGCperiod_' + "slide" + '.png')
        title = "Poly(d%s:d%s) tract length" % tuple(list(nts))
        plt.title(title, fontsize=8)
        plt.savefig('tractlength_' + nts + '.svg', format='svg', bbox_inches='tight')
        #plt.show()
        plt.close()
    

    
    sys.exit(1)

    # 4-mer analysis
    klen = 4
    MMname = 'MM' + str(klen-1)
    M1 = LinModel.LinSliderModel(key_slider1)
    M1.report (MM_orders=[klen-1], Kmer_k_b=[], PolyA_b=None, GC_b=None, Harmonic=None)
    M2 = LinModel.LinSliderModel(key_slider2)
    M2.report (MM_orders=[klen-1], Kmer_k_b=[], PolyA_b=None, GC_b=None, Harmonic=None)

    all_kmers = [kmer for GC, kmer in sorted([(GC_content(kmer), kmer) for kmer in LinModel.all_path(klen, 'ATCG')])]
    tlen = len(M1.freq[MMname])

    all_foldchange = []
    kmer_change_profile = {}
    for kmer in all_kmers:
        for i in range(tlen):
            freq1 = M1.freq[MMname][i][kmer] / float(sum(M1.freq[MMname][i].values()))
            freq2 = M2.freq[MMname][i][kmer] / float(sum(M2.freq[MMname][i].values()))
            #foldchange = (freq2 - freq1)/float(freq1)
            foldchange = (freq2 - freq1)
            #foldchange = float(freq2)/freq1
            all_foldchange.append(foldchange)
            if kmer not in kmer_change_profile:
                kmer_change_profile[kmer] = [0.0]*tlen
            kmer_change_profile[kmer][i] = foldchange

    group_num = 2
    group_size = len(all_kmers) / group_num
    limit = max(abs(min(all_foldchange)), abs(max(all_foldchange)))
    #limit = max(abs(min(all_foldchange)-1), abs(max(all_foldchange)-1))
    vmin, vmax = -limit, limit
    #vmin, vmax = min(all_foldchange), max(all_foldchange)
    #vmin, vmax = 1-limit, 1+limit
    #vmin, vmax = None, None
    cmap = 'seismic'
    for i in range(group_num):
        st, ed = group_size*i, group_size*(i+1)

        sub_kmers = all_kmers[st:ed]
        
        img = []
        for kmer in sub_kmers:
            change_profile = kmer_change_profile[kmer]
            img.append(change_profile)

        fig = plt.figure(figsize=(3.2, 11))
        #fig = plt.figure(figsize=(7.5, 11))
        im = plt.imshow(img, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
        if vmin == None and vmax == None:
            vmin, vmax = im.get_clim()
        plt.yticks(range(len(sub_kmers)), sub_kmers, horizontalalignment='right', fontname='monospace', fontsize=5)
        plt.xticks([tlen/2 + 10*k for k in range(-7, 8, 2)], [str(10*k) for k in range(-7, 8, 2)], fontsize=6)
        plt.savefig("4-mers_profile_change_g%s.svg" % (i+1), format='svg', bbox_inches='tight') 
        #plt.colorbar()
        #plt.show()
        plt.close()

    fig = plt.figure(figsize=(0.74,1.3))
    cbar = plt.colorbar(im, cax=plt.gca(), ticks=[vmin, vmax], format='%.0e')
    #cbar.ax.set_yticklabels([vmin, vmax], fontsize=5)
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(5)
    cbar.ax.set_ylabel('$\Delta$freq (Chd1-HS)', rotation=-90, va="bottom", fontsize=6, labelpad=-8)
    plt.tight_layout()
    plt.savefig('change_kmer_cbar.svg', format='svg', bbox_inches='tight')
    plt.close()

    #sys.exit(1)

    
    # Dinucleotide periodicity analysis
    M1 = LinModel.LinSliderModel(key_slider1)
    M1.report (MM_orders=[1], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None)
    M2 = LinModel.LinSliderModel(key_slider2)
    M2.report (MM_orders=[1], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None)

    AT_sig1, GC_sig1 = get_ATGC_sig(M1.freq)
    AT_sig2, GC_sig2 = get_ATGC_sig(M2.freq)


    #fig = plt.figure(figsize=(2.6, 2))
    fig = plt.figure(figsize=(3, 2.2))
    plt.plot(AT_sig1, color='hotpink', lw=1.2, label='AT-rich (HS)')
    plt.plot(GC_sig1, color='lightskyblue', lw=1.2, label='GC-rich (HS)')
    plt.plot(AT_sig2, '--', color='red', lw=1.2, label='AT-rich (Chd1)')
    plt.plot(GC_sig2, '--', color='blue', lw=1.2, label='GC-rich (Chd1)')
    line_list = [147/2 - i*20 for i in range(4)] + [147/2 + i*20 for i in range(1, 4)]
    for line in line_list:
        plt.axvline(x=line, color='k', linestyle='--', alpha=0.25, lw=1)
    plt.xticks([147/2 + 10*i for i in range(-7, 8, 2)], [str(10*i) for i in range(-7, 8, 2)], fontsize=5)
    plt.xlabel("Super Helical Location", fontsize=6)
    plt.ylabel("Relative frequency", fontsize=6)
    plt.gca().tick_params(axis='both', which='major', labelsize=5)
    plt.gca().tick_params(axis='both', which='minor', labelsize=5)
    leg = plt.legend(frameon=True, fontsize=5)

    plt.ylim([0.235, 0.275])
    #plt.ylim([0.22, 0.28])
    #plt.savefig('ATGCperiod_' + "slide" + '.png')
    plt.title("Dinucleotide periodicity", fontsize=8)
    plt.savefig('ATGCperiod_' + "slide" + '.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    sys.exit(1)


    #plt.show()
    plt.close()
    
    
    MMfreq_list1, MMfreq_list2 = [], []
    for i in range(gnum):
        M1 = LinModel.LinSliderModel(key_slider1, keys=group_keys[i])
        M1.report (MM_orders=[1], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None)
        M2 = LinModel.LinSliderModel(key_slider2, keys=group_keys[i])
        M2.report (MM_orders=[1], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None)
        MMfreq_list1.append(M1.freq)
        MMfreq_list2.append(M2.freq)

    AT_sig_list1, GC_sig_list1 = [], []
    AT_sig_list2, GC_sig_list2 = [], []
    for i in range(gnum):
        freq1 = MMfreq_list1[i]
        AT_sig1, GC_sig1 = get_ATGC_sig(freq1)
        AT_sig_list1.append(AT_sig1)
        GC_sig_list1.append(GC_sig1)
        freq2 = MMfreq_list2[i]
        AT_sig2, GC_sig2 = get_ATGC_sig(freq2)
        AT_sig_list2.append(AT_sig2)
        GC_sig_list2.append(GC_sig2)


    color_list = np.linspace(0.3, 0.8, num=gnum)
    #color_list = ["tab:red", "tab:green", "tab:blue", "tab:olive"]
    cmap1, cmap2 = cm.get_cmap("OrRd"), cm.get_cmap("GnBu")
    cmap = cm.get_cmap("rainbow")

    fig = plt.figure(figsize=(2.6, 2))
    #plt.subplot(1,2,1)
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    for i in range(gnum):
        AT_sig1, GC_sig1 = AT_sig_list1[i], GC_sig_list1[i]
        AT_sig2, GC_sig2 = AT_sig_list2[i], GC_sig_list2[i]
        ax1.plot(GC_sig1, color=cmap1(color_list[i]), lw=1)
        ax1.plot(GC_sig2, '--', color=cmap2(color_list[i]), lw=1)
        #ax2.plot(GC_sig, color=cmap2(color_list[i]), lw=1)
        #ax1.plot(AT_sig, color=color_list[i], lw=1)
        #ax2.plot(GC_sig, color=color_list[i], lw=1)
    ax1.set_ylabel('AA/AT/TA/TT freqeuncy', color='r', fontsize=8)
    ax2.set_ylabel('CC/CG/GC/GG freqeuncy', color='b', fontsize=8)
    ax1.set_xlabel("Super Helical Location", fontsize=8)
    ax1.tick_params('y', colors='r', labelsize=8)
    ax2.tick_params('y', colors='b', labelsize=8)
    line_list = [147/2 - i*20 for i in range(4)] + [147/2 + i*20 for i in range(1, 4)]
    for line in line_list:
        plt.axvline(x=line, color='k', linestyle='--', alpha=0.25)
    plt.xticks([147/2 + 10*i for i in range(-7, 8, 2)], [str(10*i) for i in range(-7, 8, 2)], fontsize=5)
    #plt.ylabel("Relative frequency")
    #plt.legend()
    #plt.ylim([0.22, 0.28])
    #plt.savefig('ATGCperiod_' + "slide" + '.png')
    #plt.savefig('ATGCperiod_' + "slide" + '.svg', format='svg', bbox_inches='tight')
    plt.title("Dinucleotide periodicity", fontsize=8)
    #plt.show()
    plt.close()

    sys.exit(1)


    group_feature = []
    for i in range(len(group_keys)):
        key_list = group_keys[i]
        feature = {'GC':[], 'CpG':[], 'TpA':[], 'polyG':[], 'polyA':[]}
        for key in key_list:
            seq = key_slider1[key].seq[NCP_len/2:ref_length-NCP_len/2]
            GC = GC_content(seq)
            CpG = get_dincount(seq, din='CG')
            TpA = get_dincount(seq, din='TA')
            Alen_pos = poly_score(seq, nts='AT', pos=True)
            polyA_count = 0
            for Alen, pos in Alen_pos.items():
                if Alen < 3:
                    continue
                polyA_count += len(Alen_pos[Alen])
            Glen_pos = poly_score(seq, nts='GC', pos=True)
            polyG_count = 0
            for Glen, pos in Glen_pos.items():
                if Glen < 3:
                    continue
                polyG_count += len(Glen_pos[Glen])
            feature['GC'].append(GC)
            feature['CpG'].append(CpG)
            feature['TpA'].append(TpA)
            feature['polyA'].append(polyA_count)
            feature['polyG'].append(polyG_count)
        group_feature.append(feature)

        fig = plt.figure()
        featname_title = {'GC':'GC content', 'CpG':'CpG count', 'TpA':'TpA count', 'polyA':'Poly-A/T count', 'polyG':'Poly-G/C count'}
        featname_list = group_feature[0].keys()
        for i in range(len(featname_list)):
            featname = featname_list[i]
            box_data = []
            for j in range(gnum):
                data = group_feature[j][featname]
                box_data.append(data)
            plt.subplot(1, len(featname_list), i+1)
            box = plt.boxplot(box_data, positions=range(gnum), widths=0.3, notch=False, patch_artist=True)
            for k in range(gnum):
                box['boxes'][k].set_facecolor(color_list[k])
            #plt.legend([box["boxes"][i] for i in range(gnum)], ['Cluster '+ str(i+1) for i in range(gnum)], loc='best')
            plt.xticks(range(gnum), ['Cluster '+ str(i+1) for i in range(gnum)])
            if featname == 'GC':
                plt.ylabel("%")
            else:
                plt.ylabel("Count")
            plt.title(featname_title[featname])
        plt.suptitle("Sequence features")
        #plt.show()
        plt.close()

#sys.exit(1)

# dinucleotide pair correlation analysis
if False:
    def din_pair_count (seq, max_dist=50):
        assert len(seq) >= max_dist
        pair_sig = {}
        for i in range(len(seq)-2):
            for j in range(i+1, min(len(seq)-1, i+max_dist+1)):
                din1 = seq[i:i+2]
                din2 = seq[j:j+2]
                dist = j - i
                pair = tuple(sorted([din1, din2]))
                if pair not in pair_sig:
                    pair_sig[pair] = [0.0]*max_dist
                pair_sig[pair][dist-1] +=1
        return pair_sig

    def correlate (sig1, sig2, max_dist=sys.maxint, circular=False):
        sig1 = np.asarray(sig1) - np.mean(sig1)
        sig2 = np.asarray(sig2) - np.mean(sig2)
        dist_products = {}
        for i in range(len(sig1)):
            for j in range(i, min(len(sig1), i+max_dist+1)):
                dist = j - i
                product = sig1[i]*sig2[j]
                if dist not in dist_products:
                    dist_products[dist] = []
                dist_products[dist].append(product)
        corr_sig = [0.0]*(max(dist_products.keys())+1)
        for dist, products in dist_products.items():
            corr_sig[dist] = np.mean(products)
        return corr_sig

    def FFT (sig):
        N = len(sig)
        sig_ft = fft(sig)[1:N/2]
        periods = [float(N)/k for k in range(1, N/2)]
        amplts = np.abs(sig_ft)/float(N)
        phases = np.arctan2(sig_ft.imag, sig_ft.real) / np.pi # unit of pi
        shifts = np.asarray([(phases[k-1] * N) / (2*np.pi*k) for k in range(1, N/2)]) /np.pi # unit of pi
        return periods, amplts, phases

    m1 = EnModel.EnergyModel(key_slider1)
    m1.report(MM_orders=[1], Kmer_k_b=None, PolyA_b=False, GC_b=False, PairCorr_k_dist=False)
    m2 = EnModel.EnergyModel(key_slider2)
    m2.report(MM_orders=[1], Kmer_k_b=None, PolyA_b=False, GC_b=False, PairCorr_k_dist=False)

    GCcontent_din = sorted([(GC_content(din), din) for din in LinModel.all_path(2, 'ATCG')], cmp=tuple_cmp)
    all_din = [din for GCcontent, din in GCcontent_din]

    din_sig_list = []
    pair_corr_list = []
    pair_best_list = []
    for pos_din_freq in [m1.freq['MM1'], m2.freq['MM1']]:
        din_sig = {}
        for i in range(len(pos_din_freq)):
            for din in pos_din_freq[i]:
                freq = pos_din_freq[i][din]
                if din not in din_sig:
                    din_sig[din] = []
                din_sig[din].append(freq)
        din_sig_list.append(din_sig)

        pair_corr = {}
        for i in range(len(all_din)):
            for j in range(i, len(all_din)):
                din1, din2 = all_din[i], all_din[j]
                sig1, sig2 = din_sig[din1], din_sig[din2]
                corr_sig = correlate(sig1, sig2, max_dist=49)
                pair = tuple(sorted([din1, din2]))
                pair_corr[pair] = corr_sig
                #fig = plt.figure()
                #plt.title(din1 + '-' + din2)
                #plt.plot(corr_sig)
                #plt.show()
                #plt.close()
        pair_corr_list.append(pair_corr)

        pair_best = {}
        for pair in pair_corr:
            corr = pair_corr[pair]
            periods, amplts, phases = FFT(corr)
            #idx = sorted([(amplts[i], i) for i in range(len(amplts))], reverse=True)[0][1]
            idx = 4
            print periods[idx], amplts[idx], phases[idx]
            pair_best[pair] = (periods[idx], amplts[idx], abs(phases[idx]))
            #fig = plt.figure()
            #plt.title(pair[0] + '-' + pair[1])
            #plt.plot(periods, phases)
            #plt.show()
            #plt.close()
        pair_best_list.append(pair_best)

    amplt_img = np.zeros((len(all_din), len(all_din)))
    amplt_img[:] = np.nan
    phase_img = np.zeros((len(all_din), len(all_din)))
    phase_img[:] = np.nan
    for i in range(len(all_din)):
        for j in range(i, len(all_din)):
            din1, din2 = all_din[i], all_din[j]
            pair = tuple(sorted([din1, din2]))
            period1, amplt1, phase1 = pair_best_list[0][pair]
            period2, amplt2, phase2 = pair_best_list[1][pair]
            amplt_img[i][j] = amplt1
            amplt_img[j][i] = amplt2
            phase_img[i][j] = phase1
            phase_img[j][i] = phase2

    fig = plt.figure(figsize=(3.4, 2.4))
    plt.imshow(amplt_img, cmap='afmhot_r', aspect='auto')
    plt.xticks(range(len(all_din)), all_din, rotation=45, fontsize=6)
    plt.yticks(range(len(all_din)), all_din, fontsize=6)
    plt.title("Amplitude", fontsize=8)
    cbar = plt.colorbar(shrink=0.6, aspect=20, format='%.0e')
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(5)
    #cbar.ax.set_yticklabels([round(min(Z), 2), round(max(Z), 2)], fontsize=5)
    #cbar.ax.set_ylabel('KL-divergence', rotation=-90, va="bottom", fontsize=5, labelpad=-5)
    plt.savefig("amplt.svg", format="svg", bbox_inches='tight')
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(3.4, 2.4))
    plt.imshow(phase_img, cmap='Spectral', vmin=0, vmax=1, aspect='auto')
    plt.xticks(range(len(all_din)), all_din, rotation=45, fontsize=6)
    plt.yticks(range(len(all_din)), all_din, fontsize=6)
    plt.title("Phase", fontsize=8)
    cbar = plt.colorbar(shrink=0.6, aspect=20)
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(5)
    plt.savefig("phase.svg", format="svg", bbox_inches='tight')
    #plt.show()
    plt.close()


    sys.exit(1)


    for din, sig in din_sig.items():
        fig = plt.figure()
        plt.title(din)
        plt.plot(sig)
        plt.show()
        plt.close()




    #pair_corr = m1.freq['PairCorr']
    #for pair in pair_corr:
    #    din1, din2 = pair
    #    corr = pair_corr[pair]
    #    fig = plt.figure()
    #    plt.plot(corr)
    #    plt.title(din1 + '-' + din2)
    #    plt.show()
    #    plt.close()


    sys.exit(1)

# make predictive model
if True:
    def read_ref (ref_fname):
        id_seq = {}
        for line in open(ref_fname):
            line = line.strip()
            if line.startswith(">"):
                id = line[1:]
                continue
            seq = line
            assert id not in id_seq
            id_seq[id] = seq
        return id_seq
    key_seq = read_ref ("plusonelib.ref")

    if False:
        # Markov Model
        M1 = LinModel.LinSliderModel(key_slider1)
        Markov1 = HMM.MarkovModel(1)
        even_sample = M1._even_sampling()
        bias_sample = M1._bias_sampling()
        Markov1.train({i:even_sample[i] for i in range(len(even_sample))}, bias_sample)

        # check Widom 601
        X, Y = [], []
        for key in key_slider1:
            psig = analysis.norm(key_slider1[key].dyadmap)
            seq = key_seq[key]
            assert len(seq) == 225
            pred_psig = Markov1.single_prob_profile (seq, background=False)

            fig = plt.figure()
            plt.plot(psig, 'k-')
            plt.plot(pred_psig, 'r-')
            plt.title(key)
            plt.show()
            plt.close()

            X += psig[NCP_len/2:len(psig)-NCP_len/2]
            Y += pred_psig[NCP_len/2:len(pred_psig)-NCP_len/2]

        print analysis.get_corr(X, Y)

        fig = plt.figure()
        plt.plot(X, Y, '.')
        plt.xlabel("Experiment")
        plt.ylabel("Prediction")
        plt.show()
        plt.close()


        sys.exit(1)


    # Exp vs Pred
    order = 0
    M1 = LinModel.LinSliderModel(key_slider1)
    M2 = LinModel.LinSliderModel(key_slider2)

    #order=1
    M1.train(MM_orders=[order], Kmer_k_b=[], PolyA_b=False, GC_b=False, Harmonic=True, graph=False)
    M2.train(MM_orders=[order], Kmer_k_b=[], PolyA_b=False, GC_b=False, Harmonic=True, graph=False)

    X, Y = [], []
    pred_key_slider1 = {}
    for key in key_slider1:
        psig = analysis.norm(key_slider1[key].dyadmap)
        seq = key_seq[key]
        assert len(seq) == 225
        pred_psig = M1._predict(seq)

        fig = plt.figure()
        plt.plot(psig, 'k-')
        plt.plot([0.0]*(NCP_len/2) + pred_psig + [0.0]*(NCP_len/2), 'r-')
        plt.show()
        plt.close()

        newslider = copy.deepcopy(key_slider1[key])
        newslider.dyadmap = [0.0]*(NCP_len/2) + pred_psig + [0.0]*(NCP_len/2)
        pred_key_slider1[key] = newslider

        X += psig[NCP_len/2:len(psig)-NCP_len/2]
        #Y += pred_psig[NCP_len/2:len(psig)-NCP_len/2]
        Y += pred_psig

    #tlen = 225
    #random.seed(0)
    #sample_list = random.sample(key_slider1.keys(), 30)

    #graph_final.plot_map(key_slider1, Slider.KDE, ids=sample_list, thickness=[8, 2, 2], cmap='YlGnBu', save=True, note='_before', xticks=[[i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)]], figscale=150, fontsize=5)

    #graph_final.plot_map(pred_key_slider1, Slider.KDE, ids=sample_list, thickness=[8, 2, 2], cmap='YlGnBu', save=True, note='_after', xticks=[[i+tlen/2 for i in range(-100, 101, 20)], [str(i) for i in range(-100, 101, 20)]], figscale=150, fontsize=5)


    

    fig = plt.figure(figsize=(2, 1.8))
    plt.plot(X, Y, 'k.', markersize=2, alpha=0.1)
    sns.kdeplot(X, Y, cmap='Reds', shade=False, shade_lowest=False, zorder=100, linewidth=0.1, alpha=0.5)
    plt.xlabel("Experiment", fontsize=6)
    plt.ylabel("Prediction", fontsize=6)
    plt.title("Chd1 data", fontsize=8)
    plt.gca().tick_params(axis='both', which='major', labelsize=5)
    plt.gca().tick_params(axis='both', which='minor', labelsize=5)
    plt.xlim([0, 0.03])
    plt.ylim([0, 0.03])
    #plt.savefig("ExVsPre.png", dpi=500, bbox_inches='tight')
    #plt.show()
    plt.close()

    sys.exit(1)

    # linear energy model
    k_fold = 10
    mname_list1, mname_list2 = [], []
    corrs_list1, corrs_list2 = [], []
    M1 = LinModel.LinSliderModel(key_slider1)
    M2 = LinModel.LinSliderModel(key_slider2)

    # 0-th order Markov
    order=0
    M1.train(MM_orders=[order], Kmer_k_b=[], PolyA_b=False, GC_b=False, Harmonic=False, graph=False)
    M2.train(MM_orders=[order], Kmer_k_b=[], PolyA_b=False, GC_b=False, Harmonic=False, graph=False)
    mname_list1.append('Mono-nt steps')
    mname_list2.append('Mono-nt steps')
    corrs_list1.append(M1.corr)
    corrs_list2.append(M2.corr)

    # 1-th order Markov
    order=1
    M1.train(MM_orders=[order], Kmer_k_b=[], PolyA_b=False, GC_b=False, Harmonic=False, graph=False)
    M2.train(MM_orders=[order], Kmer_k_b=[], PolyA_b=False, GC_b=False, Harmonic=False, graph=False)
    mname_list1.append('Di-nt steps')
    mname_list2.append('Di-nt steps')
    corrs_list1.append(M1.corr)
    corrs_list2.append(M2.corr)

    # 1-th order Markov + 5-mers
    order=1
    M1.train(MM_orders=[order], Kmer_k_b=[5, 1], PolyA_b=False, GC_b=False, Harmonic=False, graph=False)
    M2.train(MM_orders=[order], Kmer_k_b=[5, 1], PolyA_b=False, GC_b=False, Harmonic=False, graph=False)
    mname_list1.append('Di-nt + 5-mers')
    mname_list2.append('Di-nt + 5-mers')
    corrs_list1.append(M1.corr)
    corrs_list2.append(M2.corr)

    # 1-th order Markov + 5-mers + Harmonics
    order=1
    M1.train(MM_orders=[order], Kmer_k_b=[5, 1], PolyA_b=False, GC_b=False, Harmonic=True, graph=False)
    M2.train(MM_orders=[order], Kmer_k_b=[5, 1], PolyA_b=False, GC_b=False, Harmonic=True, graph=False)
    mname_list1.append('Di-nt + 5-mers + Harmo')
    mname_list2.append('Di-nt + 5-mers + Harmo')
    corrs_list1.append(M1.corr)
    corrs_list2.append(M2.corr)

    # check correlation
    color_list = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

    fig = plt.figure(figsize=(2.8, 2))
    for i in range(len(corrs_list1)):
        bp1 = plt.boxplot(corrs_list1[i], positions=[i-0.18],
                          patch_artist=True, showfliers=False, widths=[0.2])
        for patch in bp1['boxes']:
            patch.set_facecolor('tab:orange')
        for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bp1[element], color='black', lw=0.5)

        bp2 = plt.boxplot(corrs_list2[i], positions=[i+0.18],
                          patch_artist=True, showfliers=False, widths=[0.2])
        for patch in bp2['boxes']:
            patch.set_facecolor('tab:green')
        for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bp2[element], color='black', lw=0.5)

    #plt.boxplot(corrs_list, positions=range(len(mname_list)))
    plt.xticks(range(len(corrs_list1)), mname_list1, rotation=45, fontsize=8,  ha="right", rotation_mode="anchor")
    plt.gca().tick_params('y', labelsize=6)
    plt.ylabel("Pearson correlation", fontsize=8)
    plt.title("10-fold cross validation", fontsize=10)
    plt.xlim([-0.5, len(mname_list1)-0.5])
    plt.legend([bp1["boxes"][0], bp2["boxes"][0]], ['HS', 'Chd1'], frameon=True, loc='lower right', fontsize=8)
    plt.savefig('Pearson_model.svg', format='svg', bbox_inches='tight')
    #plt.show()
    plt.close()

    

    


    sys.exit(1)

    # Markov Model
    #Markov1 = HMM.MarkovModel(1)
    #even_sample = M1._even_sampling()
    #bias_sample = M1._bias_sampling()
    #Markov1.train({i:even_sample[i] for i in range(len(even_sample))}, bias_sample)

    # check Widom 601
    seq = "ATCCGACTGGCACCGGCAAGGTCGCTGTTCGCCACATGCGCAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCAGGGCGTCCTCGTATAGGGTCCATCACATAAGGGATGAACT"
    assert len(seq) == 225
    p_sig1 = analysis.norm(Control2['601'].dyadmap)
    pred_psig = M1._predict(seq)

    fig = plt.figure()
    plt.plot(p_sig1, 'k-')
    plt.plot(pred_psig, 'r-')
    #plt.show()
    plt.close()

    #exit(1)
    

    X, Y = [], []
    for key in ['97']:
        psig = analysis.norm(key_slider1[key].dyadmap)
        seq = key_seq[key]
        assert len(seq) == 225
        pred_psig = M1._predict(seq)
        #pred_psig = Markov1.single_prob_profile (seq)

        X += psig[NCP_len/2:len(psig)-NCP_len/2+1]
        Y += pred_psig[NCP_len/2:len(psig)-NCP_len/2+1]

        fig = plt.figure()
        plt.plot(psig, 'k-')
        plt.plot(pred_psig, 'r-')
        #plt.plot(X, Y,'.')
        #plt.xlabel("Experiment")
        #plt.ylabel("Prediction")
        plt.title(key)
        plt.show()
        plt.close()

    fig = plt.figure()
    plt.plot(X, Y, '.')
    plt.xlabel("Experiment")
    plt.ylabel("Prediction")
    plt.show()
    plt.close()

    print analysis.get_corr(X, Y)

    sys.exit(1)

# predictive model
#key_slider1 = name_key_slider["plusonelib_new_0"]
#key_slider2 = name_key_slider["plusonelib_new_30"]

#mask_idxs = range(NCP_len/2-50-2, NCP_len/2-50+2) + range(NCP_len/2+50-2, NCP_len/2+50+2)
mask_idxs=None

m1 = EnModel.EnergyModel(key_slider1, mask_idxs=mask_idxs)
m1.report(MM_orders=[4], Kmer_k_b=[5, 1], PolyA_b=False, GC_b=False)
#m1.train(MM_orders=[1], Kmer_k_b=[5,3], PolyA_b=False, GC_b=False, Harmonic=True, k_fold=3, adjust=True, graph=False)

m2 = EnModel.EnergyModel(key_slider2, mask_idxs=mask_idxs)
m2.report(MM_orders=[4], Kmer_k_b=[5, 1], PolyA_b=False, GC_b=False)
#m2.train(MM_orders=[1], Kmer_k_b=[5,3], PolyA_b=False, GC_b=False, Harmonic=True, k_fold=3, adjust=True, graph=False)


kmers = LinModel.all_path(5, 'ATCG')

# check all 5-mers regardless of position
if False:
    X, Y = [], []
    C = []
    picks = []
    for k in range(len(kmers)):
        kmer = kmers[k]
        freq1 = m1.freq['Kmer0'][kmer]
        freq2 = m2.freq['Kmer0'][kmer]
        GC = GC_content(kmer)
        polyA = poly_score (kmer, nts='AT', pos=False)
        poly = poly_score (kmer, nts='ATGC', pos=False)
        CpG = get_dincount(kmer, din='CG')
        TpA = get_dincount(kmer, din='TA')
        X.append(freq1)
        Y.append(freq2)
        C.append(polyA)
        #X.append(coeff2 - coeff1)
        #Y.append(polyA)
        #C.append(GC)
        if kmer == 'AAAAA':
            picks.append((kmer, freq1, freq2))
            #picks.append((kmer, coeff2-coeff1, polyA))
        elif kmer == 'GGGGG':
            picks.append((kmer, freq1, freq2))
            #picks.append((kmer, coeff2-coeff1, polyA))

    fig = plt.figure()
    plt.scatter(X, Y, c=C, s=3, alpha=0.5, cmap='jet')
    ax = plt.gca()
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    plt.plot(lims, lims, 'k--', alpha=0.5)
    for kmer, x, y in picks:
        plt.plot(x, y, 'kx')
        plt.text(x, y, kmer)
    plt.xlabel("Before")
    plt.ylabel("After")
    #plt.xlim([0.85, 1.15])
    #plt.ylim([0.85, 1.15])
    plt.title("5-mers frequency")
    plt.colorbar()
    plt.show()
    plt.close()



# check 5 mers on every location

kmer_freq_list1 = m1.freq['MM4']
kmer_freq_list2 = m2.freq['MM4']

kmer_order_list1 = []
kmer_order_list2 = []

for i in range(len(kmer_freq_list1)):
    kmer_freq1 = kmer_freq_list1[i]
    kmer_freq2 = kmer_freq_list2[i]
    freqkmer1 = sorted([[freq, kmer] for kmer, freq in kmer_freq1.items()], reverse=True)
    freqkmer2 = sorted([[freq, kmer] for kmer, freq in kmer_freq2.items()], reverse=True)
    kmer_order1 = { freqkmer1[k][1]:k+1 for k in range(len(freqkmer1)) }
    kmer_order2 = { freqkmer2[k][1]:k+1 for k in range(len(freqkmer2)) }
    kmer_order_list1.append(kmer_order1)
    kmer_order_list2.append(kmer_order2)

#for kmer in kmers:
#    freq_profile1, freq_profile2 = [], []
#    for i in range(len(kmer_freq_list1)):
#        freq_profile1.append(kmer_freq_list1[i][kmer])
#        freq_profile2.append(kmer_freq_list2[i][kmer])
#    fig = plt.figure()
#    plt.plot(freq_profile1, label='before')
#    plt.plot(freq_profile2, label='after')
#    plt.title(kmer)
#    plt.legend()
#    plt.show()
#    plt.close()

#sys.exit(1)

polyA_profile1, polyA_profile2, polyA_profile3 = [], [], []
polyG_profile1, polyG_profile2, polyG_profile3 = [], [], []

A = []
B1, B2, B3 = [], [], []
S = []

for i in range(len(kmer_freq_list1)):
    kmer_freq1 = kmer_freq_list1[i]
    kmer_freq2 = kmer_freq_list2[i]

    kmer_order1 = kmer_order_list1[i]
    kmer_order2 = kmer_order_list2[i]
    
    X, Y = [], []
    C = []
    picks = []

    mean_polyA1, mean_polyA2 = 0.0, 0.0
    mean_polyG1, mean_polyG2 = 0.0, 0.0
    total_freq1, total_freq2 = 0.0, 0.0
    
    for k in range(len(kmers)):
        kmer = kmers[k]
        freq1 = kmer_freq1[kmer]
        freq2 = kmer_freq2[kmer]
        order1 = kmer_order1[kmer]
        order2 = kmer_order2[kmer]
        GC = GC_content(kmer)
        polyA = poly_score (kmer, nts='AT', pos=False)
        polyG = poly_score (kmer, nts='GC', pos=False)
        poly = poly_score (kmer, nts='ATGC', pos=False)
        CpG = get_dincount(kmer, din='CG')
        TpA = get_dincount(kmer, din='TA')
        X.append(freq1)
        Y.append(freq2)
        mean_polyA1 += freq1*polyA
        mean_polyA2 += freq2*polyA
        mean_polyG1 += freq1*polyG
        mean_polyG2 += freq2*polyG
        total_freq1 += freq1
        total_freq2 += freq2
        #C.append(polyA)
        #X.append(freq2 - freq1)
        #Y.append(polyA)
        C.append(GC)
        if kmer == 'AAAAA':
            picks.append((kmer, freq1, freq2))
            #picks.append((kmer, freq2-freq1, polyA))
        elif kmer == 'GGGGG':
            picks.append((kmer, freq1, freq2))
            #picks.append((kmer, freq2-freq1, polyA))

        A.append(i)
        #B1.append(freq1)
        #B2.append(freq2)
        #B3.append(freq2 - freq1)

        B1.append(order1)
        B2.append(order2)
        B3.append(order2 - order1)

        S.append(polyA)
        #S.append(100-GC)
            
    mean_polyA1 = mean_polyA1 / total_freq1
    mean_polyA2 = mean_polyA2 / total_freq2

    polyA_profile1.append(mean_polyA1)
    polyA_profile2.append(mean_polyA2)
    polyA_profile3.append(mean_polyA2 - mean_polyA1)

    mean_polyG1 = mean_polyG1 / total_freq1
    mean_polyG2 = mean_polyG2 / total_freq2

    polyG_profile1.append(mean_polyG1)
    polyG_profile2.append(mean_polyG2)
    polyG_profile3.append(mean_polyG2 - mean_polyG1)

    fig = plt.figure()
    plt.scatter(X, Y, c=C, s=3, alpha=0.5, cmap='jet')
    ax = plt.gca()
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    plt.plot(lims, lims, 'k--', alpha=0.5)
    for kmer, x, y in picks:
        plt.plot(x, y, 'kx')
        plt.text(x, y, kmer)
    plt.xlabel("Before")
    plt.ylabel("After")
    #plt.xlim([0.85, 1.15])
    #plt.ylim([0.85, 1.15])
    plt.title("5-mers frequency " + str(i))
    cbar = plt.colorbar()
    cbar.set_label("Poly-A/T score")
    #plt.show()
    plt.close()

    freq_kmer1 = sorted([ [freq, kmer] for kmer, freq in kmer_freq1.items()])
    freq_kmer2 = sorted([ [freq, kmer] for kmer, freq in kmer_freq2.items()])
    freq_kmer3 = sorted([ [kmer_freq2[kmer] - kmer_freq1[kmer], kmer] for kmer in kmers ])


    f = open(str(i) + "_temp.txt", 'w')

    for value, kmer in freq_kmer3:
        print >> f, value, kmer
    f.close()
    
    X1, X2, X3 = [], [], []
    Y1, Y2, Y3 = [], [], []
    Z1, Z2, Z3 = [], [], []
    for freq, kmer in freq_kmer1:
        X1.append(kmer)
        Y1.append(freq)
        #Z1.append(GC_content(kmer))
        Z1.append(poly_score(kmer, nts='AT', pos=False))

    for freq, kmer in freq_kmer2:
        X2.append(kmer)
        Y2.append(freq)
        #Z2.append(GC_content(kmer))
        Z2.append(poly_score(kmer, nts='AT', pos=False))

    picks = []
    for k in range(len(freq_kmer3)):
        freq, kmer = freq_kmer3[k]
        if kmer == 'AAAAA':
            picks.append((kmer, k))
        elif kmer == 'GGGGG':
            picks.append((kmer, k))
        X3.append(kmer)
        Y3.append(freq)
        #Z3.append(GC_content(kmer))
        Z3.append(poly_score(kmer, nts='AT', pos=False))

    fig = plt.figure(figsize=(15,5))
    plt.scatter(range(len(Y1)), Y1, c=Z1, s=1, alpha=0.5, cmap='jet')
    plt.title("Heat Shift")
    plt.ylabel("Relative frequency" + str(i))
    plt.xlabel("5-mers")
    plt.colorbar()
    #plt.xticks(range(len(Y1)), X1, rotation=90)
    #plt.savefig('5mer_HS_' + "slide" + '.png', bbox_inches='tight')
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(15,5))
    plt.scatter(range(len(Y2)), Y2, c=Z2, s=1, alpha=0.5, cmap='jet')
    plt.title("Chd1 sliding " + str(i))
    plt.ylabel("Relative frequency")
    plt.xlabel("5-mers")
    plt.colorbar()
    #plt.savefig('5mer_Chd1_' + "slide" + '.png', bbox_inches='tight')
    #plt.xticks(range(len(Y1)), X1, rotation=90)
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(15,5))
    plt.scatter(range(len(Y3)), Y3, c=Z3, s=1, alpha=0.5, cmap='jet')
    xtick_pos, xtick_label = [], []
    for kmer, pos in picks:
        plt.axvline(x=pos, linestyle = '--', linewidth=1, alpha=0.5)
        xtick_pos.append(pos)
        xtick_label.append(kmer)
    plt.axhline(y=0, color='black', linestyle = '-', linewidth=1, alpha=1)
    plt.xticks(xtick_pos, xtick_label)
    plt.title("Chd1-HS " + str(i))
    plt.ylabel("frequency change")
    plt.xlabel("5-mers")
    plt.colorbar()
    #plt.savefig('5mer_change_' + "slide" + '.png', bbox_inches='tight')
    #plt.xticks(range(len(Y1)), X1, rotation=90)
    #plt.show()
    plt.close()

fig = plt.figure()
plt.plot(polyA_profile1, 'tab:orange', linestyle='--', label='Heat shifted')
plt.plot(polyA_profile2, 'k', label='Chd1 Slided')
plt.xticks([(147-4)/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
plt.ylabel("Poly-A/T score")
#plt.ylabel("Relative frequency")
#plt.ylabel("Coefficient")
plt.legend()
#plt.ylim([0.22, 0.28])
#plt.savefig('ATGCperiod_' + "slide" + '.png')
#plt.show()
plt.close()

fig = plt.figure()
plt.plot(polyG_profile1, 'tab:orange', linestyle='--', label='Heat shifted')
plt.plot(polyG_profile2, 'k', label='Chd1 Slided')
plt.xticks([(147-4)/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
plt.ylabel("Poly-G/C score")
#plt.ylabel("Relative frequency")
#plt.ylabel("Coefficient")
plt.legend()
#plt.ylim([0.22, 0.28])
#plt.savefig('ATGCperiod_' + "slide" + '.png')
#plt.show()
plt.close()

fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(10, 5))
g=axs[0].scatter(A, B1, c=S, s=1, alpha=0.25, cmap='jet')
axs[1].scatter(A, B2, c=S, s=1, alpha=0.25, cmap='jet')
#axs[0].set_ylim([0.5,2.5])
#axs[1].set_ylim([0.5,2.5])
plt.xticks([(147-4)/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
axs[0].set_title("Heat shift", position=(0.8, 0.8))
axs[1].set_title("Chd1 sliding", position=(0.8, 0.8))
axs[0].set_ylabel("frequency")
axs[1].set_ylabel("frequency")
axs[1].set_xlabel("Super Helical Location")
fig.colorbar(g, ax=axs, location='right', shrink=0.6, alpha=1, label='AT content(%)')
#plt.tight_layout()
fig.suptitle("5-mer frequencies along SHL")
#plt.savefig('5mer_position.png', bbox_inches='tight')
plt.show()
plt.close()

fig = plt.figure(figsize=(10,5))
plt.plot(A, B1, '.', markersize=2, alpha=0.2, label='Heat shift')
plt.plot(A, B2, 'x', markersize=2, alpha=0.2, label='Chd1')
plt.xticks([(147-4)/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.title("5-mer frequency")
plt.ylabel("Normalized frequency")
plt.xlabel("Super Helical Location")
plt.ylim([0.4, 1.7])
plt.legend()
#plt.savefig('5mer_change_position.png', bbox_inches='tight')
#plt.tight_layout()
#plt.show()
plt.close()


fig = plt.figure(figsize=(10,5))
plt.scatter(A, B3, c=S, s=1, alpha=0.2, cmap='jet')
plt.xticks([(147-4)/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.title("5-mer frequency change")
plt.ylabel("Frequency change (Chd1-HS)")
plt.xlabel("Super Helical Location")
#plt.ylim([-0.7, 0.7])
cbar = plt.colorbar()
#cbar.set_label("Poly-A/T score")
cbar.set_label("AT content(%)")
cbar.set_alpha(1)
cbar.draw_all()
plt.savefig('5mer_change_position.png', bbox_inches='tight')
plt.tight_layout()
plt.show()
plt.close()




sys.exit(1)


    



GCcontent_din = sorted([(GC_content(din), din) for din in LinModel.all_path(2, 'ATCG')], cmp=tuple_cmp)
all_din = [din for GCcontent, din in GCcontent_din]

#parameters = [m1.coeff, m2.coeff]
#cmap = 'RdBu'
#vlimit = [-0.05, 0.05]

parameters = [m1.freq, m2.freq]
cmap = 'magma'
#vlimit = [None, None]
vlimit = [0.85, 1.15]

img1, img2 = [], []
for nt in all_din:
    row1, row2 = [], []
    for i in range(NCP_len-1):
        row1.append(parameters[0]['MM1'][i][nt])
        row2.append(parameters[1]['MM1'][i][nt])
    img1.append(row1)
    img2.append(row2)

for img, title in zip([img1, img2], ['Before sliding', 'After sliding']):
    fig = plt.figure()
    plt.imshow(img, aspect='auto', cmap=cmap, vmin=vlimit[0], vmax=vlimit[1])
    plt.yticks(range(len(all_din)), all_din)
    plt.xticks([NCP_len/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
    plt.title(title)
    plt.colorbar()
    plt.show()
    plt.close()

def get_ATGC_sig (freq):
    freq_MM1 = freq['MM1']
    length = len(freq_MM1)
    nts = LinModel.all_path(2, 'ATCG')
    AT_sig, GC_sig = np.zeros(length), np.zeros(length)
    for nt in nts:
        row = [ freq_MM1[i][nt] for i in range(length)]
        if nt in ['AA', 'AT', 'TA', 'TT']:
            AT_sig += np.asarray(row)
        if nt in ['GG', 'GC', 'CG', 'CC']:
            GC_sig += np.asarray(row)
    #for i in range(length):
    #    AT_sig[i] = float(AT_sig[i]) / sum(freq_MM1[i].values())
    #    GC_sig[i] = float(GC_sig[i]) / sum(freq_MM1[i].values())
    return AT_sig, GC_sig

AT_sig1, GC_sig1 = get_ATGC_sig (parameters[0])
AT_sig2, GC_sig2 = get_ATGC_sig (parameters[1])
fig = plt.figure()
plt.plot(AT_sig1, 'hotpink', label='AT-rich HS')
plt.plot(GC_sig1, 'lightskyblue', label='GC-rich HS')
plt.plot(AT_sig2, 'r--', label='AT-rich Chd1')
plt.plot(GC_sig2, 'b--', label='GC-rich Chd1')
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
#plt.ylabel("Relative frequency")
plt.ylabel("Coefficient")
plt.legend()
#plt.ylim([0.22, 0.28])
#plt.savefig('ATGCperiod_' + "slide" + '.png')
plt.show()
plt.close()

kmers = LinModel.all_path(5, 'ATCG')
for i in range(5):
    X, Y = [], []
    C = []
    picks = []
    for k in range(len(kmers)):
        kmer = kmers[k]
        coeff1 = parameters[0]['Kmer'+str(i)][kmer]
        coeff2 = parameters[1]['Kmer'+str(i)][kmer]
        GC = GC_content(kmer)
        polyA = poly_score (kmer, nts='AT', pos=False)
        poly = poly_score (kmer, nts='ATGC', pos=False)
        CpG = get_dincount(kmer, din='CG')
        TpA = get_dincount(kmer, din='TA')
        X.append(coeff1)
        Y.append(coeff2)
        C.append(poly)
        #X.append(coeff2 - coeff1)
        #Y.append(polyA)
        #C.append(GC)
        if kmer == 'AAAAA':
            picks.append((kmer, coeff1, coeff2))
            #picks.append((kmer, coeff2-coeff1, polyA))
        elif kmer == 'GGGGG':
            picks.append((kmer, coeff1, coeff2))
            #picks.append((kmer, coeff2-coeff1, polyA))

    fig = plt.figure()
    plt.scatter(X, Y, c=C, s=3, alpha=0.5, cmap='jet')
    ax = plt.gca()
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    plt.plot(lims, lims, 'k--', alpha=0.5)
    for kmer, x, y in picks:
        plt.plot(x, y, 'kx')
        plt.text(x, y, kmer)
    plt.xlabel("Before")
    plt.ylabel("After")
    #plt.xlim([0.85, 1.15])
    #plt.ylim([0.85, 1.15])
    plt.title("5-mers frequency " + str(i))
    plt.colorbar()
    plt.show()
    plt.close()

fig = plt.figure()
for para, title in zip(parameters, ['Before sliding', 'After sliding']):
    k = para['Harmonic']
    Y = [0.5*k*(x-ref_length/2)**2 for x in range(ref_length)] 
    plt.plot(range(1, ref_length+1),Y, label=title)
plt.title("Harmonic energy")
plt.xlabel("Position (bp)")
plt.ylabel("Energy (A.U.)")
plt.legend()
plt.show()
plt.close()

sample_mode = "polyA:" + "-".join([str(size) for size in [5, 12]])
sample_list = sample.sampling(name_key_slider["polyAlib_old_0_1rep"], sample_mode)

key_pslider1 = m1.predict(name_key_slider["polyAlib_old_0_1rep"], sample_list[0]+sample_list[1], scale=3.0)
graph_edit.plot_map(key_pslider1, sample_list, True, Slider.peak_signal, draw = "key", slicing=0, note='_pred1')

key_pslider2 = m2.predict(name_key_slider["polyAlib_old_0_1rep"], sample_list[0]+sample_list[1], scale=3.0)
graph_edit.plot_map(key_pslider2, sample_list, True, Slider.peak_signal, draw = "key", slicing=0, note='_pred2')




"""
# dinucleotide cross correlation analysis
def get_din_sig (seq):
    din_sig = {}
    for i in range(len(seq)-1):
        din = seq[i:i+2]
        if din not in din_sig:
            din_sig[din] = [0.0]*(len(seq)-1)
        din_sig[din][i] += 1
    return din_sig

key_din_sig = {}
for key in keys:
    seq = key_slider1[key].seq[NCP_len/2:ref_length-NCP_len/2]
    din_sig = get_din_sig(seq)
    key_din_sig[key] = din_sig
"""
seq_list = []
score_list1, score_list2 = [], []
count_list1, count_list2 = [], []
for key in keys:
    template = key_slider1[key].seq
    dyadmap1 = analysis.norm(key_slider1[key].dyadmap)
    dyadmap2 = analysis.norm(key_slider2[key].dyadmap)
    for k in range(NCP_len/2, ref_length-NCP_len/2):
        seq = template[k - NCP_len/2 : k + NCP_len/2 + 1]
        count1, count2 = dyadmap1[k], dyadmap2[k]
        score1, score2 = -np.log(count1 + 10**-15), -np.log(count2 + 10**-15)
        #score1, score2 = key_KL[key], key_KL[key]
        seq_list.append(seq)
        score_list1.append(count1)
        score_list2.append(count2)
        count_list1.append(count1)
        count_list2.append(count2)
#m1 = LinModel.SeqLinearModel(seq_list, score_list1, count_list1)
#m1.report(MM_orders=[1], Kmer_k_b=[2, 3], PolyA_b=False, GC_b=False, Harmonic=False, sym=True)
#m1.train(MM_orders=[1], Kmer_k_b=False, PolyA_b=False, GC_b=False, Harmonic=False, sym=True, k_fold=3)
#m2 = LinModel.SeqLinearModel(seq_list, score_list2, count_list2)
#m2.report(MM_orders=[1], Kmer_k_b=[2, 3], PolyA_b=False, GC_b=False, Harmonic=False, sym=True)
#m2.train(MM_orders=[1], Kmer_k_b=False, PolyA_b=False, GC_b=False, Harmonic=False, sym=True, k_fold=3)

def get_ATGC_sig (freq):
    freq_MM1 = freq['MM1']
    length = len(freq_MM1)
    nts = LinModel.all_path(2, 'ATCG')
    AT_sig, GC_sig = np.zeros(length), np.zeros(length)
    for nt in nts:
        row = [ freq_MM1[i][nt] for i in range(length)]
        if nt in ['AA', 'AT', 'TA', 'TT']:
            AT_sig += np.asarray(row)
        if nt in ['GG', 'GC', 'CG', 'CC']:
            GC_sig += np.asarray(row)
    #for i in range(length):
    #    AT_sig[i] = float(AT_sig[i]) / sum(freq_MM1[i].values())
    #    GC_sig[i] = float(GC_sig[i]) / sum(freq_MM1[i].values())
    return AT_sig, GC_sig

AT_sig1, GC_sig1 = get_ATGC_sig (m1.coeff)
AT_sig2, GC_sig2 = get_ATGC_sig (m2.coeff)
fig = plt.figure()
plt.plot(AT_sig1, 'hotpink', label='AT-rich HS')
plt.plot(GC_sig1, 'lightskyblue', label='GC-rich HS')
plt.plot(AT_sig2, 'r--', label='AT-rich Chd1')
plt.plot(GC_sig2, 'b--', label='GC-rich Chd1')
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
plt.ylabel("Relative frequency")
plt.legend()
#plt.ylim([0.22, 0.28])
#plt.savefig('ATGCperiod_' + "slide" + '.png')
plt.show()
plt.close()

kmers = LinModel.all_path(2, 'ATCG')

for i in range(2):
    X, Y = [], []
    C = []
    picks = []
    for k in range(len(kmers)):
        kmer = kmers[k]
        freq1 = m1.freq['Kmer'+str(i)][kmer]
        freq2 = m2.freq['Kmer'+str(i)][kmer]
        GC = GC_content(kmer)
        polyA = poly_score (kmer, nts='AT', pos=False)
        poly = poly_score (kmer, nts='ATGC', pos=False)
        CpG = get_dincount(kmer, din='CG')
        TpA = get_dincount(kmer, din='TA')
        X.append(freq1)
        Y.append(freq2)
        C.append(polyA)
        #X.append(freq2 - freq1)
        #Y.append(polyA)
        #C.append(GC)
        if kmer == 'AAAAA':
            picks.append((kmer, freq1, freq2))
            #picks.append((kmer, freq2-freq1, polyA))
        elif kmer == 'GGGGG':
            picks.append((kmer, freq1, freq2))
            #picks.append((kmer, freq2-freq1, polyA))


    fig = plt.figure()
    plt.scatter(X, Y, c=C, s=3, alpha=0.5, cmap='jet')
    ax = plt.gca()
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    plt.plot(lims, lims, 'k--', alpha=0.5)
    for kmer, x, y in picks:
        plt.plot(x, y, 'kx')
        plt.text(x, y, kmer)
    plt.xlabel("Before")
    plt.ylabel("After")
    #plt.xlim([0.85, 1.15])
    #plt.ylim([0.85, 1.15])
    plt.title("5-mers frequency " + str(i))
    plt.colorbar()
    plt.show()
    plt.close()

    freq_kmer1 = sorted([ [freq, kmer] for kmer, freq in m1.freq['Kmer'+str(i)].items()])
    freq_kmer2 = sorted([ [freq, kmer] for kmer, freq in m2.freq['Kmer'+str(i)].items()])
    freq_kmer3 = sorted([ [m2.freq['Kmer'+str(i)][kmer] - m1.freq['Kmer'+str(i)][kmer], kmer] for kmer in m1.freq['Kmer'+str(i)]])

    f = open(str(i) + "_temp.txt", 'w')

    for value, kmer in freq_kmer3:
        print >> f, value, kmer
    f.close()
    
    X1, X2, X3 = [], [], []
    Y1, Y2, Y3 = [], [], []
    Z1, Z2, Z3 = [], [], []
    for freq, kmer in freq_kmer1:
        X1.append(kmer)
        Y1.append(freq)
        #Z1.append(GC_content(kmer))
        Z1.append(poly_score(kmer, nts='AT', pos=False))

    for freq, kmer in freq_kmer2:
        X2.append(kmer)
        Y2.append(freq)
        #Z2.append(GC_content(kmer))
        Z2.append(poly_score(kmer, nts='AT', pos=False))

    picks = []
    for k in range(len(freq_kmer3)):
        freq, kmer = freq_kmer3[k]
        if kmer == 'AAAAA':
            picks.append((kmer, k))
        elif kmer == 'GGGGG':
            picks.append((kmer, k))
        X3.append(kmer)
        Y3.append(freq)
        #Z3.append(GC_content(kmer))
        Z3.append(poly_score(kmer, nts='ATGC', pos=False))

    fig = plt.figure(figsize=(15,5))
    plt.scatter(range(len(Y1)), Y1, c=Z1, s=1, alpha=0.5, cmap='jet')
    plt.title("Heat Shift")
    plt.ylabel("Relative frequency")
    plt.xlabel("5-mers")
    plt.colorbar()
    #plt.xticks(range(len(Y1)), X1, rotation=90)
    #plt.savefig('5mer_HS_' + "slide" + '.png', bbox_inches='tight')
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(15,5))
    plt.scatter(range(len(Y2)), Y2, c=Z2, s=1, alpha=0.5, cmap='jet')
    plt.title("Chd1 sliding")
    plt.ylabel("Relative frequency")
    plt.xlabel("5-mers")
    plt.colorbar()
    #plt.savefig('5mer_Chd1_' + "slide" + '.png', bbox_inches='tight')
    #plt.xticks(range(len(Y1)), X1, rotation=90)
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(15,5))
    plt.scatter(range(len(Y3)), Y3, c=Z3, s=1, alpha=0.5, cmap='jet')
    xtick_pos, xtick_label = [], []
    for kmer, pos in picks:
        plt.axvline(x=pos, linestyle = '--', linewidth=1, alpha=0.5)
        xtick_pos.append(pos)
        xtick_label.append(kmer)
    plt.axhline(y=0, color='black', linestyle = '-', linewidth=1, alpha=1)
    plt.xticks(xtick_pos, xtick_label)
    plt.title("Chd1-HS")
    plt.ylabel("frequency change")
    plt.xlabel("5-mers")
    plt.colorbar()
    #plt.savefig('5mer_change_' + "slide" + '.png', bbox_inches='tight')
    #plt.xticks(range(len(Y1)), X1, rotation=90)
    plt.show()
    plt.close()

    
sys.exit(1)

group_freq = m1.spectrum(MM_orders=[1], Kmer_k_b=[5, 3], PolyA_b=False, GC_b=False, Harmonic=False, gnum=gnum, norm=True)
fig = plt.figure()
for i in range(gnum):
    freq = group_freq[i]
    AT_sig, GC_sig = get_ATGC_sig (freq)
    #fig = plt.figure()
    plt.plot(AT_sig, color=color_list[i], label='Cluster ' + str(i+1))
    #plt.plot(GC_sig, color=color_list[i], linestyle='--', label='Cluster')
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
plt.ylabel("Relative frequency")
plt.title("AT-rich Dinucleotide")
plt.legend()
plt.ylim([0.21, 0.29])
#plt.savefig('ATGCperiod_' + "slide" + '.png')
plt.show()
plt.close()

fig = plt.figure()
for i in range(gnum):
    freq = group_freq[i]
    AT_sig, GC_sig = get_ATGC_sig (freq)
    #fig = plt.figure()
    #plt.plot(AT_sig, color=color_list[i], label='Cluster ' + str(i+1))
    plt.plot(GC_sig, color=color_list[i], label='Cluster ' +str(i+1))
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
plt.ylabel("Relative frequency")
plt.title("GC-rich Dinucleotide")
plt.legend()
plt.ylim([0.21, 0.29])
#plt.savefig('ATGCperiod_' + "slide" + '.png')
plt.show()
plt.close()

for i in range(gnum):
    kmers = LinModel.all_path(5, 'ATCG')

    for i in range(2):
        X, Y = [], []
        C = []
        picks = []
        for k in range(len(kmers)):
            kmer = kmers[k]
            freq1 = m1.freq['Kmer'+str(i)][kmer]
            freq2 = m2.freq['Kmer'+str(i)][kmer]
            GC = GC_content(kmer)
            polyA = poly_score (kmer, nts='AT', pos=False)
            poly = poly_score (kmer, nts='ATGC', pos=False)
            CpG = get_dincount(kmer, din='CG')
            TpA = get_dincount(kmer, din='TA')
            X.append(freq1)
            Y.append(freq2)
            C.append(GC)
            #X.append(freq2 - freq1)
            #Y.append(polyA)
            #C.append(GC)
            if kmer == 'AAAAA':
                picks.append((kmer, freq1, freq2))
                #picks.append((kmer, freq2-freq1, polyA))
            elif kmer == 'GGGGG':
                picks.append((kmer, freq1, freq2))
                #picks.append((kmer, freq2-freq1, polyA))


        fig = plt.figure()
        plt.scatter(X, Y, c=C, s=3, alpha=0.5, cmap='jet')
        ax = plt.gca()
        lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
        plt.plot(lims, lims, 'k--', alpha=0.5)
        for kmer, x, y in picks:
            plt.plot(x, y, 'kx')
            plt.text(x, y, kmer)
        plt.xlabel("Before")
        plt.ylabel("After")
        #plt.xlim([0.85, 1.15])
        #plt.ylim([0.85, 1.15])
        plt.title("5-mers frequency " + str(i))
        plt.colorbar()
        plt.show()
        plt.close()








sys.exit(1)

group_model1 = []
group_model2 = []

for i in range(gnum):
    seq_list = []
    score_list1, score_list2 = [], []
    count_list1, count_list2 = [], []
    for key in group_keys[i]:
        template = key_slider1[key].seq
        dyadmap1 = analysis.norm(key_slider1[key].dyadmap)
        dyadmap2 = analysis.norm(key_slider2[key].dyadmap)
        for k in range(NCP_len/2, ref_length-NCP_len/2):
            seq = template[k - NCP_len/2 : k + NCP_len/2 + 1]
            count1, count2 = dyadmap1[k], dyadmap2[k]
            score1, score2 = -np.log(count1 + 10**-15), -np.log(count2 + 10**-15)
            seq_list.append(seq)
            score_list1.append(score1)
            score_list2.append(score2)
            #count_list1.append(count1)
            #count_list2.append(count2)
            count_list1.append(0.01)
            count_list2.append(0.01)
    m1 = LinModel.SeqLinearModel(seq_list, score_list1, count_list1)
    m1.report(MM_orders=[1], Kmer_k_b=[5, 1], PolyA_b=False, GC_b=False, Harmonic=False, sym=True, norm=False)
    m2 = LinModel.SeqLinearModel(seq_list, score_list2, count_list2)
    m2.report(MM_orders=[1], Kmer_k_b=[5, 1], PolyA_b=False, GC_b=False, Harmonic=False, sym=True, norm=False)

    group_model1.append(m1)
    group_model2.append(m2)

   
    freq_MM1_1 = m1.freq['MM1']
    freq_MM1_2 = m2.freq['MM1']

    length = len(freq_MM1_1)
    nts = LinModel.all_path(2, 'ATCG')
    AT_row1, GC_row1 = np.zeros(length), np.zeros(length)
    AT_row2, GC_row2 = np.zeros(length), np.zeros(length)
    for nt in nts:
        row = [ freq_MM1_1[i][nt] for i in range(length)]
        if nt in ['AA', 'AT', 'TA', 'TT']:
            AT_row1 += np.asarray(row)/len(nts)
        if nt in ['GG', 'GC', 'CG', 'CC']:
            GC_row1 += np.asarray(row)/len(nts)
    for nt in nts:
        row = [ freq_MM1_2[i][nt] for i in range(length)]
        if nt in ['AA', 'AT', 'TA', 'TT']:
            AT_row2 += np.asarray(row)/len(nts)
        if nt in ['GG', 'GC', 'CG', 'CC']:
            GC_row2 += np.asarray(row)/len(nts)

    fig = plt.figure()
    plt.plot(AT_row1, 'hotpink', label='AT-rich HS')
    plt.plot(GC_row1, 'lightskyblue', label='GC-rich HS')
    plt.plot(AT_row2, 'r--', label='AT-rich Chd1')
    plt.plot(GC_row2, 'b--', label='GC-rich Chd1')
    plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
    plt.xlabel("Super Helical Location")
    plt.ylabel("Relative frequency")
    plt.legend()
    #plt.ylim([0.21, 0.285])
    #plt.savefig('ATGCperiod_' + "slide" + '.png')
    plt.show()
    plt.close()            

    kmers = LinModel.all_path(5, 'ATCG')

    X, Y = [], []
    C = []
    picks = []
    for kmer in kmers:
        freq1 = m1.freq['Kmer0'][kmer]
        freq2 = m2.freq['Kmer0'][kmer]
        GC = GC_content(kmer)
        polyA = poly_score (kmer, nts='AT', pos=False)
        X.append(freq1)
        Y.append(freq2)
        #C.append(GC)
        C.append(polyA)
        if kmer == 'AAAAA':
            picks.append((kmer, freq1, freq2))
        elif kmer == 'GGGGG':
            picks.append((kmer, freq1, freq2))


    fig = plt.figure()
    plt.scatter(X, Y, c=C, s=2, alpha=0.5, cmap='jet')
    ax = plt.gca()
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    plt.plot(lims, lims, 'k--', alpha=0.5)
    for kmer, x, y in picks:
        plt.plot(x, y, 'kx')
        plt.text(x, y, kmer, horizontalalignment='center', verticalalignment='center')
    plt.xlabel("Before")
    plt.ylabel("After")
    #plt.xlim([0.85, 1.15])
    #plt.ylim([0.85, 1.15])
    plt.title("5-mers frequency")
    plt.colorbar()
    plt.show()
    plt.close()
    

    freq_kmer1 = sorted([ [freq, kmer] for kmer, freq in m1.freq['Kmer0'].items()])
    freq_kmer2 = sorted([ [freq, kmer] for kmer, freq in m2.freq['Kmer0'].items()])
    freq_kmer3 = sorted([ [m2.freq['Kmer0'][kmer] - m1.freq['Kmer0'][kmer], kmer] for kmer in m1.freq['Kmer0']])
    #freq_kmer3 = sorted([ [100*(m2.freq['Kmer0'][kmer] - m1.freq['Kmer0'][kmer])/m1.freq['Kmer0'][kmer], kmer] for kmer in m1.freq['Kmer0']])

    X1, X2, X3 = [], [], []
    Y1, Y2, Y3 = [], [], []
    Z1, Z2, Z3 = [], [], []
    for freq, kmer in freq_kmer1:
        X1.append(kmer)
        Y1.append(freq)
        Z1.append(GC_content(kmer))

    for freq, kmer in freq_kmer2:
        X2.append(kmer)
        Y2.append(freq)
        Z2.append(GC_content(kmer))

    picks = []
    for k in range(len(freq_kmer3)):
        freq, kmer = freq_kmer3[k]
        if kmer == 'AAAAA':
            picks.append((kmer, k))
        elif kmer == 'GGGGG':
            picks.append((kmer, k))
        X3.append(kmer)
        Y3.append(freq)
        Z3.append(GC_content(kmer))

    fig = plt.figure(figsize=(15,5))
    plt.scatter(range(len(Y1)), Y1, c=Z1, s=1, alpha=0.5, cmap='jet')
    plt.title("Heat Shift")
    plt.ylabel("Relative frequency")
    plt.xlabel("5-mers")
    plt.colorbar()
    #plt.xticks(range(len(Y1)), X1, rotation=90)
    #plt.savefig('5mer_HS_' + "slide" + '.png', bbox_inches='tight')
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(15,5))
    plt.scatter(range(len(Y2)), Y2, c=Z2, s=1, alpha=0.5, cmap='jet')
    plt.title("Chd1 sliding")
    plt.ylabel("Relative frequency")
    plt.xlabel("5-mers")
    plt.colorbar()
    #plt.savefig('5mer_Chd1_' + "slide" + '.png', bbox_inches='tight')
    #plt.xticks(range(len(Y1)), X1, rotation=90)
    #plt.show()
    plt.close()

    fig = plt.figure(figsize=(15,5))
    plt.scatter(range(len(Y3)), Y3, c=Z3, s=1, alpha=0.5, cmap='jet')
    xtick_pos, xtick_label = [], []
    for kmer, pos in picks:
        plt.axvline(x=pos, linestyle = '--', linewidth=1, alpha=0.5)
        xtick_pos.append(pos)
        xtick_label.append(kmer)
    plt.xticks(xtick_pos, xtick_label)
    plt.title("(Chd1-HS)/HS (%)")
    plt.ylabel("Fold change (%)")
    plt.xlabel("5-mers")
    plt.colorbar()
    #plt.savefig('5mer_change_' + "slide" + '.png', bbox_inches='tight')
    #plt.xticks(range(len(Y1)), X1, rotation=90)
    #plt.show()
    plt.close()

fig = plt.figure()
for i in range(gnum):
    m1 = group_model1[i]
    AT_sig, GC_sig = get_ATGC_periodic(m1)
    plt.plot(AT_sig, 'hotpink', label='AT-rich HS')
    plt.plot(GC_sig, 'lightskyblue', label='GC-rich HS')
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
plt.ylabel("Relative frequency")
plt.legend()
plt.ylim([0.21, 0.285])
#plt.savefig('ATGCperiod_' + "slide" + '.png')
plt.show()
plt.close()            





    





"""                
        
# Diffusion map
neighbor_params = {'n_jobs': -1, 'algorithm': 'ball_tree'}

mydmap = dm.DiffusionMap.from_sklearn(n_evecs=2, k=10, epsilon=1, alpha=1.0, neighbor_params=neighbor_params)
dmap = mydmap.fit_transform(X)

fig = plt.figure()
for i in range(len(X)):
    plt.plot(dmap[i][0], dmap[i][1], '.', color=color_list[num_cID[i]])
plt.title("Diffusion map")
plt.show()
plt.close()

num_cID = []
for i in range(len(dmap)):
    x = dmap[i][0]
    if x < -0.52:
        num_cID.append(0)
    elif x > 0.7:
        num_cID.append(1)
    else:
        num_cID.append(2)

cID_keys = {}
for i in range(len(num_cID)):
    cID = num_cID[i]
    if cID not in cID_keys:
        cID_keys[cID] = []
    cID_keys[cID].append(num_key[i])

fig = plt.figure()
for i in range(len(X)):
    plt.plot(dmap[i][0], dmap[i][1], '.', color=color_list[num_cID[i]])
plt.title("Diffusion map")
plt.show()
plt.close()

fig = plt.figure()
for i in range(len(cID_keys)):
    key_list = cID_keys[i]
    for key in key_list:
        plt.subplot(len(cID_keys),1,i+1)
        plt.plot(key_dyadmap7[key], alpha=0.2)
        loc, mtype, nts = key.split('-')
        st = int(loc)
        ed = st+len(nts)
        plt.axvspan(st, ed-1, alpha=0.1, color='red')
plt.show()
plt.close()

    
def plot_dendrogram(model, **kwargs):

    # Children of hierarchical clustering
    children = model.children_

    # Distances between each pair of children
    # Since we don't have this information, we can use a uniform one for plotting
    distance = np.arange(children.shape[0])

    # The number of observations contained in each cluster level
    no_of_observations = np.arange(2, children.shape[0]+2)

    # Create linkage matrix and then plot the dendrogram
    linkage_matrix = np.column_stack([children, distance, no_of_observations]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

plot_dendrogram(model, labels=model.labels_)
#plot_dendrogram(model)
plt.show()
"""
