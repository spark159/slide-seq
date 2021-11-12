import sys
import math
import copy
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


def key_cmp (key1, key2, priority='length'):
    st1, mtype1, nts1 = key1.split('-')
    st2, mtype2, nts2 = key2.split('-')
    st1, st2 = int(st1), int(st2)
    ed1, ed2 = st1 + len(nts1), st2 + len(nts2)
    if priority == 'length':
        a = (len(nts1), st1)
        b = (len(nts2), st2)
    elif priority == 'start':
        a = (st1, len(nts1))
        b = (st2, len(nts2))
    elif priority == 'end':
        a = (ed1, len(nts1))
        b = (ed2, len(nts2))
    return tuple_cmp(a, b)

def key_cmp_st (key1, key2):
    st1, mtype1, nts1 = key1.split('-')
    st2, mtype2, nts2 = key2.split('-')
    st1, st2 = int(st1), int(st2)
    ed1, ed2 = st1 + len(nts1), st2 + len(nts2)
    a = (st1, len(nts1))
    b = (st2, len(nts2))
    return tuple_cmp(a, b)

def key_cmp_ed (key1, key2):
    st1, mtype1, nts1 = key1.split('-')
    st2, mtype2, nts2 = key2.split('-')
    st1, st2 = int(st1), int(st2)
    ed1, ed2 = st1 + len(nts1), st2 + len(nts2)
    a = (ed1, len(nts1))
    b = (ed2, len(nts2))
    return tuple_cmp(a, b)


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
for mtype in ['M']:
    if mtype == 'M':
        library_type = 'mm'
    elif mtype in ['I', 'D']:
        library_type = 'ID'
    for condition in ['bubble']:
        for time in [5]:
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
#for condition in ['new:corrected']:
for condition in []:
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

# heatmap plotting
if True:
    for name in name_key_slider:
        if name.startswith('Control'):
            continue
        print name
        print
        library_type, condition, time, rep = name.split('_')

        key_slider = name_key_slider[name]
        keys = name_keys[name]
        
        size_keys = {}
        for key in keys:
            #print key
            st, mtype, nts = key.split('-')
            size = len(nts)
            if size not in size_keys:
                size_keys[size] = []
            size_keys[size].append(key)

        sample1, sample2, sample3, sample4 = [], [], [], []

        for size in sorted(size_keys.keys()):
            sample1 += size_keys[size]
            #if size >=3 and size <=5:
            #    sample1 += size_keys[size]
            #elif size >=6 and size <=8:
            #    sample2 += size_keys[size]
            #elif size >=9 and size <=11:
            #    sample3 += size_keys[size]
            #else:
            #    assert size >=12 and size <=15
            #    sample4 += size_keys[size]

        #sample_list = [sample1, sample2, sample3, sample4]
        sample_list = [sample1]

        # order by size
        graph_edit.plot_map(key_slider, [sorted(sample, cmp=key_cmp) for sample in sample_list], norm_choice=True, obs_func = Slider.peak_signal, draw = 'key', slicing=0, note='order_by_size')
        # order by start location
        graph_edit.plot_map(key_slider, [sorted(sample, cmp=key_cmp_st) for sample in sample_list], norm_choice=True, obs_func = Slider.peak_signal, draw = 'key', slicing=0, note='order_by_start')
        # order by end location
        graph_edit.plot_map(key_slider, [sorted(sample, cmp=key_cmp_ed) for sample in sample_list], norm_choice=True, obs_func = Slider.peak_signal, draw = 'key', slicing=0, note='order_by_end')
            
        

# clustering analysis
fakekey = '0-M-A'
if False:
    name_key_cID, name_cID_keys = {}, {}
    name_outlier_cID, name_cID_outliers = {}, {}
    for name in name_key_slider:
        if name.startswith('Control'):
            continue
        print name
        print
        library_type, condition, time, rep = name.split('_')

        prediction = ('predict' in condition.split('-'))

        key_slider = name_key_slider[name]
        keys = name_keys[name]

        # temporally add fake 601 key
        key_slider[fakekey] = name_key_slider['Control2']['601']
        keys.insert(0, fakekey) 
    

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

        #print num
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

        algorithm = 'Spectral'
        #algorithm = 'Serial_spectral'
        metric = 'JS'
        cnum_min, cnum_max = 9, 9
        #out_threshold = 0 #polyAlib
        out_threshold = -np.inf #mmlib/IDlib
        cnum_ms_list = {}

        color_list = ['r', 'g', 'b', 'orange', 'm', 'c', 'darksalmon', 'navy', 'limegreen', 'skyblue', 'grey']*2

        # distance metric
        if metric == 'JS':
            # load distrance matrix
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


        elif metric == 'PCA_euc':
            key1_key2_dist = {}
            for i in range(len(keys)):
                key1 = keys[i]
                for j in range(i, len(keys)):
                    key2 = keys[j]
                    vect1, vect2 = Xr[i], Xr[j]
                    dist = np.sqrt(sum([(vect2[k]-vect1[k])**2 for k in range(len(vect1))]))
                    if key1 not in key1_key2_dist:
                        key1_key2_dist[key1] = {}
                    key1_key2_dist[key1][key2] = dist
                    if key2 not in key1_key2_dist:
                        key1_key2_dist[key2] = {}
                    key1_key2_dist[key2][key1] = dist

                    
        D = np.zeros((len(keys), len(keys))) # distance matrix
        A = np.zeros((len(keys), len(keys))) # similarity matrix
        for i in range(len(keys)-1):
            for j in range(i+1, len(keys)):
                key1, key2 = keys[i], keys[j]
                dist = key1_key2_dist[key1][key2]
                D[i][j], D[j][i] = dist, dist
                score = np.exp(-10*dist)
                A[i][j], A[j][i] = score, score


        # change the number of clusters                
        for cnum in range(cnum_min, cnum_max+1):
            #Kmeans clustering
            if algorithm == 'Kmeans':

                idx_cID = sklearn.cluster.KMeans(init='k-means++', n_init=10, n_clusters = cnum, max_iter = 10000).fit_predict(np.asarray(X))

            # Bottom-up hierarchial clustering
            if algorithm == 'Hierarchy':
                # hierarichal clustering by PCA components
                #Z = linkage(Xr, 'ward')
                #Z = linkage(Xr, 'average', optimal_ordering=True)

                # hierarchial clustering with JS distance metric of raw data
                y = squareform(D)
                #Z = linkage(y, 'average', optimal_ordering=True)
                Z = linkage(y, 'single', optimal_ordering=True)

                idx_cID = [cID-1 for cID in fcluster(Z, cnum, 'maxclust')]

            
            # Spectral clustering
            if algorithm == 'Spectral':                
                idx_cID = sklearn.cluster.spectral_clustering(affinity=A, n_clusters=cnum, random_state=0)
                
            # Top-down hierarchial spectral clustering
            if algorithm == 'H_Spectral':

                # general form
                #def split(group):
                #    assert len(group) > 1
                #    left, right = split_machine(group)
                #    return left, right

                def Spectral_bisect (key_list):
                    assert len(key_list) > 1
                    nA = np.zeros((len(key_list), len(key_list))) # similarity matrix
                    for i in range(len(key_list)-1):
                        for j in range(i+1, len(key_list)):
                            key1, key2 = key_list[i], key_list[j]
                            idx1, idx2 = key_num[key1], key_num[key2]
                            score = A[idx1][idx2]
                            nA[i][j], nA[j][i] = score, score

                    idx_cID = sklearn.cluster.spectral_clustering(affinity=nA, n_clusters=2, random_state=0)
                    output = [[], []]
                    for i in range(len(key_list)):
                        key = key_list[i]
                        cID = idx_cID[i]
                        output[cID].append(key)

                    left, right = output
                    return left, right
                
                def recursive_split (group):
                    btree = []
                    if len(group) <= 1:
                        return group
                    #left, right = split(group)
                    left, right = Spectral_bisect(group)
                    btree.append(recursive_split(left))
                    btree.append(recursive_split(right))
                    return btree

                def read_tree (tree):
                    # BFS tree search
                    visited = []
                    tnode_children = {}
                    queue = [tree]
                    while len(queue) > 0 :
                        t = queue.pop(0)
                        tnode = str(t)
                        visited.append(tnode)
                        tnode_children[tnode] = set([])
                        for subt in t:
                            if type(subt) != list:
                                continue
                            queue.append(subt)
                            tnode_children[tnode].add(str(subt))

                    # index the nodes in the order of clustering
                    tnode_node = {}
                    count = 0
                    for tnode in visited[::-1]:
                        if len(tnode_children[tnode]) < 1:
                            tnode_node[tnode] = count
                            count +=1
                    for tnode in visited[::-1]:
                        if len(tnode_children[tnode]) >= 1:
                            tnode_node[tnode] = count
                            count +=1
                    node_tnode = dict_sort(tnode_node)
                    
                    # node to children nodes
                    node_children = {}
                    for tnode, children in tnode_children.items():
                        node = tnode_node[tnode]
                        node_children[node] = set([tnode_node[child] for child in list(children)])

                    # node to all leaves(data) belong to it
                    node_leaves = {}
                    for node in sorted(node_children.keys()):
                        if node not in node_leaves:
                            node_leaves[node] = set([])
                        children = list(node_children[node])
                        if len(children) < 1:
                            node_leaves[node].add(node_tnode[node][2:-2])
                        else:
                            for child in children:
                                node_leaves[node] |= node_leaves[child]
                            
                    return node_children, node_leaves

                # make linkage matrix Z for binary tree
                def make_Z (node_children, node_leaves):
                    def average_dist (keys1, keys2, key1_key2_dist):
                        dist_list = []
                        for key1 in list(keys1):
                            for key2 in list(keys2):
                                dist = key1_key2_dist[key1][key2]
                                dist_list.append(dist)
                        return np.median(dist_list)
                    Z = []
                    for node in sorted(node_children.keys()):
                        children = sorted(list(node_children[node]))
                        if len(children) <= 1:
                            continue
                        left, right = children
                        left_leaves, right_leaves = node_leaves[left], node_leaves[right]
                        dist = average_dist(left_leaves, right_leaves, key1_key2_dist)
                        row = [left, right, dist, len(left_leaves) + len(right_leaves)]
                        Z.append(row)
                    return Z

                btree = recursive_split (keys)
                node_children, node_keys = read_tree(btree)
                Z = make_Z(node_children, node_keys)

                fig = plt.figure()
                #dn = dendrogram(Z, truncate_mode='lastp')
                dn = dendrogram(Z)
                #plt.savefig("dendrogram.png")
                #plt.show()
                plt.close()

                # clustering
                key_cID, cID_keys = {}, {}
                #leafnode_cID = [cID-1 for cID in fcluster(Z, cnum, 'maxclust')]
                #leafnode_cID = [cID-1 for cID in fcluster(Z, cnum, 'maxclust_monocrit')]
                leafnode_cID = [cID-1 for cID in fcluster(Z, t=0.65, criterion='distance')]
                for i in range(len(leafnode_cID)):
                    leafnode = i
                    cID = leafnode_cID[leafnode]
                    assert len(node_keys[leafnode]) == 1
                    key = list(node_keys[leafnode])[0]
                    key_cID[key] = cID
                    if cID not in cID_keys:
                        cID_keys[cID] = []
                    cID_keys[cID].append(key)

                #sys.exit(1)

            # user defined serial spectral clustering
            elif algorithm == 'Serial_spectral':
                def Spectral_clustering (key_list, cnum):
                    assert len(key_list) > 1
                    nA = np.zeros((len(key_list), len(key_list))) # similarity matrix
                    for i in range(len(key_list)-1):
                        for j in range(i+1, len(key_list)):
                            key1, key2 = key_list[i], key_list[j]
                            idx1, idx2 = key_num[key1], key_num[key2]
                            score = A[idx1][idx2]
                            nA[i][j], nA[j][i] = score, score

                    idx_cID = sklearn.cluster.spectral_clustering(affinity=nA,
                                                                  n_clusters=cnum,
                                                                  random_state=0)
                    key_cID, cID_keys = {}, {}
                    for i in range(len(idx_cID)):
                        key = key_list[i]
                        cID = idx_cID[i]
                        key_cID[key] = cID
                        if cID not in cID_keys:
                            cID_keys[cID] = []
                        cID_keys[cID].append(key)

                    return [cID_keys[cID] for cID in sorted(cID_keys.keys())]

                def recursive_clustering (data, div_tree):
                    data_tree = []
                    if len(div_tree) <= 0:
                        return data
                    clusters = Spectral_clustering (data, len(div_tree))
                    for cluster, sub_tree in zip(clusters, div_tree):
                        data_tree.append(recursive_clustering(cluster, sub_tree))
                    return data_tree

                def all_leaves (tree):
                    def is_leaf (data_list):
                        for data in data_list:
                            if type(data) == list:
                                return False
                        return True
                    output = []
                    if is_leaf(tree):
                        return [{'path':[], 'data':tree}]
                    for i in range(len(tree)):
                        subtree = tree[i]
                        for leaf in all_leaves(subtree):
                            leaf['path'].append(i)
                            output.append(leaf)
                    return output

                div_tree = [[[],[]], [], [[],[[],[]],[]], [], [[],[]], []] #polyAlib
                #div_tree = [[],[],[],[],[],[]]

                key_tree = recursive_clustering (keys, div_tree)
                all_leaves = all_leaves (key_tree)

                cID_keys = {}
                for i in range(len(all_leaves)):
                    leaf = all_leaves[i]
                    cID = '-'.join([str(num+1) for num in leaf['path'][::-1]])
                    #cID = i
                    cID_keys[cID] = sorted(leaf['data'])

                key_cID = {}
                for cID in cID_keys:
                    for key in cID_keys[cID]:
                        key_cID[key] = cID

                #for cID in cID_keys:
                #    key_list = sorted(cID_keys[cID], cmp=key_cmp)
                #    graph_edit.plot_map(key_slider, [key_list], norm_choice=True, obs_func = Slider.peak_signal, draw='key', slicing=0, note= '_' + cID + "_clustering")
                    
                #sys.exit(1)
                    
                        
            # read clustering output
            if algorithm in ['Spectral', 'Kmeans']:
                key_cID, cID_keys = {}, {}
                for i in range(len(idx_cID)):
                    key = keys[i]
                    cID = idx_cID[i]
                    key_cID[key] = cID
                    if cID not in cID_keys:
                        cID_keys[cID] = []
                    cID_keys[cID].append(key)
                
                
            # compute Silhouette coefficients
            key_s = analysis.Silhouette (key_cID, cID_keys, key1_key2_dist)
            #key_s = {key:1 for key in keys}
            
            # take out outliers from clusters
            cID_outliers = {}
            outlier_cID = {}
            for key, s in key_s.items():
                if s < out_threshold:
                    cID = key_cID[key]
                    
                    if cID not in cID_outliers:
                        cID_outliers[cID] = []
                    cID_outliers[cID].append(key)

                    outlier_cID[key] = cID

                    cID_keys[cID].remove(key)
                    del key_cID[key]

            assert set(keys) == set(outlier_cID.keys()) | set(key_cID.keys())


            # change cluster ID in the order of the mean value of signal
            cID_meanpos = {}
            for cID in cID_keys:
                temp = []
                for key in cID_keys[cID]:
                    meanpos = key_slider[key].median_pos()
                    temp.append(meanpos)
                cID_meanpos[cID] = np.mean(temp)

            cIDs = dict_sort(cID_meanpos)

            cID_newcID = {}
            for i in range(len(cIDs)):
                cID = cIDs[i]
                cID_newcID[cID] = i

            key_newcID = {}
            for key, cID in key_cID.items():
                key_newcID[key] = cID_newcID[cID]

            newcID_keys = {}
            for cID, key_list in cID_keys.items():
                newcID_keys[cID_newcID[cID]] = key_list

            key_cID = key_newcID
            cID_keys = newcID_keys

            del key_newcID
            del newcID_keys

            cIDs = sorted(cID_keys.keys())

            # save clustering result
            name_key_cID[name] = key_cID
            name_cID_keys[name] = cID_keys
            name_outlier_cID[name] = outlier_cID
            name_cID_outliers[name] = cID_outliers

            cID_601 = key_cID[fakekey]

            
            # Silhouette plot for all data
            fig = plt.figure()
            i = 0
            outX_list, outY_list = [], []
            for cID in sorted(cID_keys):
                color = color_list[cID]
                try:
                    total = cID_keys[cID]+cID_outliers[cID]
                except:
                    total = cID_keys[cID]
                s_key = sorted([(key_s[key], key) for key in total], cmp=tuple_cmp, reverse=True)
                X_list, Y_list = [], [] 
                for s, key in s_key:
                    if s < out_threshold:
                        outX_list.append(i)
                        outY_list.append(s)
                    else:
                        X_list.append(i)
                        Y_list.append(s)
                    i +=1
                plt.plot(X_list, Y_list, '.', color=color, label='Cluster #' + str(cID+1))
            if len(outX_list) > 0:
                plt.plot(outX_list, outY_list, 'k.', label='Outliers')
                
            plt.xlabel('seq ID')
            plt.ylabel('s-value')
            plt.title("Silhouette plot")
            plt.legend()
            #plt.show()
            plt.close()


            # compute mean Silhouette coefficient per cluster
            cID_s = {}
            for cID in cID_keys:
                for key in cID_keys[cID]:
                    s = key_s[key]
                    if cID not in cID_s:
                        cID_s[cID] = []
                    cID_s[cID].append(s)

            cID_ms = {}
            for cID in cID_s:
                ms = np.mean(cID_s[cID])
                cID_ms[cID] = ms

            cnum_ms_list[cnum] = cID_ms.values()


        # Silhouette coefficients V.S. cluster number
        fig = plt.figure()
        X_list = []
        box_data = []
        for cnum in sorted(cnum_ms_list.keys()):
            X_list.append(cnum)
            box_data.append(cnum_ms_list[cnum])
        plt.boxplot(box_data, positions=X_list)
        #plt.show()
        plt.close()

        
        # PCA plot
        fig = plt.figure()
        for i in range(len(cIDs)):
            cID = cIDs[i]
            color = color_list[cID]
            X_list, Y_list = [], []
            for key in cID_keys[cID]:
                idx = key_num[key]
                X_list.append(Xr[idx][0])
                Y_list.append(Xr[idx][1])
            plt.plot(X_list, Y_list, '.', color=color, label='Cluster #' + str(cID+1))
            #plt.plot(X_list, Y_list, '.', label='Cluster #' + str(cID+1))

        outX_list, outY_list = [], []
        for key in outlier_cID:
            idx = key_num[key]
            outX_list.append(Xr[idx][0])
            outY_list.append(Xr[idx][1])
        if len(outX_list) > 0:
            plt.plot(outX_list, outY_list, 'k.', label='Outliers')
                    
        plt.title("PCA plot")
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.legend()
        #plt.show()
        plt.close()


        # Diffusion map
        neighbor_params = {'n_jobs': -1, 'algorithm': 'ball_tree'}
        mydmap = dm.DiffusionMap.from_sklearn(n_evecs=2, k=10, epsilon=1, alpha=1.0, neighbor_params=neighbor_params)
        dmap = mydmap.fit_transform(X)

        fig = plt.figure()
        for i in range(len(cIDs)):
            cID = cIDs[i]
            color = color_list[cID]
            X_list, Y_list = [], []
            for key in cID_keys[cID]:
                idx = key_num[key]
                X_list.append(dmap[idx][0])
                Y_list.append(dmap[idx][1])
            plt.plot(X_list, Y_list, '.', color=color, label='Cluster #' + str(cID+1))
            #plt.plot(X_list, Y_list, '.', label='Cluster #' + str(cID+1))

        outX_list, outY_list = [], []
        for key in outlier_cID:
            idx = key_num[key]
            outX_list.append(dmap[idx][0])
            outY_list.append(dmap[idx][1])
        if len(outX_list) > 0:
            plt.plot(outX_list, outY_list, 'k.', label='Outliers')
                    
        plt.title("Diffusion map")
        plt.legend()
        #plt.show()
        plt.close()

        
        # tSNE plot
        perp_min, perp_max, perp_step=80, 80, 1 #polyAlib
        #perp_min, perp_max, perp_step=60, 60, 1 #mmlib

        for perp in range(perp_min, perp_max+1, perp_step):
            tsne = sklearn.manifold.TSNE(n_components=2, metric='precomputed', perplexity=perp, random_state=0)
            trans_data = tsne.fit_transform(D).T

            fig = plt.figure()
            for i in range(len(cIDs)):
                cID = cIDs[i]
                color = color_list[cID]
                X_list, Y_list = [], []
                for key in cID_keys[cID]:
                    idx = key_num[key]
                    X_list.append(trans_data[0][idx])
                    Y_list.append(trans_data[1][idx])
                plt.plot(X_list, Y_list, '.', color=color, label='Cluster #' + str(cID+1))
                #plt.plot(X_list, Y_list, '.', label='Cluster #' + str(cID+1))

            outX_list, outY_list = [], []
            for key in outlier_cID:
                idx = key_num[key]
                outX_list.append(trans_data[0][idx])
                outY_list.append(trans_data[1][idx])
            if len(outX_list) > 0:
                plt.plot(outX_list, outY_list, 'k.', label='Outliers')

            plt.title("tSNE plot (perplexity=%d)" % (perp))
            plt.xlabel("tSNE 1")
            plt.ylabel("tSNE 2")
            plt.legend()
            #plt.show()
            plt.close()

        
        # plot similarity matrix
        img = np.zeros((len(keys), len(keys)))
        img[:] = np.nan
        for i in range(len(keys)-1):
            for j in range(i+1, len(keys)):
                key1, key2 = keys[i], keys[j]
                dist = key1_key2_dist[key1][key2]
                score = np.exp(-10*dist)
                img[i][j] = score

        new_keys = []
        for cID in sorted(cID_keys.keys()):
            temp = cID_keys[cID]
            random.shuffle(temp)
            new_keys += temp

        outliers = outlier_cID.keys()
        random.shuffle(outliers)
        new_keys += outliers

        for i in range(len(new_keys)-1): 
           for j in range(i+1, len(new_keys)):
                key1, key2 = new_keys[i], new_keys[j]
                dist = key1_key2_dist[key1][key2]
                score = np.exp(-10*dist)
                img[j][i] = score

        fig = plt.figure()
        #plt.imshow(img, cmap="YlOrRd")
        plt.imshow(img, cmap="jet")
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Similarity value', rotation=-90, va="bottom")
        plt.title("Similarity Matrix")
        plt.xlabel("seq ID")
        plt.ylabel("seq ID")
        #plt.show()
        plt.close()


# plot maps
color_list = ['r', 'g', 'b', 'orange', 'm', 'c', 'darksalmon', 'navy', 'limegreen', 'skyblue', 'grey']*2
#color_list = ['r', 'g', 'b', 'orange', 'm', 'c', 'darksalmon', 'navy', 'seagreen', 'skyblue']*2
if False:
    for name, key_slider in name_key_slider.items():
        if name.startswith('Control'):
            continue
        print name
        print

        keys = name_keys[name]

        key_cID = name_key_cID[name]
        cID_keys = name_cID_keys[name]
        outlier_cID = name_outlier_cID[name]
        cID_outliers = name_cID_outliers[name]

        # remove 601 data
        try:
            keys.remove(fakekey)
        except:
            pass
        try:
            del key_cID[fakekey]
            cID_keys[cID_601].remove(fakekey)
        except:
            pass
        try:
            del outlier_cID[fakekey]
            cID_outliers[cID_601].remove(fakekey)
        except:
            pass

        cIDs = sorted(cID_keys.keys())

        # plot heatmap for each cluster
        #obs_img = graph_edit.plot_map(key_slider, [sorted(cID_keys[cID], cmp=key_cmp_st) for cID in cIDs], norm_choice=True, obs_func = Slider.get_dyadmap, draw = 'key', slicing=0, note='Clustering_raw')
        #graph_edit.plot_map(key_slider, [sorted(cID_keys[cID], cmp=key_cmp_st) for cID in cIDs], norm_choice=True, obs_func = Slider.peak_signal, draw = 'key', slicing=0, note='Clustering')
        #graph_edit.plot_map(key_slider, [sorted(outlier_cID.keys(), cmp=key_cmp)], norm_choice=True, obs_func = Slider.get_dyadmap, draw = 'key', slicing=0, note='out_Clustering_raw')
        #graph_edit.plot_map(key_slider, [sorted(outlier_cID.keys(), cmp=key_cmp)], norm_choice=True, obs_func = Slider.peak_signal, draw='key', slicing=0, note='out_Clustering')

        
        # plot clustering map
        fig = plt.figure(figsize=(18,4))
        for i in range(len(cIDs)):
            cID = cIDs[i]
            color = color_list[cID]
            X_list, Y_list = [], []
            for key in cID_keys[cID]:
                loc, mtype, nts = key.split('-')
                st = int(loc)
                X_list.append(st+1)
                Y_list.append(len(nts))
            plt.plot(X_list, Y_list, 'o', color=color, label='Cluster #' + str(cID+1))
        outX_list, outY_list = [], []
        for key in outlier_cID:
            loc, mtype, nts = key.split('-')
            st = int(loc)
            outX_list.append(st+1)
            outY_list.append(len(nts))
        if len(outX_list) > 0:
            plt.plot(outX_list, outY_list, 'ko', label='Outliers')
        plt.title("Clustering map")
        plt.xlabel("Poly-A tract start position on 601 DNA")
        plt.ylabel('Size (bp)')
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.close()

        
        # plot KL divergence map
        fig = plt.figure(figsize=(18,4))
        X_list, Y_list = [], []
        C_list = []
        for key in key_cID.keys():
            loc, mtype, nts = key.split('-')
            st = int(loc)
            KL = key_KL[key]
            X_list.append(st+1)
            Y_list.append(len(nts))
            C_list.append(KL)
        plt.scatter(X_list, Y_list, c=C_list, cmap = 'RdPu', vmin=min(key_KL.values()), vmax=max(key_KL.values()))
        cbar = plt.colorbar(pad=0.01)
        cbar.ax.set_ylabel('KL-div', rotation=-90, va="bottom", labelpad=-1)
        #cbar.ax.set_label("test", fontsize=15, labelpad=10)
        plt.title("KL divergence map")
        plt.xlabel("Poly-A tract start position on 601 DNA")
        plt.ylabel('Size (bp)')
        plt.tight_layout()
        #plt.show()
        plt.close()


        # make Widom 601 perturbation map
        size_loc_KLs = {}
        for key in key_cID.keys():
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

        total_cIDs = cIDs
        group_num = min([len(cIDs), 6])
        group_num = len(cIDs)
        #group_num = 5
        page_num = int(math.ceil(float(len(cID_keys)) / group_num))
        for u in range(page_num):
            fig = plt.figure()
            cIDs = total_cIDs[group_num*u:min(group_num*(u+1), len(total_cIDs))]
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
            #plt.show()
            plt.close()



# plot all signals in pdf file
# cluster colors
color_list = ['r', 'g', 'b', 'orange', 'm', 'c', 'darksalmon', 'navy', 'limegreen', 'skyblue', 'grey']*2
#color_list = ['r', 'g', 'b', 'orange', 'm', 'c', 'darksalmon', 'navy', 'seagreen', 'skyblue']
if False:
    # draw data on pdf

    m, n = 10, 4 # number of rows and colums
    for name, key_slider in name_key_slider.items():
        if name.startswith('Control'):
            continue
        print name
        print
        pdf = matplotlib.backends.backend_pdf.PdfPages(name + "_cluster.pdf")
        key_cID = name_key_cID[name]
        cID_keys = name_cID_keys[name]
        outlier_cID = name_outlier_cID[name]
        cID_outliers = name_cID_outliers[name]

        # plot clusters
        for cID in sorted(cID_keys.keys()):
            keys = sorted(cID_keys[cID], cmp=key_cmp)
            color = color_list[cID]
            page_nums = int(math.ceil(len(keys)/float(m*n))) # number of pages per cluster
            #page_nums = 1
            for i in range(page_nums):
                fig = plt.figure(figsize=(15,20))
                j = 0
                while j < min(m*n, len(keys)-m*n*i):
                    key = keys[m*n*i + j]
                    if key == fakekey:
                        label = '601'
                    else:
                        label = key
                    dyad_map = [value for value in key_slider[key].dyadmap]
                    plt.subplot(m, n, j+1)
                    ax = plt.gca()
                    ax.spines['bottom'].set_color(color)
                    ax.spines['top'].set_color(color)
                    ax.spines['left'].set_color(color)
                    ax.spines['right'].set_color(color)
                    plt.plot(dyad_map, 'k-', label=label)
                    loc, mtype, nts = key.split('-')
                    st = int(loc)
                    ed = st+len(nts)
                    plt.axvspan(st, ed-1, alpha=0.5, color='red')
                    #plt.title(key)
                    #plt.title('Cluster ' + str(cID+1))
                    leg = plt.legend(loc='upper right', frameon=False)
                    for item in leg.legendHandles:
                        item.set_visible(False)
                    plt.xlim([0,225])
                    j +=1
                fig.suptitle('Cluster #' + str(cID+1) + ' (' + str(i+1) + '/' + str(page_nums) + ')', fontsize=20, color=color, y=0.92)
                pdf.savefig(fig)
                plt.close()

        # plot outliers
        keys = sorted(outlier_cID.keys(), cmp=key_cmp)
        page_nums = int(math.ceil(len(keys)/float(m*n))) # number of pages per cluster
        #page_nums = 1
        for i in range(page_nums):
            fig = plt.figure(figsize=(15,20))
            j = 0
            while j < min(m*n, len(keys)-m*n*i):
                key = keys[m*n*i + j]
                if key == fakekey:
                    label = '601'
                else:
                    label = key
                dyad_map = [value for value in key_slider[key].dyadmap]
                plt.subplot(m, n, j+1)
                plt.plot(dyad_map, 'k-', label=label)
                loc, mtype, nts = key.split('-')
                st = int(loc)
                ed = st+len(nts)
                plt.axvspan(st, ed-1, alpha=0.5, color='red')
                #plt.title(key)
                #plt.title('Cluster ' + str(cID+1))
                leg = plt.legend(loc='upper right', frameon=False)
                for item in leg.legendHandles:
                    item.set_visible(False)
                plt.xlim([0,225])
                j +=1
            fig.suptitle('Outliers ' +  '(' + str(i+1) + '/' + str(page_nums) + ')', fontsize=20, y=0.92)
            pdf.savefig(fig)
            plt.close()

        pdf.close()
   

