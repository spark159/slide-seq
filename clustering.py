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
for condition in ['old']:
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
    name_keys[name] = sorted(keys, cmp=key_cmp)

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
if True:
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
        pca = PCA(n_components=10).fit(X)
        variance_ratio_list = 100* pca.explained_variance_ratio_

        num = 1
        while sum(variance_ratio_list[:num]) <= 90: # can explain 90% of variance
            num +=1

        #num=10

        fig = plt.figure()
        plt.plot(range(1, len(variance_ratio_list)+1), variance_ratio_list, '.-')
        plt.axvline(x=num, color = 'red', linestyle='--', alpha=0.5)
        plt.xlabel("PCA component")
        plt.ylabel("Variance (%)")
        plt.xticks(range(1, len(variance_ratio_list)+1))
        #plt.show()
        plt.close()

        #sys.exit(1)

        pca = PCA(n_components=num).fit(X)
        Xr = pca.transform(X)

        fig = plt.figure()
        for i in range(len(pca.components_)):
            plt.subplot(pca.n_components_, 1, i+1) 
            plt.plot(pca.components_[i])
        #plt.show()
        plt.close()


        # hierarichal clustering by PCA components
        Z = linkage(Xr, 'ward')

        # hierarchial clustering with JS distance metric of raw data
        #Z = linkage(X, 'weighted', metric=analysis.JS_dist)

        node_children = {i:{} for i in range(len(num_key))}
        node_keys = {i:{num_key[i]} for i in range(len(num_key))}
        for i in range(len(Z)):
            node1, node2 = int(Z[i][0]), int(Z[i][1])
            new_node = max(node_keys.keys()) + 1
            node_children[new_node] = set([node1, node2])
            node_keys[new_node] = node_keys[node1] | node_keys[node2]

        num_cID = [cID-1 for cID in fcluster(Z, 8, 'maxclust')]
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
        color_list = np.linspace(0, 1, num=len(keys))
        cmap = cm.get_cmap('RdPu')
        key_color = {}
        key_list = dict_sort(key_KL)
        for i in range(len(key_list)):
            key = key_list[i]
            color = cmap(color_list[i])
            key_color[key] = color

        # cluster color map
        #color_list = np.linspace(0.01, 1, num=len(cID_nodes))
        #cmap = cm.get_cmap("jet")
        color_list = ['r', 'g', 'b', 'y', 'm', 'c', 'darksalmon', 'seagreen', 'skyblue', 'navy']
        node_color = ['black'] * (2 * len(Xr) - 1)
        for i in range(len(cID_nodes)):
            cID = cID_nodes.keys()[i]
            #color = cmap(color_list[i])
            color = color_list[i]
            for node in cID_nodes[cID]:
                node_color[node] = color

        # plot dendogram
        fig = plt.figure()
        dn = dendrogram(Z, link_color_func=lambda k: node_color[k])

        ax = plt.gca()
        xlbls = ax.get_xmajorticklabels()
        for lbl in xlbls:
                lbl.set_color(key_color[num_key[int(lbl.get_text())]])

        #plt.savefig("dendrogram.png")
        #plt.show()
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

        ## Other clustering approach
        ##Kmeans clustering
        #key_xr = {}
        #for i in range(len(Xr)):
        #    key = num_key[i]
        #    key_xr[key] = Xr[i]
        #key_cID, cID_keys = analysis.Kmeans(key_xr, cluster_num=5, sample_list=None, type_targets=[None, []])
        #color_list = ['r', 'g', 'b', 'y', 'm', 'c', 'darksalmon', 'seagreen', 'skyblue', 'navy']
    
        #Spectral clustering
        A = np.zeros((len(keys), len(keys))) # similarity matrix
        for i in range(len(keys)-1):
            for j in range(i+1, len(keys)):
                key1, key2 = keys[i], keys[j]
                dyadmap1 = analysis.norm(key_slider[key1].dyadmap)
                dyadmap2 = analysis.norm(key_slider[key2].dyadmap)
                JS = analysis.JS_dist(dyadmap1, dyadmap2)
                score = np.exp(-10*JS)
                A[i][j] = score
                A[j][i] = score

        idx_cID = sklearn.cluster.spectral_clustering(affinity=A, n_clusters=10)
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
        color_list = ['r', 'g', 'b', 'orange', 'm', 'c', 'darksalmon', 'navy', 'seagreen', 'skyblue']
        eigen_key = sorted([(v[i,1], keys[i]) for i in range(len(keys))], cmp=tuple_cmp)
        fig = plt.figure()
        for i in range(len(eigen_key)):
            eigen, key = eigen_key[i]
            cID = key_cID[key]
            plt.plot(i, eigen, '.', color=color_list[cID])
        plt.ylabel("Eigenvector")
        #plt.show()
        plt.close()

        # plot of similarity matrix
        img = np.zeros((len(keys), len(keys)))
        img[:] = np.nan
        for i in range(len(keys)-1):
            for j in range(i, len(keys)):
                img[i][j] = A[i][j]

        new_keys = []
        for cID in sorted(cID_keys.keys()):
            new_keys += cID_keys[cID]

        for i in range(len(new_keys)-1):
            for j in range(i+1, len(new_keys)):
                key1, key2 = new_keys[i], new_keys[j]
                idx1, idx2 = key_num[key1], key_num[key2]
                img[j][i] = A[idx1][idx2]

        fig = plt.figure()
        plt.imshow(img, cmap="YlOrRd")
        plt.colorbar()
        plt.show()
        plt.close()
        
        # mean position for each cluster
        cID_meanpos = {}
        for cID in cID_keys:
            temp = []
            for key in cID_keys[cID]:
                meanpos = key_slider[key].median_pos()
                temp.append(meanpos)
            cID_meanpos[cID] = np.mean(temp)

        # make perturbation map
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

        # clustering map
        fig = plt.figure()
        for key in keys:
            loc, mtype, nts = key.split('-')
            st = int(loc)
            ed = st+len(nts)
            length = len(nts)
            pos = (st + ed)/2
            cID = key_cID[key]
            #cID = num_cID[key_num[key]]
            alpha = 0.01 + (key_KL[key]-min(key_KL.values()))*float(1.0 - 0.01)/(max(key_KL.values())-min(key_KL.values()))
            plt.plot(range(st, ed), [len(nts)]*len(nts), '.', color=color_list[cID], alpha=alpha)
        #plt.show()
        plt.close()

        group_num = 5
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

        sys.exit(1)

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
                #plt.show()
                plt.close()

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
    min_energy, max_energy = None, None
    size_shl_energy = {}
    for size in size_dyad_shl_energy:
        if size not in size_shl_energy:
            size_shl_energy[size] = {}
        for dyad in size_dyad_shl_energy[size]:
            for shl in size_dyad_shl_energy[size][dyad]:
                if shl not in size_shl_energy[size]:
                    size_shl_energy[size][shl] = []
                size_shl_energy[size][shl] += size_dyad_shl_energy[size][dyad][shl]
                energy = np.mean(size_dyad_shl_energy[size][dyad][shl])
                if min_energy == None or energy < min_energy:
                    min_energy = energy
                if max_energy == None or energy > max_energy:
                    max_energy = energy

    print min_energy, max_energy
    extreme = min(abs(min_energy), max_energy)

    #fig = plt.figure()
    for size in sorted(size_dyad_shl_energy.keys()):
        fig = plt.figure()
        ax1 = plt.gca()
        #ax1.plot(analysis.norm(control.dyadmap)[(ref_length - NCP_len)/2:ref_length-(ref_length - NCP_len)/2], 'k')
        ax2 = ax1.twinx()
        color_list = np.linspace(0, 1, num=len(size_dyad_shl_energy[size]))
        cmap = cm.get_cmap("brg")
        mean_signal = [np.nan]*NCP_len
        for shl in sorted(size_shl_energy[size].keys()):
            mean_signal[shl] = np.mean(size_shl_energy[size][shl])
        for i in range(len(size_dyad_shl_energy[size])):
            dyad = sorted(size_dyad_shl_energy[size].keys())[i]
            signal = [np.nan]*NCP_len
            for shl in sorted(size_dyad_shl_energy[size][dyad].keys()):
                energy = np.mean(size_dyad_shl_energy[size][dyad][shl])
                signal[shl] = energy
                weight = float(control.dyadmap[dyad]-min(control.dyadmap))/(max(control.dyadmap)-min(control.dyadmap))
            ax1.plot(range(-NCP_len/2, NCP_len/2), signal, '.', color='k', alpha=0.05)
        ax2.plot(range(-NCP_len/2, NCP_len/2), mean_signal, 'r.', alpha=1)
            #plt.scatter(range(-NCP_len/2, NCP_len/2), [dyad]*NCP_len, c=signal, cmap='RdBu_r', vmin=-extreme, vmax=extreme)
            #plt.scatter([dyad]*NCP_len, range(-NCP_len/2, NCP_len/2), c=signal, cmap='RdBu_r', vmin=-extreme, vmax=extreme)
            #ax = ax2.scatter([dyad - (ref_length - NCP_len)/2]*NCP_len, range(-NCP_len/2, NCP_len/2), c=signal, cmap='RdBu_r', vmin=-extreme, vmax=extreme, alpha=0.5)
            #for tick in [-60, -40, -20, -10, 0, 10, 20, 40, 60]:
            #    ax2.axhline(y=tick, linestyle = '--', linewidth=1, color='k', alpha=0.25)
            #plt.axhline(y=dyad, linestyle = '--', color='r')
        #plt.colorbar()
        #plt.axis("equal")
        ax1.set_ylabel("$\Delta$ E")
        ax2.set_ylabel("Average $\Delta$ E" )
        plt.title("Size " + str(size))
        plt.savefig('EFmean' + str(size) + '.png', bbox_inches='tight')
        #plt.colorbar(ax)
        #plt.ylabel("SHL coordinate")
        #plt.xlabel("Dyad location")
        #plt.show()
        #plt.close()
    #plt.show()
    #plt.close()
    """

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
