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
from sklearn.decomposition import PCA

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
mtype_choice = ["M"]
size_st, size_ed = 1, 5

# load 601 control data
fnames1 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_before_.combined.sort"]
fnames2 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_after_.combined.sort"]
Control1 = load.load_files(fnames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)['601']
Control2 = load.load_files(fnames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)['601'] 


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


# data clustering
key_dyadmap3 = {}
key_dyadmap7 = {}
key_change = {}
for key in keys:
    key_dyadmap3[key] = analysis.norm(key_slider3[key].dyadmap)
    key_dyadmap7[key] = analysis.norm(key_slider7[key].dyadmap)
    key_change[key] = [key_dyadmap7[key][k] - key_dyadmap3[key][k] for k in range(ref_length)]

key_cluster, cluster_keys, model= analysis.Kmeans(key_dyadmap7, cluster_num=4, type_targets=['slicing', [(NCP_len/2, ref_length-NCP_len/2)]])

#sample_list = []
#for i in range(len(cluster_keys)):
#    sample_list.append(sorted(cluster_keys[i], cmp=key_cmp))
#graph.plot_map(key_slider3, sample_list, norm_choice=True, draw_key=True, note='_clustering')

cluster_num = len(cluster_keys)
fig = plt.figure()
for i in range(cluster_num):
    key_list = cluster_keys[i]
    plt.subplot(cluster_num, 1, i+1)
    for key in key_list:
        #plt.plot(key_change[key], alpha=0.2)
        plt.plot(key_dyadmap7[key], alpha=0.2)
        loc, mtype, nts = key.split('-')
        st = int(loc)
        ed = st+len(nts)
        plt.axvspan(st, ed-1, alpha=0.1, color='red')
    plt.title("Cluster #%d" % (i+1))
    plt.xlim([0,225])

#plt.show()
plt.close()

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

# hierarchial clustering
X = []
key_KL = {}
key_num, num_key = {}, []
for i in range(len(keys)):
    key = keys[i]
    key_num[key] = i
    num_key.append(key)
    #X.append([key_dyadmap7[key][k] - key_dyadmap3[key][k] for k in range(ref_length)])
    #key_KL[key] = analysis.KL_div(key_dyadmap7[key], key_dyadmap3[key])
    X.append(key_dyadmap7[key][NCP_len/2:ref_length-NCP_len/2])
    #X.append([key_dyadmap7[key][k] - analysis.norm(Control2.dyadmap)[k] for k in range(ref_length)])
    key_KL[key] = analysis.KL_div(key_dyadmap7[key], analysis.norm(Control2.dyadmap))


# PCA
pca = PCA(n_components=5).fit(X)
Xr = pca.transform(X)

fig = plt.figure()
plt.plot(pca.explained_variance_)
plt.show()
plt.close()

fig = plt.figure()
for i in range(len(pca.components_)):
    plt.subplot(pca.n_components_, 1, i+1) 
    plt.plot(pca.components_[i])
plt.show()
plt.close()


# hierarichal clustering
color_list = np.linspace(0.01, 1, num=len(keys))
cmap = cm.get_cmap("Reds")
key_color = {}
key_list = dict_sort(key_KL)
for i in range(len(key_list)):
    key = key_list[i]
    color = cmap(color_list[i])
    key_color[key] = color

Z = linkage(Xr, 'ward')
node_children = {i:{} for i in range(len(num_key))}
node_keys = {i:{num_key[i]} for i in range(len(num_key))}
for i in range(len(Z)):
    node1, node2 = int(Z[i][0]), int(Z[i][1])
    new_node = max(node_keys.keys()) + 1
    node_children[new_node] = set([node1, node2])
    node_keys[new_node] = node_keys[node1] | node_keys[node2]

num_cID = [cID-1 for cID in fcluster(Z, 5, 'maxclust')]
cID_keys = {}
for i in range(len(num_cID)):
    cID = num_cID[i]
    key = num_key[i]
    if cID not in cID_keys:
        cID_keys[cID] = set([])
    cID_keys[cID].add(key)

def get_descendants (node_children, root):
    output = set([])
    for node in node_children[root]:
        output.add(node)
        output |= get_descendants(node_children, node)
    return output

cID_nodes = {}
for cID, keys1 in cID_keys.items():
    for node, keys2 in node_keys.items():
        if keys1 == keys2:
            if cID not in cID_nodes:
                cID_nodes[cID] = set([])
            cID_nodes[cID].add(node)
            cID_nodes[cID] |= get_descendants(node_children, node)
            break

#color_list = np.linspace(0.01, 1, num=len(cID_nodes))
#cmap = cm.get_cmap("jet")

color_list = ['r', 'g', 'b', 'y', 'm', 'c']
node_color = ['black'] * (2 * len(Xr) - 1)
for i in range(len(cID_nodes)):
    cID = cID_nodes.keys()[i]
    #color = cmap(color_list[i])
    color = color_list[i]
    for node in cID_nodes[cID]:
        node_color[node] = color

fig = plt.figure()
dn = dendrogram(Z, link_color_func=lambda k: node_color[k])

ax = plt.gca()
xlbls = ax.get_xmajorticklabels()
for lbl in xlbls:
        lbl.set_color(key_color[num_key[int(lbl.get_text())]])

#plt.savefig("dendrogram.png")
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

fig = plt.figure()
for i in range(len(Xr)):
    plt.plot(Xr[i][0], Xr[i][1], '.', color=color_list[num_cID[i]])
plt.show()
plt.close()



"""
    
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
