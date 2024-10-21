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
fnames1 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_before_.combined.sort"]
fnames2 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_after_.combined.sort"]
Control1 = load.load_files(fnames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)
Control2 = load.load_files(fnames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None) 
name_key_slider['Control1'] = Control1
name_key_slider['Control2'] = Control2

# PolyA library
for condition in []:
    if condition == 'old':
        path = "/home/spark159/../../media/spark159/sw/AscanlibFinal/"
    elif condition == 'new':
        #path = "/home/spark159/../../media/spark159/sw/all_slide_seq_data/"
        path = "/home/spark159/../../media/spark159/sw/AscanlibFinal/"
    for time in [0]:
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
    name_keys[name] = sorted(keys, cmp=key_cmp)

# noise subtraction
#name_strand_threshold = {"Control2":{"top":0.05, "bott":0.02}, "mmlib_control_5_1rep":{"top":0.02, "bott":0.06}, "mmlib_control_5_2rep":{"top":0.04, "bott":0.085}, "polyAlib_old_5_1rep":{"top":0.02, "bott":0.06}}
name_strand_threshold = {}
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

# multiplexing check
name_key_slider = {}
path = "/home/spark159/../../media/spark159/sw/polyAplate/"
# PolyA plate
for condition in ['individual', 'mixed', 'blended']:
    if condition == 'individual':
        fill = None
    else:
        fill = 'linear'
    for time in [0, 5]:
        name = "%s_%s_%s" % ('polyAplate', condition, time)
        if time == 0:
            sort_fname = "%s_%s_%s.sort" % ('polyAplate', condition, 'before')
        elif time == 5:
            sort_fname = "%s_%s_%s.sort" % ('polyAplate', condition, 'after')
        try:
            with open(name + ".pickle", "rb") as f:
                key_slider = pickle.load(f)
        except:
            key_slider = load.load_files([path +  sort_fname], ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=fill)
            with open(name + ".pickle", "wb") as f:
                pickle.dump(key_slider, f)
        assert name not in name_key_slider
        name_key_slider[name] = key_slider

names = sorted(name_key_slider.keys())
        
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
        for k in range(len(common_keys)):
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
                library_type, condition, time = name1.split('_')
                s = condition.upper() + "\n" + time + "min "
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
plt.show()
plt.close()






sys.exit(1)

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
