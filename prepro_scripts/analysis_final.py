# -*- coding: utf-8 -*-
import os, sys, subprocess, re
import copy
import math
import random
import numpy as np
import scipy
import sklearn
import sklearn.cluster
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.decomposition import NMF

# get histogram
def get_hist (data, binnum=1000, prob=False):
    hist={};
    if prob:
        deno=float(len(data))
    else:
        deno=1.0
    binwidth=float(max(data)-min(data))/binnum
    for value in data:
        bin=int((value-min(data))/binwidth)
        bincenter=min(data)+(bin+0.5)*binwidth
        if bincenter not in hist:
            hist[bincenter]=0
        hist[bincenter]+=1/deno
    return hist


# compute GC content of sequence
def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100


# find a reverse complementary of sequence
def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq:
        nt = nt.upper()
        if nt == 'A':
            new_seq += 'T'
        elif nt == 'T':
            new_seq += 'A'
        elif nt == 'C':
            new_seq += 'G'
        elif nt == 'G':
            new_seq += 'C'
        else:
            new_seq += nt
    return new_seq


# check the sequence is palindromic
def is_pal (seq):
    if len(seq) % 2 != 0:
        return False
    for i in range(len(seq)/2):
        if nt[i] != rev_comp(nt[len(seq)-1-i]):
            return False
    return True


# find the longest length of poly-nt (pos=False)
# find all poly-nt locations and counts (pos=True)
def polynt_count (seq, nts, pos=False):
    len_pos = {}
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
            if count not in len_pos:
                len_pos[count] = []
            len_pos[count].append(i)
            i = j
        else:
            i +=1
    if pos:
        return len_pos    
    if len(len_pos) == 0:
        return 0
    return max(len_pos.keys())


# find the count of the given dinucleotide
# find the count of all existing dinucleotide (din=None)
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


# all possible k-mers with given nts
def all_kmer(k, nts='ATCG'):
    if N==1:
        return list(nts)
    output=[]
    for kmer in all_kmer(k-1):
        for nt in nts:
            output.append(kmer+nt)
    return output


# compare function for sorting tuple
# first component is higher priority
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


# compare function for sorting tuple
# second component is higher priority
def tuple_cmp_rev (a, b, priority=[1,0]):
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


# sort keys by the values of the dictionary
def dict_sort (dict, reverse=False):
    value_key = [(value, key) for key, value in dict.items()]
    value_key = sorted(value_key, cmp=tuple_cmp, reverse=reverse)
    return [value_key[i][1] for i in range(len(value_key))]


# compare function for sorting window-ids
# window length is the highest priority
def wid_cmp_len (wid1, wid2):
    loc1, mtype1, nts1 = wid1.split('-')
    loc2, mtype2, nts2 = wid2.split('-')
    len1, len2 = len(nts1), len(nts2)
    wst1, wst2 = int(loc1), int(loc2)
    a = (len1, wst1)
    b = (len2, wst2)
    return tuple_cmp(a, b)


# compare function for sorting window-ids
# window start location is the highest priority
def wid_cmp_st (wid1, wid2):
    loc1, mtype1, nts1 = wid1.split('-')
    loc2, mtype2, nts2 = wid2.split('-')
    wst1, wst2 = int(loc1), int(loc2)
    wed1, wed2 = wst1 + len(nts1), wst2 + len(nts2)
    a = (wst1, len(nts1))
    b = (wst2, len(nts2))
    return tuple_cmp(a, b)


# compare function for sorting window-ids
# window end location is the highest priority
def wid_cmp_ed (wid1, wid2):
    loc1, mtype1, nts1 = wid1.split('-')
    loc2, mtype2, nts2 = wid2.split('-')
    wst1, wst2 = int(loc1), int(loc2)
    wed1, wed2 = wst1 + len(nts1), wst2 + len(nts2)
    a = (wed1, len(nts1))
    b = (wed2, len(nts2))
    return tuple_cmp(a, b)


# return the normalized list
def normalize_list(L):
    total = float(sum(L))
    if total <=0:
        return L
    return [value/total for value in L]


# return the normalized matrix (a list of list of same length)
def normalize_matrix (M):
    total=0.0
    for i in range(len(M)):
        for j in range(len(M[i])):
            total += M[i][j]
    if total <= 0:
        return new
    new = [[0]*len(M[0]) for i in range(len(M))]
    for i in range(len(M)):
        for j in range(len(M[i])):
            new[i][j] = (M[i][j]/float(total))
    return new


# background substraction (To do)
def sub_background (map, frac=0.1):
    thres= min(map) + frac*(max(map)-min(map))
    #new = [0 for i in range(len(map))]
    #for i in range(len(map)):
    #    if map[i] > thres:
    #        new[i] = map[i]
    #return new
    return [max(value-thres, 0.0) for value in map]


# find peaks of 1-d signal
def find_peaks(sig, num=None):
    pos_value={}
    for i in range(1, len(sig)-1):
        if sig[i] > sig[i-1] and sig[i] > sig[i+1]:
            pos_value[i] = sig[i]
    if num == None:
        num = len(pos_value)
    else:
        num = min(num, len(pos_value))
    return dict_sort(pos_value, reverse=True)[:num]


# find/broden peaks and igonore others in 1-d siganl 
def highlight_peaks (sig, num=15, broden=1):
    new_sig = [0.0]*len(sig)
    for pos in find_peaks(sig, num=num):
        value = sig[pos]
        new_sig[pos] += value
        for offset in range(1, broden+1):
            new_sig[pos-offset] +=value
            new_sig[pos+offset] +=value
    return new_sig


# Kernel density estimation of signal
# return natural log of KDE
def logKDE (sig, scale=1.0, band_width=1.0):
    X, X_plot = [], []
    for k in range(len(sig)):
        for num in range(int(scale*sig[k])):
            X.append([k])
        X_plot.append([k])
    kde = sklearn.neighbors.KernelDensity(kernel="gaussian", bandwidth=band_width).fit(X)
    log_density = kde.score_samples(X_plot)
    return log_density


# get pearson correlation
def get_corr(X, Y):
    assert len(X) == len(Y)
    X_mean = np.mean(X)
    Y_mean = np.mean(Y)
    XYterm = 0.0
    Xterm, Yterm = 0.0, 0.0
    for i in range(len(X)):
        dX = X[i] - X_mean
        dY = Y[i] - Y_mean
        XYterm += dX*dY
        Xterm += dX*dX
        Yterm += dY*dY
    return XYterm / np.sqrt(Xterm * Yterm)


# check the vector is a probabilty mass function
def is_pmf (vec):
    for value in vec:
        if value < 0:
            return False
    if int(round(sum(vec))) != 1:
        return False
    return True

    
# get Kullback–Leibler divergence
# pmf1 divergence from pmf2
def KL_div (pmf1, pmf2):
    assert len(pmf1) == len(pmf2)
    assert is_pmf(pmf1)
    assert is_pmf(pmf2)
    div = 0.0
    for i in range(len(pmf1)):
        if  pmf1[i] == 0 or pmf2[i] == 0:
            continue
        div += pmf1[i] * np.log2(float(pmf1[i])/pmf2[i])
    return div


# get Jensen–Shannon divergence
def JS_div (pmf1, pmf2):
    m = [(pmf1[i] + pmf2[i])*0.5 for i in range(len(pmf1))]
    return 0.5*KL_div(pmf1, m) + 0.5*KL_div(pmf2, m)


# get Jensen–Shannon distance
def JS_dist (pmf1, pmf2):
    return np.sqrt(JS_div(pmf1, pmf2))


# PCA analysis
# input: id to feature, standardization option, component number
# output: id to PCA vector
def PCA_analysis (id_feature, stand_mean=True, stand_std=True, n_components=False):
    ids = id_feature.keys()
    X = [id_feature[id] for id in ids]

    # check the dimensionality of data
    dim = len(X[0])
    for i in range(len(X)):
        assert dim == len(X[i]) 

    # standardization of data
    scaler = preprocessing.StandardScaler(with_mean=stand_mean, with_std=stand_std).fit(X)
    X_scaled = scaler.transform(X)
    pca = PCA(n_components=10).fit(X_scaled)

    # Dimensionality reduction
    if not n_components:
        variance_ratio_list = 100* pca.explained_variance_ratio_
        num = 1 # component number
        while num < dim and sum(variance_ratio_list[:num]) <= 90: # can explain 90% of variance
            num +=1
    else:
        num = n_components
    assert num <= dim
    
    if num >= dim:
        print >> sys.stderr, "Warning:dimensionality reduction is not possible."

    # save PCA components
    Xr = [vect[:num] for vect in pca.transform(X)]
    id_pca = {ids[i]:Xr[i] for i in range(len(ids))}
    return id_pca, pca.components_[:num], variance_ratio_list[:num]

# None-negative Matrix Factorization
# input: data features, cluster number
# output: cluster-ID to basis, id to weight of each basis, approximation error
def NMF_analysis (id_feature, cnum):
    ids = id_feature.keys()
    X = [id_feature[id] for id in ids]

    model = NMF(n_components=cnum, init='random', random_state=0, verbose=False)
    W = model.fit_transform(X)
    H = model.components_

    cID_basis = []
    for i in range(cnum):
        cID_basis.append(H[i])

    id_weight = {}
    for i in range(len(ids)):
        id = ids[i]
        id_weight[id] = W[i]

    error = model.reconstruction_err_

    return cID_basis, id_weight, error

# Kmeans clustering
# input: data features, cluster number
# output: id to cluster-ID, cluster-ID to id
def Kmeans_clustering (id_feature, cnum):
    ids = id_feature.keys()
    X = [id_feature[id] for id in ids]

    # Kmeans clustering
    idx_cID = sklearn.cluster.KMeans(init='k-means++', n_init=10, n_clusters = cnum, max_iter = 10000, random_state=0).fit_predict(X)

    id_cID, cID_ids = {}, {}
    for i in range(len(idx_cID)):
        id = ids[i]
        cID = idx_cID[i]
        id_cID[id] = cID
        if cID not in cID_ids:
            cID_ids[cID] = []
        cID_ids[cID].append(id)

    return id_cID, cID_ids


# compute Euclidean distance between two vectors
def Euc_dist (vec1, vec2):
    assert len(vec1) == len(vec2)
    total = 0.0
    for i in range(len(vec1)):
        total += (vec1[i] - vec2[i])**2
    return np.sqrt(total)


# compute pair-wise distance
def get_pair_dist (id_vec, metric):
    ids = id_vec.keys()
    id1_id2_dist = {}
    for i in range(len(ids)-1):
        for j in range(i+1, len(ids)):
            id1, id2 = ids[i], ids[j]
            vec1, vec2 = id_vec[id1], id_vec[id2]
            if metric == 'Euc':
                dist = Euc_dist(vec1, vec2)
            elif metric == 'JS':
                dist = JS_dist(vec1, vec2)
            if id1 not in id1_id2_dist:
                id1_id2_dist[id1] = {}
            id1_id2_dist[id1][id2] = dist
            if id2 not in id1_id2_dist:
                id1_id2_dist[id2] = {}
            id1_id2_dist[id2][id1] = dist
    return id1_id2_dist


# compute pair-wise similarity score
def get_pair_score (id1_id2_dist, scale):
    ids = id1_id2_dist.keys()
    id1_id2_score = {}
    for i in range(len(ids)-1):
        for j in range(i+1, len(ids)):
            id1, id2 = ids[i], ids[j]
            dist = id1_id2_dist[id1][id2]
            score = np.exp(-scale*dist)
            if id1 not in id1_id2_score:
                id1_id2_score[id1] = {}
            id1_id2_score[id1][id2] = score
            if id2 not in id1_id2_score:
                id1_id2_score[id2] = {}
            id1_id2_score[id2][id1] = score
    return id1_id2_score


# Spectral clustering
# input: pair-wise similarity score, cluster number
# output: id to cluster-ID, cluster-ID to id
def Spectral_clustering (id1_id2_score, cnum):
    ids = id1_id2_score.keys()

    # get similarity matrix
    A = np.zeros((len(ids), len(ids)))
    for i in range(len(ids)-1):
        for j in range(i+1, len(ids)):
            id1, id2 = ids[i], ids[j]
            score = id1_id2_score[id1][id2]
            A[i][j], A[j][i] = score, score

    # spectral clustering
    idx_cID = sklearn.cluster.spectral_clustering(affinity=A, n_clusters=cnum, random_state=0)

    id_cID, cID_ids = {}, {}
    for i in range(len(idx_cID)):
        id = ids[i]
        cID = idx_cID[i]
        id_cID[id] = cID
        if cID not in cID_ids:
            cID_ids[cID] = []
        cID_ids[cID].append(id)

    return id_cID, cID_ids


# compute Silhouette score for each clustered data
# input: id to cluster-ID, cluster-ID to id, pair-wise distance
# output: id to Silhouette score
def Silhouette (id_cID, cID_ids, id1_id2_dist):
    id_cID_dists = {}
    for id in id_cID:
        for cID in cID_ids:
            for other_id in cID_ids[cID]:
                if id == other_id:
                    continue
                dist = id1_id2_dist[id][other_id]
                if id not in id_cID_dists:
                    id_cID_dists[id] = {}
                if cID not in id_cID_dists[id]:
                    id_cID_dists[id][cID] = []
                id_cID_dists[id][cID].append(dist)

    id_cID_mdist = {}
    for id in id_cID_dists:
        for cID in id_cID_dists[id]:
            mdist = np.mean(id_cID_dists[id][cID])
            if id not in id_cID_mdist:
                id_cID_mdist[id] = {}
            id_cID_mdist[id][cID] = mdist

    id_s = {}
    for id in id_cID_mdist:
        cID = id_cID[id]
        cID_mdist = id_cID_mdist[id]
        a = cID_mdist[cID]
        b = np.min([cID_mdist[i] for i in list(set(cID_mdist.keys())-set([cID]))])
        if a < b:
            s = 1.0 - float(a)/b
        elif a > b:
            s = float(b)/a - 1.0
        else:
            assert a == b
            s = 0
        id_s[id] = s
        
    return id_s


# tSNE dimesionality reduction
# input: pair-wise distance, component number, perplexity
# output: id to tsne vector
def tSNE_analysis (id1_id2_dist, n_components, perp):
    ids = id1_id2_dist.keys()

    # get distance matrix
    D = np.zeros((len(ids), len(ids)))
    for i in range(len(ids)-1):
        for j in range(i+1, len(ids)):
            id1, id2 = ids[i], ids[j]
            dist = id1_id2_dist[id1][id2]
            D[i][j], D[j][i] = dist, dist

    # tSNE transformation
    tsne = sklearn.manifold.TSNE(n_components=n_components, metric='precomputed', perplexity=perp, random_state=0)
    trans_data = tsne.fit_transform(D).T

    # save tSNE components
    id_tsne = {ids[i]:[trans_data[0][i], trans_data[1][i]] for i in range(len(ids))}
    return id_tsne


# Diffusion map analysis
# input: pair-wise distance, time per step, step number, component number
# output: id to diffusion vector
def Diffusion_map (id1_id2_dist, epsilon, step, n_components):
    ids = id1_id2_dist.keys()

    # get distance matrix
    D = np.zeros((len(ids), len(ids)))
    for i in range(len(ids)-1):
        for j in range(i+1, len(ids)):
            id1, id2 = ids[i], ids[j]
            dist = id1_id2_dist[id1][id2]
            D[i][j], D[j][i] = dist, dist
            
    dist_matrix = np.asarray(D)

    # make diffusion matrix
    L = np.exp(-dist_matrix**2 / (2*epsilon))
    norm = np.sum(L, axis=0)
    Di = np.diag(norm)
    inv_Di = np.diag(1.0/norm)

    M = np.matmul(inv_Di, L)
    sym_M = np.matmul(Di**0.5, np.matmul(M, inv_Di**0.5))

    # solve eigenvalue problem
    eigenValues, eigenVectors = np.linalg.eigh(sym_M)
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]

    time = epsilon*step
    eigenScaler = np.diag(eigenValues**time)

    # compute diffusion coordinates
    diff_coords = np.matmul(np.matmul(inv_Di**0.5, eigenVectors), eigenScaler)
    id_diff = {ids[i]:diff_coords[i][:n_components] for i in range(len(ids))}
    
    return id_diff


# Fast Fourier transformation
# input: 1-d siganl
# output: periods, amplts, phases(unit of pi) of basis function
def FFT (sig):
    N = len(sig)
    sig_ft = scipy.fft(sig)[1:N/2]
    periods = [float(N)/k for k in range(1, N/2)]
    amplts = np.abs(sig_ft)/float(N)
    phases = np.arctan2(sig_ft.imag, sig_ft.real) / np.pi # unit of pi
    shifts = np.asarray([(phases[k-1] * N) / (2*np.pi*k) for k in range(1, N/2)]) /np.pi # unit of pi
    return periods, amplts, phases


# rescale the data in old range (old_st, old_ed) into new range (new_st, new_ed)
def rescale (value_list, old_st, old_ed, new_st, new_ed):
    output = []
    for value in value_list:
        assert value >= old_st and value <= old_ed
        new_value = new_st + (new_ed - new_st)*float(value-old_st)/(old_ed-old_st)
        output.append(new_value)
    return output

