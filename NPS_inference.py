import math
import copy
import random
import matplotlib.pyplot as plt

def normalize(L):
    total = 0.0
    for i in range(len(L)):
        total += L[i]
    return [value/total for value in L]

def EM_algorithm (cutmap, noisemap, offset_cleav, LR_ratio, max_iternum, udate='poisson', covg=0.001):
    def distance (kmap1, kmap2):
        dyadmap1 = [ (kmap1['R'][i] + kmap1['L'][i])*0.5 for i in range(len(kmap1['R']))]
        dyadmap2 = [ (kmap2['R'][i] + kmap2['L'][i])*0.5 for i in range(len(kmap2['R']))]
        dyadmap1, dyadmap2 = normalize(dyadmap1), normalize(dyadmap2)
        dist = 0.0
        for i in range(len(dyadmap1)):
            dist += (dyadmap1[i] - dyadmap2[i])**2
        return math.sqrt(dist)
    
    def poisson_update(oldmap, cutmap, noisemap, offset_cleav, LR_ratio):
        N = len(oldmap['R'])
        t = {'R':[0.0]*N, 'L':[0.0]*N}
        for side in ['R', 'L']:
            for i in range(N):
                for offset in offset_cleav[side].keys():
                    try:
                        t[side][i] += oldmap[side][i-offset]*offset_cleav[side][offset]
                    except:
                        continue
                t[side][i] += noisemap[side][i]
        newmap = {'R':[0.0]*N, 'L':[0.0]*N}
        alpha = {'R':1.0, 'L':LR_ratio}
        for i in range(N):
            factor = 0.0
            for side in ['R', 'L']:
                for offset in offset_cleav[side].keys():
                    try:
                        factor += alpha[side]*(cutmap[side][i+offset]*offset_cleav[side][offset])/t[side][i+offset]
                    except:
                        continue
            factor /= (1 + LR_ratio)
            newmap['R'][i] = oldmap['R'][i] * factor
            newmap['L'][i] = newmap['R'][i] * LR_ratio
        return newmap
    
    N = len(cutmap['R'])
    oldmap = {'R':[1.0]*N, 'L':[1.0]*N}
    for i in range(max_iternum):
        if udate == 'poisson':
            newmap = poisson_update(oldmap, cutmap, noisemap, offset_cleav, LR_ratio)
        elif udate == 'NB':
            newmap = NB_update(oldmap, cutmap, noisemap, offset_cleav, LR_ratio)
        if distance(oldmap, newmap) < covg:
            break
        oldmap = copy.deepcopy(newmap)
    return [ (newmap['R'][i] + newmap['L'][i])*0.5 for i in range(len(newmap['R']))]


"""
noisemap = {'R':[0.0]*225, 'L':[0.0]*225}
offset_cleav = {'R':{-32:1}, 'L':{+32:1.0}}
LR_ratio = 1.0
cutmap = {'R':[ random.randint(0, 10000) for i in range(225)], 'L':[ random.randint(0, 10000) for i in range(225)]}

dyadmap = [0.0]*225
for side in ['R', 'L']:
    for i in range(len(cutmap[side])):
        if side == 'R':
            offset = 32
        else:
            offset = -32
        try:
            dyadmap[i + offset] += cutmap[side][i]
        except:
            continue

EMmap = EM_algorithm(cutmap, noisemap, offset_cleav, LR_ratio, 1)
plt.plot(normalize(EMmap), label='EM')
plt.plot(normalize(dyadmap), label='Sum')
print sum(normalize(EMmap))
print sum(normalize(dyadmap))
plt.legend()
plt.show()
"""    
    
    
    
