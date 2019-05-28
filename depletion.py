import load
import numpy as np
import EnModel
import matplotlib.pyplot as plt

# read pluse-one library
ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
#filenames1 = ["../../Illumina/plusoneHS/data/Plslib-HS_S1_L001_R.sort"]
filenames1 = ["../../Illumina/SW_Ascan_new/data/Ascan0_S1_L001_R.sort"]
#filenames2 = ["../../Illumina/plusoneHS/data/Plslib-HS-30min_S2_L001_R.sort"]
filenames2 = ["../../Illumina/SW_Ascan_new/data/Ascan-5min_S1_L001_R.sort"]
#ref_fname = "../../Illumina/plusoneHS/plusonelib.ref"
ref_fname = "../../Illumina/SW_Ascan_new/polyAscanlib.ref"
key_slider1 = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', load_ref=ref_fname)
key_slider2 = load.load_files(filenames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill='linear', load_ref=ref_fname)
#id_seq = read_library("../../Illumina/plusoneHS/plusonelib.ref")

keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
size_fold1 = {}
size_fold2 = {}
ref_key = 'AAA-46'

def norm (L):
    total = 0.0
    for value in L:
        total += value
    return [value/total for value in L]

for key in keys:
    win,loc = key.split('-')
    size, loc = len(win), int(loc)
    st, ed = loc, loc+size
    KDE1 = norm(key_slider1[key].dyadmap)
    KDE2 = norm(key_slider2[key].dyadmap)
    refKDE1 = norm(key_slider1[ref_key].dyadmap)
    refKDE2 = norm(key_slider2[ref_key].dyadmap)
    if sum(refKDE1[st:ed]) <= 0:
        fold1 = np.NaN
    else:
        fold1 = 1 - sum(KDE1[st:ed])/sum(refKDE1[st:ed])
    if sum(refKDE2[st:ed]) <= 0:
        fold2 = np.NaN
    else:
        fold2 = 1 - sum(KDE2[st:ed])/sum(refKDE2[st:ed])
    if size not in size_fold1:
        size_fold1[size] = []
    if size not in size_fold2:
        size_fold2[size] = []
    size_fold1[size].append(fold1)
    size_fold2[size].append(fold2)

X1, Y1, Z1 = [], [], []
X2, Y2, Z2 = [], [], []
for size in sorted(size_fold1.keys()):
    X1.append(size)
    Y1.append(np.mean(size_fold1[size]))
    Y2.append(np.mean(size_fold2[size]))    
    Z1.append(np.std(size_fold1[size]))
    Z2.append(np.std(size_fold2[size]))

fig = plt.figure()
plt.plot(X1,Y1, '.')
plt.errorbar(X1,Y1,yerr=Z1,fmt='o')
plt.show()
plt.close()


fig = plt.figure()
plt.plot(X1,Y2, '.')
plt.errorbar(X1,Y2,yerr=Z2,fmt='o')
plt.show()
plt.close()

