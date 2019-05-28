import load
import numpy as np
from SliderClass import Slider
import readLoopseq
import pickle
import LinModel
import matplotlib.pyplot as plt
import readLoopseq
import analysis

def GC_content (seq):
    count = 0
    for nt in seq:
        if nt in ['G', 'C']:
            count += 1
    return float(count)/len(seq)

def get_shape_list(fname, bound=2):
    id_MGW, id_HelT, id_ProT, id_Roll = load.read_DNAshape(fname)
    MGW_list, HelT_list, ProT_list, Roll_list = [], [], [], []
    for i in range(1, len(id_MGW)+1):
        MGW = id_MGW[str(i)]
        HelT = id_HelT[str(i)]
        ProT = id_ProT[str(i)]
        Roll = id_Roll[str(i)]
        MGW_list.append(MGW[bound:len(MGW)-bound])
        HelT_list.append(HelT[bound:len(HelT)-bound])
        ProT_list.append(ProT[bound:len(ProT)-bound])
        Roll_list.append(Roll[bound:len(Roll)-bound])
    return MGW_list, HelT_list, ProT_list, Roll_list

def get_SeqScoreCountShape_Slider(key_slider, fname, bound=2):
    id_MGW, id_HelT, id_ProT, id_Roll = load.read_DNAshape(fname)
    MGW_list, HelT_list, ProT_list, Roll_list = [], [], [], []
    seq_list, score_list, count_list = [], [], []
    for key in key_slider:
        slider = key_slider[key]
        seq_profile = slider.seq
        score_profile = slider.KDE()
        count_profile = score_profile
        MGW_profile = id_MGW[key]
        HelT_profile = id_HelT[key]
        ProT_profile = id_ProT[key]
        Roll_profile = id_Roll[key]
        for i in range(147/2+bound, 225-147/2-bound):
            seq = seq_profile[i-147/2:i+147/2+1]
            score = score_profile[i]
            count = count_profile[i]
            MGW = MGW_profile[i-147/2:i+147/2+1]
            HelT = HelT_profile[i-147/2:i+147/2]
            ProT = ProT_profile[i-147/2:i+147/2+1]
            Roll = Roll_profile[i-147/2:i+147/2]
            seq_list.append(seq)
            score_list.append(score)
            count_list.append(count)
            MGW_list.append(MGW)
            HelT_list.append(HelT)
            ProT_list.append(ProT)
            Roll_list.append(Roll)
    return seq_list, score_list, count_list, [MGW_list, HelT_list, ProT_list, Roll_list]
        
def plot_shape_spectrum(group_freq, bound=0, note=""):
    names = ["MGW", "HelT", "ProT", "Roll"]
    for name in names:
        fig = plt.figure()
        for i in range(len(group_freq)):
            freq = group_freq[i]
            plt.plot(freq[name][bound:len(freq[name])-bound], label="group " + str(i+1), alpha=0.8)
        plt.title(name)
        #plt.ylabel("Frequency")
        plt.xlabel("Position (bp)")
        plt.legend()
        plt.savefig('shape_spectrum_' + name + '_' + note + '.png')
        #plt.show()
        plt.close()


key_slider1 = pickle.load(open("slider1.p", "rb"))
key_slider2 = pickle.load(open("slider2.p", "rb"))
slide_seq_list1, slide_score_list1, slide_count_list1, slide_shape_list1 = get_SeqScoreCountShape_Slider(key_slider1, "Plusone_shape/phpIpl6Te", bound=2)
slide_seq_list2, slide_score_list2, slide_count_list2, slide_shape_list2 = get_SeqScoreCountShape_Slider(key_slider2, "Plusone_shape/phpIpl6Te", bound=2)

slide_seq_m1 = LinModel.SeqLinearModel(slide_seq_list1, slide_score_list1, slide_count_list1)
slide_seq_m2 = LinModel.SeqLinearModel(slide_seq_list2, slide_score_list2, slide_count_list2)

slide_seq_m1.report([1], [2,1], 1, 1, None)
slide_seq_m2.report([1], [2,1], 1, 1, None)

slide_freq1 = slide_seq_m1.freq
slide_freq2 = slide_seq_m2.freq

freq_Kmer2 = []
for Kmer, freq in slide_freq2['Kmer0'].items():
    freq_Kmer2.append([freq, Kmer])

freq_Kmer2 = sorted(freq_Kmer2)
nts = []
X = []
Y1, Y2 = [], []
Y3 = []
for freq, Kmer in freq_Kmer2:
    nts.append(Kmer)
    X.append(Kmer)
    Y1.append(slide_freq1['Kmer0'][Kmer])
    Y2.append(freq)
    Y3.append(100*float(freq - slide_freq1['Kmer0'][Kmer])/slide_freq1['Kmer0'][Kmer])
    
#fig = plt.figure()
#plt.plot(range(len(X)), Y1, 'x-', label='Heat Shift')
#plt.plot(range(len(X)), Y2, 'x-', label='Chd1 Sliding')
#plt.xticks(range(len(X)), X)
#plt.legend()
#plt.ylabel("Relative frequency")
#plt.title("Slide-seq")
#plt.savefig('dinfreq_' + "slide" + '.png')
##plt.show()
#plt.close()

fig = plt.figure()
plt.plot(range(len(X)), Y3, 'rx-', label='(Chd1-HS)/HS')
plt.axhline(y=0, color='k', linestyle='--')
plt.xticks(range(len(X)), X)
plt.legend()
plt.ylabel("Fold change (%)")
plt.title("(Chd1-HS)/HS (%)")
plt.savefig('dinfreqchange_' + "slide" + '.png')
#plt.show()
plt.close()

freq_MM1_1 = slide_freq1['MM1']
freq_MM1_2 = slide_freq2['MM1']

#for freq_MM1 in [freq_MM1_1, freq_MM1_2]:
#    fig = plt.figure()
#    for nt in nts:
#        row = []
#        for i in range(len(freq_MM1)):
#            row.append(freq_MM1[i][nt])
#        plt.plot(range(len(row)), row, label=nt)
#    plt.legend()
#    #plt.show()
#    plt.close()

nts = nts[::-1]
img = []
for nt in nts:
    row = [ 100*float(freq_MM1_2[i][nt]-freq_MM1_1[i][nt])/freq_MM1_1[i][nt] for i in range(length)]
    img.append(row)
    
fig = plt.figure()
plt.imshow(img, aspect='auto', vmin=-10, vmax=10, cmap='bwr')
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.yticks(range(len(nts)), nts)
plt.title("(Chd1-HS)/HS (%)")
plt.xlabel("Super Helical Location")
plt.colorbar()
plt.savefig('MM1change_' + "slide" + '.png')
plt.close()

length = len(freq_MM1_1)
AT_row1, GC_row1 = np.zeros(length), np.zeros(length)
AT_row2, GC_row2 = np.zeros(length), np.zeros(length)
for nt in nts:
    row = [ freq_MM1_1[i][nt] for i in range(length)]
    if nt in ['AA', 'AT', 'TA', 'TT']:
        AT_row1 += np.asarray(row)
    if nt in ['GG', 'GC', 'CG', 'CC']:
        GC_row1 += np.asarray(row)
for nt in nts:
    row = [ freq_MM1_2[i][nt] for i in range(length)]
    if nt in ['AA', 'AT', 'TA', 'TT']:
        AT_row2 += np.asarray(row)
    if nt in ['GG', 'GC', 'CG', 'CC']:
        GC_row2 += np.asarray(row)

fig = plt.figure()
plt.plot(AT_row1, 'hotpink', label='AT-rich HS')
plt.plot(GC_row1, 'lightskyblue', label='GC-rich HS')
plt.plot(AT_row2, 'r--', label='AT-rich Chd1')
plt.plot(GC_row2, 'b--', label='GC-rich Chd1')
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
plt.ylabel("Relative frequency")
plt.legend()
plt.savefig('ATGCperiod_' + "slide" + '.png')
#plt.show()
plt.close()

AT_row3 = 100*(AT_row2 - AT_row1) / AT_row1
GC_row3 = 100*(GC_row2 - GC_row1) / GC_row1

fig = plt.figure()
plt.plot(AT_row3, 'C1', label='AT-rich Chd1-HS')
plt.plot(GC_row3, 'C0', label='GC-rich Chd1-HS')
plt.xticks([147/2 + 10*i for i in range(-7, 8)], [str(10*i) for i in range(-7,8)])
plt.xlabel("Super Helical Location")
plt.ylabel("Fold change (%)")
plt.title("(Chd1-HS)/HS (%)")
plt.legend()
plt.savefig('ATGCperiodchange_' + "slide" + '.png')
#plt.show()
plt.close()

slide_seq_m1.report(None, [5,1], None, None, None)
slide_seq_m2.report(None, [5,1], None, None, None)

freq_kmer1 = sorted([ [freq, kmer] for kmer, freq in slide_seq_m1.freq['Kmer0'].items()])
freq_kmer2 = sorted([ [freq, kmer] for kmer, freq in slide_seq_m2.freq['Kmer0'].items()])
freq_kmer3 = sorted([ [100*(slide_seq_m2.freq['Kmer0'][kmer] - slide_seq_m1.freq['Kmer0'][kmer])/slide_seq_m1.freq['Kmer0'][kmer], kmer] for kmer in slide_seq_m1.freq['Kmer0']])

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

for freq, kmer in freq_kmer3:
    X3.append(kmer)
    Y3.append(freq)
    Z3.append(GC_content(kmer))

fig = plt.figure(figsize=(15,5))
plt.scatter(range(len(Y1)), Y1, c=Z1, s=2, cmap='cool')
plt.title("Heat Shift")
plt.ylabel("Relative frequency")
plt.xlabel("5-mers")
plt.colorbar()
#plt.xticks(range(len(Y1)), X1, rotation=90)
plt.savefig('5mer_HS_' + "slide" + '.png', bbox_inches='tight')
plt.close
#plt.show()

fig = plt.figure(figsize=(15,5))
plt.scatter(range(len(Y2)), Y2, c=Z2, s=2, cmap='cool')
plt.title("Chd1 sliding")
plt.ylabel("Relative frequency")
plt.xlabel("5-mers")
plt.savefig('5mer_Chd1_' + "slide" + '.png', bbox_inches='tight')
plt.close()
#plt.xticks(range(len(Y1)), X1, rotation=90)
#plt.show()

fig = plt.figure(figsize=(15,5))
plt.scatter(range(len(Y3)), Y3, c=Z3, s=2, cmap='cool')
plt.title("(Chd1-HS)/HS (%)")
plt.ylabel("Fold change (%)")
plt.xlabel("5-mers")
plt.savefig('5mer_change_' + "slide" + '.png', bbox_inches='tight')
plt.close()
#plt.xticks(range(len(Y1)), X1, rotation=90)
#plt.show()


    







        
"""
    
key_flex = readLoopseq.get_flex("plusonelib.ref", "/home/spark159/../../media/spark159/sw/all_slide_seq_data/F0_yeast.loopseq")
key_dAlen = readLoopseq.get_polyA("plusonelib.ref")
key_dGC = readLoopseq.get_GC("plusonelib.ref")
key_dCpG = readLoopseq.get_Din("plusonelib.ref", 'CG')
key_dTpA = readLoopseq.get_Din("plusonelib.ref", 'TA')

key_slider1 = pickle.load(open("slider1.p", "rb"))
key_slider2 = pickle.load(open("slider2.p", "rb"))

keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))

X1, X2, X3 = [], [], []
Y1, Y2, Y3 = [], [], []

titles = ['GC', 'Alen', 'flex', 'CpG', 'TpA']
labels = ["GC content:Right-left", "poly-A len:Right-left", "Right-Left flexibility", "CpG content:Right-left", "TpA content:Right-left"]
key_props = [key_dGC, key_dAlen, key_flex, key_dCpG, key_dTpA]

for i in range(len(titles)):
    label = labels[i]
    title = titles[i]
    print title
    key_prop = key_props[i]
    for key in keys:
        #dyadmap1 = key_slider1[key].dyadmap
        #dyadmap2 = key_slider2[key].dyadmap
        mean_pos1 = key_slider1[key].median_pos()
        mean_pos2 = key_slider2[key].median_pos()
        #left, right = key_flex[key]['left'], key_flex[key]['right']
        Y1.append(mean_pos1)
        Y2.append(mean_pos2)
        Y3.append(mean_pos2 - mean_pos1)
        #X3.append(right-left)
        #X3.append(key_dGC[key])
        #X3.append(key_dAlen[key])
        if i == 2:
            #X3.append(np.exp(key_prop[key]['right']) - np.exp(key_prop[key]['left']))
            X3.append(key_prop[key]['right'] - key_prop[key]['left'])
        else:
            X3.append(key_prop[key])

    fig = plt.figure()
    #plt.plot(X3, Y3, '.')
    plt.plot(X3, Y1, 'b.', alpha=0.5, label='Before')
    print "before: " + str(analysis.get_corr(X3,Y1))
    #plt.plot(X2, Y3, '.', label='After')
    plt.legend(loc='best', numpoints=1)
    plt.xlabel(label)
    #plt.xlabel("GC content:Right-left")
    #plt.xlabel("Right-Left flexibility")
    plt.ylabel("Mean Position")
    #plt.ylim([-15,15])
    plt.savefig("before_" + title + ".png")
    #plt.show()
    plt.close()

    fig = plt.figure()
    #plt.plot(X3, Y3, '.')
    #plt.plot(X1, Y3, '.', label='Before')
    plt.plot(X3, Y2, 'r.', alpha=0.5, label='After')
    print "after: " + str(analysis.get_corr(X3,Y2))
    plt.legend(loc='best', numpoints=1)
    #plt.xlabel("GC content:Right-left")
    plt.xlabel(label)
    #plt.xlabel("Right-Left flexibility")
    plt.ylabel("Mean Position")
    #plt.ylim([-15,15])
    plt.savefig("after_" + title + ".png")
    #plt.show()
    plt.close()

    fig = plt.figure()
    plt.plot(X3, Y3, 'g.', alpha=0.5, label='After-Before')
    print "after-before: " + str(analysis.get_corr(X3,Y3))
    #plt.plot(X1, Y3, '.', label='Before')
    #plt.plot(X2, Y3, '.', label='After')
    plt.legend(loc='best', numpoints=1)
    plt.xlabel(label)
    #plt.xlabel("GC content:Right-left")
    #plt.xlabel("Right-Left flexibility")
    plt.ylabel("Mean Position")
    #plt.ylim([-15,15])
    plt.savefig("after_before_" + title + ".png")
    #plt.show()
    plt.close()
    print


        
loop_seq_list1, loop_score_list1 = readLoopseq.read_loopseq("/home/spark159/../../media/spark159/sw/all_slide_seq_data/F0_random.loopseq", bound=2)
loop_count_list1 = [2**value for value in loop_score_list1]
loop_shape_list1 = get_shape_list("loopseq_random_shape/phptD5tBg")


loop_seq_list2, loop_score_list2 = readLoopseq.read_loopseq("/home/spark159/../../media/spark159/sw/all_slide_seq_data/F0_yeast.loopseq", bound=2)
loop_count_list2 = [2**value for value in loop_score_list2]
loop_shape_list2 = get_shape_list("loopseq_yeast_shape/phpTeWft0")


loop_seq_m1 = LinModel.SeqLinearModel(loop_seq_list1, loop_score_list1, loop_count_list1)
loop_seq_m2 = LinModel.SeqLinearModel(loop_seq_list2, loop_score_list2, loop_count_list2)


loop_seq_m1.report([1], [2,1], 1, 1, None)
loop_seq_m2.report([1], [2,1], 1, 1, None)
loop_freq1 = loop_seq_m1.freq
loop_freq2 = loop_seq_m2.freq

freq_Kmer1 = []
for Kmer, freq in loop_freq1['Kmer0'].items():
    freq_Kmer1.append([freq, Kmer])

freq_Kmer1 = sorted(freq_Kmer1)
X = []
Y1, Y2 = [], []
for freq, Kmer in freq_Kmer1:
    X.append(Kmer)
    Y1.append(freq)
    Y2.append(loop_freq2['Kmer0'][Kmer])
    
fig = plt.figure()
plt.plot(range(len(X)), Y1, 'x-', label='Random')
plt.plot(range(len(X)), Y2, 'x-', label='Yeast')
plt.xticks(range(len(X)), X)
plt.legend()
plt.ylabel("Relative frequency")
plt.title("Loop-seq")
plt.savefig('dinfreq_' + "loop" + '.png')
#plt.show()
plt.close()


group_freq1 = loop_seq_m1.spectrum([1], None, None, None, None)
group_freq2 = loop_seq_m2.spectrum([1], None, None, None, None)
#plot_shape_spectrum(group_freq1, note='loop_random')
#plot_shape_spectrum(group_freq2, note='loop_yeast')


loop_shape_m1 = LinModel.ShapeLinearModel(loop_shape_list1, loop_score_list1, loop_count_list1)
loop_shape_m2 = LinModel.ShapeLinearModel(loop_shape_list2, loop_score_list2, loop_count_list2)

group_freq1 = loop_shape_m1.spectrum(gnum=5)
group_freq2 = loop_shape_m2.spectrum(gnum=5)
plot_shape_spectrum(group_freq1, note='loop_random')
plot_shape_spectrum(group_freq2, note='loop_yeast')


key_slider1 = pickle.load(open("slider1.p", "rb"))
slide_seq_list1, slide_score_list1, slide_count_list1, slide_shape_list1 = get_SeqScoreCountShape_Slider(key_slider1, "Plusone_shape/phpIpl6Te", bound=2)
key_slider2 = pickle.load(open("slider2.p", "rb"))
slide_seq_list2, slide_score_list2, slide_count_list2, slide_shape_list2 = get_SeqScoreCountShape_Slider(key_slider2, "Plusone_shape/phpIpl6Te", bound=2)

slide_seq_m1 = LinModel.SeqLinearModel(slide_seq_list1, slide_score_list1, slide_count_list1)
slide_seq_m2 = LinModel.SeqLinearModel(slide_seq_list2, slide_score_list2, slide_count_list2)
slide_seq_m1.report([1], [2,1], 1, 1, None)
slide_seq_m2.report([1], [2,1], 1, 1, None)
slide_freq1 = slide_seq_m1.freq
slide_freq2 = slide_seq_m2.freq

freq_Kmer2 = []
for Kmer, freq in slide_freq2['Kmer0'].items():
    freq_Kmer2.append([freq, Kmer])

freq_Kmer2 = sorted(freq_Kmer2)
X = []
Y1, Y2 = [], []
for freq, Kmer in freq_Kmer2:
    X.append(Kmer)
    Y1.append(slide_freq1['Kmer0'][Kmer])
    Y2.append(freq)
    
fig = plt.figure()
plt.plot(range(len(X)), Y1, 'x-', label='Heat shift')
plt.plot(range(len(X)), Y2, 'x-', label='Chd1 sliding')
plt.xticks(range(len(X)), X)
plt.legend()
plt.ylabel("Relative frequency")
plt.title("Slide-seq")
plt.savefig('dinfreq_' + "slide" + '.png')
#plt.show()
plt.close()

slide_shape_m1 = LinModel.ShapeLinearModel(slide_shape_list1, slide_score_list1, slide_count_list1)
slide_shape_m2 = LinModel.ShapeLinearModel(slide_shape_list2, slide_score_list2, slide_count_list2)

group_freq1 = slide_shape_m1.spectrum(gnum=5)
group_freq2 = slide_shape_m2.spectrum(gnum=5)
plot_shape_spectrum(group_freq1, bound=25, note='slide_HS')
plot_shape_spectrum(group_freq2, bound=25, note='slide_Chd1')


loop_seq_list1, loop_score_list1 = readLoopseq.read_loopseq("/home/spark159/../../media/spark159/sw/all_slide_seq_data/F0_random.loopseq", bound=2)
loop_count_list2 = list(2**(np.asarray(values)))

#id_MGW, id_HelT, id_ProT, id_Roll = load.read_DNAshape("loopseq_yeast_shape/phpTeWft0")
#id_MGW, id_HelT, id_ProT, id_Roll = load.read_DNAshape("loopseq_random_shape/phptD5tBg")
#id_MGW, id_HelT, id_ProT, id_Roll = load.read_DNAshape("Plusone_shape/phpIpl6Te")

#shape_m = LinModel.ShapeLinearModel([MGW_list, HelT_list, ProT_list, Roll_list], values, counts)
#shape_m.report(scale=10)
#shape_m.train(graph=True)

m = LinModel.SeqLinearModel(seqs, values, counts)

m.train(MM_orders=[0, 1, 2], Kmer_k_b=None, PolyA_b=None, GC_b=None, Harmonic=None, graph=False, sym=True)
#m.report(MM_orders=[0,1,2], Kmer_k_b=[2,1], PolyA_b=None, GC_b=None, Harmonic=None, graph=False, sym=True, scale=10)

key_slider1 = pickle.load(open("slider1.p", "rb"))
key_slider2 = pickle.load(open("slider2.p", "rb"))

keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))

X = []
Y1, Y2 = [], []
for key in keys:
    slider1 = key_slider1[key]
    slider2 = key_slider2[key]
    Y1.append(slider1.median_pos())
    Y2.append(slider2.median_pos())
    seq = slider1.seq
    pred_flex = m.predict([seq])[0]
    pred_flex = 2**np.asarray(pred_flex)
    pred_count = pred_flex - min(pred_flex)
    length = len(pred_flex)
    #X.append(np.mean(pred_flex[:length/2]) - np.mean(pred_flex[length/2:]))

    average = 0.0
    num = 0
    for i in range(len(pred_count)):
        pos = i + 50/2
        count = pred_count[i]
        average += pos*count
        num += count
    average = float(average)/num
    X.append(average)

fig = plt.figure()
plt.plot(X, Y1, '.', alpha=0.5)
#plt.plot(X, Y2, 'x')
plt.show()
plt.close()

fig = plt.figure()
#plt.plot(X, Y1, '.')
plt.plot(X, Y2, 'x', alpha=0.5)
plt.show()
plt.close()

Y3 = np.asarray(Y2) - np.asarray(Y1)
fig = plt.figure()
#plt.plot(X, Y1, '.')
plt.plot(X, Y3, 'x', alpha=0.5)
plt.show()
plt.close()
    
freq_Kmer = []
for Kmer, freq in m.freq['Kmer0'].items():
    freq_Kmer.append([freq, Kmer])

freq_Kmer = sorted(freq_Kmer)
X, Y = [], []
for freq, Kmer in freq_Kmer:
    X.append(freq)
    Y.append(Kmer)

fig = plt.figure()
plt.plot(range(len(X)), X, 'x-')
plt.xticks(range(len(X)), Y)
plt.show()
"""
