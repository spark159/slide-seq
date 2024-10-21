import math
import random
import pickle
import EnModel
import SliderClass
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf


# physical parameters
# periodicity of DNA in nucleosome
L = 10
#L = 11.5

# JJ's stiffness coefficient matrix
K1 = {'TA':12, 'TG':6, 'TT':68, 'TC':60, 'CA':6, 'CG':12, 'CT':35, 'CC':34, 'AA':68, 'AG':35, 'AT':48, 'AC':44, 'GA':60, 'GG':34, 'GT':44, 'GC':58}
# Aakash's Loop-seq parameters
K2 = {'TA':-0.0505, 'TG':0.0169, 'TT':0.0028   , 'TC':0.0276, 'CA':0.0144, 'CG':0.0432, 'CT':-0.0242, 'CC':-0.0188, 'AA':-0.0007, 'AG':-0.0197, 'AT':0.0043, 'AC':0.0112, 'GA':0.0347, 'GG':-0.0172, 'GT':0.0099, 'GC':-0.0061}
# Sangwoo's Slide-seq parameters (Heat shift)
K3 = {'AA': -0.012876517419095658, 'AC': -0.0020720202526200043, 'GT': -0.0020720202526200043, 'AG': -0.0035479073236012355, 'CC': 0.01874051980894128, 'CA': -0.004600379011252342, 'CG': 0.018381325258072546, 'TT': -0.012876517419095663, 'GG': 0.01874051980894127, 'GC': 0.014126910581500979, 'AT': -0.010490625823870931, 'GA': -0.0018208848564576926, 'TG': -0.004600379011252344, 'CT': -0.0035479073236012364, 'TC': -0.0018208848564576937, 'TA': -0.009663231906070343}
# Sangwoo's Slide-seq parameters (Chd1 sliding) 
K4 = {'AA': -0.006563952473581738, 'AC': -0.008489309364065396, 'GT': -0.00848930936406539, 'AG': -0.0017946835269604506, 'CC': 0.01737341520832092, 'CA': -0.0033801634933784052, 'CG': 0.010003006650292771, 'TT': -0.006563952473581737, 'GG': 0.01737341520832093, 'GC': 0.01922497123078911, 'AT': -0.0057390232570493785, 'GA': -0.005617342683749458, 'TG': -0.0033801634933784074, 'CT': -0.001794683526960448, 'TC': -0.005617342683749461, 'TA': -0.006544881954959417}
K_list = [K1, K2, K3, K4]

# Nucleosomal DNA length
NCPlen = 147
#NCPlen = 81

# JJ model sanity check
sanity_check=False

def key_cmp (key1, key2):
    if float(key1) < float(key2):
        return -1
    else:
        return 1

def rev_comp(seq):
    rev = ""
    for nt in seq[::-1]:
        nt = nt.upper()
        if nt == 'A':
            rev += 'T'
        elif nt == 'T':
            rev += 'A'
        elif nt == 'G':
            rev += 'C'
        elif nt == 'C':
            rev += 'G'
    return rev

def angle_to_radian(angle):
    return math.pi*float(angle)/180.0

def randian_to_angle(radian):
    return float(180*radian)/math.pi
    
def omega (index, phi=0):
    angle = (float(math.pi) / L)*index - 0.5*(phi - 1.5)
    return angle_to_radian(9)*(math.sin(angle))**2

def get_energy (seq, K, phi=0, circular=False):
    energy = 0.0
    if circular:
        sequence = seq + seq[0]
    else:
        sequence = seq
    for i in range(len(sequence)-1):
        din = sequence[i:i+2]
        k = K[din]
        omg = omega(i, phi)
        #omg = omega(i-len(sequence)/2, phi)
        energy += k*(omg**2)
    return energy

def get_energy_NCP (seq, K, phi=math.pi+1.5, circular=False):
    energy = 0.0
    for i in range(len(seq)-1):
        din = seq[i:i+2]
        k = K[din]
        #omg = omega(i, phi)
        omg = omega(i-len(seq)/2, phi)
        energy += k*(omg**2)
    return energy

def get_JJ_energy (seq, K, min_angle=0, max_angle=360, step=5):
    energy_list = []
    for angle in range(min_angle, max_angle+1, step):
        phi = angle_to_radian(angle)
        energy_list.append(get_energy(seq, K, phi=phi))
    return min(energy_list)

def mask (seq, idxs, type="skips", space=""):
    new_seq = ""
    idxs = sorted(idxs)
    if type == "skips":
        new_seq += seq[0:idxs[0]] + space
        for i in range(len(idxs)-1):
            st, ed = idxs[i]+1, idxs[i+1]
            if st == ed:
                continue
            new_seq += seq[st:ed] + space
        new_seq += seq[idxs[-1]+1:]
    elif type == "includes":
        for i in range(len(idxs)-1):
            idx = idxs[i]
            next_idx = idxs[i+1]
            new_seq += seq[idx]
            if next_idx - idx > 1:
                new_seq += space
        new_seq += seq[next_idx] 
    return new_seq

# JJ model sanity check
if sanity_check:
    # check NGO sequence energy
    NGO_seq = 'CAGAATCCGTGCTAGTACCTCAATATAGACTCCCTCCGGTGCCGAGGCCGCTCAATTGGTCGTAGGACTATCCTCACCTCCACCGTTTCA'
    X, Y = [], []
    for i in range(60, 421):
        X.append(i)
        phi = angle_to_radian(i)
        Y.append(get_energy(NGO_seq, K1, phi=phi))

    fig = plt.figure()
    plt.plot(X, Y)
    plt.ylim([28, 38])
    #plt.show()
    plt.close()

    # check Widom 601 sequence energy
    widom_seq = 'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGC'
    widom_LH = widom_seq[:len(widom_seq)/2]
    widom_RH = widom_seq[len(widom_seq)/2+1:]

    X, Y = [], []
    YL, YR = [], []
    for i in range(-180, 181, 30):
        X.append(i)
        phi = angle_to_radian(i)
        Y.append(get_energy(widom_seq, K1, phi=phi))
        YL.append(get_energy(widom_LH, K1, phi=phi))
        YR.append(get_energy(widom_RH, K1, phi=phi))

    fig = plt.figure()
    plt.subplot(2,1,1)
    plt.plot(X, Y, 'k', label='601')
    plt.ylim([48,54])
    plt.legend()
    plt.subplot(2,1,2)
    plt.plot(X, YL, 'r', label='601LH')
    plt.plot(X, YR, 'b', label='601RH')
    plt.ylim([22,28])
    plt.legend()
    #plt.show()
    plt.close()

    # check the NCP positioning on 601 40N40 template
    template = "ATCCGACTGGCACCGGCAAGGTCGCTGTTCAATACATGCACAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCAGGGCGGCCGCGTATAGGGTCCATCACATAAGGGATGAACT"

    energy_list = []
    for i in range(NCPlen/2, len(template)-NCPlen/2):
        NCPseq = template[i-NCPlen/2:i+NCPlen/2+1]
        energy_list.append(get_energy_NCP(NCPseq, K1))

    fig = plt.figure()
    plt.plot(range(len(energy_list)), energy_list)
    plt.axvline(x=len(energy_list)/2)
    plt.show()
    plt.close()
    

    # check yeast genome nucleosome occupancy vs energy prediction
    Roman_dict = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16}
    def read_genome(fname):
        chr_seq = {}
        chr_name, sequence = "", ""
        for line in open(fname):
            if line.startswith(">"):
                if chr_name and sequence:
                    chr_seq[chr_name] = sequence
                chr_name = line.strip().split()[0][1:]
                if chr_name.startswith('chr'):
                    try:
                        chr_name = "chr" + str(Roman_dict[chr_name[3:]])
                    except:
                        None
                sequence = ""
            else:
                sequence += line.strip()
        if chr_name and sequence:
            chr_seq[chr_name] = sequence
        return chr_seq
    chr_seq = read_genome("sacCer2.fa")
    
    def read_occfile (fname, chr_seq, chr_choices=['chr1']):
        chr_occ = {}
        for line in open(fname):
            cols = line.strip().split()
            chr, st, ed, values = cols
            chr = 'chr' + chr
            if chr not in chr_choices:
                continue
            st, ed = int(st), int(ed)
            if chr not in chr_occ:
                chr_occ[chr] = [np.nan]*len(chr_seq[chr])
            values = values.strip().split(';')
            assert len(values) == ed-st+1
            for i in range(len(values)):
                chr_occ[chr][i+st] = float(values[i])
        return chr_occ
    chr_occ = read_occfile("InVitro_dMean.chv", chr_seq)

    chr = "chr1"
    st, ed = 100000, 120000

    energy_list = []
    for pos in range(st, ed):
        NCPseq = chr_seq[chr][pos-NCPlen/2:pos+NCPlen/2+1]
        JJ_energy = get_JJ_energy(NCPseq, K1)
        energy_list.append(JJ_energy)

    occ_list = chr_occ[chr][st:ed]

    fig, ax1 = plt.subplots(figsize=(10,5))
    ax1.plot(range(len(energy_list)), energy_list, 'r', label='JJ predict')
    ax1.set_xlabel("Position (bp)")
    ax1.set_ylabel("JJ predicted Energy (A.U.)", color='r')
    ax2 = ax1.twinx()
    ax2.plot(range(len(occ_list)), occ_list, 'b', label='NCP occupancy')
    ax2.set_ylabel("Nucleosome occupancy", color='b')
    plt.show()
    plt.close()

    fig = plt.figure()
    plt.plot(energy_list, occ_list, '.', alpha=0.1)
    plt.xlabel("JJ predicted Energy (A.U.)")
    plt.ylabel("Nucleosome occupancy")
    plt.show()
    plt.close()
                

# read slider data
key_slider1 = pickle.load(open("slider1.p", "rb"))
key_slider2 = pickle.load(open("slider2.p", "rb"))

keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
keys = sorted(keys, cmp=key_cmp)


# JJ model vs Linear model
mask_type, mask_idxs = None, None

#mask_type, mask_idxs = "includes", range(NCPlen/2-70-3, NCPlen/2-70+3) + range(NCPlen/2+70-3, NCPlen/2+70+3)
#space = 'G'

#mask_type, mask_idxs = "skips", range(NCPlen/2-70-3, NCPlen/2-70+3) + range(NCPlen/2+70-3, NCPlen/2+70+3)
#space = 'G'

m1 = EnModel.EnergyModel(key_slider1, NCPlen=NCPlen, mask_type=mask_type, mask_idxs=mask_idxs)
#m1.train(MM_orders=[1], Kmer_k_b=False, PolyA_b=False,  GC_b=False, Harmonic=False, graph=True)
m2 = EnModel.EnergyModel(key_slider2, NCPlen=NCPlen, mask_type=mask_type, mask_idxs=mask_idxs)
#m2.train(MM_orders=[1], Kmer_k_b=False, PolyA_b=False,  GC_b=False, Harmonic=False, graph=True)

print "done"

exp_list1, exp_list2 = [], []
JJ_pred_matrix = [ [] for i in range(len(K_list)) ] # JJ model with different parameters
pred_list1, pred_list2 = [], [] # Linear regression model
for key in keys:
    slider1 = key_slider1[key]
    slider2 = key_slider2[key]
    energy_profile1 = slider1.energy_profile()
    energy_profile2 = slider2.energy_profile()
    seq = slider1.seq
    for i in range(NCPlen/2, len(seq)-NCPlen/2):
        NCP_seq = seq[i-NCPlen/2:i+NCPlen/2+1]
        for j in range(len(K_list)):
            #JJ_pred_matrix[j].append(get_JJ_energy(NCP_seq, K_list[j]))
            JJ_pred_matrix[j].append(get_energy_NCP(NCP_seq, K_list[j]))
        #NCPseq = mask(NCPseq, mask_idxs, type=mask_type, space='G')
        #pred_energy1 = m1.energy_predict(NCP_seq)
        #pred_energy2 = m2.energy_predict(NCP_seq)
        exp_list1.append(energy_profile1[i])
        exp_list2.append(energy_profile2[i])
        #pred_list1.append(pred_energy1)
        #pred_list2.append(pred_energy2)

name_list = ['JJ MD parameters', 'Loop-seq parameters', "Slide-seq parameters (HS)", "Slide-seq parameters (Chd1)"]
for i in range(len(K_list)):
    fig = plt.figure()
    plt.title("JJ model energy prediction with " + name_list[i] )
    plt.plot(exp_list1, JJ_pred_matrix[i], '.', alpha=0.1, label="Heat Shift")
    plt.plot(exp_list2, JJ_pred_matrix[i], '.', alpha=0.1, label="Chd1 Sliding")
    plt.xlabel("Measured Energy (A.U.)")
    plt.ylabel("Predicted Energy (A.U.)")
    leg = plt.legend()
    for lh in leg.legendHandles:
        lh._legmarker.set_alpha(1)                    
    plt.savefig("JJmodel_" + name_list[i] + ".png", bbox_inches='tight')
    #plt.show()
    plt.close()


"""
fig = plt.figure()
plt.title("Linear model energy prediction")
plt.plot(exp_list1, pred_list1, '.', alpha=0.2, label="Heat Shift")
plt.plot(exp_list2, pred_list2, '.', alpha=0.2, label="Chd1 Sliding")
plt.xlabel("Measured Energy (A.U.)")
plt.ylabel("Predicted Energy (A.U.)")
leg = plt.legend()
for lh in leg.legendHandles:
    lh._legmarker.set_alpha(1)                    
plt.savefig("Linmodel.png", bbox_inches='tight')
#plt.show()
plt.close()



# check one sequence of yeast
pdf = matplotlib.backends.backend_pdf.PdfPages("JJ.pdf")

#random.seed(42)
for key in keys[:50]:
    #key = random.choice(keys)
    slider1 = key_slider1[key]
    slider2 = key_slider2[key]
    seq = slider1.seq

    JJ_profile1 = [np.nan for i in range(NCPlen/2)]
    JJ_profile2 = [np.nan for i in range(NCPlen/2)]
    for i in range(NCPlen/2, len(seq)-NCPlen/2):
        NCPseq = seq[i - NCPlen/2 : i + NCPlen/2 + 1]
        JJ_profile1.append(get_JJ_energy(NCPseq, K1))
        JJ_profile2.append(get_JJ_energy(NCPseq, K2))
    JJ_profile1 += [np.nan for i in range(NCPlen/2)]
    JJ_profile2 += [np.nan for i in range(NCPlen/2)]

    energy_profile1 = slider1.energy_profile()
    energy_profile2 = slider2.energy_profile()
    energy_profile1 = [np.nan for i in range(10)] + list(energy_profile1[10:len(energy_profile1)-10]) + [np.nan for i in range(10)]
    energy_profile2 = [np.nan for i in range(10)] + list(energy_profile2[10:len(energy_profile1)-10]) + [np.nan for i in range(10)]

    mid = 225/2
    energy_profile1 = energy_profile1[mid-40:mid+41]
    energy_profile2 = energy_profile2[mid-40:mid+41]
    JJ_profile1 = JJ_profile1[mid-40:mid+41]
    JJ_profile2 = JJ_profile2[mid-40:mid+41]

    fig = plt.figure()
    
    plt.subplot(2, 1, 1)
    ax1 = plt.gca()
    #fig, ax1 = plt.subplot(figsize=(10, 5))
    ax1 = plt.gca()
    ax1.plot(range(len(energy_profile1)), energy_profile1, 'b', label='Heat Shift')
    ax1.plot(range(len(energy_profile2)), energy_profile2, 'g', label='Chd1 Sliding')
    ax1.set_xlabel("Position (bp)")
    ax1.set_ylabel("Measured Energy (A.U.)", color='b')
    ax1.tick_params('y', colors='k')
    plt.legend(loc='upper left')
    ax2 = ax1.twinx()
    ax2.plot(range(len(JJ_profile1)), JJ_profile1, 'r', label='JJ')
    ax2.set_ylabel("JJ predicted Energy (A.U.)", color='r')
    ax2.tick_params('y', colors='r')
    #if xtick_loc_name:
    #    xtick_locs, xtick_names = xtick_loc_name
    #    ax2.set_xticks(xtick_locs)
    #    ax2.set_xticklabels(xtick_names)
    plt.title(key)
    plt.legend(loc='upper right')
    #plt.savefig("JJ_" + str(k) + ".png",bbox_inches='tight')
    #plt.savefig("AK_" + str(k) + ".png",bbox_inches='tight')
    #plt.show()

    plt.subplot(2, 1, 2)
    ax1 = plt.gca()
    #fig, ax1 = plt.subplot(figsize=(10, 5))
    ax1 = plt.gca()
    ax1.plot(range(len(energy_profile1)), energy_profile1, 'b', label='Heat Shift')
    ax1.plot(range(len(energy_profile2)), energy_profile2, 'g', label='Chd1 Sliding')
    ax1.set_xlabel("Position (bp)")
    ax1.set_ylabel("Measured Energy (A.U.)", color='b')
    ax1.tick_params('y', colors='k')
    plt.legend(loc='upper left')
    ax2 = ax1.twinx()
    ax2.plot(range(len(JJ_profile2)), JJ_profile2, 'r', label='Aakash')
    ax2.set_ylabel("JJ predicted Energy (A.U.)", color='r')
    ax2.tick_params('y', colors='r')
    #if xtick_loc_name:
    #    xtick_locs, xtick_names = xtick_loc_name
    #    ax2.set_xticks(xtick_locs)
    #    ax2.set_xticklabels(xtick_names)
    plt.title(key)
    plt.legend(loc='upper right')

    fig.subplots_adjust(hspace=.5)
    pdf.savefig(fig)
    plt.close()
pdf.close()
"""
