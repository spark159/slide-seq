import load
import pickle
import EnModel
import matplotlib.pyplot as plt

ref_length = 225
dyad_axis = 225/2
dyad_offset = 52

# load files
ref_name = "/home/spark159/scripts/slide_seq/plusonelib.ref"
filenames1 = ["/home/spark159/../../media/spark159/sw/all_slide_seq_data/Plslib-HS_S1_L001_R.sort"]
filenames2 = ["/home/spark159/../../media/spark159/sw/all_slide_seq_data/Plslib-HS-30min_S2_L001_R.sort"]
#key_slider1 = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, key_choice='mapped_id', choice=['valid'], filter_num=0, fill=False, load_ref=ref_name, shape_fname=False)
#key_slider2 = load.load_files(filenames2, ref_length, dyad_axis, dyad_offset, key_choice='mapped_id', choice=['valid'], filter_num=0, fill=False, load_ref=ref_name, shape_fname=False)
#pickle.dump(key_slider1, open("slider1.p", "wb"))
#pickle.dump(key_slider2, open("slider2.p", "wb"))

#key_slider1 = pickle.load(open("slider1.p", "rb"))
#key_slider2 = pickle.load(open("slider2.p", "rb"))

#model1 = EnModel.EnergyModel(key_slider1)
#model2 = EnModel.EnergyModel(key_slider2)

#pickle.dump(model1, open("model1.p", "wb"))
#pickle.dump(model2, open("model2.p", "wb"))

model1 = pickle.load(open("model1.p", "rb"))
model2 = pickle.load(open("model2.p", "rb"))

model1.report(MM_orders=False, Kmer_k_b=[2, 1], PolyA_b=False, GC_b=False, Harmonic=False)
model2.report(MM_orders=False, Kmer_k_b=[2, 1], PolyA_b=False, GC_b=False, Harmonic=False)

# position independent dinucleotide frequency
freq1 = model1.freq['Kmer0']
freq2 = model2.freq['Kmer0']

freq_Kmer1, freq_Kmer2 = [], []
for Kmer, freq in freq1.items():
    freq_Kmer1.append([freq, Kmer])
for Kmer, freq in freq2.items():
    freq_Kmer2.append([freq, Kmer])

freq_Kmer1 = sorted(freq_Kmer1)
freq_Kmer2 = sorted(freq_Kmer2)

Kmers2, freqs2 = [], []
freqs1 = []
for freq, kmer in freq_Kmer2:
    Kmers2.append(kmer)
    freqs2.append(freq)
    freqs1.append(freq1[kmer])

fig = plt.figure()
plt.plot(range(len(freqs2)), freqs2, 'rx-', label='Chd1 Sliding')
plt.plot(range(len(freqs1)), freqs1, 'bx-', label='Heat Shift')
plt.xticks(range(len(freqs2)), Kmers2)
plt.ylabel("Relative frequency")
plt.legend()
plt.show()
