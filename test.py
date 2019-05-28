import load
import numpy as np
import EnModel
import matplotlib.pyplot as plt
import graph_edit
from SliderClass import Slider

# read pluse-one library
ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
filenames1 = ["../../Illumina/plusoneHS/data/Plslib-HS_S1_L001_R.sort"]
filenames2 = ["../../Illumina/plusoneHS/data/Plslib-HS-30min_S2_L001_R.sort"]
ref_fname = "../../Illumina/plusoneHS/plusonelib.ref"
key_slider1 = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, filter_num = 50, fill=None, load_ref=ref_fname, shape_fname='phpoXp6Za')
key_slider2 = load.load_files(filenames2, ref_length, dyad_axis, dyad_offset, filter_num = 50, fill=None, load_ref=ref_fname, shape_fname='phpoXp6Za')
filenames1 = ["../../Illumina/SW_Ascan_new/data/Ascan0_S1_L001_R.sort"]
filenames2 = ["../../Illumina/SW_Ascan_new/data/Ascan-5min_S1_L001_R.sort"]
ref_fname = "../../Illumina/SW_Ascan_new/polyAscanlib.ref"
#key_slider3 = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, filter_num = 50, fill='linear', load_ref=ref_fname)
#key_slider4 = load.load_files(filenames2, ref_length, dyad_axis, dyad_offset, filter_num = 50, fill='linear', load_ref=ref_fname)


#temp = {}
#for i in range(100):
#    temp[i] = key_slider1[str(i)]

m1 = EnModel.EnergyModel(key_slider1, bound=2, shape=True)
#m1 = EnModel.EnergyModel(key_slider1)
#m1.train(MM_orders=[0], Kmer_k_b=False, PolyA_b=False,  GC_b=False, Harmonic=False, k_fold=10)
#m1.report(MM_orders=[0], Kmer_k_b=False, PolyA_b=False,  GC_b=False, Harmonic=False)

m2 = EnModel.EnergyModel(key_slider2, bound=2, shape=True)
#m2.train(MM_orders=[0], Kmer_k_b=False, PolyA_b=False,  GC_b=False, Harmonic=True, k_fold=10)
#m3 = EnModel.EnergyModel(key_slider3)
#m4 = EnModel.EnergyModel(key_slider4)
