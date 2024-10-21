import matplotlib.pyplot as plt
from SliderClass import Slider
import graph
import sample
import load


directory = "/home/spark159/../../media/spark159/sw/mmlibIDlibFinal/"
#sample_mode = "polyA:1"
mtype_choice = ["I"]

for polyA_len in [1]:
    sample_mode = "polyA:" + str(polyA_len)
    for library in ["IDlib"]:
        for rep in ["1"]:
            for condition in ['control', "bubble"]:
                for time in ["0", "5"]:
                    filename = directory + "%s_%s_%s_%srep_.combined.sort" % (library, condition, time, rep)
                    filenames1 = [filename]
                    key_slider = load.load_files(filenames1, ref_length=225, dyad_axis=225/2, dyad_offset=52, filter_num = 0,  fill="linear", load_ref = None, mtype_choice=mtype_choice)

                    sample_list = sample.sampling(key_slider, sample_mode)
                    #print sample_list
                    #sample_list[0].remove('A-164')
                    #note = "%s_%s_%s_%srep_%s" % (library, condition, time, rep, sample_mode)
                    note = "%s_%s_%s_%srep_%s" % ("".join(mtype_choice), condition, time, rep, sample_mode)
                    print note
                    graph.plot_map(key_slider, sample_list, norm_choice=True, note=note, draw_key=True, draw_vert=False)
                    graph.plot_signal(key_slider, sample_list, note=note)
                    print


"""
directory = "/home/spark159/../../media/spark159/sw/mmlibIDlibQC/IDlib/"

filenames1 = [directory + "IDlib-5min_S5_L001_R.combined.sort"]
key_slider = load.load_files(filenames1, ref_length=225, dyad_axis=225/2, dyad_offset=52, no_spikes=False, filter_num = 0,  fill=None, load_ref = None)
note = "IDlib_5"
sample_list = sample.sampling(key_slider, "polyA:1")
graph.plot_map(key_slider, sample_list, norm_choice=True, note=note, draw_key=True, draw_vert=False)
graph.plot_signal(key_slider, sample_list, note=note)

filenames1 = [directory + "IDlib-bubble-5min_S7_L001_R.combined.sort"]
key_slider = load.load_files(filenames1, ref_length=225, dyad_axis=225/2, dyad_offset=52, no_spikes=False, filter_num = 0,  fill=None, load_ref = None)
note = "IDlib_5_bubble"
sample_list = sample.sampling(key_slider, "polyA:1")
graph.plot_map(key_slider, sample_list, norm_choice=True, note=note, draw_key=True, draw_vert=False)
graph.plot_signal(key_slider, sample_list, note=note)
"""
