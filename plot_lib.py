import matplotlib.pyplot as plt
import graph
import sample
from SliderClass import Slider

def open_ref (fname):
    key_slider = {}
    for line in open(fname):
        line = line.strip()
        if line.startswith('>'):
            key = line[1:]
            loc, type, nts = key.split('-')
            new_key = nts + '-' + loc
            continue
        if type != 'I':
            continue
        key_slider[new_key] = Slider (new_key,
                                      225,
                                      225/2,
                                      52,
                                      52,
                                      'N',
                                      [1]*225,
                                      [1]*225,
                                      [1]*225,
                                      None,
                                      None,
                                      None,
                                      None)
    return key_slider



#key_slider = open_ref("polyABubble_NoBB.ref")
#sample_mode = 'polyA:1-2-3-4-5-6-7-8-9-10'
key_slider = open_ref("singleInDel_NoBB.fa")
sample_mode = 'polyA:1'
sample_list = sample.sampling(key_slider, sample_mode)
graph.plot_map(key_slider, sample_list, norm_choice=True, note='Insert', draw_key=True, draw_vert=True)
                                      
