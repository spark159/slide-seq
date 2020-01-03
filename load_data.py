from SliderClass import Slider
import sample
import graph

def read_datafile (fname):
    para_value = {'ref_length':None,
                  'dyad_axis':None,
                  'left_offset':None,
                  'right_offset':None}
    key_slider = {}
    First = True
    for line in open(fname):
        if not line.strip():
            continue
        if line.startswith('@'):
            cols = line[1:].split()
            if len(cols) <= 1:
                info = cols[0]
                continue
            para, value = cols
            try:
                value = int(value)
            except:
                pass
            para_value[para] = value
            continue
        if line.startswith('>'):
            if not First:
                key_slider[key] = Slider(id = key,
                                         ref_length = para_value['ref_length'],
                                         dyad_axis = para_value['dyad_axis'],
                                         left_offset = para_value['left_offset'],
                                         right_offset = para_value['right_offset'],
                                         seq = info_value['Sequence'],
                                         dyadmap = info_value['PositioningSignal'],
                                         left_cutmap = info_value['BottomCleavageCounts'],
                                         right_cutmap = info_value['TopCleavageCounts'],
                                         MGW = info_value['MinorGrooveWidth'],
                                         HelT = info_value['HelixTwist'],
                                         ProT = info_value['PropellerTwist'],
                                         Roll = info_value['Roll'])
            key = line[1:].strip()
            info_value = {'Sequence':None,
                          'PositioningSignal':None,
                          'TopCleavageCounts':None,
                          'BottomCleavageCounts':None,
                          'MinorGrooveWidth':None,
                          'HelixTwist':None,
                          'PropellerTwist':None,
                          'Roll':None}
            First = False
            continue
        if line.startswith('['):
            line = line.strip()
            value_array = [float(value) for value in line[1:-1].split(',')]
            info_value[info] = value_array
            continue
        info_value[info] = line.strip()

    key_slider[key] = Slider(id = key,
                             ref_length = para_value['ref_length'],
                             dyad_axis = para_value['dyad_axis'],
                             left_offset = para_value['left_offset'],
                             right_offset = para_value['right_offset'],
                             seq = info_value['Sequence'],
                             dyadmap = info_value['PositioningSignal'],
                             left_cutmap = info_value['BottomCleavageCounts'],
                             right_cutmap = info_value['TopCleavageCounts'],
                             MGW = info_value['MinorGrooveWidth'],
                             HelT = info_value['HelixTwist'],
                             ProT = info_value['PropellerTwist'],
                             Roll = info_value['Roll'])

    return key_slider

#key_slider = read_datafile("IDlib_bubble_0_1rep_.combined.data")
key_slider = read_datafile("IDlib_bubble_0_.data")
sample_mode = 'polyA:1'
sample_list = sample.sampling(key_slider, sample_mode)
graph.plot_map(key_slider, sample_list, norm_choice=True, note='check', draw_key=True, draw_vert=False)
