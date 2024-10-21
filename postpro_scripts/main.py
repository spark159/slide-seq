import sys
import code
from argparse import ArgumentParser, FileType
import numpy as np
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from SliderClass import Slider
import math
import sample
import analysis
import graph_edit
import graph
import load
import random
import readLoopseq
from sklearn import linear_model
import Nucleosome_positioning as nuc

def display_graph (filenames1,
                   filenames2,
                   filenames3,
                   ref_fname,
                   ref_length,
                   dyad_offset,
                   stat_mode,
                   graph_mode,
                   norm_choice,
                   sample_mode):
    
    dyad_axis = ref_length/2  

    # read files and load data
    key_slider = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, filter_num = 0,  fill=None, load_ref = ref_fname)
    #key_slider1 = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None, load_ref = ref_fname)
    #key_slider2 = load.load_files(filenames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None, load_ref = ref_fname)

    """
    # save positioning data
    def tuple_cmp(a,b):
        win1, st1 = a.split('-')[0], a.split('-')[1]
        win2, st2 = b.split('-')[0], b.split('-')[1]
        if len(win1) < len(win2):
            return -1
        elif len(win1) > len(win2):
            return 1
        else:
            if int(st1) < int(st2):
                return -1
            elif int(st1) > int(st2):
                return 1
            else:
                return 0
        
    f = open("AscanKinetics_60min.txt", 'w')
    
    for key in sorted(key_slider.keys(), cmp=tuple_cmp):
        slider = key_slider[key]
        dyadmap = slider.dyadmap
        print >> f, '>%s' % (key)
        s = ""
        for i in range(len(dyadmap)):
            value = dyadmap[i]
            if i == 0:
                s += str(value)
                continue
            s += "," + str(value)
        print >> f, s
    f.close()

    
    # left/right scattering plot
    
    left, right = [], []
    for key in key_slider:
        win, loc = key.split('-')
        size, loc = len(win), int(loc)
        if size % 2 != 0:
            loc = loc + size/2
        else:
            loc = loc + size/2 - 0.5
        if loc < dyad_axis:
            left.append(key)
        elif loc > dyad_axis:
            right.append(key)
    sample_list = [left, right]
    
    graph.plot_corr2(key_slider, Slider.Amer_len, Slider.median_pos, xlabel='Poly-A length', ylabel='Mean position', sample_labels = ['PolyA in Left','PolyA in Right'], sample_list=sample_list)
   
    
    # mean position heat map
    keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))

    map1 = [[np.nan for i in range(ref_length)] for i in range(3,26)]
    map2 = [[np.nan for i in range(ref_length)] for i in range(3,26)]
    map3 = [[np.nan for i in range(ref_length)] for i in range(3,26)]

    for key in keys:
        slider1, slider2 = key_slider1[key], key_slider2[key]
        win, loc = key.split('-')
        size, loc = len(win), int(loc)
        mean_pos1 = slider1.mean_pos()
        mean_pos2 = slider2.mean_pos()
        mean_pos3 = slider2.mean_pos() - slider1.mean_pos()
        #mean_pos = slider2.find_peaks(choice='dyad', num=1)[0][0] - slider1.find_peaks(choice='dyad', num=1)[0][0]
        #mean_pos = slider.find_peaks(choice='dyad', num=1)[0][0] - dyad_axis
        map1[size-3][loc] = mean_pos1
        map2[size-3][loc] = mean_pos2
        map3[size-3][loc] = mean_pos3

    maps = [map1, map2, map3]


    #for i in range(len(maps)):
    #    map = maps[i]
    #    fig = plt.figure()
    #    ax1 = fig.add_subplot(111)
    #    ax1.set_yticks(range(0,23))
    #    ax1.set_yticklabels([str(i) for i in range(3,26)])
    #    ax1.set_xticks([i+0.5 for i in range(0,225,16)])
    #    ax1.set_xticklabels([str(i) for i in range(1,226,16)])
    #    ax1.set_xlabel('Poly-A start location (bp)')
    #    ax1.set_ylabel('Poly-A length (bp)')
    #    ax2 = ax1.twiny()
    #    cmap = cm.seismic
    #    cmap.set_bad('white',1.)
    #    ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-20, vmax=20)
    #    #ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-5, vmax=5)
    #    #ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-100, vmax=100)
    #    #ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-200, vmax=200)
    #    for x in [dyad_axis-30, dyad_axis, dyad_axis+30]:
    #        ax2.axvline(x, color='k', linestyle='--', linewidth=3)
    #    ax2.set_xticks([dyad_axis-30, dyad_axis, dyad_axis+30])
    #    ax2.set_xticklabels(['-30','Center','+30'])
    #    #plt.savefig('mean_position_heatmap_' + str(i) + '.png')
    #    plt.show()
    #    plt.close()

    
    # sensitivity map (mean position)
    smaps = []
    for u in range(len(maps)):
        smap = np.zeros((23, ref_length))
        smap.fill(np.nan)
        map = maps[u]
        for i in range(3,26):
            for j in range(ref_length):
                if np.isnan(map[i-3][j]):
                    continue
                for k in range(i):
                    if np.isnan(smap[i-3][j+k]):
                        smap[i-3][j+k] = 0.0
                    #smap[i-3][j+k] += map[i-3][j]/float(i)
                    smap[i-3][j+k] += map[i-3][j]
        smaps.append(smap)


    for u in range(len(smaps)):
        smap = smaps[u]
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_yticks(range(0,23))
        ax1.set_yticklabels([str(i) for i in range(3,26)])
        ax1.set_xticks([i+0.5 for i in range(0,225,16)])
        ax1.set_xticklabels([str(i) for i in range(1,226,16)])
        ax1.set_xlabel('DNA template (bp)')
        ax1.set_ylabel('Poly-A length (bp)')
        ax1.xaxis.label.set_size(20)
        ax1.yaxis.label.set_size(20)
        ax1.tick_params(axis='both', which='major', labelsize=15)
        ax2 = ax1.twiny()
        cmap = cm.seismic
        #cmap = cm.get_cmap('Spectral_r')
        cmap.set_bad('white',1.)
        #ax2.imshow(smap, cmap=cmap, interpolation='none', aspect='auto', vmin=-15, vmax=15)
        ax2.imshow(smap, cmap=cmap, interpolation='none', aspect='auto', vmin=-200, vmax=200)
        for x in [dyad_axis-30, dyad_axis, dyad_axis+30]:
            ax2.axvline(x, color='k', linestyle='--', linewidth=3)
        ax2.set_xticks([dyad_axis-30, dyad_axis, dyad_axis+30])
        ax2.tick_params(axis='both', which='major', labelsize=15)
        ax2.set_xticklabels(['-30','Center','+30'])
        plt.tight_layout()
        plt.savefig('mean_position_smap_' + str(u+1) + '.png')
        #plt.show()
        plt.close()
    
    
    # collapse the heat map (mean position)
    #lines = []
    #for u in range(len(maps)):
    #    line = np.zeros((1, ref_length))
    #    line.fill(np.nan)
    #    map = maps[u]
    #    for i in range(3,26):
    #        for j in range(ref_length):
    #            if np.isnan(map[i-3][j]):
    #                continue
    #            for k in range(i):
    #                if np.isnan(line[0][j+k]):
    #                    line[0][j+k] = 0.0
    #                line[0][j+k] += map[i-3][j]
    #    lines.append(line)

    #for i in range(len(lines)):
    #    line = lines[i]
    #    fig = plt.figure()
    #    ax1 = fig.add_subplot(111)
    #    ax2 = ax1.twiny()
    #    cmap = cm.seismic
    #    cmap.set_bad('white',1.)
    #    ax2.imshow(line, cmap=cmap, interpolation='none', aspect='auto', vmin=-3000, vmax=3000)
    #    for x in [dyad_axis-30, dyad_axis, dyad_axis+30]:
    #        ax2.axvline(x, color='k', linestyle='--', linewidth=3)
    #    ax2.set_xticks([dyad_axis-30, dyad_axis, dyad_axis+30])
    #    ax2.set_xticklabels(['-30','Center','+30'])
    #    #plt.imshow([line], interpolation='none', aspect='auto')
    #    plt.show()
    #    plt.close()
    
    
    
    # calculate effective force map
    force_map1, force_map2, force_map3 = {}, {}, {}
    energy_map1, energy_map2, energy_map3 = {}, {}, {}

    keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
    for key in keys:
        win, loc = key.split('-')
        size, loc = len(win), int(loc)

        slider1 = key_slider1[key]
        slider2 = key_slider2[key]
        dforce1 = slider1.force_profile()
        dforce2 = slider2.force_profile()
        dforce3 = dforce2 - dforce1
        denergy1 = slider1.energy_profile()
        denergy2 = slider2.energy_profile()
        denergy3 = denergy2 - denergy1
        
        for i in range(66, len(dforce1)-66):
            pos = loc - i
            new_key = (size, pos)
            if new_key not in force_map1:
                force_map1[new_key] = []
            if new_key not in force_map2:
                force_map2[new_key] = []
            if new_key not in force_map3:
                force_map3[new_key] = []
            if new_key not in energy_map1:
                energy_map1[new_key] = []
            if new_key not in energy_map2:
                energy_map2[new_key] = []
            if new_key not in energy_map3:
                energy_map3[new_key] = []
            force_map1[new_key].append(dforce1[i])
            force_map2[new_key].append(dforce2[i])
            force_map3[new_key].append(dforce3[i])
            energy_map1[new_key].append(denergy1[i])
            energy_map2[new_key].append(denergy2[i])
            energy_map3[new_key].append(denergy3[i])

    # mean force heat map
    fm_list = [force_map1, force_map2, force_map3]
    #map = [[np.nan for i in range(ref_length)] for i in range(3,26)]
    maps = []
    for k in range(len(fm_list)):
        #map_list = []
        map = np.zeros((23, ref_length))
        map.fill(np.nan)
        force_map = fm_list[k]
        for key in force_map:
            mean = np.mean(force_map[key])
            size, pos = key
            map[size - 3, pos + dyad_axis] = mean
        maps.append(map)
    
        #graph.plot_heatmap3(map_list, cmap, vmin = -30, vmax = 30, yticks=yticks, vlines=vlines, upticks=upticks, xlabels = ["", "Poly-A start location (bp)", ""], ylabels = ["Poly-A length (bp)"]*3, titles = ["Before", "After", "After-Before"], aspect=10, note="MDP")

    #for u in range(len(maps)):
    #    map = maps[u]
    #    fig = plt.figure()
    #    ax1 = fig.add_subplot(111)
    #    ax1.set_yticks(range(0,23))
    #    ax1.set_yticklabels([str(i) for i in range(3,26)])
    #    ax1.set_xticks([i+0.5 for i in range(0,225,16)])
    #    ax1.set_xticklabels([str(i) for i in range(-112, 113, 16)])
    #    ax1.set_xlabel('Location w.r.t. Dyad (bp)')
    #    ax1.set_ylabel('Poly-A length (bp)')
    #    ax2 = ax1.twiny()
    #    cmap = cm.seismic
    #    cmap.set_bad('white',1.)
    #    ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-10, vmax=10)
    #    #ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-0.4, vmax=0.4)
    #    for x in [dyad_axis-30, dyad_axis, dyad_axis+30]:
    #        ax2.axvline(x, color='k', linestyle='--', linewidth=3)
    #    ax2.set_xticks([dyad_axis-30, dyad_axis, dyad_axis+30])
    #    ax2.set_xticklabels(['-SHL3','Dyad','+SHL3'])
    #    #plt.savefig('mean_force_heatmap_' + str(i) + '.png')
    #    plt.show()
    #    plt.close()

    
    # sensitivity map (mean force)
    smaps = []
    for u in range(len(maps)):
        smap = np.zeros((23, ref_length))
        smap.fill(np.nan)
        map = maps[u]
        for i in range(3,26):
            for j in range(ref_length):
                if np.isnan(map[i-3][j]):
                    continue
                for k in range(i):
                    if np.isnan(smap[i-3][j+k]):
                        smap[i-3][j+k] = 0.0
                    #smap[i-3][j+k] += map[i-3][j]/float(i)
                    smap[i-3][j+k] += map[i-3][j]
        smaps.append(smap)

    for u in range(len(smaps)):
        map = smaps[u]
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_yticks(range(0,23))
        ax1.set_yticklabels([str(i) for i in range(3,26)])
        ax1.set_xticks([i+0.5 for i in range(0,225,16)])
        ax1.set_xticklabels([str(i) for i in range(-112, 113, 16)])
        ax1.set_xlabel('Location w.r.t. Dyad (bp)')
        ax1.set_ylabel('Poly-A length (bp)')
        ax2 = ax1.twiny()
        cmap = cm.seismic
        cmap.set_bad('white',1.)
        #cmap.set_under('lightcyan')
        #cmap.set_over('lightpink')
        ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-1, vmax=1)
        #ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-0.4, vmax=0.4)
        for x in [dyad_axis-30, dyad_axis, dyad_axis+30]:
            ax2.axvline(x, color='k', linestyle='--', linewidth=3)
        ax2.set_xticks([dyad_axis-30, dyad_axis, dyad_axis+30])
        ax2.set_xticklabels(['-SHL3','Dyad','+SHL3'])
        plt.savefig('mean_force_smap_' + str(u) + '.png')
        plt.show()
        plt.close()

    
    # collapes the heat map (mean force)
    map_list = []
    for k in range(len(fm_list)):
        position_map = np.zeros((1,ref_length))
        position_map.fill(np.nan)
        force_map = fm_list[k]
        for key in force_map:
            size, pos = key
            if size > 12:
                continue
            mean = np.mean(force_map[key])
            for i in range(pos + dyad_axis, pos + dyad_axis + size):
                if np.isnan(position_map[0, i]):
                    position_map[0, i] = 0.0
                position_map[0, i] += mean
        map_list.append(position_map)

    cmap = cm.seismic
    cmap.set_bad('black',1.)
    xticks = [range(0, ref_length, 16), [str(i) for i in range(-112, 113, 16)]]
    vlines = [dyad_axis - 65, dyad_axis - 20, dyad_axis, dyad_axis + 20, dyad_axis + 65]
    upticks = [ [dyad_axis - 65, dyad_axis - 20, dyad_axis, dyad_axis + 20, dyad_axis + 65], ['SHL6.5', 'SHL2', 'dyad', 'SHL2','SHL6.5'] ]

    graph.plot_heatmap4(map_list, dim = [3,1], cmap_list = [cmap]*3, vmin_list = [-10,-10,-5], vmax_list = [10,10,5], yticks_list = [[[],[]], [[],[]], [[],[]]], xticks_list = [[[],['']*len(xticks[1])],[[],['']*len(xticks[1])],xticks], vlines_list=[vlines]*3, upticks_list=[upticks,[[],[]], [[],[]]], ylabels = ["Before", "After", "After-Before"], aspect_list=[25]*3, note="MDP", xlabels=['','','Location w.r.t. Dyad (bp)'])
            
    #fig = plt.figure()
    #img = plt.imshow(position_map, cmap=cmap, interpolation='none', aspect='auto', vmin=-5, vmax=5)
    #img = plt.imshow(position_map, cmap=cmap, interpolation='none', aspect='auto', vmin=-10, vmax=10)
    #for x in [dyad_axis-30, dyad_axis, dyad_axis+30]:
    #    plt.axvline(x, color='k', linestyle='--', linewidth=2)
    #xlabels = [str(i) for i in range(-112, 113, 16)]
    #plt.xticks(range(0, ref_length, 16), xlabels)
    #plt.xlabel('Location w.r.t. Dyad (bp)')
    #plt.savefig("sensitivitymap.png")
    #plt.close()
    #plt.show()
    
    
    #linear regression
    def read_ref (ref_fname):
        key_seq = {}
        for line in open(ref_fname):
            line = line.strip()
            if line.startswith(">"):
                key = line[1:]
                continue
            else:
                assert key not in key_seq
                key_seq[key] = line
        return key_seq

    def Amer_len(seq, pos=False):
        num = []
        num_pos = {}
        i = 0
        while i < len(seq):
            if seq[i] in 'AT':
                nt = seq[i]
                count = 1
                j = i + 1
                while j < len(seq):
                    if seq[j] != nt:
                        break
                    count +=1
                    j +=1
                num.append(count)
                if count not in num_pos:
                    num_pos[count] = []
                num_pos[count].append(i)
                i = j
            else:
                i +=1
        if pos:
            return num_pos
        if len(num) == 0:
            return 0
        return max(num)

    key_seq = read_ref("polyAscanlib.ref")
    keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))

    y1, y2, y3 = [], [], []
    X = []
    smX = []
    
    for key in keys:
        slider1 = key_slider1[key]
        slider2 = key_slider2[key]
        mpos1 = slider1.mean_pos()
        mpos2 = slider2.mean_pos()
        mpos3 = mpos2 - mpos1
        y1.append(mpos1)
        y2.append(mpos2)
        y3.append(mpos3)
        seq = key_seq[key]
        Anum_pos = Amer_len(seq, pos=True)
        x = []
        for i in range(3, 26):
            temp = [0]* ref_length
            try:
                Apos = Anum_pos[i]
            except:
                x += temp
                continue
            for j in range(ref_length):
                if j not in Apos:
                    continue
                for u in range(j, j+i):
                    temp[u] = 1
            x += temp
        X.append(x)

    reg1 = linear_model.Ridge (alpha = .5)
    reg2 = linear_model.Ridge (alpha = .5)
    reg3 = linear_model.Ridge (alpha = .5)
    reg1.fit(X, y1)
    reg2.fit(X, y2)
    reg3.fit(X, y3)

    prey1 = reg1.predict(X)
    prey2 = reg2.predict(X)
    prey3 = reg3.predict(X)

    ys = [y1,y2,y3]
    preys = [prey1, prey2, prey3]
    titles = ["before", "after", "after-before"]

    for i in range(3):
        fig = plt.figure()
        plt.plot(ys[i], preys[i], 'k.')
        plt.xlabel("Experiment")
        plt.ylabel("Prediction")
        plt.title("Poly-A linear model: " + titles[i])
        plt.xlim([-20,30])
        plt.ylim([-20,30])
        plt.show()
        plt.close()
        

    regs = [reg1, reg2, reg3]

    maps = []
    for i in range(len(regs)):
        reg = regs[i]
        map = reg.coef_.reshape((23, ref_length))
        maps.append(map)


    for u in range(len(maps)):
        map = maps[u]
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_yticks(range(0,23))
        ax1.set_yticklabels([str(i) for i in range(3,26)])
        ax1.set_xticks([i+0.5 for i in range(0,225,16)])
        ax1.set_xticklabels([str(i) for i in range(1,226,16)])
        ax1.set_xlabel('DNA template (bp)')
        ax1.set_ylabel('Poly-A length (bp)')
        ax2 = ax1.twiny()
        cmap = cm.seismic
        cmap.set_bad('white',1.)
        ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-5, vmax=5)
        #ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=-200, vmax=200)
        for x in [dyad_axis-30, dyad_axis, dyad_axis+30]:
            ax2.axvline(x, color='k', linestyle='--', linewidth=3)
        ax2.set_xticks([dyad_axis-30, dyad_axis, dyad_axis+30])
        ax2.set_xticklabels(['-30','Center','+30'])
        #plt.savefig('mean_position_smap_' + str(u+1) + '.png')
        plt.show()
        plt.close()

    # linear model with fewer variables
    y1, y2, y3 = [], [], []
    smX = []
    
    for key in keys:
        slider1 = key_slider1[key]
        slider2 = key_slider2[key]
        mpos1 = slider1.mean_pos()
        mpos2 = slider2.mean_pos()
        mpos3 = mpos2 - mpos1
        y1.append(mpos1)
        y2.append(mpos2)
        y3.append(mpos3)
        seq = key_seq[key]
        Anum_pos = Amer_len(seq, pos=True)
        temp = [0] * ref_length
        for num, poslist in Anum_pos.items():
            for k in range(len(poslist)):
                pos = poslist[k]
                for i in range(pos, pos+num):
                    temp[i] = 1
        smX.append(temp)

    reg1 = linear_model.Ridge (alpha = .5)
    reg2 = linear_model.Ridge (alpha = .5)
    reg3 = linear_model.Ridge (alpha = .5)
    reg1.fit(smX, y1)
    reg2.fit(smX, y2)
    reg3.fit(smX, y3)

    prey1 = reg1.predict(smX)
    prey2 = reg2.predict(smX)
    prey3 = reg3.predict(smX)

    ys = [y1,y2,y3]
    preys = [prey1, prey2, prey3]
    titles = ["before", "after", "after-before"]

    for i in range(3):
        fig = plt.figure()
        plt.plot(ys[i], preys[i], 'k.')
        plt.xlabel("Experiment")
        plt.ylabel("Prediction")
        plt.title("Poly-A linear model: " + titles[i])
        plt.xlim([-20,20])
        plt.ylim([-20,20])
        plt.show()
        plt.close()
        
    regs = [reg1, reg2, reg3]

    lines = []
    for i in range(len(regs)):
        reg = regs[i]
        line = reg.coef_
        lines.append([line])

    cmap = cm.seismic
    cmap.set_bad('white',1.)
    xticks = [range(0, ref_length, 16), [str(i) for i in range(-112, 113, 16)]]
    vlines = [dyad_axis - 65, dyad_axis - 20, dyad_axis, dyad_axis + 20, dyad_axis + 65]
    upticks = [ [dyad_axis - 65, dyad_axis - 20, dyad_axis, dyad_axis + 20, dyad_axis + 65], ['SHL6.5', 'SHL2', 'dyad', 'SHL2','SHL6.5'] ]

    graph.plot_heatmap4(lines, dim = [3,1], cmap_list = [cmap]*3, vmin_list = [-10,-10,-5], vmax_list = [10,10,5], yticks_list = [[[],[]], [[],[]], [[],[]]], xticks_list = [[[],['']*len(xticks[1])],[[],['']*len(xticks[1])],xticks], vlines_list=[vlines]*3, upticks_list=[upticks,[[],[]], [[],[]]], ylabels = ["Before", "After", "After-Before"], aspect_list=[25]*3, note="MDP", xlabels=['','','Location w.r.t. Dyad (bp)'])
    """

    # plot colorbar
    #fig = plt.figure()
    #ax1 = fig.add_subplot(111)
    #ax1 = fig.add_axes([0.05, 0.80, 0.8, 0.15])
    #cmap = cm.seismic
    #norm = mpl.colors.Normalize(vmin=-200, vmax=200)
    #norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    #bounds = np.linspace(-1,1,5)
    #bounds = [-1., -.5, 0., .5, 1.]
    #cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal', extend='both', ticks=bounds)
    #cb1.set_label('Some Units')
    #note = 'mean_position'
    #note = 'mean_force'
    #plt.savefig('colorbar_' + note + '.png', bbox_inches='tight')
    #plt.show()
    #plt.close
        
    
    # select the subset of inserts
    #sample_list = sample.sampling(key_slider, sample_mode)
    #print sample_mode
    sample_list = [['601']]
    #sample_list1 = sample.sampling(key_slider1, sample_mode)
    #sample_list2 = sample.sampling(key_slider2, sample_mode)

    #sample_list = []
    #assert len(sample_list1) == len(sample_list2)
    #for i in range(len(sample_list1)):
    #    sample1 = sample_list1[i]
    #    sample2 = sample_list2[i]
    #    samplec = []
    #    for key in sample1:
    #        if key in sample2:
    #            samplec.append(key)
    #    sample_list.append(samplec)

    #keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
    #sample_list = [random.sample(keys, 30)]

    #key_sub = {}
    #for key in key_slider1.keys():
    #    slider1 = key_slider1[key]
    #    slider2 = key_slider2[key]
    #    key_sub[key] = slider2 - slider1
        
    
    #graph.plot_rmap(key_slider1, key_slider2, sample_list, norm_choice=False, note='_ratio_log', draw_key=True, draw_vert=False)
    # plot cut/dyad map

    #graph_edit.plot_map(key_slider1, sample_list, norm_choice=False, obs_func = Slider.eqm_flux, draw = False, slicing=0, note='Before')
    #graph_edit.plot_map(key_slider2, sample_list, norm_choice=False, obs_func = Slider.eqm_flux, draw = False, slicing=0, note='After')
    #graph_edit.plot_map(key_sub, sample_list, norm_choice=False, obs_func = Slider.get_dyadmap, draw ='polyA', slicing=0, note='Sub')

    #graph.plot_map(key_slider, sample_list, norm_choice=True, note='check', draw_key=False, draw_vert=False)
    #graph.plot_map(key_slider_r, sample_list, norm_choice=True, note='norm_raw', draw_key=True, draw_vert=False)
    #graph.plot_map(key_slider1, sample_list, norm_choice=True, note='_before', draw_key=True, draw_vert=False)
    #graph.plot_map(key_slider2, sample_list, norm_choice=True, note='_After', draw_key=True, draw_vert=False)
    
    # plot average cut/dyad signal
    graph.plot_signal(key_slider, sample_list, note='check')
    #graph.plot_signal(key_slider, sample_list, show_key=True)
    #graph.plot_signal(key_slider, note='test')
    #graph.plot_signal(key_slider_r, note='all_raw')
    #graph.plot_signal(key_slider1, note='before')
    #graph.plot_signal(key_slider2, note='after')
    
    # plot SeqID vs cut peaks
    #graph.plot_cpeaks(key_slider, left_peak_num = 2, right_peak_num = 2)

    # plot SeqID vs dyad peaks
    #graph.plot_dpeaks(key_slider, peak_num = 3, st_rank = 1, sample_list=sample_list)
    #graph.plot_dpeaks(key_slider, peak_num = 3, st_rank = 1, note='all')

    # plot correlation
    #graph.plot_corr2(key_slider2, Slider.Amer_len, Slider.mean_pos, xlabel='Poly-A length', ylabel='Mean position', sample_labels = ['Left','Right'], sample_list=sample_list)
    #graph.plot_corr3(key_slider1, key_slider2, Slider.Amer_len, Slider.mean_pos, xlabel='Poly-A length', ylabel='Mean position', sample_labels = ['Left','Right'], sample_list=sample_list)
    #graph.plot_corr(key_slider, Slider.Amer_len, Slider.mean_pos, xlabel='Poly-A length', ylabel='Mean position')
    #graph.plot_corr(key_slider, Slider.Amer_len, Slider.max_dis, xlabel='Poly-A length', ylabel='Maximum displacement')

    #graph.plot_energy(key_slider1, key_slider2)


    """
    # Entropy difference
    X, Y = [], []
    Z = []
    diff = []
    diff2 = []
    keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
    for key in keys:
        slider1 = key_slider1[key]
        slider2 = key_slider2[key]
        kde1 = slider1.KDE()
        kde2 = slider2.KDE()
        entropy1, entropy2 = 0.0, 0.0
        entropy3 = 0.0
        for i in range(len(kde1)):
            entropy1 += -kde1[i]*np.log(kde1[i])
            entropy2 += -kde2[i]*np.log(kde2[i])
            entropy3 += -kde2[i]*np.log(kde1[i])
        X.append(entropy1)
        Y.append(entropy2)
        Z.append(entropy3)
        diff.append(entropy2 - entropy1)
        diff2.append(entropy3 - entropy1)
    diff = sorted(diff)
    diff2 = sorted(diff2)    

    fig = plt.figure()
    plt.plot(X,Y,'.')
    plt.xlabel("Entropy (Before)")
    plt.ylabel("Entropy (After)")
    plt.xlim([4,5])
    plt.ylim([4,5])
    plt.plot([4,5],[4,5],'--')
    plt.savefig("entropy.png")
    plt.show()
    plt.close()

    fig = plt.figure()
    plt.plot(range(len(diff)), diff, '.')
    plt.xlabel("sequences")
    plt.ylabel("Entropy change")
    plt.savefig("entropy_change.png")
    plt.show()
    plt.close()

    
    fig = plt.figure()
    plt.plot(X,Z,'.')
    plt.xlabel("Entropy (Before)")
    plt.ylabel("Entropy (After)")
    plt.xlim([4,5])
    plt.ylim([4,5])
    plt.plot([4,5],[4,5],'--')
    plt.savefig("entropy_gen.png")
    plt.show()
    plt.close()

    fig = plt.figure()
    plt.plot(range(len(diff2)), diff2, '.')
    plt.xlabel("sequences")
    plt.ylabel("Entropy change")
    plt.savefig("entropy_change_gen.png")
    plt.show()
    plt.close()
    

    # Entropy vs flexibility
    key_flex = readLoopseq.get_flex("plusonelib.ref", "DataReady2.txt")
    #fig = plt.figure()
    #plt.plot(key_slider['0'].dyadmap)
    #plt.show()
    #fig = plt.figure()
    #plt.plot(key_slider['0'].KDE())
    #plt.show()
    
    X, Y = [], []
    Z = []
    L, R, E = [], [], []
    keyEntropy = []
    keyAsym = []
    for key in key_slider:
        slider = key_slider[key]
        entropy = slider.entropy(band_width=1)
        X.append(entropy)
        left, right = key_flex[key]['left'], key_flex[key]['right']
        Z.append(left-right)
        asym = abs(left-right)
        L.append(left)
        R.append(right)
        E.append(entropy)
        Y.append(asym)
        keyEntropy.append([key, entropy])
        keyAsym.append([key, asym])

    def temp_cmp (a, b):
        if a[1] < b[1]:
            return -1
        else:
            return 0

    keyEntropy = sorted(keyEntropy, cmp=temp_cmp, reverse=True)
    keyAsym = sorted(keyAsym, cmp=temp_cmp, reverse=True)

    key_Erank, key_Arank = {}, {}
    for i in range(len(keyEntropy)):
        key, _ = keyEntropy[i]
        print _
        key_Erank[key] = i+1
    for i in range(len(keyAsym)):
        key, _ = keyAsym[i]
        key_Arank[key] = i+1

    #X, Y = [], []
    #for key in key_Erank:
    #    X.append(key_Erank[key])
    #    Y.append(key_Arank[key])

    #fig = plt.figure()
    #plt.plot(X, Y, '.')
    #plt.plot(key_slider['97'].entropy(), abs(key_flex['97']['left'] - key_flex['97']['right']), 'x', label='Park-97')
    #plt.plot(key_Erank['97'], key_Arank['97'], 'x', label='Park-97')
    #plt.legend()
    #plt.title('Chd1 Sliding')
    #plt.xlabel('Entropy')
    #plt.ylabel('Asymmetricity')
    #plt.savefig('Chd1_entropyVSasym.png')
    #plt.show()
    #plt.close()

    fig = plt.figure()
    plt.scatter(L, R, c=E, s=3)
    #plt.legend()
    plt.title('Salt Gradient Dialysis')
    plt.xlabel('Left flexibility')
    plt.ylabel('Right flexibility')
    #plt.savefig('Chd1_entropyVSasym.png')
    plt.show()
    plt.close()

    code.interact(local=locals())
    
    
    
    keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))

    X1, Y1 = [], []
    X2, Y2 = [], []
    for key in keys:
        X1.append(key_slider1[key].weighted_GC())
        Y1.append(key_slider2[key].weighted_GC())
        X2.append(key_slider1[key].weighted_Amer_len())
        Y2.append(key_slider2[key].weighted_Amer_len())

    fig = plt.plot()
    plt.plot(X1, Y1, '.')
    plt.xlabel("GC mean (Before)")
    plt.xlabel("GC mean (After)")
    plt.show()
    plt.close()

    fig = plt.plot()
    plt.plot(X2, Y2, '.')
    plt.xlabel("poly-A length mean (Before)")
    plt.xlabel("poly-A length mean (After)")
    plt.show()
    plt.close()
    
    

    
    key_flex = readLoopseq.get_flex("plusonelib.ref", "DataReady2.txt")
    key_dAlen = readLoopseq.get_polyA("plusonelib.ref")
    key_dGC = readLoopseq.get_GC("plusonelib.ref")
    
    keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))

    X1, X2, X3 = [], [], []
    Y1, Y2, Y3 = [], [], []

    titles = ['GC', 'Alen', 'flex']
    labels = ["GC content:Right-left", "poly-A len:Right-left", "Right-Left flexibility"]
    key_props = [key_dGC, key_dAlen, key_flex]
    
    for i in range(len(titles)):
        label = labels[i]
        title = titles[i]
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
    
    """
    
    
    
if __name__ == '__main__':
    parser = ArgumentParser(description='graphical analysis of slide-seq')
    parser.add_argument(metavar='-f1',
                        dest="filenames1",
                        type=str,
                        nargs='+',
                        help='sort filenames for one condition')
    parser.add_argument('-f2',
                        dest="filenames2",
                        type=str,
                        nargs='+',
                        help='sort filenames for another condition')
    parser.add_argument('-f3',
                        dest="filenames3",
                        type=str,
                        nargs='+',
                        help='sort filenames for another condition')
    parser.add_argument('-x',
                        dest="ref_fname",
                        type=str,
                        help='ref filename')   
    parser.add_argument('--ref-length',
                        dest="ref_length",
                        default=225,
                        type=int,
                        help='length of reference template')
    parser.add_argument('--dyad-offset',
                        dest="dyad_offset",
                        default=52,
                        type=int,
                        help='off-set length from cut site to dyad position')
    parser.add_argument('--stat',
                        dest="stat_mode",
                        default=False,
                        type=bool,
                        help='stat mode')
    parser.add_argument('--graph-mode',
                        dest="graph_mode",
                        default=1,
                        type=int,
                        help='graph mode')
    parser.add_argument('-n',
                        dest="norm_choice",
                        default=False,
                        type=bool,
                        help='normalize intesnity for each sequence')
    parser.add_argument('--sample',
                        dest="sample_mode",
                        default='r',
                        type=str,
                        help='sampling mode \nRandom r:num \nUserInput k:key1,key2.. \nReadcounts counts:order/div_num/num \nGCcontents GC:order/div_num/num PolyAlen polyA:order/div_num/num \nMeanpos mpos:order/div_num/num \nMaxDis mdis:order/div_num/num \nBinEnrich bin#:order/div_num/num \nCluster cluster:obs/div_num/num'
                        )
    
    args = parser.parse_args()

    if not args.ref_fname:
        ref_fname = False
    else:
        ref_fname = args.ref_fname
    
    display_graph(args.filenames1,
                  args.filenames2,
                  args.filenames3,
                  ref_fname,
                  args.ref_length,
                  args.dyad_offset,
                  args.stat_mode,
                  args.graph_mode,
                  args.norm_choice,
                  args.sample_mode)
