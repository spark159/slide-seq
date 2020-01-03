from termcolor import colored
import random
import math
from SliderClass import Slider
import analysis
import sys

def color_A (seq):
    text = ''
    for nt in seq:
        if nt == 'A':
            text += colored(nt, 'red')
        else:
            text += nt
    return text

def ID_cmp(a, b):
    loc1, mtype1, nts1 = a.split('-')
    loc2, mtype2, nts2 = b.split('-')
    loc1, loc2 = int(loc1), int(loc2)
    if nts1 < nts2:
        return -1
    elif nts1 == nts2:
        if loc1 < loc2:
            return -1
        else:
            return 1
    else:
        return 1

def sampling (key_slider, sample_mode):
    allkey_list = key_slider.keys()
    
    def sort (key_slider, order, obs_func, *args):
        def key_cmp(a, b):
            if a[0] < b[0]:
                return -1
            elif a[0] == b[0]:
                
                alist = a[1].split('-')
                blist = b[1].split('-')
                alist = [alist[2], int(alist[0])]
                blist = [blist[2], int(blist[0])]
                if alist[0] < blist[0]:
                    return -1
                elif alist[0] == blist[0]:
                    if alist[1] <= blist[1]:
                        return -1
                    else:
                        return 1
                else:
                    return 1
                """
                if a[1] <= b[1]:
                    return -1
                else:
                    return 1
                """
            else:
                return 1
        if order == 't':
            choice = True
        elif order == 'b':
            choice = False    
        obs_key = [[obs_func(slider, *args), key] for key, slider in key_slider.items()]
        obs_key = sorted(obs_key, cmp=key_cmp, reverse=choice)
        return obs_key
    
    def quantile (obs_key, div_num, num):
        sample_list = []
        block_size = max([len(obs_key)/div_num, num])
        block_num = int(math.ceil(float(len(obs_key))/float(block_size)))
        win = (block_size-num+1)/block_num
        print block_size
        for k in range(block_num):
            key_list = []
            for obs, key in obs_key[k*block_size+k*win:k*block_size+k*win+num]:
                key_list.append(key)
            sample_list.append(key_list)
        for u in range(len(sample_list)):
            key_list = sample_list[u]
            print
            print "cond" + str(u+1) + " selected: " + str(len(key_list))
            for key in key_list[0:min([50,len(key_list)])]:
                print color_A(key)
        return sample_list

    def pick_keys (key_slider, filter_list, obs_func, *args):
        def key_cmp(a, b):
            if a[0] < b[0]:
                return -1
            elif a[0] == b[0]:
                
                alist = a[1].split('-')
                blist = b[1].split('-')
                alist = [alist[2], int(alist[0])]
                blist = [blist[2], int(blist[0])]
                if alist[0] < blist[0]:
                    return -1
                elif alist[0] == blist[0]:
                    if alist[1] <= blist[1]:
                        return -1
                    else:
                        return 1
                else:
                    return 1
                """
                if a[1] <= b[1]:
                    return -1
                else:
                    return 1
                """
            else:
                return 1

        mmrange, flist = [] , []
        for i in range(len(filter_list)):
            if filter_list[i] == '..':
                if i == 0:
                    mmrange.append([-sys.maxint, float(filter_list[i+1])])
                elif i == len(filter_list) - 1:
                    mmrange.append([float(filter_list[i-1]),sys.maxint])
                else:
                    mmrange.append([float(filter_list[i-1]),filter_list[i+1]])
            else:
                flist.append(float(filter_list[i]))
        key_list = []
        for key, slider in key_slider.items():
            obs = obs_func(slider, *args)
            if obs in flist and key not in key_list:
                    key_list.append(key)
            else:
                for min, max in mmrange:
                    if obs >= min and obs <= max and key not in key_list:
                        key_list.append(key)
        sort_list = [["",key] for key in key_list]
        sort_list = sorted(sort_list, cmp=key_cmp)
        key_list = [element[1] for element in sort_list]
        return key_list

    mode, detail = sample_mode.split(':')
    sample_list = []
    # random choice
    if mode == 'r':
        num = int(detail.strip())
        key_list = random.sample(allkey_list, num)
        sample_list.append(key_list)
        return sample_list
    # select by refID
    if mode == 'ID':
        input_list = detail.strip().split(',')
        #key_list = []
        for input in input_list:
            input = input.upper()
            key_list = []
            for key in allkey_list:
                win, size = key.split('-')
                if win == input:
                    key_list.append(key)            
            key_list = sorted(key_list, cmp=ID_cmp)
            sample_list.append(key_list)
        return sample_list
    # select by position of insertion
    if mode == 'POS':
        input_list = detail.strip().split(',')
        #key_list = []
        for input in input_list:
            key_list = []
            for key in allkey_list:
                win, pos = key.split('-')
                pos = int(pos)
                if pos == int(input):
                    key_list.append(key)            
            key_list = sorted(key_list, cmp=ID_cmp)
            sample_list.append(key_list)
        return sample_list

    # select input keys
    if mode == 'k':
        input_list = detail.strip().split(',')
        key_list = []
        for key in input_list:
            key = key.upper()
            if key in allkey_list:
                key_list.append(key)
            else:
                print key + " is missing"
        sample_list.append(key_list)
        return sample_list
    # quantile analysis
    if '/' in detail:
        order, div_num, num = detail.strip().split('/')
        div_num, num = int(div_num), int(num)
        if mode == 'counts':    
            obs_key = sort(key_slider, order, Slider.read_counts, 'dyad')
            sample_list = quantile(obs_key, div_num, num)
        elif mode == 'GC':
            obs_key = sort(key_slider, order, Slider.GC_content)
            sample_list = quantile(obs_key, div_num, num)
        elif mode == 'polyA':
            obs_key = sort(key_slider, order, Slider.Amer_len)
            sample_list = quantile(obs_key, div_num, num)
        elif mode == 'mpos':
            obs_key = sort(key_slider, order, Slider.mean_pos)
            sample_list = quantile(obs_key, div_num, num)
        elif mode == 'mdis':
            obs_key = sort(key_slider, order, Slider.max_dis)
            sample_list = quantile(obs_key, div_num, num)
        elif mode.startswith('bin'):
            bin_num = int(mode[3:].strip())
            obs_key = sort(key_slider, order, Slider.rel_bin_counts, 'dyad', bin_num)
            sample_list = quantile(obs_key, div_num, num)
        return sample_list
    # pick keys based on observables
    else:
        sample_list = []
        filter_list = detail.strip().split('-')
        for i in range(len(filter_list)):
            filter_list[i] = filter_list[i].strip().split(',')
        if mode == 'counts':
            for i in range(len(filter_list)):
                key_list = pick_keys(key_slider, filter_list[i], Slider.read_counts, 'dyad')
                sample_list.append(key_list)
        elif mode == 'GC':
            for i in range(len(filter_list)):
                key_list = pick_keys(key_slider, filter_list[i], Slider.GC_content)
                sample_list.append(key_list)
        elif mode == 'polyA':
            for i in range(len(filter_list)):
                key_list = pick_keys(key_slider, filter_list[i], Slider.Amer_len)
                sample_list.append(key_list)
        elif mode == 'mpos':
            for i in range(len(filter_list)):
                key_list = pick_keys(key_slider, filter_list[i], Slider.mean_pos)
                sample_list.append(key_list)
        elif mode == 'mdis':
            for i in range(len(filter_list)):
                key_list = pick_keys(key_slider, filter_list[i], Slider.max_dis)
                sample_list.append(key_list)
        #elif mode.startswith('bin'):
            #bin_num = int(mode[3:].strip())
            #obs_key = sort(key_slider, order, Slider.rel_bin_counts, 'dyad', bin_num)
            #sample_list = quantile(obs_key, div_num, num)
        return sample_list
    # cluster analysis
    if mode == 'cluster':
        key_obs = {}
        for key in key_slider:
            slider = key_slider[key]
            key_obs[key] = eval("slider.%s" % (order))
        key_cdx, cdx_key = analysis.Kmeans(key_obs, div_num)
        for cdx in cdx_key:
            temp = cdx_key[cdx]
            sample_list.append(random.sample(temp, num))
    return sample_list
        
    
