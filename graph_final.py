import os, sys, subprocess, re
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors
import numpy as np
import seaborn as sns
import random  
import math
import pandas as pd
import copy
import analysis_final as analysis
from sklearn.neighbors import KernelDensity
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.ticker as tck
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap

pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
                                             [('darkblue'),
                                              ('blue'),
                                              ('cyan'),
                                              ('yellow'),
                                              ('orange'),
                                              ('red'),
                                              ('darkred')],
                                               N=256)



#pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
#                                             [(0, '#ffffff'),
#                                              (0.6, 'tab:blue'),
#                                              (0.7, 'tab:green'),
#                                              (0.8, 'yellow'),
#                                              (0.9, 'tab:orange'),
#                                              (1, 'darkred')
#                                             ], N=256)

#pastel_jet = LinearSegmentedColormap.from_list('white_viridis',
#                                             [(0, '#ffffff'),
#                                              (0.6, 'tab:cyan'),
#                                              (0.7, 'tab:blue'),
#                                              (0.8, 'yellow'),
#                                              (0.9, 'tab:orange'),
#                                              (1, 'darkred')
#                                             ], N=256)


#pastel_jet = LinearSegmentedColormap.from_list('pastel_jet',
#                                              [(0, 'tab:blue'),
#                                               (0.2, 'tab:cyan'),
#                                               (0.4, 'tab:green'),
#                                               (0.6, 'yellow'),
#                                               (0.8, 'tab:orange'),
#                                               (1, 'tab:red'),],
#                                               N=256)


#pastel_jet_r = ListedColormap(pastel_jet.colors[::-1])
#pastel_jet_r = pastel_jet[::-1]

# plot (stacked) heatmap for 1-d signal observables of sliders
def plot_map(id_slider,       # dictionary of id:Slider
             obs_func,        # set observables to plot
             ids=[],          # id selection
             xticks=None,     # [xticks, xticklabels]
             norm=True,       # normalization option
             thickness=None,  # [strip, space, marker]
             mark=None,       # mark id or poly-N in the heatmap
             cmap='jet',      # set colormap 
             mcolor='red',    # set marker color
             slicing=0,       # slice out the both ends of data
             note="",         # add note to the name
             save=False,      # save heatmap
             figscale=100,    # set the figure size len(data)/figscale (inches)
             figsize=None,    # set figure size [height, width] (inches)
             fontsize=6,      # set xick label font size
             *args):

    if len(ids) <= 0:
        ids = sorted(id_slider.keys(), cmp=analysis.wid_cmp_len)

    if type(obs_func) != list:
        obs_funcs = [obs_func]
    else:
        obs_funcs = obs_func

    if type(cmap) != list:
        cmaps = [cmap]*len(obs_funcs)
    else:
        cmaps = cmap

    if thickness == None:
        if len(ids) <= 50:
            strip, space, marker = 10, 5, 2
        else:
            if mark:
                strip, space, marker = 2, 0, 1
            else:
                strip, space, marker = 1, 0, 0
    else:
        strip, space, marker = thickness

    substrip = int(math.ceil(float(strip)/len(obs_funcs)))

    # make stack of images
    img_stack = [[] for i in range(len(obs_funcs))]

    for i in range(len(img_stack)):
        img = img_stack[i]
        
        for j in range(len(ids)):
            id = ids[j]
            slider = id_slider[id]
            obs_sig = obs_funcs[i](slider, *args)
            obs_sig = obs_sig[slicing:len(obs_sig)-slicing]

            if norm:
                obs_sig = analysis.normalize_list(obs_sig)

            for k in range(len(img_stack)*substrip):
                if k >= i*substrip and k < (i+1)*substrip:
                    img.append(obs_sig)
                else:
                    img.append([np.nan]*len(obs_sig))

            if j < len(ids)-1:
                for k in range(space):
                    img.append([np.nan]*len(obs_sig))

    # make marker image
    if mark:
        mark_img = []
        
        for j in range(len(ids)):
            id = ids[j]
            slider = id_slider[id]
            obs_sig = obs_funcs[0](slider, *args)
            mark_sig = np.zeros(len(obs_sig))
            mark_sig[:] = np.nan

            if mark == 'wid':
                loc, mtype, nts = id.split('-')
                size, loc = len(nts), int(loc)
                st, ed = loc, loc+size
                mark_sig[st:ed] = 1

            elif mark.startswith('poly'):
                try:
                    nts = mark.split('-')[1].split(':')[0]
                    nts += analysis.rev_comp(nts)
                except:
                    nts = 'AT'
                try:
                    _, minlen = mark.split(':')
                    minlen = int(minlen)
                except:
                    minlen = 5

                len_pos = slider.polynt_count(nts='AT', pos=True)
                for length, pos in len_pos.items():
                    if length >= minlen:
                        for loc in pos:
                            mark_sig[loc:loc+length] = 1

            mark_sig = mark_sig[slicing:len(mark_sig)-slicing]

            for k in range(len(img_stack)*substrip):
                if k < marker:
                    mark_img.append(mark_sig)
                else:
                    mark_img.append([np.nan]*len(mark_sig))

            if j < len(ids)-1:
                for k in range(space):
                    mark_img.append([np.nan]*len(mark_sig))

    fig = plt.figure()
    if figsize == None:
        ncols, nrows = np.shape(img_stack[0])
        height, width = float(ncols)/figscale, float(nrows)/figscale
        aspect='equal'
    else:
        height, width = figsize
        aspect='auto'
    fig.set_size_inches(width, height, forward = False)
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0, bottom=0.3/height, right=1, top=1, wspace=0, hspace=0)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tick_params(top='off', left='off', right='off', labelleft='off', labelbottom='on')

    for img, cmap in zip(img_stack, cmaps):
        colormap = cm.get_cmap(cmap)
        colormap.set_bad(alpha=0)
        ax.imshow(img, cmap=colormap, interpolation='none', aspect=aspect)

    if mark:
        colormap = colors.ListedColormap([mcolor])
        colormap.set_bad(alpha=0)
        ax.imshow(mark_img, cmap=colormap, interpolation='none', aspect=aspect)

    if xticks:
        ax.set_xticks(xticks[0])
        ax.set_xticklabels(xticks[1], fontsize=fontsize)
        ax.xaxis.set_minor_locator(tck.AutoMinorLocator(2))

    #plt.tight_layout()

    if save:
        #plt.savefig('Heatmap' + note + '.png', dpi=1000)
        plt.savefig('Heatmap' + note + '.svg', format='svg')
    else:
        plt.show()

    plt.close()

    if len(img_stack) <= 1:
        return img_stack[0]
    return img_stack

# plot (stacked) 1-d signal observables of sliders
def plot_sig(id_slider,
             obs_func,
             ids = [],
             row=8,
             col=4,
             norm=False,
             mark=None,
             color=None,
             label=None,
             alpha=None,
             xticks=None,
             yticks=None,
             xlims=[None, None],
             ylims=[None, None],
             frame_color='k',
             note="",
             save=False,
             *args):

    if len(ids) <= 0:
        ids = sorted(id_slider.keys(), cmp=analysis.wid_cmp_len)

    if type(obs_func) != list:
        obs_funcs = [obs_func]
    else:
        obs_funcs = obs_func

    if type(color) != list:
        colors = [color]*len(obs_funcs)
    else:
        colors = color

    if type(alpha) != list:
        alphas = [alpha]*len(obs_funcs)
    else:
        alphas = alpha

    if type(label) != list:
        labels = [label]*len(obs_funcs)
    else:
        labels = label

    page_nums = int(math.ceil(len(ids)/float(row*col))) # number of pages

    if save == True:
        if page_nums > 1:
            save = 'pdf'
        else:
            save = 'png'
    
    if save == 'pdf':
        pdf = matplotlib.backends.backend_pdf.PdfPages('Signals' + note + ".pdf")

    #page_nums = 1
    for i in range(page_nums):
        fig = plt.figure(figsize=(8.5,11))
        j = 0
        while j < min(row*col, len(ids)-row*col*i):
            id = ids[row*col*i + j]
            slider = id_slider[id]
            
            plt.subplot(row, col, j+1)
            ax = plt.gca()
            ax.spines['bottom'].set_color(frame_color)
            ax.spines['top'].set_color(frame_color)
            ax.spines['left'].set_color(frame_color)
            ax.spines['right'].set_color(frame_color)

            for obs_func, color, alpha, label in zip(obs_funcs, colors, alphas, labels):
                sig = obs_func(slider, *args)
                if norm:
                    sig = analysis.normalize_list(sig)
                plt.plot(sig, color=color, alpha=alpha, label=label, lw=1)

            if mark == 'wid':
                loc, mtype, nts = id.split('-')
                size, loc = len(nts), int(loc)
                st, ed = loc, loc+size
                if size > 1:
                    plt.axvspan(st, ed-1, alpha=0.6, facecolor='red', lw=1)
                else:
                    plt.axvline(x=st, color='red', alpha=0.6, lw=1)

            elif mark != None and mark.startswith('poly'):
                try:
                    nts = mark.split('-')[1].split(':')[0]
                    nts += analysis.rev_comp(nts)
                except:
                    nts = 'AT'
                try:
                    _, minlen = mark.split(':')
                    minlen = int(minlen)
                except:
                    minlen = 5
                    
                len_pos = slider.polynt_count(nts='AT', pos=True)
                for length, pos_list in len_pos.items():
                    if length < minlen:
                        continue
                    for pos in pos_list:
                        st, ed = pos, pos+length
                        plt.axvspan(st, ed-1, alpha=0.5, color='red')
            
            if mark == 'wid':
                loc, mtype, nts = id.split('-')
                size, loc = len(nts), int(loc)
                st, ed = loc-225/2, loc-225/2+size
                #mtype = 'A' #temporal
                if mtype == 'I':
                    new_id = '%s$_{%d}$[%d^%d]' % (mtype, size, st, ed)
                else:
                    if size <= 1:
                        new_id = '%s$_{%d}$[%d]' % (mtype, size, st)
                    else:
                        new_id = '%s$_{%d}$[%d:%d]' % (mtype, size, st, ed-1) 
                plt.title(new_id, fontsize=8)
                #plt.title(id + '(%dbp)' % (size))
            else:
                plt.title(id, fontsize=8)

            if xticks:
                plt.xticks(xticks[0], xticks[1])
                ax.xaxis.set_minor_locator(tck.AutoMinorLocator(2))
            if yticks:
                plt.yticks(yticks[0], yticks[1])
            plt.tick_params(axis='x', labelsize=5)
            plt.tick_params(axis='y', labelsize=5)

            leg = plt.legend(loc='upper right', frameon=False, fontsize=5)
            for lh in leg.legendHandles: 
                #lh.set_alpha(1)
                lh.set_linewidth(1.0)
                
            plt.xlim(xlims)
            plt.ylim(ylims)
            j +=1

        plt.subplots_adjust(wspace=0.25, hspace=0.8)
        #fig.tight_layout()
        if save == 'pdf':
            fig.text(4.25/8.5, 0.5/11., str(i+1), ha='center', fontsize=8)
            pdf.savefig(fig)
        elif save == 'png':
            plt.savefig('Signals' + note + '_' + str(i+1) + '.png', dpi=1000)
        else:
            plt.show()
        plt.close()

    if save == 'pdf':
        pdf.close()

    


"""
# plot heatmap for 1-d signal observables of sliders
def plot_map(id_slider,
             obs_func,
             ids=[],
             norm=True,
             thickness=None,
             mark=None,
             cmap='jet', 
             slicing=0,
             note="",
             save=False,
             *args):    

    if len(ids) <= 0:
        ids = analysis.sorted(id_slider.keys(), cmp=analysis.wid_cmp_len)

    if thickness == None:
        if len(ids) <= 50:
            strip, space, marker = 10, 5, 2
        else:
            if mark:
                strip, space, marker = 2, 0, 1
            else:
                strip, space, marker = 1, 0, 0
    else:
        strip, space, marker = thickness

        
    img = []

    for i in range(len(ids)):
        id = ids[i]
        slider = id_slider[id]
        tmp_img = []
        obs_map = obs_func(slider, *args)

        if mark == 'wid':
            loc, mtype, nts = id.split('-')
            size, loc = len(nts), int(loc)
            st, ed = loc, loc+size
            mark_map = np.zeros(len(obs_map))
            mark_map[st:ed] = [np.nan] * size

        if mark.startswith('poly'):
            try:
                nts = mark.split('-')[1].split(':')[0]
                nts += analysis.rev_comp(nts)
            except:
                nts = 'AT'
            try:
                _, minlen = mark.split(':')
                minlen = int(minlen)
            except:
                minlen = 5

            len_pos = slider.polynt_count(nts='AT', pos=True)
            mark_map = np.zeros(len(obs_map))
            for length, pos in len_pos.items():
                if length >= minlen:
                    for loc in pos:
                        mark_map[loc:loc+num] = [np.nan] * length

        obs_map = obs_map[slicing:len(obs_map)-slicing]
        if mark:
            mark_map = mark_map[slicing:len(mark_map)-slicing]

        for k in range(strip):
            tmp_img.append(obs_map)

        if norm:
            tmp_img = normalize_matrix(tmp_img)

        for k in range(len(tmp_img)):
            row = tmp_img[k]
            if mark and k < marker:
                img.append(np.asarray(row) + mark_map)
            else:
                img.append(row)

        if i < len(id_list)-1:
            for k in range(space):
                img.append([0]*len(obs_map))            

    fig = plt.figure()
    height, width = np.shape(img)
    fig.set_size_inches(float(width)/100, float(height)/100, forward = False)
    ax = plt.subplot(111, aspect = 'equal')
    plt.subplots_adjust(left=0, bottom=0.3/(float(height)/100), right=1, top=1, wspace=0, hspace=0)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tick_params(top='off', left='off', right='off', labelleft='off', labelbottom='on')
    #ax.set_xticks(range(0, width, 5))
    #ax.set_xticklabels([str(i+1) for range(0, width, 5)])

    colormap = cm.get_cmap(cmap)
    if mark:
        colormap.set_bad((1, 0, 0, 1)) # red marker

    ax.imshow(img, cmap=colormap, interpolation='none')

    if save:
        plt.savefig('Heatmap' + note + '.png', dpi=1000)
    else:
        plt.show()

    plt.close()
    return img


# plot cut/dyad map
def plot_map(key_slider, sample_list, norm_choice, note = "", draw_key=False, draw_vert=False):    
    for i in range(len(sample_list)):
        key_list = sample_list[i]
        cut_Rimg, cut_Limg = [], []
        dyad_img = []

        if len(key_list) <= 50:
        #if len(key_list) <= 1:
            A,B,C,D = 5, 10, 5, 2
        else:
            A,B,C,D = 1, 1, 0, 0
            if draw_key:
                A,B,C,D = 4, 8, 0, 2
                #A,B,C,D = 1, 2, 0, 1
                
        for j in range(len(key_list)):
            key = key_list[j]
            slider = key_slider[key]
            cut_Rtmp_img, cut_Ltmp_img = [], []
            dyad_tmp_img = []

            if draw_key:
                #win, loc = key.split('-')
                #size, loc = len(win), int(loc)
                try:
                    loc, mtype, nts = key.split('-')
                    size, loc = len(nts), int(loc)
                    st, ed = loc, loc+size
                except:
                    st, ed = 0, 1

            for side in ['R','L']:
                for k in range(A):
                    if side == 'R':
                        cut_map = slider.right_cutmap
                        cut_Rtmp_img.append(cut_map)
                        cut_Ltmp_img.append([0.0]*len(cut_map))
                    else:
                        assert side == 'L'
                        cut_map = slider.left_cutmap
                        cut_Ltmp_img.append(cut_map)
                        cut_Rtmp_img.append([0.0]*len(cut_map))

            dyad_map = slider.dyadmap    
            for k in range(B):
                dyad_tmp_img.append(dyad_map)
            
            if norm_choice:
                cut_Rtmp_img = normalize(cut_Rtmp_img)
                cut_Ltmp_img = normalize(cut_Ltmp_img)
                dyad_tmp_img = normalize(dyad_tmp_img)
            
            for k in range(len(cut_Rtmp_img)):
                row = cut_Rtmp_img[k]
                if draw_key and np.mean([st,ed]) < len(cut_Rtmp_img[0])/2 and k < D:
                    cut_Rimg.append(row[:st] + [np.nan]*(ed-st) + row[ed:])
                else:
                    cut_Rimg.append(row)

            for k in range(len(cut_Ltmp_img)):
                row = cut_Ltmp_img[k]
                if draw_key and np.mean([st,ed]) >= len(cut_Ltmp_img[0])/2 and k >= len(cut_Ltmp_img) - D:
                    cut_Limg.append(row[:st] + [np.nan]*(ed-st) + row[ed:])
                else:
                    cut_Limg.append(row)

            for k in range(len(dyad_tmp_img)):
                row = dyad_tmp_img[k]
                if draw_key and k < D:
                    dyad_img.append(row[:st] + [np.nan]*(ed-st) + row[ed:])
                else:
                    dyad_img.append(row)

            if j < len(key_list)-1:
                for k in range(C):
                    cut_Rimg.append([0]*len(cut_map))
                    cut_Limg.append([0]*len(cut_map))
                    dyad_img.append([0]*len(dyad_map))

        assert np.shape(cut_Limg) == np.shape(cut_Rimg)
        #cut_img = np.matrix(cut_Limg) + np.matrix(cut_Rimg)
            
        size = np.shape(cut_Limg)
        height = float(size[0])
        width = float(size[1])
        fig = plt.figure()
        fig.set_size_inches(width/100, height/100, forward = False)
        #ax = plt.Axes(fig, [10,10,0.5,0.5])
        #ax.axis('off')
        #fig.add_axes(ax)
        ax = plt.subplot(111, aspect = 'equal')
        plt.subplots_adjust(left=0, bottom=0.3/(height/100), right=1, top=1, wspace=0, hspace=0)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tick_params(top='off', left='off', right='off', labelleft='off', labelbottom='on')
        cut_Rimg = np.asarray(cut_Rimg)
        if draw_key:
            R_coordinate = np.argwhere(np.isnan(cut_Rimg))
        norm = plt.Normalize(np.nanmin(cut_Rimg), np.nanmax(cut_Rimg))
        with np.errstate(invalid='ignore'):
            cut_Rimg = np.ma.masked_where(cut_Rimg <= 0.0, cut_Rimg)
            cmap = plt.cm.Blues
            cut_Rimg = cmap(norm(cut_Rimg))
        if draw_key:
            for x, y in R_coordinate:
                cut_Rimg[x, y, :3] = 1, 0, 0
        ax.imshow(cut_Rimg, interpolation='none')
        cut_Limg = np.asarray(cut_Limg)
        if draw_key:
            L_coordinate = np.argwhere(np.isnan(cut_Limg))
        norm = plt.Normalize(np.nanmin(cut_Limg), np.nanmax(cut_Limg))
        with np.errstate(invalid='ignore'):
            cut_Limg = np.ma.masked_where(cut_Limg <= 0.0, cut_Limg)
            cmap = plt.cm.Reds
            cut_Limg = cmap(norm(cut_Limg))
        if draw_key:
            for x, y in L_coordinate:
                cut_Limg[x, y, :3] = 1, 0, 0
        ax.imshow(cut_Limg, interpolation='none')
        if draw_vert:
            center = len(cut_Limg[0])/2
            for x in [center-30, center-20, center, center+20, center+30]:
                ax.axvline(x, color='k', linestyle='--', linewidth=0.5)
        #ax.imshow(cut_img, cmap=cmap, interpolation='none')
        #plt.show()
        #plt.savefig('cut_cond' + str(i+1) + win + note + '.png', dpi=1500)
        plt.savefig('cut_cond' + str(i+1) + note + '.png', dpi=1000)
        #plt.savefig('cut_cond' + str(i+1) + '.svg', format = 'svg')
        #plt.show()
        plt.close()
        

        fig=plt.figure()
        size = np.shape(dyad_img)
        height = float(size[0])
        width = float(size[1])
        #fig.set_size_inches(width/height, 1, forward = False)
        fig.set_size_inches(width/100, height/100, forward = False)
        ax = plt.subplot(111, aspect = 'equal')
        plt.subplots_adjust(left=0, bottom=0.3/(height/100), right=1, top=1, wspace=0, hspace=0)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tick_params(top='off', left='off', right='off', labelleft='off', labelbottom='on')
        #ax=plt.Axes(fig, [0.,0.,1.,1.])
        #ax.axis('off')
        #fig.add_axes(ax)
        #cmap = plt.cm.YlGnBu
        cmap = plt.cm.Greens
        #cmap = plt.cm.Purples
        #cmap = plt.cm.binary_r
        #cmap = plt.cm.jet
        #cmap = plt.cm.binary
        if draw_key:
            cmap.set_bad((1, 0, 0, 1))
            #cmap.set_bad('y', 1)
        ax.imshow(dyad_img, cmap=cmap, interpolation='none')
        if draw_vert:
            center = len(dyad_img[0])/2
            for x in [center-30, center-20, center, center+20, center+30]:
                ax.axvline(x, color='k', linestyle='--', linewidth=0.5)
        plt.savefig('dyad_cond' + str(i+1) + note + '.png', dpi=1200)
        #plt.savefig('dyad_cond' + str(i+1) + win + note + '.png', dpi=1500)
        #plt.savefig('dyad_cond' + str(i+1) + '.svg', format = 'svg')
        #plt.show()
        plt.close()


        
# plot observable signal
def plot_signal (key_slider, sample_list, norm_choice, obs_func, mean_choice=False, draw = None, slicing = 0, note = "", *args ):
    if sample_list == None:
        sample_list = [key_slider.keys()]
    for i in range(len(sample_list)):
        key_list = sample_list[i]
        obs_list = []
        for j in range(len(key_list)):
            key = key_list[j]
            slider = key_slider[key]
            obs = obs_func(slider, *args)
            if norm_choice:
                obs = normalize([obs], percentage=False)[0]
            obs = np.asarray(obs)
            obs_list.append(obs)
        if mean_choice:
            obs_list = [np.sum(obs_list, axis=0)/len(key_list)]
        global dyad_axis        
        fig = plt.figure()
        for k in range(len(obs_list)):
            if mean_choice:
                label = 'mean:'+str(key_list[0])+'-'+str(key_list[-1])
            else:
                label = key_list[k]
            plt.plot(obs, label=label)
        if draw == 'key':
            loc, mtype, nts = key.split('-')
            size, loc = len(nts), int(loc)
            st, ed = loc, loc+size
            plt.axvspan(st, ed-1, alpha=0.5, color='red')
        if draw == 'polyA':
            num_pos = slider.Amer_len_detail()
            for num, pos in num_pos.items():
                if num >= 3:
                    for loc in pos:
                        st, ed = loc, loc+num
                        plt.axvspan(st, ed-1, alpha=0.5, color='red')
        plt.axvline(dyad_axis, color='k', linestyle='--')
        plt.xlim([0,225])
        plt.xlabel('Location (bp)')
        plt.ylabel('Signal')
        plt.legend(loc='best')
        #plt.show()
        plt.savefig('sig_fig' + str(i+1) + note + '.png', bbox_inches = 'tight')
        plt.close()

# plot SeqID vs cut peaks
def plot_cpeaks (key_slider, left_peak_num = 1, right_peak_num = 1, sample_list=None, note=""):
    if sample_list == None:
        sample_list = [key_slider.keys()]

    for i in range(len(sample_list)):
        key_list = sample_list[i]
        key_lpeaks, key_rpeaks = [], []

        for j in range(len(key_list)):
            key = key_list[j]
            slider = key_slider[key]
            left_cutmap = slider.left_cutmap
            right_cutmap = slider.right_cutmap
            if min(sum(left_cutmap), sum(right_cutmap)) < 10:
                continue
            lpeaks = find_peaks(left_cutmap, num = left_peak_num)
            rpeaks = find_peaks(right_cutmap, num = right_peak_num)
            global dyad_axis
            lpeaks = [peak - dyad_axis for peak,sig in lpeaks]
            rpeaks = [peak - dyad_axis for peak,sig in rpeaks]
            if len(lpeaks) < left_peak_num or len(rpeaks) < right_peak_num:
                continue
            key_lpeaks.append([key,lpeaks])
            key_rpeaks.append([key,rpeaks])
        
        key_lpeaks = sorted(key_lpeaks, cmp=key_cmp)
        key_rpeaks = sorted(key_rpeaks, cmp=key_cmp)        
        ly, ry = [[] for k in range(left_peak_num)], [[] for k in range(right_peak_num)]
        assert len(key_lpeaks) == len(key_rpeaks)
        for u in range(len(key_lpeaks)):
            lkey, lpeaks = key_lpeaks[u]
            rkey, rpeaks = key_rpeaks[u]
            assert lkey == rkey
            for k in range(left_peak_num):
                ly[k].append(lpeaks[k])
            for k in range(right_peak_num):
                ry[k].append(rpeaks[k])
        fig = plt.figure()
        for k in range(left_peak_num):
            plt.plot(range(len(key_lpeaks)),ly[k],'.')
        for k in range(right_peak_num):
            plt.plot(range(len(key_rpeaks)),ry[k],'.')
        plt.xlabel('Seq ID')
        plt.ylabel('Cut Peak location')
        plt.savefig('SeqID_cpeak_location_' + 'cond' + str(i+1) + '.png')
        plt.close()
        
        
        Amer_lpeaks = [[Amer_len(key), lpeak] for key, lpeak in key_lpeaks]
        Amer_rpeaks = [[Amer_len(key), rpeak] for key, rpeak in key_rpeaks]
        Amer_lpeaks = sorted(Amer_lpeaks, cmp=key_cmp)
        Amer_rpeaks = sorted(Amer_rpeaks, cmp=key_cmp)
        ly, ry = [[] for k in range(left_peak_num)], [[] for k in range(right_peak_num)]
        assert len(Amer_lpeaks) == len(Amer_rpeaks)
        for u in range(len(Amer_lpeaks)):
            lAmer, lpeaks = key_lpeaks[u]
            rAmer, rpeaks = key_rpeaks[u]
            assert lAmer == rAmer
            for k in range(left_peak_num):
                ly[k].append(lpeaks[k])
            for k in range(right_peak_num):
                ry[k].append(rpeaks[k])
        fig = plt.figure()
        for k in range(left_peak_num):
            plt.plot(range(len(Amer_lpeaks)),ly[k],'.')
        for k in range(right_peak_num):
            plt.plot(range(len(Amer_rpeaks)),ry[k],'.')
        plt.xlabel('Seq ID')
        plt.ylabel('Cut Peak location')
        plt.savefig('SeqID_cpeak_location_' + 'cond' + str(i+1) + note + '.png')
        plt.close()
        

# plot SeqID vs dyad peaks
def plot_dpeaks (key_slider, peak_num = 1, st_rank = 1, sample_list=None, note=""):
    if sample_list == None:
        sample_list=[key_slider.keys()]
    for i in range(len(sample_list)):
        key_list = sample_list[i]
        key_peaks = []
        for j in range(len(key_list)):
            key = key_list[j]
            slider = key_slider[key]
            dyadmap = slider.dyadmap
            if sum(dyadmap) < 10:
                continue
            dpeaks = find_peaks(dyadmap, num = peak_num + st_rank -1)
            global dyad_axis
            dpeaks = [peak - dyad_axis for peak,sig in dpeaks]
            if len(dpeaks) < st_rank:
                continue
            dpeaks = dpeaks[st_rank-1:]
            dpeaks = dpeaks + [np.nan]*(len(dpeaks) - peak_num)
            key_peaks.append([key,dpeaks])
        
        key_peaks = sorted(key_peaks, cmp=key_cmp)
        y = [[] for k in range(peak_num)]
        for key, peaks in key_peaks:
            for k in range(peak_num):
                y[k].append(peaks[k])
        fig = plt.figure()
        for k in range(peak_num):
            plt.plot(range(len(key_peaks)),y[k],'.')
        plt.xlabel('Seq ID')
        plt.ylabel('Dyad Peak location')
        #plt.show()
        plt.savefig('SeqID_dpeak_location_' + 'cond' + str(i+1) + '.png')
        plt.close()

        
        Amer_peaks = [[Amer_len(key), peaks] for key, peaks in key_peaks]
        Amer_peaks = sorted(Amer_peaks, cmp=key_cmp)
        y = [[] for k in range(peak_num)]
        for Amer, peaks in Amer_peaks:
            for k in range(peak_num):
                y[k].append(peaks[k])
        fig = plt.figure()
        for k in range(peak_num):
            plt.plot(range(len(Amer_peaks)),y[k],'.')
        plt.xlabel('Seq ID')
        plt.ylabel('Dyad Peak location')
        plt.savefig('SeqID_dpeak_location_' + 'cond' + str(i+1) + note + '.png')
        plt.close()
        
        
# plot correlation between seq observable and sig observable        
def plot_corr (key_slider, obs_func1, obs_func2, sample_list=None, xlabel="", ylabel=""):    
    def get_corr(x, y):
        assert len(x) == len(y)
        n = len(x)
        assert n > 0
        avg_x = np.average(x)
        avg_y = np.average(y)
        diffprod = 0
        xdiff2 = 0
        ydiff2 = 0
        for idx in range(n):
            xdiff = x[idx] - avg_x
            ydiff = y[idx] - avg_y
            diffprod += xdiff * ydiff
            xdiff2 += xdiff * xdiff
            ydiff2 += ydiff * ydiff
        return diffprod / np.sqrt(xdiff2 * ydiff2)
    
    if sample_list == None:
        sample_list = [key_slider.keys()]
    for i in range(len(sample_list)):
        key_list = sample_list[i]
        X, Y = [], []
        data_frame = pd.DataFrame()
        obs1_obs2 = {}
        for key in key_list:
            slider = key_slider[key]
            obs1 = obs_func1(slider)
            obs2 = obs_func2(slider)
            X.append(obs1); Y.append(obs2)
            if obs1 not in obs1_obs2:
                obs1_obs2[obs1] = []
            obs1_obs2[obs1].append(obs2)
        data_frame[xlabel]=X
        data_frame[ylabel]=Y

        print "%s VS %s correlation: %f" % (xlabel, ylabel, get_corr(X,Y))
        fig = plt.figure()
        plt.plot(X,Y, '.b', markersize=7, alpha=0.3)
        sns.kdeplot(X,Y, shade=False, shade_lowest=False)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.savefig(xlabel + 'VS' + ylabel + '_scatter_' + 'cond' + str(i+1) + '.png')
        plt.close()

        fig = plt.figure()
        sns.lmplot(x=xlabel, y=ylabel, data=data_frame)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.savefig(xlabel + 'VS' + ylabel + '_scatterfit_' + 'cond' + str(i+1) + '.png')
        plt.close()                             

        x, y, z = [], [], []
        for obs1 in obs1_obs2:
            x.append(obs1)
            y.append(np.mean(obs1_obs2[obs1]))
            #z.append(np.std(obs1_obs2[obs1])/np.sqrt(len(obs1_obs2[obs1])))
            z.append(np.std(obs1_obs2[obs1]))
        fig = plt.figure()
        plt.plot(x,y,'.')
        plt.errorbar(x,y,yerr=z,fmt='o')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.savefig(xlabel + 'VS' + ylabel + 'cond' + str(i+1) + '.png')
        plt.close()

# plot correlation between seq observable and sig observable (combined all)
def plot_corr2 (key_slider, obs_func1, obs_func2, sample_list=None, sample_labels=None, xlabel="", ylabel=""):    
    def get_corr(x, y):
        assert len(x) == len(y)
        n = len(x)
        assert n > 0
        avg_x = np.average(x)
        avg_y = np.average(y)
        diffprod = 0
        xdiff2 = 0
        ydiff2 = 0
        for idx in range(n):
            xdiff = x[idx] - avg_x
            ydiff = y[idx] - avg_y
            diffprod += xdiff * ydiff
            xdiff2 += xdiff * xdiff
            ydiff2 += ydiff * ydiff
        return diffprod / np.sqrt(xdiff2 * ydiff2)
    
    if sample_list == None:
        sample_list = [key_slider.keys()]
    if sample_labels == None:
        sample_labels = [None] * len(sample_list)

    color_list = ['r','b']
        
    for i in range(len(sample_list)):
        key_list = sample_list[i]
        X, Y = [], []
        data_frame = pd.DataFrame()
        obs1_obs2 = {}
        for key in key_list:
            slider = key_slider[key]
            obs1 = obs_func1(slider)
            obs2 = obs_func2(slider)
            X.append(obs1); Y.append(obs2)
            if obs1 not in obs1_obs2:
                obs1_obs2[obs1] = []
            obs1_obs2[obs1].append(obs2)
        data_frame[xlabel]=X
        data_frame[ylabel]=Y

        print "%s VS %s correlation: %f" % (xlabel, ylabel, get_corr(X,Y))
        plt.figure(1)
        plt.plot(X,Y, '.', color=color_list[i], markersize=7, alpha=0.3, label=sample_labels[i])
        sns.kdeplot(X,Y, shade=False, shade_lowest=False)

        #fig = plt.figure()
        #sns.lmplot(x=xlabel, y=ylabel, data=data_frame)
        #plt.xlabel(xlabel)
        #plt.ylabel(ylabel)
        #plt.savefig(xlabel + 'VS' + ylabel + '_scatterfit_' + 'cond' + str(i+1) + '.png')
        #plt.close()                             

        x, y, z = [], [], []
        for obs1 in obs1_obs2:
            x.append(obs1)
            y.append(np.mean(obs1_obs2[obs1]))
            #z.append(np.std(obs1_obs2[obs1])/np.sqrt(len(obs1_obs2[obs1])))
            z.append(np.std(obs1_obs2[obs1]))
        plt.figure(2)
        plt.plot(x,y,'.-', color=color_list[i], label=sample_labels[i])
        plt.errorbar(x,y,yerr=z,fmt='o', color=color_list[i])
        
    plt.figure(1)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.tight_layout()
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim([-5,30])
    plt.ylim([-20,20])
    leg = plt.legend(loc='best', numpoints=1, prop={'size': 15})
    for lh in leg.legendHandles:
        lh._legmarker.set_markersize(15)
        lh._legmarker.set_alpha(1)
    plt.savefig(xlabel + 'VS' + ylabel + '_scatter_' + '.png')
    plt.close()
    
    plt.figure(2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(xlabel + 'VS' + ylabel + '.png')
    plt.close()


# plot correlation between seq observable and sig observable (combined all)
def plot_corr3 (key_slider1, key_slider2, obs_func1, obs_func2, sample_list=None, sample_labels=None, xlabel="", ylabel=""):    
    def get_corr(x, y):
        assert len(x) == len(y)
        n = len(x)
        assert n > 0
        avg_x = np.average(x)
        avg_y = np.average(y)
        diffprod = 0
        xdiff2 = 0
        ydiff2 = 0
        for idx in range(n):
            xdiff = x[idx] - avg_x
            ydiff = y[idx] - avg_y
            diffprod += xdiff * ydiff
            xdiff2 += xdiff * xdiff
            ydiff2 += ydiff * ydiff
        return diffprod / np.sqrt(xdiff2 * ydiff2)
    
    if sample_list == None:
        sample_list = [list(set(key_slider1.keys()) & set(key_slider2.keys()))]
    if sample_labels == None:
        sample_labels = [None] * len(sample_list)

    color_list = ['b','r']
        
    for i in range(len(sample_list)):
        key_list = sample_list[i]
        X, Y = [], []
        data_frame = pd.DataFrame()
        obs1_obs2 = {}
        for key in key_list:
            slider1, slider2 = key_slider1[key], key_slider2[key]
            obs1 = obs_func1(slider1)
            obs2 = obs_func2(slider2) - obs_func2(slider1)
            X.append(obs1); Y.append(obs2)
            if obs1 not in obs1_obs2:
                obs1_obs2[obs1] = []
            obs1_obs2[obs1].append(obs2)
        data_frame[xlabel]=X
        data_frame[ylabel]=Y

        print "%s VS %s correlation: %f" % (xlabel, ylabel, get_corr(X,Y))
        plt.figure(1)
        plt.plot(X,Y, '.', color=color_list[i], markersize=7, alpha=0.3, label=sample_labels[i])
        sns.kdeplot(X,Y, shade=False, shade_lowest=False)

        #fig = plt.figure()
        #sns.lmplot(x=xlabel, y=ylabel, data=data_frame)
        #plt.xlabel(xlabel)
        #plt.ylabel(ylabel)
        #plt.savefig(xlabel + 'VS' + ylabel + '_scatterfit_' + 'cond' + str(i+1) + '.png')
        #plt.close()                             

        x, y, z = [], [], []
        for obs1 in obs1_obs2:
            x.append(obs1)
            y.append(np.mean(obs1_obs2[obs1]))
            #z.append(np.std(obs1_obs2[obs1])/np.sqrt(len(obs1_obs2[obs1])))
            z.append(np.std(obs1_obs2[obs1]))
        plt.figure(2)
        plt.plot(x,y,'.-', color=color_list[i], label=sample_labels[i])
        plt.errorbar(x,y,yerr=z,fmt='o', color=color_list[i])
        
    plt.figure(1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(xlabel + 'VS' + ylabel + '_scatter_' + '.png')
    plt.close()
    
    plt.figure(2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(xlabel + 'VS' + ylabel + '.png')
    plt.close()

def plot_energy (key_slider1, key_slider2,  sample_list=None, note=""):
    
    if sample_list == None:
        sample_list = [['8-77']]
    for i in range(len(sample_list)):
        key_list = sample_list[i]
        
        for j in range(len(key_list)):
            key = key_list[j]
            #size, loc = key.split('-')
            win, loc = key.split('-')
            size, loc = len(win), int(loc)
            st, ed = loc, loc+size
            
            slider1, slider2 = key_slider1[key], key_slider2[key]
            dyadmap1, dyadmap2 = slider1.dyadmap, slider2.dyadmap
            KDE1 = slider1.KDE()
            KDE2 = slider2.KDE()
            energy1 = slider1.energy_profile()
            force1 = slider1.force_profile()
            energy2 = slider2.energy_profile()
            force2 = slider2.force_profile()
            #denergy = energy_profile2 - energy_profile1
            dforce = force2 - force1
            
            fig = plt.figure()
            #plt.plot(range(len(dyadmap1)), dyadmap1, label="Before")
            #plt.plot(range(len(dyadmap2)), dyadmap2, label="After")
            #plt.plot(range(len(KDE1)), KDE1, label = 'Before')
            #plt.plot(range(len(KDE2)), KDE2, label = 'After')
            #plt.plot(range(len(energy1)), energy1, label='Before')
            #plt.plot(range(len(energy2)), energy2, label='After')
            #plt.plot(range(len(force1)), force1, label='Before')
            #plt.plot(range(len(force2)), force2, label='After')
            plt.plot(range(len(dforce)), dforce, label='After-Before')
            plt.axvline(dyad_axis, color='k', linestyle='--')
            plt.axhline(0, color='k', linestyle='--')
            plt.axvspan(st, ed-1, alpha=0.5, color='red')
            #for peak in av_dyadpeaks:
            #    plt.plot(peak[0], peak[1], color='r', marker='o')
            #plt.xlim([0,225])
            plt.xlim([60,225-60])
            plt.ylim([-5,5])
            plt.xlabel('Dyad position (bp)')
            #plt.ylabel('Counts')
            plt.ylabel('Force (A.U.)')
            plt.legend(loc='best')  
            plt.savefig('dForce.png', bbox_inches = 'tight')
            plt.show()
            #plt.savefig('dyadsig_cond' + str(i+1) + note + '.png', bbox_inches = 'tight')
            plt.close()


def plot_heatmap (map, cmap, vmin, vmax, xlabel="", ylabel="", xticks=None, yticks=None, vlines=None, hlines=None, upticks=None, note=""):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    if xticks:
        ax1.set_xticks(xticks[0])
        ax1.set_xticklabels(xticks[1])
    if yticks:
        ax1.set_yticks(yticks[0])
        ax1.set_yticklabels(yticks[1])
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax2 = ax1.twiny()
    #cmap = cm.seismic
    #cmap.set_bad(set_bad,1.)
    ax2.imshow(map, cmap=cmap, interpolation='none', aspect='auto', vmin=vmin, vmax=vmax)
    if vlines:
        for x in vlines:
            ax2.axvline(x, color='k', linestyle='--', linewidth=2)
    if hlines:
        for x in hlines:
            ax2.axhline(x, color='k', linestyle='--', linewidth=2)
    if upticks:
        ax2.set_xticks(upticks[0])
        ax2.set_xticklabels(upticks[1])
    plt.show()
    plt.savefig('heatmap_' + note + '.png')
    plt.close()

    
def plot_heatmap2 (map_list, cmap, vmin, vmax, aspect='auto', xlabel="", ylabels=[], xticks=None, yticks=None, vlines=None, hlines=None, upticks=None, note=""):
    n = len(map_list)
    fig, axes = plt.subplots(1, n, sharex=False)
    for i in range(len(axes)):
        ax = axes[i]
        map = map_list[i]
        if xticks and i == n-1:
            ax.set_xticks(xticks[0])
            ax.set_xticklabels(xticks[1])
        else:
            ax.set_xticks(xticks[0])
            ax.set_xticklabels([])
        if upticks and i == 0:
            ax.xaxis.set_tick_params(labeltop='on', labelbottom='off')
            ax.set_xticks(xticks[0] + upticks[0])
            ax.set_xticklabels(['']*len(xticks[1]) + upticks[1])
        if yticks:
            ax.set_yticks(yticks[0])
            ax.set_yticklabels(yticks[1])
            #ax.set_yticks([])
            #ax.set_yticklabels([])
        if vlines:
            for x in vlines:
                ax.axvline(x, color='k', linestyle='--', linewidth=2)
        if hlines:
            for x in hlines:
                ax.axhline(x, color='k', linestyle='--', linewidth=2)
        im = ax.imshow(map, cmap=cmap, interpolation='none', aspect=3, vmin=vmin, vmax=vmax)
        if ylabels:
            ax.set_ylabel(ylabels[i], rotation=0, fontsize=15, labelpad=30)        
    ax.set_xlabel(xlabel)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist())
    #cbar.ax.set_yticks([-10, 0 ,10])
    #cbar.ax.set_yticklabels(['< -10', '0', '> 10'])
    #plt.tight_layout()
    plt.savefig('heatmap_' + note + '.png', bbox_inches='tight')
    plt.show()
    plt.close()

def plot_heatmap3 (map_list, cmap, vmin, vmax, aspect = 'auto', xlabels=[], ylabels=[], titles=[], xticks=None, yticks=None, vlines=None, hlines=None, upticks=None, note=""):
    n = len(map_list)
    fig = plt.figure()
    fig.set_size_inches([20,20])
    grid = AxesGrid(fig,
                    111,
                    nrows_ncols=(3, 1),
                    axes_pad=0.05,
                    share_all=True,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="single")
    for i in range(n):
        map = map_list[i]
        ax = grid[i]
        im = ax.imshow(map, cmap=cmap, interpolation='none', aspect=aspect, vmin=vmin, vmax=vmax)
        if xticks:
            ax.set_xticks(xticks[0])
            ax.set_xticklabels(xticks[1])
        if yticks:
            ax.set_yticks(yticks[0])
            ax.set_yticklabels(yticks[1])
        if xlabels:
            ax.set_xlabel(xlabels[i])
        if ylabels:
            ax.set_ylabel(ylabels[i])
        if titles:
            ax.set_title(titles[i], y=1.10)
        if vlines:
            for x in vlines:
                ax.axvline(x, color='k', linestyle='--', linewidth=2)
        if hlines:
            for x in hlines:
                ax.axhline(x, color='k', linestyle='--', linewidth=2)
        if upticks:
                ax2 = ax.twiny()
                ax2.set_ylim(ax.get_ylim())
                ax2.set_xlim(ax.get_xlim())
                ax2.set_xticks(upticks[0])
                ax2.set_xticklabels(upticks[1])
    grid.cbar_axes[0].colorbar(im)
    #plt.savefig('heatmap_' + note + '.png', bbox_inches='tight')
    plt.savefig('test.png', bbox_inches='tight')
    plt.show()
    plt.close()

def plot_heatmap4 (map_list, dim, cmap_list, vmin_list, vmax_list, aspect_list = [], xlabels=[], ylabels=[], titles=[], xticks_list=[], yticks_list=[], vlines_list=[], hlines_list=[], upticks_list=[], note=""):
    assert dim[0]*dim[1] == len(map_list)
    fig = plt.figure()
    #fig.set_size_inches([20,20])
    grid = AxesGrid(fig,
                    111,
                    nrows_ncols=(dim[0], dim[1]),
                    axes_pad=0.2,
                    #share_all=False,
                    label_mode="L",
                    cbar_location="right",
                    cbar_size="3%",
                    cbar_pad="2%",
                    cbar_mode="each")
    
    for i in range(len(map_list)):
        map = map_list[i]
        ax = grid[i]
        if aspect_list:
            aspect = aspect_list[i]
        else:
            aspect = "auto"
        im = ax.imshow(map, cmap=cmap_list[i], interpolation='none', aspect=aspect, vmin=vmin_list[i], vmax=vmax_list[i])
        if xticks_list:
            xticks = xticks_list[i]
            ax.set_xticks(xticks[0])
            ax.set_xticklabels(xticks[1])
        if yticks_list:
            yticks = yticks_list[i]
            ax.set_yticks(yticks[0])
            ax.set_yticklabels(yticks[1])
        if xlabels:
            ax.set_xlabel(xlabels[i])
        if ylabels:
            ax.set_ylabel(ylabels[i], rotation=0, fontsize=15, labelpad=50)        
        if titles:
            ax.set_title(titles[i], y=1.10)
        if vlines_list:
            vlines = vlines_list[i]
            for x in vlines:
                ax.axvline(x, color='k', linestyle='--', linewidth=2)
        if hlines_list:
            hlines = hlines_list[i]
            for x in hlines:
                ax.axhline(x, color='k', linestyle='--', linewidth=2)

        
        if upticks_list and i == 0:
            upticks = upticks_list[i]
            ax2 = ax.twiny()
            xticks = xticks_list[i]
            ax.set_xticklabels(['']*23)
            ax2.set_xticklabels(['']*23)
            ax2.set_ylim(ax.get_ylim())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(upticks[0])
            ax2.set_xticklabels(upticks[1])
        
        #grid.cbar_axes[i].colorbar(im)
    #ax.set_xticks(xticks[0])
    #ax.set_xticklabels(xticks[1])
    #grid.cbar_axes[0].colorbar(im)
    #plt.savefig('heatmap_' + note + '.png', bbox_inches='tight')
    plt.savefig('test.png', bbox_inches='tight')
    plt.show()
    plt.close()


# plot ratio dyad map
def plot_rmap(key_slider1, key_slider2, sample_list, norm_choice, note = "", draw_key=False, draw_vert=False):    
    for i in range(len(sample_list)):
        key_list = sample_list[i]
        cut_Rimg, cut_Limg = [], []
        dyad_img = []

        if len(key_list) <= 50:
            A,B,C,D = 5, 10, 5, 2
        else:
            A,B,C,D = 1, 1, 0, 0
            if draw_key:
                A,B,C,D = 4, 8, 0, 2

        for j in range(len(key_list)):
            key = key_list[j]
            slider = key_slider1[key]
            slider2 = key_slider2[key]
            cut_Rtmp_img, cut_Ltmp_img = [], []
            dyad_tmp_img = []

            if draw_key:
                #size, loc = key.split('-')
                win, loc = key.split('-')
                size, loc = len(win), int(loc)
                st, ed = loc, loc+size

            for side in ['R','L']:
                for k in range(A):
                    if side == 'R':
                        cut_map = slider.right_cutmap
                        cut_Rtmp_img.append(cut_map)
                        cut_Ltmp_img.append([0.0]*len(cut_map))
                    else:
                        assert side == 'L'
                        cut_map = slider.left_cutmap
                        cut_Ltmp_img.append(cut_map)
                        cut_Rtmp_img.append([0.0]*len(cut_map))

            dyadmap1 = slider.dyadmap
            dyadmap2 = slider2.dyadmap
            #dyad_map = analysis.take_ratio(dyadmap2, dyadmap1, NA=1.0)
            dyad_map = analysis.take_ratio_log(dyadmap2, dyadmap1, NA=1.0)
            dyad_map = list(dyad_map)
            for k in range(B):
                dyad_tmp_img.append(dyad_map)

            #print dyad_tmp_img
            
            if norm_choice:
                cut_Rtmp_img = normalize(cut_Rtmp_img)
                cut_Ltmp_img = normalize(cut_Ltmp_img)
                dyad_tmp_img = normalize(dyad_tmp_img)
            
            for k in range(len(cut_Rtmp_img)):
                row = cut_Rtmp_img[k]
                if draw_key and np.mean([st,ed]) < len(cut_Rtmp_img[0])/2 and k < D:
                    cut_Rimg.append(row[:st] + [np.nan]*(ed-st) + row[ed:])
                else:
                    cut_Rimg.append(row)

            for k in range(len(cut_Ltmp_img)):
                row = cut_Ltmp_img[k]
                if draw_key and np.mean([st,ed]) >= len(cut_Ltmp_img[0])/2 and k >= len(cut_Ltmp_img) - D:
                    cut_Limg.append(row[:st] + [np.nan]*(ed-st) + row[ed:])
                else: 
                    cut_Limg.append(row)

            for k in range(len(dyad_tmp_img)):
                row = dyad_tmp_img[k]
                if draw_key and k < D:
                    dyad_img.append(row[:st] + [np.nan]*(ed-st) + row[ed:])
                else:
                    dyad_img.append(row)

            if j < len(key_list)-1:
                for k in range(C):
                    cut_Rimg.append([0]*len(cut_map))
                    cut_Limg.append([0]*len(cut_map))
                    dyad_img.append([0]*len(dyad_map))

        assert np.shape(cut_Limg) == np.shape(cut_Rimg)
        #cut_img = np.matrix(cut_Limg) + np.matrix(cut_Rimg)
            
        size = np.shape(cut_Limg)
        height = float(size[0])
        width = float(size[1])
        fig = plt.figure()
        fig.set_size_inches(width/100, height/100, forward = False)
        #ax = plt.Axes(fig, [10,10,0.5,0.5])
        #ax.axis('off')
        #fig.add_axes(ax)
        ax = plt.subplot(111, aspect = 'equal')
        plt.subplots_adjust(left=0, bottom=0.3/(height/100), right=1, top=1, wspace=0, hspace=0)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tick_params(top='off', left='off', right='off', labelleft='off', labelbottom='on')
        cut_Rimg = np.asarray(cut_Rimg)
        if draw_key:
            R_coordinate = np.argwhere(np.isnan(cut_Rimg))
        norm = plt.Normalize(np.nanmin(cut_Rimg), np.nanmax(cut_Rimg))
        with np.errstate(invalid='ignore'):
            cut_Rimg = np.ma.masked_where(cut_Rimg <= 0.0, cut_Rimg)
            cmap = plt.cm.Blues
            cut_Rimg = cmap(norm(cut_Rimg))
        if draw_key:
            for x, y in R_coordinate:
                cut_Rimg[x, y, :3] = 1, 0, 0
        ax.imshow(cut_Rimg, interpolation='none')
        cut_Limg = np.asarray(cut_Limg)
        if draw_key:
            L_coordinate = np.argwhere(np.isnan(cut_Limg))
        norm = plt.Normalize(np.nanmin(cut_Limg), np.nanmax(cut_Limg))
        with np.errstate(invalid='ignore'):
            cut_Limg = np.ma.masked_where(cut_Limg <= 0.0, cut_Limg)
            cmap = plt.cm.Reds
            cut_Limg = cmap(norm(cut_Limg))
        if draw_key:
            for x, y in L_coordinate:
                cut_Limg[x, y, :3] = 1, 0, 0
        ax.imshow(cut_Limg, interpolation='none')
        if draw_vert:
            center = len(cut_Limg[0])/2
            for x in [center-30, center-20, center, center+20, center+30]:
                ax.axvline(x, color='k', linestyle='--', linewidth=0.5)
        #ax.imshow(cut_img, cmap=cmap, interpolation='none')
        #plt.show()
        plt.savefig('cut_cond' + str(i+1) + note + '.png', dpi=1500)
        #plt.savefig('cut_cond' + str(i+1) + '.svg', format = 'svg')
        #plt.show()
        plt.close()
        

        
        fig=plt.figure()
        size = np.shape(dyad_img)
        height = float(size[0])
        width = float(size[1])
        #fig.set_size_inches(width/height, 1, forward = False)
        fig.set_size_inches(width/100, height/100, forward = False)
        ax = plt.subplot(111, aspect = 'equal')
        plt.subplots_adjust(left=0, bottom=0.3/(height/100), right=1, top=1, wspace=0, hspace=0)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tick_params(top='off', left='off', right='off', labelleft='off', labelbottom='on')
        #ax=plt.Axes(fig, [0.,0.,1.,1.])
        #ax.axis('off')
        #fig.add_axes(ax)
        cmap = plt.cm.Greens
        cmap = plt.cm.get_cmap('jet')
        if draw_key:
            cmap.set_bad((1, 0, 0, 1))
        ax.imshow(dyad_img, cmap=cmap, interpolation='none', vmin=-3, vmax=3)
        if draw_vert:
            center = len(dyad_img[0])/2
            for x in [center-30, center-20, center, center+20, center+30]:
                ax.axvline(x, color='k', linestyle='--', linewidth=0.5)
        plt.savefig('dyad_cond' + str(i+1) + note + '.png', dpi=1500)
        #plt.savefig('dyad_cond' + str(i+1) + '.svg', format = 'svg')
        #plt.show()
        plt.close()
        

    
#plot the sequence logo by weblogo
def plot_weblogo(seq_list, note=''):
    write_FASTA(seq_list, 'temp')
    weblogo_cmd = ["weblogo", "-f", 'temp.fasta', '-o', 'weblogo_' + note + '.png', '-F', 'png']
    subprocess.call(weblogo_cmd)
    os.system("rm %s" % ('temp.fasta'))
"""
