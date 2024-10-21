import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import load
import numpy as np

def norm (L):
    total = 0.0
    for value in L:
        total += value
    return [ value/total for value in L]

def key_cmp(k1,k2):
    win1, loc1 = k1.split('-')
    st1 = int(loc1)
    win2, loc2 = k2.split('-')
    st2 = int(loc2)
    if len(win1) < len(win2):
        return -1
    elif len(win1) > len(win2):
        return 1
    else:
        if st1 < st2:
            return -1
        else:
            return 1
        
ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
filenames1 = ["Plslib-HS_S1_L001_R.sort"]
filenames2 = ["Plslib-HS-30min_S2_L001_R.sort"]
ref_fname = "plusonelib.ref"

#filenames2 = ["../../Illumina/temp/plusone-0_S1_L001_R.sort"]
#filenames1 = ["../../Illumina/temp/plusone-30_S2_L001_R.sort"]
#ref_fname = "../../Illumina/temp/plusonelib.ref"


key_slider1 = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, filter_num = 50, fill=None, load_ref=ref_fname)
key_slider2 = load.load_files(filenames2, ref_length, dyad_axis, dyad_offset, filter_num = 50, fill=None, load_ref=ref_fname)

#filenames1 = ["../../Illumina/SW_Ascan_new/data/Ascan0_S1_L001_R.sort"]
#filenames2 = ["../../Illumina/SW_Ascan_new/data/Ascan-5min_S1_L001_R.sort"]
#ref_fname = "../../Illumina/SW_Ascan_new/polyAscanlib.ref"
#key_slider3 = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, filter_num = 50, fill='linear', load_ref=ref_fname)
#key_slider4 = load.load_files(filenames2, ref_length, dyad_axis, dyad_offset, filter_num = 50, fill='linear', load_ref=ref_fname)


# key for plusone
pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
keys = list(set(key_slider1.keys()) & set(key_slider2.keys()))
keys = [int(key) for key in keys]
keys = sorted(keys)
keys = [str(key) for key in keys]
labels = ['Before', 'After']



# side by side
i = 0
while i < len(keys):
#while i < 10:
    j = 0
    #fig = plt.figure(figsize=(15,20))
    fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    #fig = plt.figure(figsize=(11.69,8.27))
    while j < min(5, len(keys) - i):
        key = keys[i]
        #dyadmap1 = norm(key_slider1[key].dyadmap)
        dyadmap1 = key_slider1[key].dyadmap
        #dyadmap2 = norm(key_slider2[key].dyadmap)
        dyadmap2 = key_slider2[key].dyadmap
        dyadmap2 = dyadmap1
        #cutmaps1= norm(key_slider1[key].left_cutmap + key_slider1[key].right_cutmap)
        cutmaps1= key_slider1[key].left_cutmap + key_slider1[key].right_cutmap
        bcutmap1, tcutmap1 = cutmaps1[:225], cutmaps1[225:]
        #cutmaps2= norm(key_slider2[key].left_cutmap + key_slider2[key].right_cutmap)
        cutmaps2= key_slider2[key].left_cutmap + key_slider2[key].right_cutmap
        cutmaps2= cutmaps1
        bcutmap2, tcutmap2 = cutmaps2[:225], cutmaps2[225:]
        #maps = [[2*(np.asarray(tcutmap1)+np.asarray(bcutmap1)),2*(np.asarray(tcutmap2)+np.asarray(bcutmap2))],[dyadmap1,dyadmap2] ]
        maps = [[2*np.asarray(tcutmap1), 2*np.asarray(bcutmap1), 2*np.asarray(tcutmap2), 2*np.asarray(bcutmap2)],[dyadmap1, dyadmap2] ]
        labels = [["Before(t)", "Before(b)", "After(t)", "After(b)"], ["Before", "After"]] 
        for k in range(2):
            plt.subplot(5,2,j*2+k+1)
            for u in range(len(maps[k])):
                map = maps[k][u]
                label = labels[k][u]
                plt.plot(map, label=key +":"+label, alpha=1)
            leg = plt.legend(loc='upper right', frameon=False)
            #for lh in leg.legendHandles:
            #    lh.set_visible(False)
            #plt.title(key+" " +labels[k])
            plt.xlim([0,225])
            #plt.ylim([0,1])
        i +=1
        j +=1
    pdf.savefig( fig )
    plt.close()
pdf.close()

"""
# same figure
i = 0
while i < len(keys):
#while i < 10:
    j = 0
    #fig = plt.figure(figsize=(15,20))
    fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    #fig = plt.figure(figsize=(11.69,8.27))
    while j < 5:
        for k in range(min(2, len(keys) - i)):
            key = keys[i]
            dyadmap1 = norm(key_slider1[key].dyadmap)
            dyadmap2 = norm(key_slider2[key].dyadmap)
            plt.subplot(5,2,j*2+k+1)            
            plt.plot(dyadmap1, label=key +":Before")
            plt.plot(dyadmap2, label=key +":After")
            leg = plt.legend(loc='upper right', frameon=False)
            #for lh in leg.legendHandles:
            #    lh.set_visible(False)
            #plt.title(key+" " +labels[k])
            plt.xlim([0,225])
            #plt.ylim([0,1])
            i +=1
        j +=1
    pdf.savefig( fig )
    plt.close()
pdf.close()


# key for Ascan
pdf = matplotlib.backends.backend_pdf.PdfPages("polyAscan_overlap.pdf")
keys = list(set(key_slider3.keys()) & set(key_slider4.keys()))
#keys = [int(key) for key in keys]
keys = sorted(keys, cmp=key_cmp)
#keys = [str(key) for key in keys]
labels = ['Before', 'After']

# side by side 
i = 0
while i < len(keys):
#while i < 10:
    j = 0
    #fig = plt.figure(figsize=(15,20))
    fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    #fig = plt.figure(figsize=(11.69,8.27))
    while j < min(5, len(keys) - i):
        key = keys[i]
        dyadmap1 = norm(key_slider3[key].dyadmap)
        dyadmap2 = norm(key_slider4[key].dyadmap)
        cutmaps1= norm(key_slider3[key].left_cutmap + key_slider3[key].right_cutmap)
        bcutmap1, tcutmap1 = cutmaps1[:225], cutmaps1[225:]
        cutmaps2= norm(key_slider4[key].left_cutmap + key_slider4[key].right_cutmap)
        bcutmap2, tcutmap2 = cutmaps2[:225], cutmaps2[225:]
        #maps = [[2*(np.asarray(tcutmap1)+np.asarray(bcutmap1)),2*(np.asarray(tcutmap2)+np.asarray(bcutmap2))],[dyadmap1,dyadmap2] ]
        maps = [[2*np.asarray(tcutmap1), 2*np.asarray(bcutmap1), 2*np.asarray(tcutmap2), 2*np.asarray(bcutmap2)],[dyadmap1, dyadmap2] ]
        labels = [["Before(t)", "Before(b)", "After(t)", "After(b)"], ["Before", "After"]]
        
        for k in range(2):
            plt.subplot(5,2,j*2+k+1)
            # for Ascan
            win, loc = key.split('-')
            st = int(loc)
            ed = st + len(win)
            name = str(len(win))+"-"+str(loc)
            plt.axvspan(st, ed-1, alpha=0.2, color='red')
            for u in range(len(maps[k])):
                map = maps[k][u]
                label = labels[k][u]
                plt.plot(map, label=name +":"+label, alpha=0.5)
            leg = plt.legend(loc='upper right', frameon=False)
            #for lh in leg.legendHandles:
            #    lh.set_visible(False)
            #plt.title(key+" " +labels[k])
            plt.xlim([0,225])
            #plt.ylim([0,1])
        i +=1
        j +=1
    pdf.savefig( fig )
    plt.close()
pdf.close()



# same figure
i = 0
while i < len(keys):
#while i < 10:
    j = 0
    #fig = plt.figure(figsize=(15,20))
    fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    #fig = plt.figure(figsize=(11.69,8.27))
    while j < 5:
        for k in range(min(2, len(keys) - i)):
            key = keys[i]
            # for Ascan
            win, loc = key.split('-')
            st = int(loc)
            ed = st + len(win)
            plt.axvspan(st, ed-1, alpha=0.2, color='red')

            dyadmap1 = norm(key_slider3[key].dyadmap)
            dyadmap2 = norm(key_slider4[key].dyadmap)
            plt.subplot(5,2,j*2+k+1)            
            plt.plot(dyadmap1, label=str(len(win))+'-'+str(loc) +":" +":Before")
            plt.plot(dyadmap2, label=str(len(win))+'-'+str(loc) +":" +":After")
            leg = plt.legend(loc='upper right', frameon=False)
            #for lh in leg.legendHandles:
            #    lh.set_visible(False)
            #plt.title(key+" " +labels[k])
            plt.xlim([0,225])
            #plt.ylim([0,1])
            i +=1
        j +=1
    pdf.savefig( fig )
    plt.close()
pdf.close()

# side by side
i = 0
while i < len(keys):
#while i < 10:
    j = 0
    #fig = plt.figure(figsize=(15,20))
    fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    #fig = plt.figure(figsize=(11.69,8.27))
    while j < min(5, len(keys) - i):
        key = keys[i]
        dyadmap1 = norm(key_slider3[key].dyadmap)
        dyadmap2 = norm(key_slider4[key].dyadmap)
        dyadmaps = [dyadmap1, dyadmap2]
        for k in range(2):
            plt.subplot(5,2,j*2+k+1)

            # for Ascan
            win, loc = key.split('-')
            st = int(loc)
            ed = st + len(win)
            plt.axvspan(st, ed-1, alpha=0.2, color='red')
            plt.plot(dyadmaps[k], label=str(len(win))+'-'+str(loc) +":"+labels[k], color='blue')
            # for Ascan
            
            #plt.plot(dyadmaps[k], label=key +":"+labels[k])
            leg = plt.legend(loc='upper right', frameon=False)
            for lh in leg.legendHandles:
                lh.set_visible(False)
            #plt.title(key+" " +labels[k])
            plt.xlim([0,225])
            #plt.ylim([0,1])
        i +=1
        j +=1
    pdf.savefig( fig )
    plt.close()
pdf.close()
"""
