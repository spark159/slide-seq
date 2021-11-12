import numpy as np
import analysis_final as analysis
from sklearn.neighbors import KernelDensity

"""
def read_DNAshape(fname, id_seq):
    names = ['MGW', 'HelT', 'ProT', 'Roll']
    dic_list = [{} for i in range(len(names))]
    for i in range(len(names)):
        data = []
        for line in open(fname+"."+names[i]):
            line = line.strip()
            if line.startswith('>'):
                if data:
                    assert seq not in dic_list[i]
                    dic_list[i][seq] = data
                id = int(line[1:].strip())
                seq = id_seq[id]
                data =[]
                continue
            if line:
                temp = line.split(',')
                for k in range(len(temp)):
                    try:
                        temp[k] = float(temp[k])
                    except Exception:
                        pass
                data += temp
        assert seq not in dic_list[i]
        dic_list[i][seq] = data
    return dic_list
#seq_MGW, seq_HelT, seq_ProT, seq_Roll = read_DNAshape('phpoXp6Za', analysis.all_path(16))    
"""
class Slider:
    
    def __init__(self,
                 id,
                 ref_length,
                 tlen,
                 left_offset,
                 right_offset,
                 seq,
                 win,
                 dyadmap,
                 left_cutmap,
                 right_cutmap,
                 MGW,
                 HelT,
                 ProT,
                 Roll):
        
        self.id = id
        self.ref_length = ref_length
        self.tlen = tlen
        self.left_offset = left_offset
        self.right_offset = right_offset
        self.seq = seq
        self.win = win
        self.dyadmap = dyadmap
        self.left_cutmap = left_cutmap
        self.right_cutmap = right_cutmap
        self.MGW = MGW
        self.HelT = HelT
        self.ProT = ProT
        self.Roll = Roll

        
    def __add__ (self, other, norm_choice=True):
        assert self.tlen == other.tlen
        assert self.left_offset == other.left_offset
        assert self.right_offset == other.right_offset
        if self.id == other.id:
            new_id = self.id
        else:
            new_id = (self.id, other.id)
        if self.seq == other.seq:
            new_seq = self.seq
        else:
            new_seq = (self.seq, other.seq)

        if norm_choice:
            dyadmap1 = analysis.norm(self.dyadmap)
            dyadmap2 = analysis.norm(other.dyadmap)
            left_cutmap1 = analysis.norm(self.left_cutmap)
            left_cutmap2 = analysis.norm(other.left_cutmap)
            right_cutmap1 = analysis.norm(self.right_cutmap)
            right_cutmap2 = analysis.norm(other.right_cutmap)
            
        new_dyadmap = [dyadmap1[i] + dyadmap2[i] for i in range(self.ref_length)]
        new_left_cutmap = [left_cutmap1[i] + left_cutmap2[i] for i in range(self.ref_length)]
        new_right_cutmap = [right_cutmap1[i] + right_cutmap2[i] for i in range(self.ref_length)]
        return Slider(new_id,
                      self.ref_length,
                      self.dyad_axis,
                      self.left_offset,
                      self.right_offset,
                      new_seq,
                      new_dyadmap,
                      new_left_cutmap,
                      new_right_cutmap)

    def __sub__ (self, other, norm_choice=True):
        assert self.tlen == other.tlen
        assert self.left_offset == other.left_offset
        assert self.right_offset == other.right_offset
        if self.id == other.id:
            new_id = self.id
        else:
            new_id = (self.id, other.id)
        if self.seq == other.seq:
            new_seq = self.seq
        else:
            new_seq = (self.seq, other.seq)

        if norm_choice:
            #dyadmap1 = analysis.norm(self.dyadmap)
            #dyadmap2 = analysis.norm(other.dyadmap)
            dyadmap1 = self.KDE(band_width=0.5)
            dyadmap2 = other.KDE(band_width=0.5)
            left_cutmap1 = analysis.norm(self.left_cutmap)
            left_cutmap2 = analysis.norm(other.left_cutmap)
            right_cutmap1 = analysis.norm(self.right_cutmap)
            right_cutmap2 = analysis.norm(other.right_cutmap)
            
        new_dyadmap = [dyadmap1[i] / dyadmap2[i] for i in range(self.ref_length)]
        new_left_cutmap = [left_cutmap1[i] - left_cutmap2[i] for i in range(self.ref_length)]
        new_right_cutmap = [right_cutmap1[i] - right_cutmap2[i] for i in range(self.ref_length)]
        return Slider(new_id,
                      self.ref_length,
                      self.dyad_axis,
                      self.left_offset,
                      self.right_offset,
                      new_seq,
                      new_dyadmap,
                      new_left_cutmap,
                      new_right_cutmap)

    # get top strand cutmap
    def get_top_cutmap(self):
        return self.right_cutmap

    # get bottom strand cutmap
    def get_bottom_cutmap(self):
        return self.left_cutmap

    # get raw dyadmap
    def get_dyadmap(self):
        return self.dyadmap

    # get nucleosome positioning signal
    def get_psig(self, left_bound=0, right_bound=0):
        dyadmap = self.dyadmap[left_bound:len(self.dyadmap)-right_bound]
        return np.asarray(analysis.normalize_list(dyadmap))

    # find GC content of the reference sequence
    def GC_content(self):
        return analysis.GC_content(self.seq)

    # find the longest length of poly-nt (pos=False)
    # find all poly-nt locations and counts (pos=True)
    def polynt_count(self, nts, pos=False):
        return analysis.polynt_count(self.seq, nts, pos=False)

    # find total read counts
    def read_counts (self, choice):
        if choice == 'dyad':
            return sum(self.dyadmap)
        elif choice == 'left':
            return sum(self.left_cutmap)
        elif choice == 'right':
            return sum(self.right_cutmap)

    # find peaks of chosen data
    def find_peaks (self, choice, num=None):
        if choice == 'dyad':
            return analysis.find_peaks(self.dyadmap, num=num)
        elif choice == 'left':
            return analysis.find_peaks(self.left_cutmap, num=num)
        elif choice == 'right':
            return analysis.find_peaks(self.right_cutmap, num=num)

    # median value of nucleosome positioning
    def median_pos (self, scale=1.0, excluded=[], selected=None):
        data = []
        if selected == None:
            selected = range(len(self.dyadmap))
        for i in selected:
            if i in excluded:
                continue
            for k in range(int(scale*self.dyadmap[i])):
                data.append(i)
        return np.median(data)

    # mean value of nucleosome positioning
    def mean_pos (self, excluded=[], selected=None):
        mean_pos, total = 0.0, 0.0
        if selected == None:
            selected = range(len(self.dyadmap))
        for i in selected:
            if i in excluded:
                continue
            mean_pos += i*self.dyadmap[i]
            total += self.dyadmap[i]
        return mean_pos/total

    # find/broden peaks of signal and igonore others 
    def highlight_peaks (self, choice='dyad', num=15, broden=1):
        if choice == 'dyad':
            return analysis.highlight_peaks(self.dyadmap, num=num, broden=broden)
        elif choice == 'left':
            return analysis.highlight_peaks(self.left_cutmap, num=num, broden=broden)
        elif choice == 'right':
            return analysis.highlight_peaks(self.right_cutmap, num=num, broden=broden)

    # Kernel density estimation of nucleosome positioning signal
    def KDE (self, scale=1.0, band_width=0.5):
        log_density = analysis.logKDE(self.dyadmap, scale=scale, band_width=band_width)
        return np.exp(log_density)

    # compute entropy of nucleosome positioning signal
    def entropy (self, scale=1.0, band_width=0.5):
        log_density = analysis.logKDE(self.dyadmap, scale=scale, band_width=band_width)
        kde = np.exp(log_density)
        entropy = 0.0
        for i in range(len(log_density)):
            entropy += -kde[i]*log_density[i]
        return entropy

    # energy landscape of nucleosome positioning 
    def energy_profile (self, kT=1, scale=1.0, band_width=0.5):
        log_density = analysis.logKDE(self.dyadmap, scale=scale, band_width=band_width)
        return -kT * log_density

    # force profile of nucleosome positioning (To do)
    def force_profile (self, kT=1, scale=1.0, band_width=0.5):
        X = []
        for k in range(len(self.dyadmap)):
            for num in range(int(self.dyadmap[k])):
                X.append([k])
        #X_plot = np.linspace(0,len(self.dyadmap), num = (len(self.dyadmap)+1)*5)[:,np.newaxis]
        X_plot = np.arange(0,len(self.dyadmap), 0.1)[:,np.newaxis]
        kde = KernelDensity(kernel="gaussian", bandwidth=band_width).fit(X)
        energy_profile = -kT* kde.score_samples(X_plot)    
        X_axis = list(X_plot[:,0])
        
        energy_profile = self.energy_profile(self.dyadmap, kT=kT, scale=scale, band_width=band_width)
        force_profile = []
        for i in range(len(self.dyadmap)):
            idx = X_axis.index(i)
            dX = X_axis[idx+1] - X_axis[idx]
            denergy = energy_profile[idx+1] - energy_profile[idx]
            force = - denergy / dX
            force_profile.append(force)
        return np.asarray(force_profile)
