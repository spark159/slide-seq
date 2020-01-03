import sys
from argparse import ArgumentParser, FileType
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import random  
import math

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def library_check (sort_fnames,
                   ref_fname,
                   type_choice,
                   out_note):
    
    def read_ref (ref_fname):
        id_seq = {}
        for line in open(ref_fname):
            line = line.strip()
            if line.startswith('>'):
                id = line[1:].strip()
                continue
            if line:
                assert id not in id_seq
                id_seq[id] = line
        return id_seq
        
    def read_sort (sort_fname, choice = type_choice):
        id_seq = {}
        type_count = {}
        chosen_ids, chosen_seqs = [], []
        readlens = []
        intErrs = []
        for line in open(sort_fname):
            if line.strip():
                read_id, type, mapped_id, cut_loc, read_seq = line.strip().split()
                if type == choice:
                    chosen_ids.append(mapped_id)
                    chosen_seqs.append(read_seq)
                if type.startswith('invalid:instErr'):
                    type = 'invalid:instErr'
                    intErrs.append(type[15:])
                if type.startswith('invalid:multHit'):
                    type = 'invalid:multHit'
                if type not in type_count:
                    type_count[type] = 0
                type_count[type] +=1
                readlens.append(len(read_seq))
        return type_count, chosen_ids, chosen_seqs, readlens

    def key_counting(key_list):
        output = {}
        for key in key_list:
            if key not in output:
                output[key] = 0
            output[key] +=1
        return output
    
    def GC_content(seq):
        num=0.0
        for nt in seq:
            if nt in 'GC':
                num+=1
        return (num/float(len(seq)))*100

    def GC_counting(seq_list):
        GC_count={}
        for seq in seq_list:
            GC=GC_content(seq)
            if GC not in GC_count:
                GC_count[GC] = 0
            GC_count[GC] +=1
        return GC_count
    
    def normalize (prob):
        total=0.0
        for key in prob:
            total += prob[key]
        for key in prob:
            prob[key]=prob[key]/float(total)
        return

    def get_hist (data, binnum=1000, prob=True):
        hist={};
        if prob:
            deno=float(len(data))
        else:
            deno=1
        binwidth=float(max(data)-min(data))/binnum
        for value in data:
            bin=int((value-min(data))/binwidth)
            bincenter=min(data)+(bin+0.5)*binwidth
            if bincenter not in hist:
                hist[bincenter]=0
            hist[bincenter]+=1/deno
        return hist
    
    def all_path(N, states='AC'):
        if N==1:
            return list(states)
        output=[]
        for path in all_path(N-1):
            for state in states:
                output.append(path+state)
        return output

    def random_pick (data, N, bias = 0):
        output = []; i = 0
        while i < N:
            pick = data[random.randint(0,len(data)-1)]
            GC = GC_content(pick)/100.0
            AT = 1.0 - GC
            if random.random() < math.exp(-bias*AT):
                output.append(pick)
                i += 1
        return output

    def nt_freq (seq_list):
        #print [len(seq) for seq in seq_list]
        maxlen = max([len(seq) for seq in seq_list])
        freq = [ {} for i in range(maxlen) ]
        for i in range(maxlen):
            N = 0.0
            for seq in seq_list:
                if i > len(seq) - 1:
                    continue
                nt = seq[i]
                if nt not in 'ATGC':
                    continue
                N += 1
                if nt not in freq[i]:
                    freq[i][nt] = 0.0
                freq[i][nt] += 1.0
            for key in freq[i]:
                freq[i][key] = freq[i][key] / N
        return freq

    def GC_profile(freq):
        AT=[]; GC=[]
        for pos in freq:
            count1=0.0; count2=0.0
            for nt in pos.keys():
                if nt in 'AT':
                    count1 +=pos[nt]
                else:
                    count2 +=pos[nt]
            AT.append(count1)
            GC.append(count2)
        return AT, GC

    def count_cmp(a, b):
        if a[1] >= b[1]:
            return -1
        else:
            return 1

    def key_cmp(a, b):
        if a[0] <= b[0]:
            return -1
        else:
            return 1


    def pick_best(key_count, N=1000):
        keycount = [[key,count] for key, count in key_count.items()]
        keycount = sorted(keycount, cmp=count_cmp)
        return keycount[:N]

    
    id_seq = read_ref(ref_fname)
    type_count_list, chosen_ids_list, chosen_seqs_list, readlens_list = [], [], [], []
    for i in range(len(sort_fnames)):
        sort_fname = sort_fnames[i]
        type_count, chosen_ids, chosen_seqs, readlens = read_sort(sort_fname)
        type_count_list.append(type_count)
        chosen_ids_list.append(chosen_ids)
        chosen_seqs_list.append(chosen_seqs)
        readlens_list.append(readlens)
        print "condition " + str(i+1)
        print "total reads:" + str(sum([type_count[type] for type in type_count]))
        print type_count
        print
    
    # pie chart for classify
    for i in range(len(sort_fnames)):
        type_count = type_count_list[i]
        types = type_count.keys(); counts = type_count.values()
        explode = []
        for key in types:
            if key == type_choice:
                explode.append(0.05)
            else:
                explode.append(0)
        colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'lightpink', 'lime']
        plt.pie(counts,colors=colors, labels=types, shadow=True, startangle=90, autopct='%1.1f%%', explode=explode)
        plt.axis('equal')
        plt.savefig('pie'+'_cond'+str(i+1) + out_note + '.png')
        #plt.show()
        plt.close()

    # read length distribution
    for i in range(len(sort_fnames)):
        chosen_seqs = chosen_seqs_list[i]
        length_list = [len(seq) for seq in chosen_seqs]
        #length_list = readlens_list[i]
        length_count = key_counting(length_list)
        lengthcount = [[length, count] for length, count in length_count.items()]
        lengthcount = sorted(lengthcount, cmp=count_cmp)
        toplen = lengthcount[0][0]
        lengthcount = sorted(lengthcount, cmp=key_cmp)
        X, Y = [], []
        for length, count in lengthcount:
            X.append(length)
            Y.append(count)
        fig = plt.figure()
        bar_width = 0.5
        plt.bar(X,Y, width=bar_width)
        X = np.asarray(X)
        plt.xticks(X + bar_width/2, X, rotation='vertical')
        plt.xlim([toplen-10, toplen+10])
        plt.xlabel('Read length (bp)')
        plt.ylabel('Read counts')
        plt.savefig("bar_" + str(i+1) + out_note + '.png')
        #plt.show()
        plt.close()

    # use reference seq instead of read seq
    chosen_seqs_list = []
    for i in range(len(sort_fnames)):
        chosen_ids = chosen_ids_list[i]
        chosen_seqs = [id_seq[id] for id in chosen_ids]
        chosen_seqs_list.append(chosen_seqs)

    
    # make simulated data
    all_ids = id_seq.keys()
    N = sys.maxint
    for i in range(len(sort_fnames)):
        chosen_seqs = chosen_seqs_list[i]
        n = len(chosen_seqs)
        if n < N:
            N = n
    # print N
    #ideal_ids = random_pick(all_ids, N)
    #ideal_seqs = []
    #for id in ideal_ids:
    #    ideal_seqs.append(id_seq[id])
    #bias_ids = random_pick(all_ids, N, bias = 1)
    #bias_seqs = []
    #for id in bias_ids:
    #    bias_seqs.append(id_seq[id])
    
    # plot GC content vs. frequency
    colors = cm.rainbow(np.linspace(0, 1, len(sort_fnames)))
    GC_count_list = []
    for i in range(len(sort_fnames)):
        chosen_seqs = chosen_seqs_list[i]
        GC_count = GC_counting(chosen_seqs); normalize(GC_count)
        GC_count_list.append(GC_count)
        plt.scatter(GC_count.keys(), GC_count.values(), color = colors[i], label='Data'+str(i+1))
        
    #idealGC_count=GC_counting(ideal_seqs); normalize(idealGC_count)
    #plt.scatter(idealGC_count.keys(), idealGC_count.values(),color='green',marker='x', label='Ideal')

    #biasGC_count=GC_counting(bias_seqs); normalize(biasGC_count)
    #plt.scatter(biasGC_count.keys(), biasGC_count.values(),color='red',marker='*', label='Biased')
    plt.xlabel('GC contents (%)')
    plt.ylabel('Frequency')
    plt.legend(loc='best')
    plt.savefig('GCfreq' + out_note + '.png')
    #plt.show()
    plt.close()

    # plot seqID vs. coverage
    id_count_list = []
    for i in range(len(sort_fnames)):
        chosen_ids = chosen_ids_list[i]
        id_count = key_counting(chosen_ids)
        id_count_list.append(id_count)
        p = plt.plot(range(len(id_count)), sorted(id_count.values()), label='Data'+str(i+1))
        plt.axvline(len(id_count), color = p[len(p)-1].get_color(), linestyle='dotted')

    #idealid_count = key_counting(ideal_ids)
    #p = plt.plot(range(len(idealid_count)), sorted(idealid_count.values()), label='Ideal')
    #plt.axvline(len(idealid_count), color = p[len(p)-1].get_color(), linestyle='dotted')

    #biasid_count = key_counting(bias_ids)
    #p = plt.plot(range(len(biasid_count)), sorted(biasid_count.values()),  label='Biased')
    #plt.axvline(len(biasid_count), color = p[len(p)-1].get_color(), linestyle='dotted')
    plt.ylim([0,max(id_count.values())*0.95])
    plt.yticks(range(0,int(max(id_count.values())*0.95),int(max(id_count.values())*0.95/15.0)))
    plt.xlabel('Seq ID')
    plt.ylabel('Counts')
    plt.legend(loc='best')
    plt.savefig('coverage' + out_note + '.png')
    #plt.show()
    plt.close()

    """
    # GC/AT contents profile
    ATprofile_list, GCprofile_list = [], []
    for i in range(len(sort_fnames)):
        chosen_seqs = chosen_seqs_list[i]
        freq = nt_freq(chosen_seqs)
        ATprofile, GCprofile = GC_profile(freq)
        ATprofile_list.append(ATprofile)
        GCprofile_list.append(GCprofile)
        plt.plot(range(len(ATprofile)), ATprofile, 'r', label='AT frequency'+str(i+1))
        plt.plot(range(len(GCprofile)), GCprofile, 'b', label='GC frequency'+str(i+1))   
    ideal_AT, ideal_GC = GC_profile(nt_freq(ideal_seqs))
    #bias_AT, bias_GC = GC_profile(nt_freq(bias_seqs))
    plt.plot(range(len(ideal_AT)), ideal_AT, 'r--')
    plt.plot(range(len(ideal_GC)), ideal_GC, 'b--')
    #plt.plot(range(len(bias_AT)), bias_AT, 'r:')
    #plt.plot(range(len(bias_GC)), bias_GC, 'b:')
    plt.xlabel('Position (bp)')
    plt.ylabel('frequency')
    plt.legend(loc='upper right')
    plt.savefig('posbias.png')
    #plt.show()
    plt.close()
    """

    # print 10 best abundance sequences
    for i in range(len(sort_fnames)):
        id_count = id_count_list[i]
        idcount = pick_best(id_count)
        print "condition: " + str(i+1)
        for j in range(10):
            #print "%s\t%d\t%s" % (idcount[j][0], idcount[j][1], id_seq[idcount[j][0]])
            print "%s\t%d" % (idcount[j][0], idcount[j][1])
        print

        size_counts = {}
        dtype_counts = {}
        for id, count in id_count.items():
            #win_size = int(id.split('-')[0])
            loc, dtype, nts = id.split('-')
            win_size = len(nts)
            #if win_size == 0:
            #    continue
            if win_size not in size_counts:
                size_counts[win_size] = []
            size_counts[win_size].append(count)
            if dtype not in dtype_counts:
                dtype_counts[dtype] = []
            dtype_counts[dtype].append(count)

        X, Y, Z = [], [], []
        X0,Y0,Z0 = [], [], []
        for size, counts in size_counts.items():
            if size == 0:
                big = np.mean(counts)
            #    X0.append(size)
            #    Y0.append(np.mean(counts))
            #    Z0.append(np.std(counts))
            #    continue
            X.append(size)
            Y.append(np.mean(counts))
            Z.append(np.std(counts))
        plt.plot(X,Y,'.')
        plt.errorbar(X,Y,yerr=Z,fmt='o')
        plt.xlabel('Poly-A len')
        plt.ylabel('Counts')
        plt.savefig('polyAVScounts_' + str(i+1) + out_note + '.png')
        #plt.show()
        plt.close()

        dtype_list = dtype_counts.keys()
        mean_list = [np.mean(dtype_counts[dtype]) for dtype in dtype_list]
        err_list = [np.std(dtype_counts[dtype]) for dtype in dtype_list]
        plt.bar(range(len(mean_list)), mean_list, yerr=err_list, width=0.25)
        plt.xticks(range(len(mean_list)), dtype_list)
        plt.xlabel('Type')
        plt.ylabel('Counts')
        plt.savefig('typeVScounts_' + str(i+1) + out_note + '.png')
        #plt.show()
        plt.close()
            

        """
        f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

        # plot the same data on both axes
        ax.plot(X,Y,'.')
        ax.errorbar(X,Y,yerr=Z,fmt='o')
        ax2.plot(X,Y,'.')
        ax2.errorbar(X,Y,yerr=Z,fmt='o')

        # zoom-in / limit the view to different portions of the data
        ax.set_ylim(big-100, big+100)  # outliers only
        ax2.set_ylim(0, 20000)  # most of the data
        ax.set_xlim([-5, 30])
        ax2.set_xlim([-5, 30])

        # hide the spines between ax and ax2
        ax.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax.xaxis.tick_top()
        ax.tick_params(labeltop='off')  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()

        # This looks pretty good, and was fairly painless, but you can get that
        # cut-out diagonal lines look with just a bit more work. The important
        # thing to know here is that in axes coordinates, which are always
        # between 0-1, spine endpoints are at these locations (0,0), (0,1),
        # (1,0), and (1,1).  Thus, we just need to put the diagonals in the
        # appropriate corners of each of our axes, and so long as we use the
        # right transform and disable clipping.

        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        # What's cool about this is that now if we vary the distance between
        # ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
        # the diagonal lines will move accordingly, and stay right at the tips
        # of the spines they are 'breaking'

        plt.xlabel('Poly-A length(bp)')
        ax.set_ylabel('Counts')
        plt.savefig('polyAVScounts_' + str(i+1) + '.png')
        plt.close()        
        """

if __name__ == '__main__':
    parser = ArgumentParser(description='Check library quality and diversity')
    parser.add_argument(metavar='-f',
                        dest="sort_fnames",
                        type=str,
                        nargs='+',
                        help='input sort filenames')
    parser.add_argument('-x',
                        dest = 'ref_fname',
                        type=str,
                        help='input reference filename')
    parser.add_argument('-c',
                        type=str,
                        dest='choice',
                        help='DNA type choice for analysis')
    parser.add_argument('-o',
                        dest='out_note',
                        type=str,
                        help='output figure note')
    
    args = parser.parse_args()

    if args.out_note == None:
        out_note = ""
    else:
        out_note = '_' + args.out_note

    
    library_check (args.sort_fnames,
                   args.ref_fname,
                   args.choice,
                   out_note)
