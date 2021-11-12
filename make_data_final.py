import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import imputation_final as imput
from SliderClass_final import Slider
import numpy as np

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def read_DNAshape(fname):
    names = ['MGW', 'HelT', 'ProT', 'Roll']
    dic_list = [{} for i in range(len(names))]
    for i in range(len(names)):
        data = []
        for line in open(fname+"."+names[i]):
            line = line.strip()
            if line.startswith('>'):
                if data:
                    assert id not in dic_list[i]
                    dic_list[i][id] = data
                id = line[1:].strip()
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
        assert id not in dic_list[i]
        dic_list[i][id] = data
    return dic_list

def group_cmp (a, b):
    mtype1, nts1 = a.split('-')
    mtype2, nts2 = b.split('-')
    if mtype1 < mtype2:
        return -1
    elif mtype1 > mtype2:
        return 1
    else:
        if nts1 < nts2:
            return -1
        else:
            return 1

def Bayesian_wrapper (arg):
    args, kwargs = arg
    return imput.Bayesian(*args, **kwargs)

def pick_median (data_list):
    median_data = []
    n, m = len(data_list[0]), len(data_list[0][0])
    for j in range(n):
        row = []
        for k in range(m):
            row.append(np.median([data_list[i][j][k] for i in range(len(data_list))]))
        median_data.append(row)
    return median_data
    
def make_data (sort_filenames,
               ref_fname,
               shape_fname,
               out_fname,
               ref_length,
               NCP_len,
               dyad_offset,
               data_choice,
               mtype_choice,
               filter_num,
               cut_pad,
               fill):

    # collect the cleavage counts from sort files
    cutmaps, dyadmaps = {}, {} # read counts on the sliding template
    id_win = {} # window on the sliding template
    id_tlen = {} # sliding template length
    for sort_filename in sort_filenames:
        for line in open(sort_filename, 'r'):
            if line.strip():
                read_id, type, mapped_id, cutloc, seq = line.strip().split()
                if type not in data_choice:
                    continue
                if mapped_id == 'BACKBONE':
                    continue

                cols = cutloc.split(':')
                side, cut = cols[0], int(cols[1])

                tlen = ref_length
                
                # backbone-based library case
                try:
                    loc, mtype, nts = mapped_id.split('-')
                    loc = int(loc)
                    if mtype_choice != None:
                        if mtype not in mtype_choice:
                            continue

                    # identification window on the reference
                    # inclusive end, 1nt padding
                    if mtype == 'M':
                        wst = loc - 1
                        wed = loc + len(nts)
                    elif mtype == 'I':
                        wst = loc
                        wed = loc + len(nts) + 1
                    elif mtype == 'D':
                        wst = loc - 1
                        wed = loc 

                    # giving additional padding both sides of windows
                    wst -= cut_pad
                    wed += cut_pad

                    # filter out reads not including the window
                    if side == 'R' and cut >= wst:
                        continue
                    if side == 'L' and cut <= wed:
                        continue

                    # transform data from the ref sequence to the sliding template
                    # shorten the template length by deletions
                    if mtype == 'D':
                        tlen -= len(nts)

                    # adjust cut locations based on the template
                    if side == 'L' and mtype == 'I':
                        cut -= len(nts)

                    # adjust window locations based on the template
                    if mtype == 'I':
                        wed -= len(nts)

                    win = (wst, wed)
                                            
                except:
                    win = None
                    pass
                    
                if side == 'L':
                    offset = -dyad_offset
                else:
                    assert side == 'R'
                    offset = dyad_offset
                if cut+offset < 0  or cut+offset >= tlen:
                    continue
                if mapped_id not in id_tlen:
                    id_tlen[mapped_id] = tlen
                if mapped_id not in id_win:
                    id_win[mapped_id] = win
                if mapped_id not in dyadmaps:
                    dyadmaps[mapped_id] = [0.0]*tlen
                dyadmaps[mapped_id][cut+offset] += 1.0
                if mapped_id not in cutmaps:
                    cutmaps[mapped_id] = {}
                    cutmaps[mapped_id]['L'], cutmaps[mapped_id]['R'] = [0.0]*tlen, [0.0]*tlen
                cutmaps[mapped_id][side][cut] += 1


    # read reference file
    id_seq = {}
    if ref_fname:
        for line in open(ref_fname):
            line = line.strip()
            if line.startswith('>'):
                id = line[1:]
            else:
                seq = line
                assert id not in id_seq
                id_seq[id] = seq
    else:
        for id in dyadmaps.keys():
            id_seq[id] = None

    # read shape file
    if shape_fname:
        id_MGW, id_HelT, id_ProT, id_Roll = read_DNAshape(shape_fname)
        #print id_ProT
    else:
        id_MGW, id_HelT, id_ProT, id_Roll = {}, {}, {}, {}
        for id in dyadmaps.keys():
            id_MGW[id] = None; id_HelT[id] = None; id_ProT[id] = None; id_Roll[id] = None
            
        
    # imputation of missing data
    if fill != 'NONE':
        target_ids = []
        for id in id_win:
            win = id_win[id]
            if win != None:
                target_ids.append(id)
        
    if fill == "NAIVE":
        if len(target_ids) <= 0:
            print >> sys.stderr, "Warning:imputation is not possible."

        else:            
            wins, tops, botts = [], [], []
            for id in target_ids:
                wins.append(id_win[id])
                tops.append(cutmaps[id]['R'])
                botts.append(cutmaps[id]['L'])

            dyads, mtops, mbotts = imput.Naive (wins,
                                                tops,
                                                botts,
                                                dyad_offset)        

            for i in range(len(target_ids)):
                id = target_ids[i]
                dyadmaps[id] = dyads[i]
                wst, wed = id_win[id]
                tlen = len(dyads[i])
                for j in range(tlen):
                    if j <= wed:
                        cutmaps[id]['L'][j] = mbotts[i][j]
                    if j >= wst:
                        cutmaps[id]['R'][j] = mtops[i][j]
        
    elif fill == "LINEAR":
        if len(target_ids) <= 0:
            print >> sys.stderr, "Warning:imputation is not possible."

        else:            
            # group by mtype and size
            group_ids = {}
            for id in target_ids:
                loc, mtype, nts = id.split('-')
                group = (mtype, len(nts))
                if group not in group_ids:
                    group_ids[group] = []
                group_ids[group].append(id)

            groups = sorted(group_ids.keys())

            for group in groups:
                ids = group_ids[group]
                wins = [id_win[id] for id in ids]
                tops = [cutmaps[id]['R'] for id in ids]
                botts = [cutmaps[id]['L'] for id in ids]

                tlen = len(tops[0])
                graph = False
                
                ## temporal
                #if group[1] == 15:
                #    graph = True
                #else:
                #    graph = False

                #print >> sys.stderr, "Data for mtype:%s size:%d" % group
                dyads, mtops, mbotts = imput.Linear_reg (wins,
                                                         tops,
                                                         botts,
                                                         tlen,
                                                         dyad_offset,
                                                         graph=graph,
                                                         silent=True)

                for i in range(len(ids)):
                    id = ids[i]
                    dyadmaps[id] = dyads[i]
                    wst, wed = id_win[id]
                    tlen = len(dyads[i])
                    for j in range(tlen):
                        if j <= wed:
                            cutmaps[id]['L'][j] = mbotts[i][j]
                        if j >= wst:
                            cutmaps[id]['R'][j] = mtops[i][j]
                            
            print >> sys.stderr, "Linear regression imputation is done"

                
    elif fill.startswith("BAYESIAN"):
        if len(target_ids) <= 0:
            print >> sys.stderr, "Warning:imputation is not possible."
        else:            
            try:
                fill, cycle = fill.split(':')
                cycle = int(cycle)
            except:
                cycle = 1000

            # group by mtype and size
            group_ids = {}
            for id in target_ids:
                loc, mtype, nts = id.split('-')
                group = (mtype, len(nts))
                if group not in group_ids:
                    group_ids[group] = []
                group_ids[group].append(id)

            groups = sorted(group_ids.keys())

            arg_list = []
            left_bound, right_bound = NCP_len/2, NCP_len/2
            for group in groups:
                ids = group_ids[group]
                wins = [id_win[id] for id in ids]
                tops = [cutmaps[id]['R'] for id in ids]
                botts = [cutmaps[id]['L'] for id in ids]
                tlen = len(tops[0])
                args = (cycle,
                        wins,
                        tops,
                        botts,
                        tlen,
                        dyad_offset,
                        left_bound,
                        right_bound)
                kwargs = {"dyad_alphas":None,
                          "dyad_betas":None,
                          "r_alpha":None,
                          "r_beta":None,
                          "l_mu":None,
                          "l_sigma":None,
                          "initial":"linear",
                          "note":"_bayesian_" + group}
                arg_list.append((args, kwargs))

            # parallel computation for each group
            import multiprocessing as mp
            pool = mp.Pool()
            pool.map(Bayesian_wrapper, arg_list)
            
            for group in groups:
                dyads = pick_median(imput.load_data("dyads_bayesian_" + group + ".txt"))
                mtops = pick_median(imput.load_data("mtops_bayesian_" + group + ".txt"))
                mbotts = pick_median(imput.load_data("mbotts_bayesian_" + group + ".txt"))
                ids = group_ids[group]
                for i in range(len(ids)):
                    id = ids[i]
                    dyadmaps[id] = dyads[i]
                    wst, wed = id_win[id]
                    tlen = len(dyads[i])
                    for j in range(tlen):
                        if j <= wed:
                            cutmaps[id]['L'][j] = mbotts[i][j]
                        if j >= wst:
                            cutmaps[id]['R'][j] = mtops[i][j]

    else:
        assert fill == 'NONE'


    # filter by data counts
    if filter_num:
        filtered = []
        for mapped_id in cutmaps:
            if sum(dyadmaps[mapped_id]) > filter_num:
                filtered.append(mapped_id)
        cutmaps = {k:cutmaps[k] for k in filtered}
        dyadmaps = {k:dyadmaps[k] for k in filtered}        

    assert len(cutmaps) == len(dyadmaps)

    # write output file
    extension = out_fname.rsplit('.', 1)[1]

    if extension == 'data':
        f = open(out_fname, 'w')
        print >> f, '@ref_length\t%s' % (ref_length)
        print >> f, '@left_offset\t%s' % (dyad_offset)
        print >> f, '@right_offset\t%s' % (dyad_offset)
        print >> f, ''

    elif extension == 'pickle':
        import pickle
        key_slider = {}

    count = 0
    for mapped_id in sorted(cutmaps.keys()):
        if extension == 'data':
            print >> f, '>%s' % (mapped_id)        

            if ref_fname:
                print >> f, '@Sequence'
                print >> f, id_seq[mapped_id]

            if shape_fname:
                print >> f, '@MinorGrooveWidth'
                print >> f, ','.join([str(value) for value in id_MGW[mapped_id]])
                print >> f, '@HelixTwist'
                print >> f, ','.join([str(value) for value in id_HelT[mapped_id]])
                print >> f, '@PropellerTwist'
                print >> f, ','.join([str(value) for value in id_ProT[mapped_id]])
                print >> f, '@Roll'
                print >> f, ','.join([str(value) for value in id_Roll[mapped_id]])

            if id_win[mapped_id] != None:
                wst, wed = id_win[mapped_id]
                print >> f, '@IdentificationWindow'
                print >> f, wst, wed
                print >> f, '@ImputationMethod'
                print >> f, fill

            print >> f, '@TopCleavageCounts'
            print >> f, ','.join([str(value) for value in cutmaps[mapped_id]['R']])
            print >> f, '@BottomCleavageCounts'
            print >> f, ','.join([str(value) for value in cutmaps[mapped_id]['L']])
            print >> f, '@PositioningSignal'
            print >> f, ','.join([str(value) for value in dyadmaps[mapped_id]])
            print >> f, ''

        elif extension == 'pickle':
            key_slider[mapped_id] = Slider(id = mapped_id,
                                           ref_length = ref_length,
                                           tlen = id_tlen[mapped_id],
                                           left_offset = dyad_offset,
                                           right_offset = dyad_offset,
                                           seq = id_seq[mapped_id],
                                           win = id_win[mapped_id],
                                           dyadmap = dyadmaps[mapped_id],
                                           left_cutmap = cutmaps[mapped_id]['L'],
                                           right_cutmap = cutmaps[mapped_id]['R'],
                                           MGW = id_MGW[mapped_id],
                                           HelT = id_HelT[mapped_id],
                                           ProT = id_ProT[mapped_id],
                                           Roll = id_Roll[mapped_id])
        count +=1

    if extension == 'pickle':
        f = open(out_fname, "wb")
        pickle.dump(key_slider, f)

    f.close()
    
    if ref_fname:
        report = str(count) + " / " + str(len(id_seq))
    else:
        report = str(count)
    print >> sys.stderr, "Data sequence number: " + report
    print >> sys.stderr, "Done"


if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')
                
    parser = ArgumentParser(description='process sort-files to make data-files')
    parser.add_argument(dest="sort_filenames",
                        type=str,
                        nargs='+',
                        help='sort-file names')
    parser.add_argument('-x',
                        dest='ref_fname',
                        type=str,
                        help='reference filename')
    parser.add_argument('-s',
                        dest='shape_fname',
                        type=str,
                        help='DNA shape filename')
    parser.add_argument('-o',
                        dest='out_fname',
                        type=str,
                        help='output filename prefix')
    parser.add_argument('--oe',
                        dest='out_extension',
                        type=str,
                        default='data',
                        help='output file extension (.data or .pickle)')    
    parser.add_argument('--reflen',
                        dest="ref_length",
                        type=int,
                        default=225,
                        help='reference\backbone DNA length (bp)')
    parser.add_argument('--nlen',
                        dest="NCP_len",
                        type=int,
                        default=147,
                        help='Nucleosomal DNA length (bp)')
    parser.add_argument('--offset',
                        dest="dyad_offset",
                        type=int,
                        default=53,
                        help='dyad offset (bp)')
    parser.add_argument('--choice',
                        dest="data_choice",
                        type=str,
                        nargs='+',
                        help='data type choice')
    parser.add_argument('--mchoice',
                        dest="mtype_choice",
                        type=str,
                        nargs='+',
                        help='mutation type choice for backbone-based library')
    parser.add_argument('--filter',
                        dest="filter_num",
                        type=int,
                        help='filtering threshold')
    parser.add_argument('--pad',
                        dest='cut_pad',
                        type=int,
                        default=0,
                        help='cut padding size (bp)')
    parser.add_argument('--fill',
                        dest='fill',
                        type=str,
                        help='imputation method (naive, linear, bayesian:cycle)')
    
    args = parser.parse_args()

    if not args.data_choice:
        data_choice = ['valid']
    else:
        data_choice = args.data_choice

    if not args.mtype_choice:
        mtype_choice = None
    else:
        mtype_choice = args.mtype_choice

    if not args.out_extension:
        extension = 'data'
    else:
        extension = args.out_extension    

    if not args.out_fname:
        def common (strings):
            end = min([len(string) for string in strings])
            if end <= 0:
                return None
            while end > 0:
                check = strings[0][:end]
                find = True
                for i in range(1, len(strings)):
                    part = strings[i][:end]
                    if check != part:
                        find = False
                        break
                if find:
                    break
                end -= 1
            if find:
                return check
            else:
                return None            
        common_name = common([fname.rsplit('.',1)[0] for fname in args.sort_filenames])

        if common_name != None:
            out_fname = common_name + '.' + extension
        else:
            out_fname = 'out' + '.' + extension
    else:
        out_fname = args.out_fname + '.' + extension

    if args.fill:
        fill = args.fill.upper()
    else:
        fill = 'NONE'

    if args.filter_num:
        filter_num = args.filter_num
    else:
        filter_num = None
    
    make_data (args.sort_filenames,
               args.ref_fname,
               args.shape_fname,
               out_fname,
               args.ref_length,
               args.NCP_len,
               args.dyad_offset,
               data_choice,
               mtype_choice,
               filter_num,
               args.cut_pad,
               fill)
