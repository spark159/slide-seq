import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import imputation_edit as imput
from SliderClass import Slider
import numpy as np

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def size_loc_cmp (a, b):
    loc1, mtype1, nts1 = a.split('-')
    loc2, mtype2, nts2 = b.split('-')
    loc1, loc2 = int(loc1), int(loc2)
    size1, size2 = len(nts1), len(nts2)
    if size1 < size2:
        return -1
    elif size1 == size2:
        if loc1 < loc2:
            return -1
        else:
            return 1
    else:
        return 1

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
    cutmaps, dyadmaps = {}, {}
    for i in range(len(sort_filenames)):
        sort_filename = sort_filenames[i]
        for line in open(sort_filename, 'r'):
            if line.strip():
                read_id, type, mapped_id, cutloc, seq = line.strip().split()
                if type not in data_choice:
                    continue

                cols = cutloc.split(':')
                if int(cols[1]) < 0 or int(cols[1]) >= ref_length:
                    continue
                side, cut = cols[0], int(cols[1])

                # backbone-based library case
                try:
                    loc, mtype, nts = mapped_id.split('-')
                    loc = int(loc)
                    if mtype_choice != None:
                        if mtype not in mtype_choice:
                            continue

                    if side == 'R' and loc - cut <= cut_pad:
                        continue
                    if side == 'L' and cut - (loc + len(nts) - 1) <= cut_pad:
                        continue

                    ## temporal deletion case
                    #if mtype == 'D':
                    #    nts = len(nts)*"N"
                    #    mapped_id = '-'.join([str(loc), mtype, nts])

                    ## temporal 1bp case with cut_pad = 3
                    #if mtype in 'MID' and len(nts) == 1:
                    #    if side == 'R' and loc - cut <= 3:
                    #        continue
                    #    if side == 'L' and cut - loc <= 3:
                    #        continue
                    
                except:
                    pass
                    
                if side == 'L':
                    offset = -dyad_offset
                else:
                    assert side == 'R'
                    offset = dyad_offset
                if cut+offset < 0  or cut+offset >= ref_length:
                    continue
                if mapped_id not in dyadmaps:
                    dyadmaps[mapped_id] = [0.0]*ref_length
                dyadmaps[mapped_id][cut+offset] += 1.0
                if mapped_id not in cutmaps:
                    cutmaps[mapped_id] = {}
                    cutmaps[mapped_id]['L'], cutmaps[mapped_id]['R'] = [0.0]*ref_length, [0.0]*ref_length
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
    if fill == "NAIVE":
        keys = sorted(cutmaps.keys(), cmp=size_loc_cmp)

        sub_keys = []
        temp_id_seq = {}
        top, bott = [], []
        for key in keys:
            try:
                loc, mtype, nts = key.split('-')
            except:
                continue
            sub_keys.append(key)
            temp_loc = int(loc) - cut_pad
            temp_nts = 'N'*(len(nts) + 2*cut_pad)
            #temp_nts = 'N'*(cut_pad) + nts + 'N'*(cut_pad)
            temp_id = temp_nts + '-' + str(temp_loc)
            if temp_id not in temp_id_seq:
                temp_id_seq[temp_id] = id_seq[key]
            top.append(cutmaps[key]['R'])
            bott.append(cutmaps[key]['L'])

        if len(sub_keys) <= 0:
            print >> sys.stderr, "Warning:imputation is not possible."
        else:
            
            dyad, mtop, mbott = imput.Naive (temp_id_seq,
                                             top,
                                             bott,
                                             ref_length,
                                             dyad_offset)        
            dyad, mtop, mbott = dyad[0], mtop[0], mbott[0]

            for i in range(len(sub_keys)):
                key = sub_keys[i]
                dyadmaps[key] = dyad[i]
                for j in range(len(mtop[i])):
                    if not np.isnan(mtop[i][j]):
                        cutmaps[key]['R'][j] = mtop[i][j]
                    if not np.isnan(mbott[i][j]):
                        cutmaps[key]['L'][j] = mbott[i][j]
        
    elif fill == "LINEAR":
        keys = sorted(cutmaps.keys(), cmp=size_loc_cmp)

        sub_keys = []
        temp_id_seq = {}
        top, bott = [], []
        for key in keys:
            try:
                loc, mtype, nts = key.split('-')
            except:
                continue
            sub_keys.append(key)
            temp_loc = int(loc) - cut_pad
            temp_nts = 'N'*(len(nts) + 2*cut_pad)
            #temp_nts = 'N'*(cut_pad) + nts + 'N'*(cut_pad)
            temp_id = temp_nts + '-' + str(temp_loc)
            if temp_id not in temp_id_seq:
                temp_id_seq[temp_id] = id_seq[key]
            top.append(cutmaps[key]['R'])
            bott.append(cutmaps[key]['L'])

        if len(sub_keys) <= 0:
            print >> sys.stderr, "Warning:imputation is not possible."
        else:
        
            dyad, mtop, mbott = imput.Linear_reg (temp_id_seq,
                                                  top,
                                                  bott,
                                                  ref_length,
                                                  dyad_offset)
            dyad, mtop, mbott = dyad[0], mtop[0], mbott[0]

            for i in range(len(sub_keys)):
                key = sub_keys[i]
                dyadmaps[key] = dyad[i]
                for j in range(len(mtop[i])):
                    if not np.isnan(mtop[i][j]):
                        cutmaps[key]['R'][j] = mtop[i][j]
                    if not np.isnan(mbott[i][j]):
                        cutmaps[key]['L'][j] = mbott[i][j]
        
    elif fill.startswith("BAYESIAN"):
        try:
            fill, cycle = fill.split(':')
            cycle = int(cycle)
        except:
            cycle = 1000
            
        group_keys = {}
        for key in cutmaps:
            try:
                loc, mtype, nts = key.split('-')
            except:
                continue
            group = mtype + '-' + 'N'*len(nts)
            if group not in group_keys:
                group_keys[group] = []
            group_keys[group].append(key)

        groups = sorted(group_keys.keys(), cmp=group_cmp)        
        for group in groups:
            group_keys[group] = sorted(group_keys[group], cmp=size_loc_cmp)

        if len(group_keys) <= 0:
            print >> sys.stderr, "Warning:imputation is not possible."
            
        else:
            arg_list = []
            left_bound, right_bound = NCP_len/2, NCP_len/2
            for group in groups:
                keys = group_keys[group]
                top = [cutmaps[key]['R'] for key in keys]
                bott = [cutmaps[key]['L'] for key in keys]
                temp_id_seq = {}
                for key in keys:
                    loc, mtype, nts = key.split('-')
                    temp_loc = int(loc) - cut_pad
                    temp_nts = 'N'*(len(nts) + 2*cut_pad)
                    #temp_nts = 'N'*(cut_pad) + nts + 'N'*(cut_pad)
                    temp_id = temp_nts + '-' + str(temp_loc)
                    temp_id_seq[temp_id] = id_seq[key]
                args = (cycle, temp_id_seq, top, bott, ref_length, dyad_offset, left_bound, right_bound)
                kwargs = {"dyad_alphas":None, "dyad_betas":None, "r_alphas":None, "r_betas":None, "l_sigmas":None, "initial":"linear", "note":"_bayesian_" + group}
                arg_list.append((args, kwargs))

            import multiprocessing as mp
            pool = mp.Pool()
            pool.map(Bayesian_wrapper, arg_list)
            
            for group in groups:
                dyad = pick_median(imput.load_data("dyad_bayesian_" + group + ".txt"))
                mtop = pick_median(imput.load_data("mtop_bayesian_" + group + ".txt"))
                mbott = pick_median(imput.load_data("mbott_bayesian_" + group + ".txt"))
                keys = group_keys[group]
                for i in range(len(keys)):
                    key = keys[i]
                    dyadmaps[key] = dyad[i]
                    for j in range(len(mtop[i])):
                        if not np.isnan(mtop[i][j]):
                            cutmaps[key]['R'][j] = mtop[i][j]
                        if not np.isnan(mbott[i][j]):
                            cutmaps[key]['L'][j] = mbott[i][j]


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
    dyad_axis = ref_length/2

    if extension == 'data':
        f = open(out_fname, 'w')
        print >> f, '@ref_length\t%s' % (ref_length)
        print >> f, '@dyad_axis\t%s' % (dyad_axis)
        print >> f, '@left_offset\t%s' % (dyad_offset)
        print >> f, '@right_offset\t%s' % (dyad_offset)
        print >> f, ''

    elif extension == 'pickle':
        import pickle
        key_slider = {}

    count = 0
    for mapped_id in sorted(cutmaps.keys(), cmp=size_loc_cmp):
        if extension == 'data':
            print >> f, '>%s' % (mapped_id)        

            if ref_fname:
                print >> f, '@Sequence'
                print >> f, id_seq[mapped_id]

            if shape_fname:
                print >> f, '@MinorGrooveWidth'
                print >> f, id_MGW[mapped_id]
                print >> f, '@HelixTwist'
                print >> f, id_HelT[mapped_id]
                print >> f, '@PropellerTwist'
                print >> f, id_ProT[mapped_id]
                print >> f, '@Roll'
                print >> f, id_Roll[mapped_id]    

            print >> f, '@TopCleavageCounts'
            print >> f, cutmaps[mapped_id]['R']
            print >> f, '@BottomCleavageCounts'
            print >> f, cutmaps[mapped_id]['L']
            print >> f, '@PositioningSignal'
            print >> f, dyadmaps[mapped_id]
            print >> f, ''

        elif extension == 'pickle':
            key_slider[mapped_id] = Slider(id = mapped_id,
                                           ref_length = ref_length,
                                           dyad_axis = dyad_axis,
                                           left_offset = dyad_offset,
                                           right_offset = dyad_offset,
                                           seq = id_seq[mapped_id],
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
                        help='output filename (.data or .pickle)')
    parser.add_argument('--tlen',
                        dest="ref_length",
                        type=int,
                        default=225,
                        help='reference DNA length (bp)')
    parser.add_argument('--nlen',
                        dest="NCP_len",
                        type=int,
                        default=147,
                        help='Nucleosomal DNA length (bp)')
    parser.add_argument('--offset',
                        dest="dyad_offset",
                        type=int,
                        default=52,
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
            out_fname = common_name + '.data'
        else:
            out_fname = 'out.data'
    else:
        out_fname = args.out_fname

    if args.fill:
        fill = args.fill.upper()
    else:
        fill = 'None'

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
