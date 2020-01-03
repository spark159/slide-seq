import sys
from SliderClass import Slider
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model

def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for nt in rev_seq:
        nt = nt.upper()
        if nt == 'A':
            new_seq += 'T'
        if nt == 'T':
            new_seq += 'A'
        if nt == 'C':
            new_seq += 'G'
        if nt == 'G':
            new_seq += 'C'
    return new_seq

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

# read sort files and return Slider objects
def load_files(filenames,
               ref_length,
               dyad_axis,
               dyad_offset,
               key_choice='mapped_id',
               choice=['valid'],
               mtype_choice=None,
               filter_num=0,
               fill=False,
               cut_pad=0,
               load_ref=False,
               shape_fname=False):

    # collect the cleavage counts from sort file
    cutmaps, dyadmaps = {}, {}
    for i in range(len(filenames)):
        filename = filenames[i]
        for line in open(filename, 'r'):
            if line.strip():
                read_id, type, mapped_id, cutloc, seq = line.strip().split()
                if type not in choice:
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

                    # temporal delesion case
                    if mtype == 'D':
                        nts = len(nts)*"N"
                        mapped_id = '-'.join([str(loc), mtype, nts])

                    # temporal M 1bp case
                    if mtype in 'MID' and len(nts) == 1:
                        cut_pad = 3
                        if side == 'R' and loc - cut <= cut_pad:
                            continue
                        if side == 'L' and cut - loc <= cut_pad:
                            continue
                except:
                    pass
                    
                if side == 'L':
                    offset = -dyad_offset
                    #nt = id_seq[id][cut+1]
                    #weight = weight_dict[side][nt]
                else:
                    assert side == 'R'
                    offset = dyad_offset
                    #nt = id_seq[id][cut-1]
                    #weight = weight_dict[side][nt]
                if cut+offset < 0  or cut+offset >= ref_length:
                    continue
                if mapped_id not in dyadmaps:
                    dyadmaps[mapped_id] = [0.0]*ref_length
                dyadmaps[mapped_id][cut+offset] += 1.0 #/ weight
                if mapped_id not in cutmaps:
                    cutmaps[mapped_id] = {}
                    cutmaps[mapped_id]['L'], cutmaps[mapped_id]['R'] = [0.0]*ref_length, [0.0]*ref_length
                cutmaps[mapped_id][side][cut] += 1

    # simple copy and paste correction
    if fill == "simple":
        for mapped_id in cutmaps:
            #win, loc = mapped_id.split('-')
            #size, loc = len(win), int(loc)           
            try:
                loc, mtype, nts = mapped_id.split('-')
            except:
                print >> sys.stderr, "Error: there is no imputation position guide."
                sys.exit(1)
            size, loc = len(nts), int(loc)
            if size <= 0:
                continue
            st, ed = loc, loc + size
            Lrange = range(min(st + 2*dyad_offset, ref_length), ref_length)
            Rrange = range(0, max(ed -1 - 2*dyad_offset, -1) + 1)
            for i in Lrange:
                counts = cutmaps[mapped_id]['L'][i]
                cutmaps[mapped_id]['R'][i - 2*dyad_offset] += counts
                dyadmaps[mapped_id][i - dyad_offset] += counts
            for i in Rrange:
                counts = cutmaps[mapped_id]['R'][i]
                cutmaps[mapped_id]['L'][i + 2*dyad_offset] += counts
                dyadmaps[mapped_id][i + dyad_offset] += counts

    # correction by linear regression 
    if fill == "linear":
        def robust_fit(X, Y):
            mX = [[x] for x in X]
            mX = np.asarray(mX)
            Y = np.asarray(Y)
            model = linear_model.HuberRegressor()
            model.fit(mX, Y)
            inmask = np.array([b == False for b in model.outliers_])
            mX = mX[inmask]
            Y = Y[inmask]
            rsquare = model.score(mX, Y)
            coefs = [model.coef_[0], model.intercept_]
            return coefs, rsquare

        def zero_frac (L):
            num = 0.0
            for e in L:
                if e == 0:
                    num +=1
            return num/len(L)

        # collect available T/B data pair
        win_pos_TB = {}
        for mapped_id in cutmaps:
            #win, loc = mapped_id.split('-')
            #size, loc = len(win), int(loc)
            try:
                loc, mtype, nts = mapped_id.split('-')
            except:
                print >> sys.stderr, "Error: there is no imputation position guide."
                sys.exit(1)
            win = len(nts)*'N'
            #win = nts
            size, loc = len(nts), int(loc)
            if size <= 0:
                continue
            if win not in win_pos_TB:
                win_pos_TB[win] = {}
            st, ed = loc - cut_pad, loc + size + cut_pad
            Rrange = range(max(ed - 2*dyad_offset, 0), min(ref_length - 2*dyad_offset, st))
            for i in Rrange:
                if i not in win_pos_TB[win]:
                    win_pos_TB[win][i] = {}
                    win_pos_TB[win][i]['T'] = []
                    win_pos_TB[win][i]['B'] = []
                win_pos_TB[win][i]['T'].append(cutmaps[mapped_id]['R'][i])
                win_pos_TB[win][i]['B'].append(cutmaps[mapped_id]['L'][i + 2*dyad_offset])


        # build linear model of T/B cleavage on each position
        win_pos_model = {}
        for win in win_pos_TB:
            for pos in win_pos_TB[win]:
                X, Y = win_pos_TB[win][pos]['T'], win_pos_TB[win][pos]['B']
                if len(X) < 3: # too small data set
                    continue
                if zero_frac(X) > 0.7: # X domain is empty
                    continue
                if zero_frac(Y) > 0.7: # Y domain is empty
                    continue

                try:
                    coefs1, rsquare1 = robust_fit(X, Y)  # Top to bottom
                except:
                    coefs1, rsquare1 = None, -sys.maxint # fitting error
                    
                try:
                    coefs2, rsquare2 = robust_fit(Y, X)  # Bottom to Top
                except:
                    coefs2, rsquare2 = None, -sys.maxint # fitting error

                # pick the best fitting
                if rsquare1 >= 0 and rsquare2 >= 0:
                    if rsquare1 < rsquare2:
                        coefs1 = [1/coefs2[0], -1.0 * coefs2[1]/coefs2[0]]
                    if rsquare2 < rsquare1:
                        coefs2 = [1/coefs1[0], -1.0 * coefs1[1]/coefs1[0]]
                else:
                    # fitting failures
                    if rsquare1 < 0:
                        coefs1 = [0, np.median(Y)]
                    if rsquare2 < 0:
                        coefs2 = [0, np.median(X)]

                # too dangerous fitting
                if abs(coefs1[0]) > 200:
                    #print pos
                    #print coefs1[0]
                    coefs1 = [0, np.median(Y)]
                if abs(coefs2[0]) > 200:
                    #print pos
                    #print coefs2[0]
                    coefs2 = [0, np.median(X)]
                    
                if False:
                    fig = plt.figure()
                    plt.plot(X, Y, '.k')
                    bound = max(max(X),max(Y))
                    Z = [np.poly1d(coefs1)(x) for x in X]
                    plt.plot(X, Z, 'b', alpha=0.5, label="TopToBott")
                    print pos
                    print "TopToBott", coefs1[0]
                    Z = [np.poly1d(coefs2)(y) for y in Y]
                    plt.plot(Z, Y, 'r', alpha=0.5, label="BottToTop")
                    print "BottToTop", coefs2[0]
                    print 
                    plt.xlim([-2, bound])
                    plt.ylim([-2, bound])
                    plt.title(str(pos))
                    plt.xlabel("Top strand")
                    plt.ylabel("Bott strand")
                    plt.legend()
                    #plt.show()
                    plt.close()
    
                if win not in win_pos_model:
                    win_pos_model[win] = {}
                if pos not in win_pos_model[win]:
                    win_pos_model[win][pos] = {}
                win_pos_model[win][pos]['T'] = np.poly1d(coefs1)
                pos = pos + 2*dyad_offset
                if pos not in win_pos_model[win]:
                    win_pos_model[win][pos] = {}
                win_pos_model[win][pos]['B'] = np.poly1d(coefs2)
                                
        # fill the lost data based on linear model
        for mapped_id in cutmaps:
            #win, loc = mapped_id.split('-')
            #size, loc = len(win), int(loc)           
            try:
                loc, mtype, nts = mapped_id.split('-')
            except:
                print >> sys.stderr, "Error: there is no imputation guide."
                sys.exit(1)
            win = len(nts)*'N'
            #win = nts
            size, loc = len(nts), int(loc)
            if size <= 0:
                continue
            st, ed = loc, loc + size
            Lrange = range(min(st-cut_pad + 2*dyad_offset, ref_length), ref_length)
            Rrange = range(0, max(ed+cut_pad -1 - 2*dyad_offset, -1) + 1)
            for i in Lrange:
                try:
                    model = win_pos_model[win][i]['B']
                    x = cutmaps[mapped_id]['L'][i]
                    y = model(x)
                except:
                    continue
                if np.isnan(y) or y < 0:
                    continue
                cutmaps[mapped_id]['R'][i - 2*dyad_offset] += y
                dyadmaps[mapped_id][i - dyad_offset] += y
            for i in Rrange:
                try:
                    model = win_pos_model[win][i]['T']
                    x = cutmaps[mapped_id]['R'][i]
                    y = model(x)
                except:
                    continue
                if np.isnan(y) or y < 0:
                    continue
                cutmaps[mapped_id]['L'][i + 2*dyad_offset] += y
                dyadmaps[mapped_id][i + dyad_offset] += y
    
    if filter_num:
        filtered = []
        for mapped_id in cutmaps:
            if sum(dyadmaps[mapped_id]) > filter_num:
                filtered.append(mapped_id)
        cutmaps = {k:cutmaps[k] for k in filtered}
        dyadmaps = {k:dyadmaps[k] for k in filtered}        
    
    assert len(cutmaps) == len(dyadmaps)
    
    id_seq = {}
    if load_ref:
        for line in open(load_ref):
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

    if shape_fname:
        id_MGW, id_HelT, id_ProT, id_Roll = read_DNAshape(shape_fname)
        #print id_ProT
    else:
        id_MGW, id_HelT, id_ProT, id_Roll = {}, {}, {}, {}
        for id in dyadmaps.keys():
            id_MGW[id] = None; id_HelT[id] = None; id_ProT[id] = None; id_Roll[id] = None
        
    key_slider = {}
    i = 0
    for mapped_id in cutmaps:
        if key_choice == 'mapped_id':
            key = mapped_id
        elif key_choice == 'num':
            key = i
        key_slider[key] = Slider(id = mapped_id,
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
        i +=1

    if load_ref:
        report = str(i) + " / " + str(len(id_seq))
    else:
        report = str(i)
    print "data size: " + report 

    return key_slider
    
