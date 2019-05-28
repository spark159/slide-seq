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

def all_path(N, states='AC'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

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
def load_files(filenames, ref_length, dyad_axis, dyad_offset, key_choice='mapped_id', choice=['valid'], filter_num=0, fill=False, deconvolute = False, load_ref=False, shape_fname=False):

    #weight_dict = {'R':{'A':1.308, 'C':0.879, 'T':0.642, 'G':1.152}, 'L':{'T':1.308, 'G':0.879, 'A':0.642, 'C':1.152}}
    #weight_dict = {'R':{'A':10, 'C':0.5, 'T':0.01, 'G':2}, 'L':{'T':10, 'G':0.5, 'A':0.01, 'C':2}}

    cutmaps, dyadmaps = {}, {}
    for i in range(len(filenames)):
        filename = filenames[i]
        for line in open(filename, 'r'):
            if line.strip():
                read_id, type, mapped_id, cutloc, seq = line.strip().split()
                if type not in choice:
                    continue            
                cols=cutloc.split(':')
                if int(cols[1]) < 0 or int(cols[1]) >= ref_length:
                    continue
                side, cut = cols[0], int(cols[1])
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
            #size, loc = mapped_id.split('-')
            win, loc = mapped_id.split('-')
            size, loc = len(win), int(loc)           
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

        # collect valid cleavage frequencies
        win_pos_TB = {}
        for mapped_id in cutmaps:
            win, loc = mapped_id.split('-')
            size, loc = len(win), int(loc)
            if size <= 0:
                continue
            if win not in win_pos_TB:
                win_pos_TB[win] = {}
            st, ed = loc, loc + size
            Rrange = range(max(ed - 2*dyad_offset, 0), min(ref_length - 2*dyad_offset, st))
            for i in Rrange:
                if i not in win_pos_TB[win]:
                    win_pos_TB[win][i] = {}
                    win_pos_TB[win][i]['T'] = []
                    win_pos_TB[win][i]['B'] = []
                win_pos_TB[win][i]['T'].append(cutmaps[mapped_id]['R'][i])
                win_pos_TB[win][i]['B'].append(cutmaps[mapped_id]['L'][i + 2*dyad_offset])


        # get linear model of T/B cleavage on each position
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
                #fig = plt.figure()
                #plt.plot(X, Y, '.k')
                #X = [[x] for x in X]
                #Z = model.predict(X)
                #X = [x[0] for x in X]
                #if rsquare < 0.03: # fitting failure
                    #continue

                #else:
                    #plt.plot(X, Z, 'k')
                coefs1, rsquare1 = robust_fit(X, Y)  # Top to bottom
                coefs2, rsquare2 = robust_fit(Y, X)  # Bottom to Top

                if rsquare1 >= 0 and rsquare2 >= 0:
                    if rsquare1 < rsquare2:
                        coefs1 = [1/coefs2[0], -1.0 * coefs2[1]/coefs2[0]]
                    if rsquare2 < rsquare1:
                        coefs2 = [1/coefs1[0], -1.0 * coefs1[1]/coefs1[0]]
                else:
                    if rsquare1 < 0:
                        coefs1 = [0, np.median(Y)]
                    if rsquare2 < 0:
                        coefs2 = [0, np.median(X)]
                    
                #bound = max(max(X),max(Y))
                #Z = [np.poly1d(coefs1)(x) for x in X]
                #plt.plot(X, Z, 'k')
                #Z = [np.poly1d(coefs2)(y) for y in Y]
                #plt.plot(Z, Y, 'r')
                #plt.xlim([-2, bound])
                #plt.ylim([-2, bound])
                #plt.show()
                #plt.close()
    
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
            #size, loc = mapped_id.split('-')
            win, loc = mapped_id.split('-')
            size, loc = len(win), int(loc)           
            if size <= 0:
                continue
            st, ed = loc, loc + size
            Lrange = range(min(st + 2*dyad_offset, ref_length), ref_length)
            Rrange = range(0, max(ed -1 - 2*dyad_offset, -1) + 1)
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

    """
    if sym:
        for mapped_id in cutmaps:
            cutmap = cutmaps[mapped_id]
            lcutmap, rcutmap = cutmap['L'], cutmap['R']
            for i in range(len(rcutmap) - 2*dyad_offset):
                maxsig = max(lcutmap[i + 2*dyad_offset], rcutmap[i])
                lcutmap[i + 2*dyad_offset], rcutmap[i] = maxsig, maxsig
            tempmap = [0.0]*ref_length
            for i in range(ref_length):
                lsig, rsig = lcutmap[i], rcutmap[i]
                if i+dyad_offset < ref_length:
                    tempmap[i+dyad_offset] += rsig
                if i-dyad_offset >= 0:
                    tempmap[i-dyad_offset] += lsig
            dyadmaps[mapped_id] = tempmap
    """
    
    assert len(cutmaps) == len(dyadmaps)
    all_len = len(all_path(16))
    print "id counts: " + str(len(cutmaps)) + '/' + str(all_len)

    
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

    return key_slider
    
