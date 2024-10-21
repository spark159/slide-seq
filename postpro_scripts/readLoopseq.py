def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

def din_count(seq, din):
    count = 0
    for i in range(len(seq)-1):
        if seq[i:i+2] == din:
            count +=1
    return count

def Amer_score(seq, weight=1):
    def Amer_len(seq):
        num = []
        i = 0
        while i < len(seq):
            if seq[i] == 'A' or 'T':
                nt = seq[i]
                count = 1
                j = i + 1
                while j < len(seq):
                    if seq[j] != nt:
                        break
                    count +=1
                    j +=1
                num.append(count)
                i = j + 1
            else:
                i +=1
        return num
    num = Amer_len(seq)
    score = 0.0
    for i in range(len(num)):
        if num[i] < 4:
            continue
        score += (num[i]-1) * weight
    return score

def read_plusone(fname):
    id_seq = {}
    for line in open(fname):
        line = line.strip()
        if line.startswith('>'):
            id = line[1:]
            continue
        assert id not in id_seq
        id_seq[id] = line
    return id_seq

def read_loopseq(fname, bound=0):
    seqs = []
    values = []
    for line in open(fname):
        seq, score = line.strip().split()
        seq = seq[bound:len(seq)-bound]
        score = float(score)
        seqs.append(seq)
        values.append(score)
    return seqs, values

def get_flex(ref_fname, loopseq_fname):
    id_seq = read_plusone(ref_fname)
    seqs, values = read_loopseq(loopseq_fname)
    
    id_values = {}
    for i in range(len(seqs)):
        seq = seqs[i]
        seq = seq[2:len(seq)-2]
        for id in id_seq:
            fullseq = id_seq[id]
            axis = len(fullseq) / 2
            left = fullseq[axis-50:axis]
            right = fullseq[axis+1:axis+50+1]
            if seq == left:
                if id not in id_values:
                    id_values[id] = {}
                try:
                    if id_values[id]['left']:
                        assert False
                except:
                    pass
                id_values[id]['left'] = values[i]
            elif seq == right:
                if id not in id_values:
                    id_values[id] = {}
                try:
                    if id_values[id]['right']:
                        assert False
                except:
                    pass
                id_values[id]['right'] = values[i]

    return id_values

def get_polyA (ref_fname):
    id_seq = read_plusone(ref_fname)
    id_dAlen = {}
    for id in id_seq:
        seq = id_seq[id]
        seq = seq[40:len(seq)-40]
        axis = len(seq) / 2
        left, right = seq[axis-72:axis], seq[axis+1:axis+73]
        dAlen = Amer_score(right) - Amer_score(left)
        id_dAlen[id] = dAlen
    return id_dAlen

def get_GC (ref_fname):
    id_seq = read_plusone(ref_fname)
    id_dGC = {}
    for id in id_seq:
        seq = id_seq[id]
        seq = seq[40:len(seq)-40]
        axis = len(seq) / 2
        left, right = seq[axis-72:axis], seq[axis+1:axis+73]
        #if id == '0':
        #    print left
        #    print right
        dGC = GC_content(right) - GC_content(left)
        id_dGC[id] = dGC
    return id_dGC

def get_Din (ref_fname, din):
    id_seq = read_plusone(ref_fname)
    id_dDin = {}
    for id in id_seq:
        seq = id_seq[id]
        seq = seq[40:len(seq)-40]
        axis = len(seq) / 2
        left, right = seq[axis-72:axis], seq[axis+1:axis+73]
        dDin = din_count(right, din) - din_count(left, din)
        id_dDin[id] = dDin
    return id_dDin

    

    
