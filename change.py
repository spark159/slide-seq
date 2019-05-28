import math
import numpy as np
import matplotlib.pyplot as plt

def read_data (fname):
    id_count, id_control = [], []
    for line in open(fname):
        cols = line.strip().split(',')
        count, control = float(cols[0]), float(cols[1])
        if count <= 0 :
            count = 1.0
        if control <= 0:
            control = 1.0
        id_count.append(count)
        id_control.append(control)
    assert len(id_count) == len(id_control)
    id_count = [ value/sum(id_count) for value in id_count ]
    id_control = [ value/sum(id_control) for value in id_control ]
    #id_score = [ math.log(id_count[i]/id_control[i]) for i in range(len(id_count)) ]
    id_score = [ id_count[i]/id_control[i] for i in range(len(id_count)) ]
    return id_score

def read_ref(fname):
    id_seq = []
    for line in open(fname):
        cols = line.strip().split()
        seq, _ = cols
        id_seq.append(seq)
    return id_seq

def write(id_score, id_seq):
    assert len(id_score) == len(id_seq)
    f = open("data/DataReady2.txt", 'w')
    for i in range(len(id_score)):
        print >> f, "%s\t%f" % (id_seq[i], id_score[i])
    f.close()
    return None

def GC_content(seq):
    num=0.0
    for nt in seq:
        if nt in 'GC':
            num+=1
    return (num/float(len(seq)))*100

def Amer_len(seq, pos=False):
    num = []
    num_pos = {}
    i = 0
    while i < len(seq):
        if seq[i] in 'AT':
            nt = seq[i]
            count = 1
            j = i + 1
            while j < len(seq):
                if seq[j] != nt:
                    break
                count +=1
                j +=1
            num.append(count)
            if count not in num_pos:
                num_pos[count] = []
            num_pos[count].append(i)
            i = j
        else:
            i +=1
    if pos:
        return num_pos
    if len(num) == 0:
        return 0
    return max(num)

id_score = read_data("data/loopseqdata.csv")
id_seq = read_ref("data/DataReady.txt")
write(id_score, id_seq)

id_GC = [ GC_content(seq) for seq in id_seq ]
id_Amer = [ Amer_len(seq) for seq in id_seq ]

"""
fig = plt.figure()
plt.plot(id_GC, id_score, '.', alpha=0.3)
plt.show()
plt.close()
"""


