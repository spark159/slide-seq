import numpy as np
import matplotlib.pyplot as plt

def read_ref (ref_fname):
    id_seq = {}
    for line in open(ref_fname):
        line = line.strip()
        if line.startswith(">"):
            id = line[1:]
            continue
        seq = line
        assert id not in id_seq
        id_seq[id] = seq
    return id_seq
id_seq = read_ref ("/home/spark159/../../media/spark159/sw/plusonelibFinal/plusonelib.ref")

def Amer_len(seq, pos=True):
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

size_count = {}
for id, seq in id_seq.items():
    size_pos = Amer_len(seq, pos=True)
    for size in size_pos:
        if size < 5:
            continue
        if size not in size_count:
            size_count[size] = 0
        size_count[size] += len(size_pos[size])
    
X, Y = [], []
for size in sorted(size_count.keys()):
    X.append(size)
    Y.append(size_count[size])

fig = plt.figure()
plt.plot(X, Y, '.')
plt.show()
plt.close()
