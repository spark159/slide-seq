import subprocess
import matplotlib.pyplot as plt

def read_sort (sort_fname, choice):
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

def find_crosstalk (fname, choices):
    id_xcount = {}
    total, count = 0, 0
    for line in open(fname):
        read_id, type, mapped_id, cut_loc, read_seq = line.strip().split()
        read_id = read_id.split(':')[-1]
        if type in choices:
            total += 1
            if read_id != mapped_id:
                count += 1
                if read_id not in id_xcount:
                    id_xcount[read_id] = 0
                id_xcount[read_id] += 1
    return count, total, id_xcount

for i in range(1):
    print i
    #fname = "polyABubble" + ".sort"
    fname = "singleInDel" + ".sort"
    count, total, id_xcount = find_crosstalk(fname, ["freeDNA", "valid"])
    xcount_ids = {}
    for id, xcount in id_xcount.items():
        if xcount not in xcount_ids:
            xcount_ids[xcount] = []
        xcount_ids[xcount].append(id)
    xcount_list = sorted(xcount_ids.keys(), reverse=True)
    for k in range(5):
        xcount = xcount_list[k]
        print xcount
        print xcount, xcount_ids[xcount]
    fig = plt.figure()
    plt.plot(range(len(xcount_list)), xcount_list)
    plt.xlabel("Seq ID")
    plt.ylabel("Cross-talk Counts")
    plt.savefig("test_" + str(i) + ".png")
    plt.close()
    print
    print count, total
    print float(count)*100/total
    print
