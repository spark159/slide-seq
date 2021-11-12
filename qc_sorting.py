import os, sys, subprocess, re
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

def read_simulated_sort (fname, outfname):
    id_info_list = {}
    f = open(outfname, 'w')
    for line in open(fname):
        read_id, mapped_type, mapped_id, mapped_cutloc, read_seq = line.strip().split()
        read_id, cutloc = read_id.split('/')[1:]
        if read_id not in id_info_list:
            id_info_list[read_id] = []
        info = {}
        info['cutloc'] = cutloc
        info['mapped_id'] = mapped_id
        info['mapped_cutloc'] = mapped_cutloc
        info['mapped_type'] = mapped_type
        #id_info['read_seq'] = read_seq
        id_info_list[read_id].append(info)

        # record failed mapping
        if mapped_type == 'valid' and (read_id != mapped_id or cutloc != mapped_cutloc):
            print >> f, line.strip()

    f.close()

    return id_info_list



#fname_list = ['polyAscanlib_reindexed']
#fname_list = ['singleInDel']
fname_list = ['polyAMismatch']
for fname in fname_list:
    # make simulated fastq file
    subprocess.call(['python', 'simulator_edit.py', fname, '--left_prob', str(0.8), '--right_prob', str(0.8), '--mist_prob', str(0.02), '--indel_prob', str(0.002)])
    #sys.exit(1)
    # make .sort file
    subprocess.call(['python', 'sort_final.py', fname+'_simulated.fastq', fname, '-a'])

    # read sort file and check failed ones
    type_ids, fail_ids = {}, {}
    f = open(fname + '_simulated_fail.sort', 'w')
    for line in open(fname+'_simulated.sort'):
        read_id, type, mapped_id, mapped_cutloc, read_seq = line.strip().split()
        ref_id, cutloc = read_id.split('/')[1:]

        if type not in type_ids:
            type_ids[type] = []
        type_ids[type].append(read_id)

        if type == 'valid':
            fail = []
            if ref_id != mapped_id:
                fail.append('id')
            if cutloc != mapped_cutloc:
                fail.append('cutloc')

            if len(fail) > 0:
                fail = 'miss ' + '&'.join(fail)
            else:
                fail = 'success'

            if fail not in fail_ids:
                fail_ids[fail] = []
            fail_ids[fail].append(id)

        if type == 'freeDNA' and ref_id != mapped_id:
            print >> f, line
        if type == 'valid' and (ref_id != mapped_id or cutloc != mapped_cutloc):
            print >> f, line
            
    f.close()

    # plot pie chart for types
    labels, counts = [], []
    for type in type_ids:
        print type, len(type_ids[type])
        labels.append(type)
        counts.append(len(type_ids[type]))
    print
        
    fig = plt.figure()
    plt.pie(counts, labels=labels, shadow=True, startangle=90, autopct='%1.1f%%')
    plt.axis('equal')
    plt.title("Reads types", pad=20)
    plt.show()
    plt.close()

    # plot pie chart for fails
    labels, counts = [], []
    for fail in fail_ids:
        print fail, len(fail_ids[fail])
        labels.append(fail)
        counts.append(len(fail_ids[fail]))
    print
        
    fig = plt.figure()
    plt.pie(counts, labels=labels, shadow=True, startangle=90, autopct='%1.1f%%')
    plt.axis('equal')
    plt.title("Within valid reads", pad=20)
    plt.show()
    plt.close()

    






        
"""
    # read .sort file
    id_info_list = read_simulated_sort(fname+'_simulated.sort', fname + '_simulated_fail.sort')

    # make statistics
    id_stat = {}
    for id, info_list in id_info_list.items():
        if id == 'BACKBONE':
            continue
        if id not in id_stat:
            id_stat[id] = {}
        for info in info_list:
            mapped_type = info['mapped_type']
            mapped_id = info['mapped_id']
            mapped_cutloc = info['mapped_cutloc']
            if mapped_type != 'valid':
                cate = 'invalid'
                if cate not in id_stat[id]:
                    id_stat[id][cate] = {}
                if mapped_type not in id_stat[id][cate]:
                    id_stat[id][cate][mapped_type] = 0
                id_stat[id][cate][mapped_type] +=1
            else:
                cate = 'valid'
                if cate not in id_stat[id]:
                    id_stat[id][cate] = {}
                fails = []
                if mapped_id != id:
                    fails.append('id')
                if mapped_cutloc != info['cutloc']:
                    fails.append('loc')
                if len(fails) > 0:
                    subcate = 'fail:' + '-'.join(fails)
                else:
                    subcate = 'success'
                if subcate not in id_stat[id][cate]:
                    id_stat[id][cate][subcate] = 0
                id_stat[id][cate][subcate] +=1

    # make graph
    fig = plt.figure()
    ids = id_stat.keys()
    for i in range(len(ids)):
        id = ids[i]
        if id == 'BACKBONE':
            continue
        stat = id_stat[id]
        data_list = [stat['valid']['success'], sum(stat['invalid'].values())]
        color_list = ['red', 'blue']
        bottom = 0
        for data, color in zip(data_list, color_list):
            plt.bar([i], [data], color = color, bottom=bottom)
            bottom += data
    plt.show()
    plt.close()

    

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
"""
