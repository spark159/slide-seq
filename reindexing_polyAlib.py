import copy
import pickle

def read_polyAlib (fname):
    id_seq = {}
    for line in open(fname):
        line = line.strip()
        if line.startswith('>'):
            id = line[1:]
            continue
        seq = line.strip()
        assert id not in id_seq
        id_seq[id] = seq
    return id_seq

path = '/home/spark159/../../media/spark159/sw/polyAlibFinal/'
id_seq = read_polyAlib(path+'polyAscanlib.ref')
bb_seq = id_seq['BACKBONE']
del id_seq['BACKBONE'] 

# reindexing by checking neighbors 
id_newid = {}
newid_id = {}
for id, seq in id_seq.items():
    pos, mtype, tract = id.split('-')
    pos = int(pos)

    st, ed = pos, pos + len(tract)

    new_st = copy.deepcopy(pos)
    while new_st >= 0:
        nt = seq[new_st]
        if nt != 'A':
            new_st +=1
            break
        new_st -=1

    if new_st < 0:
        new_st = 0

    new_ed = copy.deepcopy(ed)
    while new_ed <= len(seq):
        nt = seq[new_ed]
        if nt != 'A':
            break
        new_ed +=1
        
    if new_ed > len(seq):
        new_ed = len(seq)
        
    newid = '-'.join([str(new_st), mtype, 'A'*(new_ed-new_st)])

    #if id != newid:
        #print id, newid
        #print len(tract), new_ed-new_st

    id_newid[id] = newid

    assert newid not in newid_id # check 1-to-1 mapping 
    newid_id[newid] = id


# save reindexing information
pickle.dump(id_newid, open("polyAlib_id_newid.p", "wb"))
pickle.dump(newid_id, open("polyAlib_newid_id.p", "wb"))

# write new ref
def id_cmp (a,b):
    pos1, mtype1, tract1 = a.split('-')
    pos2, mtype2, tract2 = b.split('-')
    pos1, pos2 = int(pos1), int(pos2)
    if len(tract1) < len(tract2):
        return -1
    elif len(tract1) > len(tract2):
        return 1
    else:
        if pos1 < pos2:
            return -1
        else:
            return 1

newids = sorted(newid_id.keys(), cmp=id_cmp)

f = open('polyAscanlib_reindexed.ref', 'w')
for newid in newids:
    seq = id_seq[newid_id[newid]]
    print >> f, '>%s' % (newid)
    print >> f, seq
print >> f, '>BACKBONE'
print >> f, bb_seq
f.close()

    
