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

def make_revcomp (fname):
    f = open("rev_" + fname, 'w')
    for line in open(fname):
        line = line.strip()
        if line.startswith('>'):
            id = line[1:]
            print >> f, '>%s' % id
            continue
        else:
            seq = line
            rev_seq = rev_comp(seq)
            print >> f, rev_seq
    f.close()

make_revcomp ("plusonelib.ref")
            
