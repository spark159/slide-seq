#12bp
#left = "ATCCGACTGGCACCGGCAAGGTCGCTGTTCGCCACATGCG"
#right = "GGGCGTCCTCGTATAGGGTCCATCACATAAGGGATGAACT"
#ref = "CAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCA"

#18bp
#left = "ATCCGACTGGCACCGGCAAGGTCGCTGTTCGCCACATGCGCAGGAT"
#right = "TCTCCAGGGCGTCCTCGTATAGGGTCCATCACATAAGGGATGAACT"
left = "TCGCCACATGCGCAGGAT"
right = "TCTCCAGGGCGTCCTCGT"
ref = "GTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGAT"

#18bp
left = "ATCCGACTGGCACCGGCAAGGTCGCTGTTCGCCACATGCGCAGGAT"
right = "TCTCCAGGGCGTCCTCGTATAGGGTCCATCACATAAGGGATGAACT"
#left = "TCGCCACATGCGCAGGAT"
#right = "TCTCCAGGGCGTCCTCGT"
ref = "GTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGAT"
print len(ref)

def rev_comp (seq):
    rev_seq=seq[::-1]; new_seq=''
    for base in rev_seq:
        if base == 'A':
            new_seq += 'T'
        if base == 'T':
            new_seq += 'A'
        if base == 'C':
            new_seq += 'G'
        if base == 'G':
            new_seq += 'C'
    return new_seq

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

def reduce1 (kmer_list):
    output = []
    for kmer in kmer_list:
        if kmer not in output and rev_comp(kmer) not in output:
            output.append(kmer)
    return output

def reduce2 (kmer_list):
    dinu = {"AA":1, "TT":1, "AC":2, "GT":2, "AG":3, "CT":3, "AT":4, "CA":5, "TG":5, "CC":6, "GG":6, "CG":7, "GA":8, "TC":8, "GC":9, "TA":10}
    id_list, output = [], []
    for kmer in kmer_list:
        id = ""
        for i in range(len(kmer)-1):
            id += str(dinu[kmer[i:i+2]])
        if id not in id_list:
            id_list.append(id)
            output.append(kmer)
    return output

print all_path(3)
print len(all_path(3))
insert_list = reduce2(all_path(3))
print insert_list
print len(insert_list)

print range(0, len(ref)-3+1, 5)

def write_FASTA (filename, seq_list, id_list=None):
    if id_list:
        assert len(id_list) == len(seq_list)
    f = open(filename + '', 'w')
    for i in range(len(seq_list)):
        if id_list:
            id = id_list[i]
        else:
            id = str(i)
        print >> f, '>%s' % (id)
        print >> f, seq_list[i]
    f.close()

seq_list = []
id_list =[]
for i in range(len(insert_list)):
    insert = insert_list[i]
    for j in range(0, len(ref)-len(insert)+1, 5):
        seq = left + ref[0:j] + insert + ref[j+len(insert):len(ref)] + right
        #assert len(seq) == 169
        assert len(seq) == 225
        if seq not in seq_list:
            seq_list.append(seq)
            id = insert + '-' + str(j+len(left))
            #id = str(i) + '-' + str(j+len(left))
            id_list.append(id)

print len(seq_list)
write_FASTA("3merscanlib.ref", seq_list, id_list)
