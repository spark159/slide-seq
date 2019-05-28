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

seq_list = [left+ref+right]
id_list =['0-0']
for i in range(3,26):
    for j in range(len(ref)-i+1):
        seq = left + ref[0:j] + 'A'*i + ref[j+i:len(ref)] + right
        assert len(seq) == 169
        if seq not in seq_list:
            seq_list.append(seq)
            id = str(i) + '-' + str(j+len(left))
            id_list.append(id)

print len(seq_list)
write_FASTA("polyAscanlib.ref", seq_list, id_list)
