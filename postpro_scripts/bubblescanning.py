import subprocess
import random


#18bp
left =  "ATCCGACTGGCACCGGCAAGGTCGCTGTTCGCCACATGCGCAGGAT"
right = "TCTCCAGGGCGTCCTCGTATAGGGTCCATCACATAAGGGATGAACT"
ref = "GTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGAT"

##for synthesis reference (shorter left/right, total 178bp, No Backbone)
#left_s = "CTGTTCGCCACATGCGCAGGAT"
#right_s= "TCTCCAGGGCGTCCTCGTATAGG"
#ref = "GTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGAT"

#for synthesis reference (shorter left/right, total 169bp, No Backbone)
left_s = "TCGCCACATGCGCAGGAT"
right_s= "TCTCCAGGGCGTCCTCGT"
ref = "GTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGAT"


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

def all_bubbles (kmer, states='ATCG'):
    def all_path (states_list):
        N = len(states_list)
        if N==1:
            return list(states_list[0])
        output = []
        for path in all_path(states_list[:-1]):
            for state in states_list[-1]:
                output.append(path+state)
        return output
    options_list = []
    for i in range(len(kmer)):
        nt = kmer[i]
        options = list(set(list(states)) - set(nt))
        if not options:
            return []
        options_list.append(options)
    output = all_path(options_list)
    # sanity check
    for word in output:
        assert len(word) == len(kmer)
        for i in range(len(word)):
            assert word[i] != kmer[i]
    return output

def min_bubbles (kmer, states='AC'):
    output = ""
    for i in range(len(kmer)):
        nt = kmer[i]
        alt = None
        for j in range(len(states)):
            st = states[j]
            if st != nt:
                alt = st
                break
        if not alt:
            return []
        output += alt
    return [output]
        
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

def mismatch_library (left, right, ref, st_num, ed_num, step_size=1, states='ATCG', ID_offset=len(left)):
    seq_list = []
    id_list = []
    for k in range(st_num, ed_num+1):
        for i in range(0, len(ref)-k+1, step_size):
            nts = ref[i:i+k]
            #options = all_bubbles(nts, states=states)
            options = min_bubbles(nts, states=states)
            #if k == 1:
            #    options = reduce1(options)
            #else:
            #    options = reduce1(reduce2(options))
            for kmer in options:
                seq = left + ref[:i] + kmer + ref[i+k:] + right
                #seq = left_s + ref[:i] + kmer + ref[i+k:] + right_s
                #assert len(seq) == len(left + ref + right) == 225
                #assert len(seq) == len(left_s + ref + right_s) == 178
                #ID = kmer + '-' + str(len(left) + i)
                #ID = str(len(left) + i) +'-' + 'M-' + kmer
                ID = str(ID_offset + i) +'-' + 'M-' + kmer
                seq_list.append(seq)
                id_list.append(ID)
    return seq_list, id_list

def deletion_library (left, right, ref, st_num, ed_num, step_size=1, states='ATCG', ID_offset=len(left)):
    seq_list = []
    id_list = []
    for k in range(st_num, ed_num+1):
        for i in range(1, len(ref)-k, step_size):
            if ref[i] not in states:
                continue
            if ref[i] == ref[i-1] or ref[i] == ref[i+1]: # no confusion (1bp case)
                continue
            seq = left + ref[:i] + ref[i+k:] + right
            #ID = str(len(left) + i) +'-' + 'D-' + ref[i:i+k]
            ID = str(ID_offset + i) +'-' + 'D-' + ref[i:i+k]
            seq_list.append(seq)
            id_list.append(ID)
    return seq_list, id_list

def insertion_library (left, right, ref, st_num, ed_num, step_size=1, states='ATCG', ID_offset=len(left)):
    seq_list = []
    id_list = []
    for k in range(st_num, ed_num+1):
        kmers = all_path(k, states=states)
        #if k == 1:
        #    kmers = reduce1(kmers)
        #else:
        #    kmers = reduce1(reduce2(kmers))
        for i in range(0, len(ref)-1, step_size):
            for kmer in kmers:
                if kmer == ref[i] or kmer == ref[i+1]: # no confusion (1bp case)
                    continue
                seq = left + ref[:i+1] + kmer + ref[i+1:] + right
                #ID = str(len(left) + i) +'-' + 'I-' + kmer
                ID = str(ID_offset + i) +'-' + 'I-' + kmer
                seq_list.append(seq)
                id_list.append(ID)
    return seq_list, id_list

#seq_list1, id_list1 = mismatch_library(1,1, step_size=1, states='A')
#seq_list2, id_list2 = mismatch_library(2,2, step_size=1, states='A')
#seq_list3, id_list3 = mismatch_library(3,3, step_size=1, states='A')
#seq_list4, id_list4 = mismatch_library(4,4, step_size=1, states='A')
#seq_list5, id_list5 = mismatch_library(5,5, step_size=1, states='A')
#seq_list6, id_list6 = insertion_library(1,1, states='A')
#seq_list7, id_list7 = deletion_library(1,1, states='A')
#seq_list8 = seq_list1 + seq_list2 + seq_list3 + seq_list4 + seq_list5 + seq_list6 + seq_list7
#id_list8 = id_list1 + id_list2 + id_list3 + id_list4 + id_list5 + id_list6 + id_list7
#seq_list_list = [seq_list1, seq_list2, seq_list3, seq_list4, seq_list5, seq_list6, seq_list7, seq_list8]
#id_list_list = [id_list1, id_list2, id_list3, id_list4, id_list5, id_list6, id_list7, id_list8]
#seq_list_list = [seq_list6, seq_list7]
#id_list_list = [id_list6, id_list7]

# polyA mismatch library
seq_list, id_list = mismatch_library(left_s, right_s, ref, 1, 7, step_size=1, states='A')
print len(seq_list)
write_FASTA ("polyAMismatch_NoBB" + ".syn", seq_list, id_list)

seq_list, id_list = mismatch_library(left, right, ref, 1, 7, step_size=1, states='A')
seq_list.append(left+ref+right)
id_list.append("BACKBONE")
write_FASTA ("polyAMismatch" + ".ref", seq_list, id_list)

# single-base indel library
seq_list1, id_list1 = insertion_library(left_s, right_s, ref, 1,1, states='A')
seq_list2, id_list2 = deletion_library(left_s, right_s, ref, 1,1, states='ATCG')
seq_list = seq_list1 + seq_list2
id_list = id_list1 + id_list2
print len(seq_list)
write_FASTA ("singleInDel_NoBB" + ".syn", seq_list, id_list)

seq_list1, id_list1 = insertion_library(left, right, ref, 1,1, states='A')
seq_list2, id_list2 = deletion_library(left, right, ref, 1,1, states='ATCG')
seq_list = seq_list1 + seq_list2
id_list = id_list1 + id_list2
seq_list.append(left+ref+right)
id_list.append("BACKBONE")
write_FASTA ("singleInDel" + ".ref", seq_list, id_list)
