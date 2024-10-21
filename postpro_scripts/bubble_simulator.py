import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import random
import numpy as np

def rev_comp(seq):
    rev = ""
    for nt in seq[::-1]:
        nt = nt.upper()
        if nt == 'A':
            rev += 'T'
        elif nt == 'T':
            rev += 'A'
        elif nt == 'G':
            rev += 'C'
        elif nt == 'C':
            rev += 'G'
    return rev

def simulation (ref_seq,
                out_fname,
                num,
                left_bound,
                right_bound,
                left_prob,
                right_prob,
                mist_prob,
                indel_prob):

    def dyad_sampling(template, left_bound, right_bound, N):
        dyad_list = []
        st, ed = left_bound, len(template) - 1 - right_bound
        for i in range(N):
            dyad_pos = random.randint(st,ed)
            dyad_list.append(dyad_pos)
        return dyad_list

    def make_cleavage(template, dyad_list, left_prob, right_prob, left_offset=52, right_offset=52):
        #left_prob = [ np.exp(-0.5*((i-110)**2)/((20)**2)) for i in range(len(template))]
        #right_prob = [ np.exp(-0.5*((i-(len(template)-1-110))**2)/((20)**2)) for i in range(len(template))]
        #left_prob = [left_prob] * len(template)
        #right_prob = [right_prob] * len(template)
        left_prob = [ (0.2*10**-13)*((x)**(10.0-1))*np.exp(-(x)/10.0) for x in range(len(template))]
        right_prob = [ (0.2*10**-13)*((-x+225)**(10.0-1))*np.exp(-(-x+225)/10.0) for x in range(len(template))]

        
        frags_list = []
        for i in range(len(dyad_list)):
            pos = dyad_list[i]
            if random.random() < left_prob[pos]:
                right_frag = template[pos-left_offset:]
            else:
                right_frag = template[:]
            if random.random() < right_prob[pos]:
                left_frag = rev_comp(template[:pos+right_offset+1])
            else:
                left_frag = rev_comp(template)
            frags_list.append(right_frag)
            frags_list.append(left_frag)
        return frags_list

    def mutations (seq_list, mist_prob, indel_prob):
        def mismatch(seq, prob):
            nts = set(['A','T','C','G'])
            new_seq = ""
            for i in range(len(seq)):
                if random.random() < prob:
                    subs = nts - set(seq[i])
                    new_seq += random.choice(list(subs))
                else:
                    new_seq += seq[i]
            return new_seq
        def indel(seq, prob):
            nts = ['A','T','C','G']
            new_seq = ""
            for i in range(len(seq)):
                if random.random() < prob:
                    if random.random() < 0.5 or i == len(seq):
                        new_seq += random.choice(nts)
                    else:
                        continue
                new_seq += seq[i]
            return new_seq
        new_list = []
        for seq in seq_list:
            new_seq = indel(mismatch(seq, prob=mist_prob), prob=indel_prob)
            new_list.append(new_seq)
        return new_list

    f = open(out_fname + '.fastq', 'w')

    for ref_id, template in ref_seq.items():
        #dyad_list = dyad_sampling(template, left_bound, right_bound, N=num)
        #frags_list = make_cleavage(template, dyad_list, left_prob, right_prob)
        frags_list = [template]*100
        seqs_list = mutations(frags_list, mist_prob=0.01, indel_prob=0.001)
        for i in range(len(seqs_list)):
            seq = seqs_list[i]
            print >> f, "@M01556:71:000000000-BHYK4:1:1101:11804:" + ref_id + " 1:N:0:CTTGTA"         # read ID
            print >> f, seq          # read seq
            print >> f, '+'          # optional
            print >> f, 'G'*len(seq) # quality score

    f.close()
        
            
if __name__ == '__main__':
    parser = ArgumentParser(description='make simulated data set for slide-seq')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence prefix filename')
    parser.add_argument('-o',
                        dest='out_fname',
                        type=str,
                        help='output prefix filename')
    parser.add_argument('-n',
                        dest='num',
                        type=int,
                        default = 10000, 
                        help='iteration number for each sequence')
    parser.add_argument('--left_bound',
                        dest='left_bound',
                        type=int,
                        default = 147/2, 
                        help='left bound length for dyad sampling')
    parser.add_argument('--right_bound',
                        dest='right_bound',
                        type=int,
                        default = 147/2, 
                        help='right bound length for dyad sampling')
    parser.add_argument('--left_prob',
                        dest='left_prob',
                        type=float,
                        default = 0.8, 
                        help='left cleavage probabilty')
    parser.add_argument('--right_prob',
                        dest='right_prob',
                        type=float,
                        default = 0.8, 
                        help='right cleavage probabilty')
    parser.add_argument('--mist_prob',
                        dest='mist_prob',
                        type=float,
                        default = 0.02, 
                        help='mismatch probabilty')
    parser.add_argument('--indel_prob',
                        dest='indel_prob',
                        type=float,
                        default = 0.03, 
                        help='insertion/deletion probabilty')

    args = parser.parse_args()

    ref_seq = {}
    for line in open(args.ref_fname + '.ref'):
        line = line.strip()
        if line.startswith('>'):
            ref_id = line[1:].strip()
            continue
        if line:
            assert ref_id not in ref_seq
            line = line.strip()
            ref_seq[ref_id] = line

    if not args.out_fname:
        out_fname = args.ref_fname
    else:
        out_fname = args.out_fname

    simulation (ref_seq,
                out_fname,
                args.num,
                args.left_bound,
                args.right_bound,
                args.left_prob,
                args.right_prob,
                args.mist_prob,
                args.indel_prob
                )
    
        
