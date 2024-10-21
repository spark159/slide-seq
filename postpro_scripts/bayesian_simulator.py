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
                dyad_offset,
                mist_prob,
                indel_prob):

    def dyad_sampling(template, left_bound, right_bound, N):
        dyadmap = [0.0]*len(template)
        st, ed = left_bound, len(template) - 1 - right_bound
        for i in range(N):
            pos = random.randint(st,ed)
            dyadmap[pos] += 1
        return dyadmap

    def get_lambda(template, dyadmap, tlist, blist, eplist, delist, offset):
        tlamb, blamb = [0.0]*len(template), [0.0]*len(template)
        for i in range(len(template)):
            if i + offset < len(template):
                tlamb[i] += dyadmap[i+offset]*tlist[i]
            if i -  offset >= 0:
                blamb[i] += dyadmap[i-offset]*blist[i]
            total1, total2 = 0.0, 0.0
            for j in range(len(template)):
                if j != i + offset:
                    total1 += dyadmap[j]
                if j != i - offset:
                    total2 += dyadmap[j]
            tlamb[i] += total1*eplist[i]
            blamb[i] += total2*delist[i]
        return tlamb, blamb

    def make_cleavage(template, tlamb, blamb):
        tnum_list, bnum_list = [], []
        frags_list = []
        for i in range(len(template)):
            tnum = np.random.poisson(tlamb[i],1)
            bnum = np.random.poisson(blamb[i],1)
            tnum_list.append(tnum)
            bnum_list.append(bnum)
            right_frag = template[i:]
            left_frag = rev_comp(template[:i+1])
            for k in range(tnum):
                frags_list.append(right_frag)
            for k in range(bnum):
                frags_list.append(left_frag)
        return tnum_list, bnum_list, frags_list
    
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
        #if not (ref_id.startswith('AAA-46') or ref_id.startswith("AAAAAAAAAAAA-157")):
        #    continue
        print ref_id
        dyadmap = dyad_sampling(template, left_bound, right_bound, N=num)
        #dyadmap = [0.0]*(225/2) + [100] + [0.0]*(225/2)
        tlist, blist = [1.0]*len(template), [0.5]*len(template)
        eplist, delist = [0.001]*len(template), [0.0005]*len(template)
        tlamb, blamb = get_lambda(template, dyadmap, tlist, blist, eplist, delist, offset=dyad_offset)
        tnum_list, bnum_list, frags_list = make_cleavage(template, tlamb, blamb)
        seqs_list = mutations(frags_list, mist_prob, indel_prob)
        for i in range(len(seqs_list)):
            seq = seqs_list[i]
            print >> f, "@M01556:71:000000000-BHYK4:1:1101:11804:1000 1:N:0:CTTGTA"         # read ID
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
                        help='dyad sampling number for each sequence')
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
    parser.add_argument('--dyad-offset',
                        dest="dyad_offset",
                        default=52,
                        type=int,
                        help='off-set length from cut site to dyad position')
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
                args.dyad_offset,
                args.mist_prob,
                args.indel_prob
                )
    
        
