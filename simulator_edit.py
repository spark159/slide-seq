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

    def dyad_sampling(ref_id, left_bound, right_bound, N):
        tlen = len(ref_seq[ref_id])
        try:
            loc, mtype, nts = ref_id.split('-')
            loc = int(loc)
            if mtype == 'I':
                tlen -=len(nts)
        except:
            pass
            
        dyad_list = []
        st, ed = left_bound, tlen - right_bound - 1
        for i in range(N):
            dyad_pos = random.randint(st,ed)
            dyad_list.append(dyad_pos)
        return dyad_list

    def make_cleavage(ref_id, dyad_list, left_prob, right_prob, left_offset=53, right_offset=53):
        seq = ref_seq[ref_id]
        tlen = len(seq)
        try:
            loc, mtype, nts = ref_id.split('-')
            loc = int(loc)
            if mtype == 'I':
                tlen -=len(nts)
        except:
            mtype = None
            pass

        left_prob = [left_prob] * tlen
        right_prob = [right_prob] * tlen
        frags_list = []
        cutloc_list = []
        for i in range(len(dyad_list)):
            pos = dyad_list[i]
            left_cut = pos - left_offset
            right_cut = pos + right_offset
            if mtype == 'I':
                if left_cut > loc:
                    left_cut += len(nts)
                if right_cut > loc:
                    right_cut += len(nts)
            # top strand cleavage
            if random.random() < left_prob[pos]:
                right_frag = seq[left_cut+1:]
                cutloc = 'R:' + str(left_cut)
            else:
                right_frag = seq[:]
                cutloc = 'NA'
            frags_list.append(right_frag)
            cutloc_list.append(cutloc)
            # Bottom strand cleavage
            if random.random() < right_prob[pos]:
                left_frag = rev_comp(seq[:right_cut])
                cutloc = 'L:' + str(right_cut)
            else:
                left_frag = rev_comp(seq)
                cutloc = 'NA'
            frags_list.append(left_frag)
            cutloc_list.append(cutloc)
        return frags_list, cutloc_list

    def mutations_simple (seq_list, mist_prob, indel_prob):
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


    def mutations (seq_list, mist_prob, indel_prob):
        def mismatch(seq, prob):
            mtype_pos_subs = {'M':{}}
            nts = set(['A','T','C','G'])
            new_seq = ""
            start = False
            for i in range(len(seq)):
                if random.random() < prob:
                    if not start:
                        pos = i
                        start = True
                    nt = random.choice(list(nts - set(seq[i])))
                    new_seq += nt
                    if pos not in mtype_pos_subs['M']:
                        mtype_pos_subs['M'][pos] = ""
                    mtype_pos_subs['M'][pos] += nt
                else:
                    new_seq += seq[i]
                    start = False
            return new_seq, mtype_pos_subs

        def indel(seq, prob):
            mtype_pos_indels = {'I':{}, 'D':{}}
            nts = ['A','T','C','G']
            new_seq = ""
            prev_s = 'M'
            
            # insert at the most left side
            indels = ""
            while random.random() < prob:
                indels += random.choice(nts)
            if len(indels) > 0:
                new_seq += indels
                mtype_pos_indels['I'][-1] = indels
                prev_s = 'I'

            # deletion on i-th base or insert right next to the i-th base
            i = 0
            while i < len(seq):
                if random.random() < prob:
                    if prev_s in 'ID':
                        s = 'I'
                    else:
                        s = random.choice(['I', 'D'])
                    pos = i
                    indels = ""
                    while i < len(seq):
                        if s == 'I':
                            indels += random.choice(nts)
                        elif s =='D':
                            indels += seq[i]
                            i +=1
                        if random.random() >= prob:
                            break
                    if s == 'I':
                        new_seq += seq[i] + indels
                        i +=1
                    mtype_pos_indels[s][pos] = indels
                    prev_s = s
                else:
                    new_seq += seq[i]
                    prev_s = 'M'
                    i +=1                    
            return new_seq, mtype_pos_indels

        newseq_list = []
        minfo_list = []
        for seq in seq_list:
            mis_seq, mtype_pos_subs = mismatch(seq, prob=mist_prob)
            new_seq, mtype_pos_indels = indel(mis_seq, prob=indel_prob)
            mtype_pos_subs.update(mtype_pos_indels)
            newseq_list.append(new_seq)
            minfo_list.append(mtype_pos_subs)
        return newseq_list, minfo_list

    
    f = open(out_fname + '_simulated.fastq', 'w')    
    
    for ref_id in ref_seq.keys():
        dyad_list = dyad_sampling(ref_id, left_bound, right_bound, N=num)
        frags_list, cutloc_list = make_cleavage(ref_id, dyad_list, left_prob, right_prob)
        reads_list, minfo_list = mutations(frags_list, mist_prob, indel_prob)
        #seqs_list = mutations_simple (frags_list, mist_prob, indel_prob)

        # make simulated fastq file
        for i in range(len(reads_list)):
            read = reads_list[i]
            cutloc = cutloc_list[i]
            print >> f, "@SIM:1:FCX:1:15:6329:%s:/%s 1:N:0:ATCCGA" % (i, ref_id + str('/') + cutloc)
            print >> f, read          # read seq
            print >> f, '+'          # optional
            print >> f, 'G'*len(read) # quality score
            

            # probably too much
            """
            # make expected sort file
            type = 'NA'

            candidates = []
            
            # check read include true ref-id
            mtype, pos, nts = red_id.split('-')
            start = int(pos)
            end = start + len(nts)
            
            cutloc = cutloc_list[i]
            if cutloc == 'NA':
                candidates.append(ref_id)
                type = 'freeDNA'
            else:
                side, loc = cutloc.split(':')
                loc = int(loc)
                if side == 'R' and loc < pos:
                    candidates.append(ref_id)
                if side == 'L' and pos >= end:
                    candidates.append(ref_id)
                    
            # include all possible ids genereated by mutations
            mtype_pos_nts = minfo_list[i]
            for mtype in mtype_pos_nts:
                for pos in mtype_pos_nts[mtype]:
                    nts = mtype_pos_nts[mtype][pos]
                    if cutloc != 'NA' and side == 'R':
                        pos += loc
                    elif cutloc != 'NA' and side == 'L':
                        pos -= loc
                        pos *= -1
                        if mtype == 'I':
                            pos -=1
                    key = '-'.join([mtype, pos, nts])
                    candidates.append(key)

            # check edge mutations 
            for key in candidates:
                mtype, pos, nts = key.split('-')
                pos = int(pos)
                if pos == -1 or pos == len(frags_list[i])-1:
                    candidates = []
                    break

            hits, nonhits = [], []
            edit_dist = 0
            for key in candidates:
                if key in ref_seq:
                    hits.append(key)
                else:
                    nonhits.append(key)
                    _, _, nts = key.split('-')
                    edit_dist += len(nts)

            if edit_dist > mm_cutoff:
                type = 'mutant'
            else:
                if len(hits) < 1:
                    type = 'unIden'
                elif len(hits) > 1:
                    type = 'multHit'
                    hit = ':'.join(hits)
                else:
                    hit = hits[0]
            """        
        

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
                        default = 1000, 
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
                        default = 0.05, 
                        help='mismatch probabilty')
    parser.add_argument('--indel_prob',
                        dest='indel_prob',
                        type=float,
                        default = 0.01, 
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

    # set output file name
    if not args.out_fname:
        out_fname = args.ref_fname.rsplit('.', 1)[0]
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
    
        
