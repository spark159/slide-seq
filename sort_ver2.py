import os, sys, subprocess, re
from argparse import ArgumentParser, FileType

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def sort_seqs (read_fname,
               bt_fname,
               out_fname,
               mm_cutoff,
               len_cutoff,
               bab_seq,
               ref_seq,
               small_sort):

    """
    # build bowtie2 index
    index_cmd = ["bowtie2-build", ref_fname + '.ref', ref_fname]
    index_proc = subprocess.Popen(index_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))
    index_proc.communicate()
    """
    
    # seq alignment
    aligner_cmd=["bowtie2", '-x', bt_fname, '-U', read_fname ]
    #aligner_cmd += ['--n-ceil', 'L,'+str(insert_len)+','+str(0.15)]
    align_proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE,stderr=open("/dev/null", 'w'))

    # start sort the reads
    seq_sort_fname=open(out_fname+".sort",'w')
    #seq_sort={}
    
    read_count=0

    for line in align_proc.stdout:
        if line.startswith('@'):
            continue
        #print line
        type, hit, cut_loc, read_seq = 'NA', 'NA', 'NA', 'NA'

        cols = line.strip().split()
        read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
        read_id=":".join(read_id.split(':')[3:7])
        read_count +=1
        flag, pos = int(flag), int(pos)
        pos-=1

        # invalid: mapping failure
        if pos < 0:
            type = 'mutant'
            #continue
        if flag & 0x4 != 0:
            type = 'mutant'
            #continue

        # invalid: ambiguous mapping
        #mapQ = float(mapQ)
        #if mapQ < 10:
        #    type = 'invalid:multimap'
        
        read_seq, qual =cols[9:11]
        ref_id = ref_id.strip()    

        AS,NM,MD = None, None, None
        for i in range(11, len(cols)):
            col = cols[i]
            if col.startswith('AS'):
                AS = int(col[5:])
            elif col.startswith('NM'):
                NM = int(col[5:])
            elif col.startswith('MD'):
                MD = col[5:]
        
        # find mismatch/indel information
        #MD = re.findall('\d+|[A-Z]|\^[A-Z]+', MD)
        #pt, mismatch = 0, {}
        #for tag in MD:
        #    if re.search('\d+', tag):
        #        pt += int(tag)
        #    elif tag.startswith('^'):
        #        pt += len(tag[1:])
        #    else:
        #        mismatch[pt] = tag
        #        pt += len(tag)
        

        # check mismatch/indel information
        bab_pt, read_pt = pos-1, -1
        mut_list = []
        edge_mut = False # check read has mutation on either ends
        
        cigar_str=re.split('(\d+)',cigar_str)[1:]
        for i in range(len(cigar_str)/2):
            s = cigar_str[2*i+1]
            num = int(cigar_str[2*i])
            for j in range(num):
                if s=='M':
                    bab_pt += 1
                    read_pt += 1
                    if bab_seq[bab_pt] != read_seq[read_pt]:
                        if read_pt == 0 or read_pt == len(read_seq) - 1:
                            edge_mut = True
                        mut = ('M', bab_pt, read_seq[read_pt])
                        mut_list.append(mut)
                elif s=='I':
                    read_pt += 1
                    mut = ('I', bab_pt, read_seq[read_pt])
                    mut_list.append(mut)
                elif s=='D':
                    bab_pt += 1
                    mut = ('D', bab_pt, bab_seq[bab_pt])
                    mut_list.append(mut)

        if edge_mut:
            mut_list = [] # consider as unIden

        mut_list = sorted(mut_list)

        # make possible seq ID by using mismatch/indel information
        candidates = []
        M_key, I_key, D_key = '', '', ''
        M_prev, I_prev, D_prev = -2, -2, -2
        for mut in mut_list:
            mut_type, pt, nt = mut
            if mut_type == 'M':
                if pt == M_prev + 1:
                    M_key += nt
                else:
                    if M_key:
                        candidates.append(M_key)
                    M_key = str(pt) + '-M-' + nt
                M_prev = pt
            elif mut_type == 'I':
                if pt == I_prev:
                    I_key += nt
                else:
                    if I_key:
                        candidates.append(I_key)
                    I_key = str(pt) + '-I-' + nt
                I_prev = pt
            elif mut_type == 'D':
                if pt == D_prev + 1:
                    D_key += nt
                else:
                    if D_key:
                        candidates.append(D_key)
                    D_key = str(pt) + '-D-' + nt
                D_prev = pt
        if M_key:
            candidates.append(M_key)
        if I_key:
            candidates.append(I_key)
        if D_key:
            candidates.append(D_key)

        # sort by hit information
        hits, nonhits = [], []
        for key in candidates:
            if key in ref_seq:
                hits.append(key)
            else:
                nonhits.append(key)

        if len(nonhits) > mm_cutoff:
            type = 'mutant'
        else:
            if len(hits) < 1:
                type = 'unIden'
            elif len(hits) > 1:
                type = 'multHit'
                hit = ':'.join(hits)
            else:
                hit = hits[0]
                    
        # sort further by read length and alignment position
        if type == 'NA':
            ref_len = len(ref_seq[hit])
            
            temp = hit.split('-')
            m, l = temp[1], len(temp[2])
            if m == 'M':
                end_pos = min(bab_pt, ref_len - 1)
            elif m == 'I':
                end_pos = min(bab_pt + l, ref_len - 1)
            elif m == 'D':
                end_pos = min(bab_pt - l, ref_len - 1)

            if pos < len_cutoff and end_pos > ref_len - len_cutoff - 1:
                type = 'freeDNA'
            elif pos < len_cutoff and end_pos <= ref_len - len_cutoff - 1: 
                cut_loc = 'L:' + str(end_pos)
            elif pos >= len_cutoff and end_pos > ref_len - len_cutoff - 1:
                cut_loc = 'R:' + str(pos)
            else:
                type = 'frag'

        # otherwise, valid data
        if type == 'NA':
            type = 'valid'

        if small_sort:
            if type !='valid':
                continue
            else:
                read_seq = 'N'
                
        print >> seq_sort_fname, "%s\t%s\t%s\t%s\t%s" % (read_id, type, hit, cut_loc, read_seq)
             
    seq_sort_fname.close()
    subprocess.call(['rm', 'temporal.ref'])
    subprocess.call('rm ' + bt_fname + '*bt2', shell=True)

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')
                
    parser = ArgumentParser(description='Sort the slide_seq data and extract valid reads')
    parser.add_argument(dest="filenames",
                        type=str,
                        nargs='+',
                        help='read file names (single or pair1,pair2)')
    parser.add_argument(dest='ref_fname',
                        type=str,
                        help='reference prefix filename')
    parser.add_argument('-o',
                        dest='out_fname',
                        type=str,
                        help='output prefix filename')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='minimun cutoff mismatch or indel length')
    parser.add_argument('-c',
                        dest="len_cutoff",
                        type=int,
                        default=10,
                        help='length cutoff for sorting free DNA')
    parser.add_argument('-s',
                        dest="small_sort",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='record only valid cleavage data in sort file')
    
    args = parser.parse_args()

    # identify single or pair-reads input
    read_type=None
    if len(args.filenames) == 1:
        read_type='single'
        read_fname = args.filenames[0]
    elif len(args.filenames) == 2:
        read_type='pair'
        read_fname1, read_fname2 = args.filenames
    elif len(args.filenames) > 2:
        print >> sys.stderr, "Error: too much file names."
        sys.exit(1)

    assert read_type

    # if pair-reads input, combine it by using FLASh
    if read_type == 'pair':
        # set combined fastq filename
        def common (str1, str2):
            N=min(len(str1), len(str2))
            name=''
            for i in range(N):
                if str1[i] == str2[i]:
                    name += str1[i]
                else:
                    break
            name=name.strip()
            if not name:
                return 'out'
            return name

        common_name = common(read_fname1.rsplit('.')[0], read_fname2.rsplit('.')[0])
        read_fname = common_name  + '.combined.fastq'
        
        # combine pair-end reads
        flash_cmd = ["flash", '-c', read_fname1, read_fname2]
        flash_proc = subprocess.Popen(flash_cmd, stdout=open(read_fname,'w'), stderr=open("/dev/null", 'w'))
        flash_proc.communicate()

    # read reference file
    if not args.ref_fname:
        print >> sys.stderr, "Error: there is no reference file input."
        sys.exit(1)

    ref_seq = {}
    bab_seq = None
    for line in open(args.ref_fname + '.ref'):
        line = line.strip()
        if line.startswith('>'):
            ref_id = line[1:].strip()
            continue
        if ref_id == 'BACKBONE':
            bab_seq = line.strip()
            continue
        assert ref_id not in ref_seq
        ref_seq[ref_id] = line.strip()

    assert ref_seq
    assert bab_seq

    # make a temporal backbone index file for bowtie2
    f = open("temporal.ref", "w")
    print >> f, ">BACKBONE"
    print >> f, bab_seq
    f.close()
    subprocess.call(["bowtie2-build", "temporal.ref", "temporal"])
    bt_fname = "temporal"
            
    # set output file name
    if not args.out_fname:
        out_fname = read_fname.rsplit('.')[0]
    else:
        out_fname = args.out_fname
    
    sort_seqs (read_fname,
               bt_fname,
               out_fname,
               args.mm_cutoff,
               args.len_cutoff,
               bab_seq,
               ref_seq,
               args.small_sort)
