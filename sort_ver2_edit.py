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
    aligner_cmd += ['--mp',  '4,4', '--rdg', '5,4', '--rfg', '6,5'] # focus on insertion
    aligner_cmd += ['--mp',  '4,4', '--rdg', '6,5', '--rfg', '5,4'] # focus on deletion
    #aligner_cmd += ['-a']
    #aligner_cmd += ['--very-sensitive']
    #aligner_cmd += ['--end-to-end']
    #aligner_cmd += ['-D', '50', '-R', '10', '-N', '1', '-L', '5', '-i', 'S,1,0.50']
    #aligner_cmd += ['--skip', str(0)]
    #aligner_cmd += ['-N', str(1)]
    aligner_cmd += ["--gbar", str(1)]
    
    #aligner_cmd += ['--very-sensitive']
    #aligner_cmd += ['-N', str(1), '-L', str(1), '-i',  'S,1,0'] # turn off heuristic seeding
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
        read_id=":".join(read_id.split(':')[3:8])
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

        
        read_seq, qual =cols[9:11]
        ref_id = ref_id.strip()    

        AS,XS,NM,MD = None, None, None, None
        for i in range(11, len(cols)):
            col = cols[i]
            if col.startswith('AS'):
                AS = int(col[5:])
            elif col.startswith('XS'):
                XS = int(cols[5:])
            elif col.startswith('NM'):
                NM = int(col[5:])
            elif col.startswith('MD'):
                MD = col[5:]

        ## invalid: ambiguous mapping
        #mapQ = float(mapQ)
        #if mapQ < 42:
        #    print mapQ
        #    print read_seq
        #    print
        #    type = 'invalid:multimap'

                
        ## ambiguous mapping (multiple best alignments)
        #if XS !=None and XS >= AS:
        #    print read_seq
        #    type = 'unIden'
        #    print >> seq_sort_fname, "%s\t%s\t%s\t%s\t%s" % (read_id, type, hit, cut_loc, read_seq)
        #    continue
        
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
        bab_pt, read_pt = pos, 0
        mtype_pt_nts = {}
        edge_mut = False 
        
        cigar_str=re.split('(\d+)',cigar_str)[1:]
        for i in range(len(cigar_str)/2):
            s = cigar_str[2*i+1]
            num = int(cigar_str[2*i])

            if i in [0, len(cigar_str)/2-1]:
                if s in 'ID':
                    edge_mut = True

            if s not in mtype_pt_nts:
                mtype_pt_nts[s] = {}

            if s == 'M':
                start = False
                for j in range(num):
                    if bab_seq[bab_pt] != read_seq[read_pt]:
                        if i in [0, len(cigar_str)/2-1]:
                            if read_pt == 0 or read_pt == len(read_seq)-1:
                                edge_mut = True
                        if not start:
                            start = True
                            pt = bab_pt
                        if pt not in mtype_pt_nts[s]:
                            mtype_pt_nts[s][pt] = ""
                        nt = read_seq[read_pt]
                        mtype_pt_nts[s][pt] +=nt
                    else:
                        start = False
                    bab_pt +=1
                    read_pt +=1
                
            elif s == 'I':
                nts = read_seq[read_pt:read_pt+num]
                mtype_pt_nts[s][bab_pt-1] = nts # insert right next to the position
                read_pt +=num

            elif s == 'D':
                nts = bab_seq[bab_pt:bab_pt+num]
                mtype_pt_nts[s][bab_pt] = nts
                bab_pt +=num

            else:
                bab_pt +=num
                read_pt +=num

                
        # if mutation on the very edge of reads, it can't be indentified
        if edge_mut:
            mtype_pt_nts = {}

            
        # make all possible ref IDs by using mismatch/indel information
        candidates = []
        for mtype in mtype_pt_nts:
            for pt in mtype_pt_nts[mtype]:
                nts = mtype_pt_nts[mtype][pt]
                candidates.append('-'.join([str(pt), mtype, nts]))

        # try to identify unique ref ID
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

                
        # try to locate cleavage sites on the reference 
        if type == 'NA':
            # alignment location w.r.t. backbone seq
            end_pos = min(bab_pt, len(bab_seq)) - 1

            # adjust alignment location w.r.t. each ref seq
            ref_len = len(ref_seq[hit])
            _, mtype, nts = hit.split('-')
            if mtype == 'I':
                end_pos += len(nts)
            elif mtype == 'D':
                end_pos -= len(nts)

            # record cleavage locations on the reference if possible
            if pos < len_cutoff and end_pos > ref_len - len_cutoff - 1:
                type = 'freeDNA'
            elif pos < len_cutoff and end_pos <= ref_len - len_cutoff - 1: 
                cut_loc = 'L:' + str(end_pos+1)
            elif pos >= len_cutoff and end_pos > ref_len - len_cutoff - 1:
                cut_loc = 'R:' + str(pos-1)
            else:
                type = 'frag'

        # if all success, valid data
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

        common_name = common(read_fname1.rsplit('.', 1)[0], read_fname2.rsplit('.', 1)[0])
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
        out_fname = read_fname.rsplit('.', 1)[0]
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
