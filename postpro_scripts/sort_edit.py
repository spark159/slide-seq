import os, sys, subprocess, re
from argparse import ArgumentParser, FileType

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def sort_seqs (read_fname,
               ref_fname,
               out_fname,
               mm_cutoff,
               len_cutoff,
               ref_length,
               ref_win,
               size_sted_IDs,
               check_window,
               discard_multHit,
               small_sort):

    """
    # build bowtie2 index
    index_cmd = ["bowtie2-build", ref_fname + '.ref', ref_fname]
    index_proc = subprocess.Popen(index_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))
    index_proc.communicate()
    """
    
    # seq alignment
    aligner_cmd=["bowtie2", '-x', ref_fname, '-U', read_fname ]
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
        type, cut_loc, read_seq = 'NA', 'NA', 'NA'

        cols = line.strip().split()
        read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
        read_id=":".join(read_id.split(':')[3:7])
        read_count +=1
        flag, pos = int(flag), int(pos)
        pos-=1

        # invalid: mapping failure
        if pos < 0:
            type = 'invalid:mutant'
            #continue
        if flag & 0x4 != 0:
            type = 'invalid:mutant'
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

        # invalid: too large edit distance 
        if NM > mm_cutoff:
            type = 'invalid:mutant'
            #continue

        # collect invalid data 
        if type != 'NA':
            assert type.startswith('invalid')
            #assert read_id not in seq_sort
            #seq_sort[read_id]=[type, insert, cut_loc, read_seq]
            if small_sort:
                continue
            print >> seq_sort_fname, "%s\t%s\t%s\t%s\t%s" % (read_id, type, ref_id, cut_loc, read_seq)
            continue

        if check_window:
            win, win_st = ref_id.split('-')
            win_st, win_size = int(win_st), len(win)
            win_ed = win_st + win_size
            if discard_multHit: # look all windows with same size
                sted_IDs = size_sted_IDs[win_size]
            else:
                sted_IDs = {} # look only aligned window
                sted_IDs[win_st] = [ref_id]
                sted_IDs[win_ed] = [ref_id]
        else: # no screen window
            sted_IDs = {}
                    
        
        # find mismatch information
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

        # get read points in the sensitive window
        ref_pt, read_pt = pos-1, -1
        ID_readpos = {}
        cigar_str=re.split('(\d+)',cigar_str)[1:]
       
        for i in range(len(cigar_str)/2):
            s = cigar_str[2*i+1]
            num = int(cigar_str[2*i])
            for j in range(num):
                if read_pt > 0  and ref_pt in sted_IDs: # not include the first nt of reads
                    for ID in sted_IDs[ref_pt]:
                        if ID not in ID_readpos:
                            ID_readpos[ID] = []
                        ID_readpos[ID].append(read_pt)
                if s=='M':
                    ref_pt += 1
                    read_pt += 1
                elif s=='I':
                    read_pt += 1
                elif s=='D':
                    ref_pt += 1
        # check the last aligned base-pair
        if read_pt < len(read_seq)  and ref_pt in sted_IDs: # not include the last nt of reads
            for ID in sted_IDs[ref_pt]:
                if ID not in ID_readpos:
                    ID_readpos[ID] = []
                ID_readpos[ID].append(read_pt)

        # sort by read length and alignment position
        ref_len = ref_length[ref_id]
        end_pos = min(ref_pt, ref_len - 1)        

        if pos < len_cutoff and end_pos > ref_len - len_cutoff - 1:
            type = 'freeDNA'
        elif pos < len_cutoff and end_pos <= ref_len - len_cutoff - 1: 
            cut_loc = 'L:' + str(end_pos)
        elif pos >= len_cutoff and end_pos > ref_len - len_cutoff - 1:
            cut_loc = 'R:' + str(pos)
        else:
            type = 'invalid:frag'

        #print check_window
        #print ID_readpos.keys()
        # screen the read sequences in windows
        if check_window:
            candidates = []
            for ID in ID_readpos:
                try:
                    st, ed = sorted(ID_readpos[ID])
                    insert = read_seq[st:ed]
                    if len(insert) != win_size:
                        continue
                    key = insert + '-' + ID.split('-')[1]
                    if key in ref_win and key not in candidates:
                        candidates.append(key)
                    #if len(candidates) > 1:
                    #    break
                except:
                    pass
            #print ref_id
            #print candidates
            #print 
            if len(candidates) <=0:
                type = 'invalid:instErr:' + 'NA'
            else:                
                if ref_id not in candidates: # check bowtie alignment
                    type = 'invalid:instErr:'
                    for insert in candidates:
                        type += ':' + insert
                elif ref_id in candidates and len(candidates) > 1: # check uniqueness
                    type = 'invalid:multHit:'
                    for insert in candidates:
                        type += ':' + insert

        # otherwise, valid data
        if type == 'NA':
            type = 'valid'

        #assert read_id not in seq_sort
        #seq_sort[read_id]=[type, insert, cut_loc, read_seq]
        if small_sort:
            if type !='valid':
                continue
            else:
                read_seq = 'N'
                
        print >> seq_sort_fname, "%s\t%s\t%s\t%s\t%s" % (read_id, type, ref_id, cut_loc, read_seq)
             
    seq_sort_fname.close()

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
    parser.add_argument('-w',
                        dest="check_window",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=True,
                        help='check the sensitive windows if possible')
    parser.add_argument('-d',
                        dest="discard_multHit",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='discard multi-hits reads')
    parser.add_argument('-s',
                        dest="small_sort",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='record only valid cleavage data in sort file')

    
    args = parser.parse_args()

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
    
    if not args.ref_fname:
        print >> sys.stderr, "Error: there is no reference file input."
        sys.exit(1)

    check_window = args.check_window
    discard_multHit = args.discard_multHit

    ref_seq, ref_length = {}, {}
    ref_win = {}
    size_sted_IDs = {}
    for line in open(args.ref_fname + '.ref'):
        line = line.strip()
        if line.startswith('>'):
            ref_id = line[1:].strip()
            try: # if library has sensitive windows
                win, win_st = ref_id.split('-')
                win_st, win_size = int(win_st), len(win)
                win_ed = win_st + win_size

                if check_window and discard_multHit: # gather all window sites if concerning multi-hits
                    if win_size not in size_sted_IDs:
                        size_sted_IDs[win_size] = {}
                    if win_st not in size_sted_IDs[win_size]:
                        size_sted_IDs[win_size][win_st] = []
                    size_sted_IDs[win_size][win_st].append(ref_id)
                    if win_ed not in size_sted_IDs[win_size]:
                        size_sted_IDs[win_size][win_ed] = []
                    size_sted_IDs[win_size][win_ed].append(ref_id)
                    
            except: # if library doesn't have sensitive windows
                win = ""
                check_window = False
            continue
        if line:
            assert ref_id not in ref_seq
            assert ref_id not in ref_length
            line = line.strip()
            ref_seq[ref_id] = line
            ref_length[ref_id] = len(line)
            ref_win[ref_id] = win

    if not check_window: # if no window, no screen multiHit
        discard_multHit = False
    
    if not args.out_fname:
        out_fname = read_fname.rsplit('.', 1)[0]
    else:
        out_fname = args.out_fname
    
    sort_seqs (read_fname,
               args.ref_fname,
               out_fname,
               args.mm_cutoff,
               args.len_cutoff,
               ref_length,
               ref_win,
               size_sted_IDs,
               check_window,
               discard_multHit,
               args.small_sort)
