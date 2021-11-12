import os, sys, subprocess, re
from argparse import ArgumentParser, FileType, ArgumentTypeError
import glob

def sort_reads (read_fname,
                ref_fname,
                out_fname,
                mm_cutoff,
                len_cutoff,
                ID_refseq,
                ID_wID_sted,
                bab_seq,
                include_all,
                direct_mode):
    
    # sequnece alignment by Bowtie2
    aligner_cmd=["bowtie2", '-x', ref_fname, '-U', read_fname ]
    #aligner_cmd += ['--mp',  '4,4', '--rdg', '5,4', '--rfg', '5,4'] # set score table
    #aligner_cmd += ['--mp',  '4,4', '--rdg', '0,5', '--rfg', '0,5'] # set score table
    aligner_cmd += ['--gbar', str(1)] # don't consider 1nt deletion on the edge of reads

    align_proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE,stderr=open("/dev/null", 'w'))

    # start sort the reads
    print >> sys.stderr, "Start sorting the reads ..."

    seq_sort_fname = open(out_fname+".sort",'w')

    type_count = {}
    for line in align_proc.stdout:
        if line.startswith('@'):
            continue
        #print line
        type, hit, cut_loc, read_seq = 'NA', 'NA', 'NA', 'NA'

        cols = line.strip().split()
        read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
        read_id=":".join(read_id.split(':')[3:])
        flag, pos = int(flag), int(pos)
        pos-=1

        # mapping failure: too low alignment score
        if pos < 0 or flag & 0x4 != 0:
            type = 'mutant'
        
        read_seq, qual =cols[9:11]
        ref_id = ref_id.strip()

        # mapping failure: too low alignment score
        if ref_id == '*':
            type = 'mutant'

        # mapping failure: mapping to backbone (non-direct mode)
        if not direct_mode and ref_id == 'BACKBONE':
            type = 'unIden'
            
        AS,XS,NM,MD = None, None, None, None
        for i in range(11, len(cols)):
            col = cols[i]
            if col.startswith('AS'):
                AS = int(col[5:])
            elif col.startswith('XS'):
                XS = int(col[5:])
            elif col.startswith('NM'):
                NM = int(col[5:])
            elif col.startswith('MD'):
                MD = col[5:]
               
        # mapping failure: ambiguous mapping (multiple best alignments)
        #if XS !=None and XS >= AS:
        #    type = 'multHit'
            
        # if mapping success, check the alignment in detail
        if type == 'NA':

            # direct alignment mode
            if direct_mode:
                # check all mismatch/indel information
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
                    if key in ID_refseq:
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
                    ref_len = len(ID_refseq[hit])
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


            # non-direct mode (default)
            else:

                hit = ref_id
                
                # invalid: too large edit distance 
                if NM > mm_cutoff:
                    type = 'mutant'

                try:
                    wID_sted = ID_wID_sted[hit]
                except:
                    wID_sted = {}

                ref_seq = ID_refseq[hit]

                # if possible, check the identification windows in read
                if len(wID_sted) > 0:
                    # read identification windows in reference
                    sted_wIDs = {}
                    for wID, sted in wID_sted.items():
                        for loc in list(sted):
                            if loc not in sted_wIDs:
                                sted_wIDs[loc] = []
                            sted_wIDs[loc].append(wID)

                    # locate the identification windows in read
                    ref_pt, read_pt = pos, 0
                    wID_readpos = {}
                    cigar_str=re.split('(\d+)',cigar_str)[1:]
                    for i in range(len(cigar_str)/2):
                        s = cigar_str[2*i+1]
                        num = int(cigar_str[2*i])
                        for j in range(num):                
                            if ref_pt in sted_wIDs:
                                for wID in sted_wIDs[ref_pt]:
                                    if wID not in wID_readpos:
                                        wID_readpos[wID] = []
                                    wID_readpos[wID].append(read_pt)
                            if s=='M':
                                ref_pt +=1
                                read_pt +=1
                            elif s=='I':
                                read_pt +=1
                            elif s=='D':
                                ref_pt +=1
                            else:
                                ref_pt +=1
                                read_pt +=1

                    # check the identification windows in read
                    for wID in wID_sted:
                        try:
                            read_st, read_ed = sorted(wID_readpos[wID])
                            # if mismatch window hit the very ends of reads, it could be truncated one.
                            if wID[0] == 'M':
                                if read_st == 0 or read_ed == len(read_seq)-1:
                                    type = 'unIden'
                                    break
                        except:
                            type = 'unIden'
                            break
                        ref_st, ref_ed = wID_sted[wID]
                        if ref_seq[ref_st:ref_ed+1] != read_seq[read_st:read_ed+1]:
                            type = 'unIden'
                            break

                else:
                    # find the right most end of alignment
                    ref_pt = pos
                    cigar_str=re.split('(\d+)',cigar_str)[1:]
                    for i in range(len(cigar_str)/2):
                        s = cigar_str[2*i+1]
                        num = int(cigar_str[2*i])
                        if s == 'I':
                            continue
                        ref_pt +=num


                # try to locate cleavage sites on the reference 
                if type == 'NA':
                    # alignment location on the reference
                    ref_len = len(ref_seq)
                    end_pos = min(ref_pt, ref_len) - 1

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

        # record read type and counts
        if type not in type_count:
            type_count[type] = 0
        type_count[type] +=1

        # record only valid data unless include_all is on
        if not include_all:
            if type !='valid':
                continue
            else:
                read_seq = 'N'

        print >> seq_sort_fname, "%s\t%s\t%s\t%s\t%s" % (read_id, type, hit, cut_loc, read_seq)

    seq_sort_fname.close()

    # remove index files
    subprocess.call('rm ' + ref_fname + '*bt2', shell=True)

    # display a short report
    total_count = sum(type_count.values())
    print >> sys.stderr, "%d reads are sorted" % (total_count)

    count_type = sorted([(count, type) for type, count in type_count.items()], reverse=True)
    for count, type in count_type:
        print >> sys.stderr, "%d reads are %s (%.2f%%)" % (count, type, 100*float(count)/total_count)
    print >> sys.stderr, "Done"
                        

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise ArgumentTypeError('Boolean value expected.')
                
    parser = ArgumentParser(description='Sort the slide-seq data and extract valid reads')
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
                        default=20,
                        help='minimun cutoff mismatch or indel count')
    parser.add_argument('-c',
                        dest="len_cutoff",
                        type=int,
                        default=10,
                        help='length cutoff for sorting free DNA')
    parser.add_argument('-a',
                        dest="include_all",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='include all reads information')
    parser.add_argument('--dw',
                        dest="no_check_window",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='disable checking the identification windows in reads')
    parser.add_argument('--direct',
                        dest="direct_mode",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='direct alignment only with backbone sequence to sort')
    
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
    else:
        assert len(args.filenames) == 0
        print >> sys.stderr, "Error: there is no read file input."
        sys.exit(1)
        
    assert read_type != None

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
        read_fname = common_name  + '.combined.fastq.gz'
        
        # combine pair-end reads
        print >> sys.stderr, "Combining pair-end reads"
        flash_cmd = ["flash", '-c', '-M', str(100), '-z', read_fname1, read_fname2]
        flash_proc = subprocess.Popen(flash_cmd, stdout=open(read_fname,'w'), stderr=open("/dev/null", 'w'))
        flash_proc.communicate()

    # read reference file
    if not args.ref_fname:
        print >> sys.stderr, "Error: there is no reference file input."
        sys.exit(1)

    ID_refseq = {}
    bab_seq = None
    for line in open(args.ref_fname + '.ref'):
        line = line.strip()
        if line.startswith('>'):
            ref_id = line[1:].strip()
            continue
        if ref_id == 'BACKBONE':
            bab_seq = line.strip()
            continue
        assert ref_id not in ID_refseq
        ID_refseq[ref_id] = line.strip()

    # no backbone sequence, no direct mode
    if bab_seq == None:
        direct_mode = False
    else:
        direct_mode = args.direct_mode

    # find the identification windows in the reference file if possible
    ID_wID_sted = { ID:{} for ID in ID_refseq }
    if not direct_mode and not args.no_check_window:
        for ID in ID_refseq:
            try:
                cols = ID.strip().split(',')
                w_list = []
                for i in range(len(cols)):
                    col = cols[i]
                    pos, mtype, nts = col.split('-')
                    pos = int(pos)
                    w_list.append((pos, mtype, nts))
                w_list = sorted(w_list)
                offset = 0
                for i in range(len(w_list)):
                    pos, mtype, nts = w_list[i]
                    pos += offset
                    # set up identification window in the reference
                    # inclusive end, 1nt padding for indel type
                    if mtype == 'M':
                        st, ed = pos, pos+len(nts)-1 
                    elif mtype == 'I':
                        st, ed = pos, pos+len(nts)+1
                        offset +=len(nts)
                    elif mtype == 'D':
                        st, ed = pos-1, pos
                        offset -=len(nts)
                    wID = (mtype, i)
                    ID_wID_sted[ID][wID] = [st, ed] 
            except:
                pass

        
    # make a reference index file for bowtie2
    if direct_mode:
        build_cmd = ["bowtie2-build", '-c', bab_seq, args.ref_fname]
    else:
        build_cmd = ["bowtie2-build", args.ref_fname + ".ref", args.ref_fname]
    print >> sys.stderr, "Building index files for Bowtie2 alignment"
    subprocess.call(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
            
    # set output file name
    if not args.out_fname:
        out_fname = read_fname.rsplit('.', 1)[0]
    else:
        out_fname = args.out_fname

    if direct_mode:
        print >> sys.stderr, "direct mode is on"
    if args.no_check_window:
        print >> sys.stderr, "disabled checking window"
    
    sort_reads (read_fname,
                args.ref_fname,
                out_fname,
                args.mm_cutoff,
                args.len_cutoff,
                ID_refseq,
                ID_wID_sted,
                bab_seq,
                args.include_all,
                direct_mode)
