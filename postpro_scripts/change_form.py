def change (fname, output_fname):
    f = open(output_fname, 'w')
    for line in open(fname):
        seq, score = line.strip().split()
        seq = seq.strip("'")
        seq = seq[23:len(seq)-23]
        assert len(seq) == 54
        s = seq + '\t' + score
        print >> f, s
    f.close()

change("/home/spark159/../../media/spark159/sw/all_slide_seq_data/F0_yeast.loopseq", "F0_yeast.loopseq")
change("/home/spark159/../../media/spark159/sw/all_slide_seq_data/F0_random.loopseq", "F0_random.loopseq")

