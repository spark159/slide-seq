
def check_bubble (fname, template):
    id_bubbles = {}
    for line in open(fname):
        if line.startswith('>'):
            id = line[1:].strip()
            continue
        seq = line.strip()
        assert len(seq) == len(template)
        i, j = 0, 1 
        while i < len(seq):
            if seq[i] != template[i]:
                j = 1
                while i+j < len(seq) and seq[i+j] != template[i+j]:
                    j += 1
                if id not in id_bubbles:
                    id_bubbles[id] = []
                id_bubbles[id].append((j, i))
                i += j 
            i += 1
    return id_bubbles

def tuple_cmp (a, b):
    if a[0] < b[0]:
        return -1
    elif a[0] > b[0]:
        return 1
    else:
        if a[1] < b[1]:
            return -1
        else:
            return 1


template = "ATCCGACTGGCACCGGCAAGGTCGCTGTTCGCCACATGCGCAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCAGGGCGTCCTCGTATAGGGTCCATCACATAAGGGATGAACT"

id_bubbles = check_bubble("polyAscanlib.ref", template)

num_ids = {}
for id, bubbles in id_bubbles.items():
    win, loc = id.split('-')
    size, loc = len(win), int(loc)
    if size > 5:
        continue
    if len(bubbles) > 1:
        continue
    for bubble in bubbles:
        num, pos = bubble
        if num not in num_ids:
            num_ids[num] = []
        num_ids[num].append((id, pos))

keys1, keys2 = [], []
for info in num_ids[1]:
    keys1.append(info[0])
for info in num_ids[2]:
    keys2.append(info[0])
s1 = ','.join(keys1)
s2 = ','.join(keys2)

print s1
print
print s2

