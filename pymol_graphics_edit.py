import os, sys, subprocess, re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import random  
import math
import copy
import pickle
import matplotlib.cm as cm
import random
import analysis_final as analysis
#import pymol
#from dssr_block import dssr_block

aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
           'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def read_FASTA_file (fname):
    id_seq = {}
    for line in open(fname):
        if not line.strip():
            continue
        if line.startswith('>'):
            id = line[1:].strip()
        else:
            seq = line.strip()
            id_seq[id] = seq
    return id_seq

def read_FASTA_string (s):
    id_seq = {}
    for line in open(fname):
        if not line.strip():
            continue
        if line.startswith('>'):
            id = line[1:].strip()
        else:
            seq = line.strip()
            id_seq[id] = seq
    return id_seq


"""
class Molecules:
    def __init__(self, code):
        self.code = code
        self.chain_resi_resn = {}
        self.chain_resi_index_atom = {}

        # read data from pymol loaded struture
        #pymol.pymol_argv = ['pymol','-qc']
        #pymol.finish_launching()
        #pymol.cmd.log_open()
        pymol.cmd.fetch(code)
        
        for atom in pymol.cmd.get_model(code).atom:
            chain = atom.chain
            name = atom.name
            resn = atom.resn.upper()
            resi = int(atom.resi)
            x, y, z = atom.coord
            x, y, z = float(x), float(y), float(z)
            index = int(atom.index)
            b = float(atom.b)
            
            if chain not in self.chain_resi_resn:
                self.chain_resi_resn[chain] = {}
            if resi not in self.chain_resi_resn[chain]:
                self.chain_resi_resn[chain][resi] = resn

            if chain not in self.chain_resi_index_atom:
                self.chain_resi_index_atom[chain] = {}
            if resi not in self.chain_resi_index_atom[chain]:
                self.chain_resi_index_atom[chain][resi] = {}
            assert index not in self.chain_resi_index_atom[chain][resi]
            self.chain_resi_index_atom[chain][resi][index] = {"name":name,
                                                             "resn":resn,
                                                             "coord":(x, y, z),
                                                             "b-factor":b}
            
        self.chain_seq = {}
        for chain in self.chain_resi_resn:
            seq = ""
            for resi in sorted(self.chain_resi_resn[chain].keys()):
                resn = self.chain_resi_resn[chain][resi]
                if len(resn) == 3:
                    try:
                        resn = AA_dict[resn]
                    except:
                        continue
                        #resn = '[' + resn + ']'
                elif len(resn) == 2 and resn.startswith('D'):
                    resn = resn[1:]
                else:
                    continue
                    #print (resn)
                    #pass
                    #assert len(resn) == 1
                seq += resn
            self.chain_seq[chain] = seq

        print("%s is loaded" % (code), file=sys.stderr)


    def print_seq(self):
        for chain in sorted(self.chain_seq.keys()):
            print("chain %s" % (chain), file=sys.stderr)
            print(self.chain_seq[chain], file=sys.stderr)

    def remove_ions(self):
        # remove waters and ions
        pymol.cmd.remove('resn hoh')
        pymol.cmd.remove('resn mn')
        pymol.cmd.remove('resn cl')
        return 
        

    def stylish(self):
        # stylish DNA
        #pymol.cmd.cartoon('oval')
        #pymol.cmd.set('cartoon_oval_length', '1')
        #pymol.cmd.set('cartoon_oval_width', '0.2')
        pymol.cmd.set('cartoon_ring_mode', '3')
        pymol.cmd.set('cartoon_ring_finder', '2')
        pymol.cmd.set('cartoon_ring_width', '0.2')
        pymol.cmd.dssr_block()
        
        # turn off all reflections
        #pymol.cmd.set('reflect', '0')
        #pymol.cmd.set('light_count', '1')
        #pymol.cmd.set('ambient', '1')

        # cartoonish ray setting
        pymol.cmd.hide('lines')
        pymol.cmd.set('ray_trace_mode', '3')
        pymol.cmd.bg_color('white')
        pymol.cmd.set('antialias', '5')
        ##pymol.cmd.set('ray_trace_fog', '0')
        pymol.cmd.set('ray_shadows', '0')

        
        #pymol.cmd.space('pymol')
        #pymol.cmd.cartoon('oval', string selection)        
        return


    def make_sphere (self, chain_resi):

        selections = []
        for chain in chain_resi:
            select = []
            for resi in chain_resi[chain].keys():
                if resi < 0:
                    select.append('\\' + str(resi))
                else:
                    select.append(str(resi))
            select = ','.join(select)
            selections.append("(" + "chain " + chain + " and " + "resi " + select + ")")
        selections = " or ".join(selections)

        pymol.cmd.show('sphere', selections)
        pymol.cmd.set("sphere_scale", '0.8')
        return

    def coloring (self, chain_resi, color):
        for chain in chain_resi:
            for resi in chain_resi[chain].keys():
                if resi < 0:
                    select = '\\' + str(resi)
                else:
                    select = str(resi)
                pymol.cmd.color(color, "(" + "chain " + chain + " and " + "resi " + select + ")")
        return
        


    def spectrum(self,
                 chain_resi_value,
                 color_list=[],
                 min=None,
                 max=None):

        selections = []
        for chain in chain_resi_value:
            select = []
            for resi in chain_resi_value[chain].keys():
                if resi < 0:
                    select.append('\\' + str(resi))
                else:
                    select.append(str(resi))
            select = ','.join(select)
            selections.append("(" + "chain " + chain + " and " + "resi " + select + ")")
        selections = " or ".join(selections)

        myspace = {"chain_resi_value":chain_resi_value}
        pymol.cmd.alter(selections, 'b=chain_resi_value[chain][resv]', space=myspace)
        
        if len(color_list) > 0:
            palette = " ".join(color_list)
        else:
            palette = "rainbow"
            
        pymol.cmd.spectrum('b', palette, selections, min, max)


    def save_session (self, fname):
        pymol.cmd.save("%s.pse" % (fname))
        return

    def clear_up (self):
        pymol.cmd.delete('all')
        return

    def done (self):
        pymol.cmd.quit()
        return

    

# load structure
#NCP = Molecules("1kx5")
#histone_chains = {'H2A':['C', 'G'], 'H2B':['D', 'H'], 'H3':['A', 'E'], 'H4':['B','F']}


# load energy profile
path = "/home/spark159/Projects/slide-seq/"
fname = "mmlib_bubble_5_1rep_energy_wrt_601"
with open(path+fname + ".pickle", "rb") as f:
    size_dyad_shl_values = pickle.load(f,encoding='latin1')

# load structure
NCPandChd1 = Molecules("5O9G")
NCPandChd1.print_seq()
#NCPandChd1.spectrum({"A":{38:1, 39:2}}, ["blue", "white", "red"], min=1, max=2)

# assign value on Widom 601 sequence
Widom_seq = "ACAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCTGCAGCACCGGGATTCTCCAG"
top_strand = NCPandChd1.chain_seq['I']
bott_strand = NCPandChd1.chain_seq['J']

size = 3
dyad = 92

chain_resi_value = {}
for shl in size_dyad_shl_values[size][dyad]:
    value = np.mean(size_dyad_shl_values[size][dyad][shl])

    if 'I' not in chain_resi_value:
        chain_resi_value['I'] = {}
    chain_resi_value['I'][shl] = value

    if 'J' not in chain_resi_value:
        chain_resi_value['J'] = {}
    chain_resi_value['J'][-shl] = value

# make spectrum plot
#pymol.cmd.log_open("test.pml")
NCPandChd1.stylish()
NCPandChd1.spectrum(chain_resi_value, color_list=['white', 'red'], min=-0.5, max=3)
#NCPandChd1.spectrum(chain_resi_value, min=-0.5, max=3)
#pymol.cmd.log_close()

pymol.cmd.hide("cartoon", "chain A")
pymol.cmd.hide("cartoon", "chain B")
pymol.cmd.hide("cartoon", "chain C")
pymol.cmd.hide("cartoon", "chain D")
pymol.cmd.hide("cartoon", "chain E")
pymol.cmd.hide("cartoon", "chain F")
pymol.cmd.hide("cartoon", "chain G")
pymol.cmd.hide("cartoon", "chain H")
pymol.cmd.hide("cartoon", "chain W")



pymol.cmd.save("test.pse")
pymol.cmd.quit()
#pymol.finish_launching()
#pymol.cmd.hide("all")
#pymol.cmd.show("")



# Load Structures
#pymol.cmd.fetch("5O9G")

# Display the information of structures
#for x in pymol.cmd.get_names():
#    print ('Protein', x)
#    for ch in pymol.cmd.get_chains(x):
#        print (x, "has chain ", ch)





#pymol.cmd.hide("all")
#pymol.cmd.show("")
#pymol.cmd.disable("all")
#pymol.cmd.enable(sname)
#pymol.cmd.png("my_image.png")

# Get out!
#pymol.cmd.quit()
"""

# read pdb file 
def read_pdb (fname):
    chain_resi_resn = {}
    chain_resi_index_atom = {}
    for line in open(fname):
        if not line.startswith('ATOM'):
            continue
        line = line.strip()

        index = line[6:11].strip()
        atomn = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21].strip()
        resi = line[22:26].strip()
        x, y, z = line[30:38].strip(), line[38:46].strip(), line[46:54].strip()
        b_factor = line[60:66].strip()

        #_, index, atomn, resn, chain, resi, x, y, z, _, b_factor, _ = cols

        index = int(index)
        resi = int(resi)
        x, y, z = float(x), float(y), float(z)
        b_factor = float(b_factor)
        
        if chain not in chain_resi_resn:
            chain_resi_resn[chain] = {}
        chain_resi_resn[chain][resi] = resn

        if chain not in chain_resi_index_atom:
            chain_resi_index_atom[chain] = {}
        if resi not in chain_resi_index_atom[chain]:
            chain_resi_index_atom[chain][resi] = {}
        chain_resi_index_atom[chain][resi][index] = {"name":atomn,
                                                     "resn":resn,
                                                     "coord":(x, y, z),
                                                     "b-factor":b_factor}

    chain_seq = {}
    chain_type = {}
    for chain in chain_resi_resn:
        seq = ""
        rtype_count = {}
        for resi, resn in sorted(chain_resi_resn[chain].items()):
            if resn in aa_dict:
                new_resn = aa_dict[resn]
                rtype = 'protein'
            elif resn.startswith('D'):
                assert resn in ['DA', 'DT', 'DC', 'DG', 'DN']
                new_resn = resn[1:]
                rtype = 'DNA'
            else:
                assert resn in 'ATCGN'
                new_resn = resn
                rtype = 'RNA'
                
            seq += new_resn

            if rtype not in rtype_count:
                rtype_count[rtype] = 0
            rtype_count[rtype] +=1

        chain_seq[chain] = seq
        
        assert len(seq) == sum(rtype_count.values())
        ctype = sorted([(count, rtype) for rtype, count in rtype_count.items()], reverse=True)[0][1]
        chain_type[chain] = ctype

    return chain_resi_resn, chain_resi_index_atom, chain_seq, chain_type

def spectrum(chain_resi_value,
             cmap='jet',
             vmin=None,
             vmax=None,
             cbar=False):

    if vmin==None or vmax==None or cbar==True:
        values = []
        for chain in chain_resi_value:
            values += chain_resi_value[chain].values()
        if vmin==None:
            vmin = min(values)
        if vmax==None:
            vmax = max(values)

    chain_resi_RGB = {}
    cmap = cm.get_cmap(cmap)

    for chain in chain_resi_value:
        for resi in chain_resi_value[chain]:
            value = chain_resi_value[chain][resi]
            if value < vmin:
                value = vmin
            if value > vmax:
                value = vmax
                
            rescaled_value = analysis.rescale([value], vmin, vmax, 0.0, 1.0)[0]
            RGB = list(cmap(rescaled_value))[:3]

            if chain not in chain_resi_RGB:
                chain_resi_RGB[chain] = {}
            chain_resi_RGB[chain][resi] = RGB


    if cbar:
        fig = plt.figure()
        plt.imshow([values], vmin=vmin, vmax=vmax, cmap=cmap)
        plt.gca().set_visible(False)
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=15)
        plt.savefig("colorbar.png", bbox_inches='tight', dpi=1000)
        plt.close()
    
    return chain_resi_RGB

# set values for coloring
# for energy plot
if True:
    path = "/home/spark159/../../media/spark159/sw/slide_seq_data(2021.07.14)/"
    fname = 'insertionlib_dyad_size_shl_values2'
    note = "insertionlibtest"

    with open(path+fname + ".pickle", "rb") as f:
        dyad_size_shl_values = pickle.load(f)

    dyad = 112
    size = 1
    chain_resi_value = {}
    base_pairs = []
    values = []
    #for shl in range(5):
    for shl in dyad_size_shl_values[dyad][size]:
        value = np.mean(dyad_size_shl_values[dyad][size][shl])
        if 'I' not in chain_resi_value:
            chain_resi_value['I'] = {}
        chain_resi_value['I'][shl] = value
        if 'J' not in chain_resi_value:
            chain_resi_value['J'] = {}
        chain_resi_value['J'][-shl] = value

        base_pair = [('I', shl), ('J', -shl)]
        base_pairs.append(base_pair)
        values.append(value)

    vmin = min(values)
    vmax = max(values)


# for NMF plot
if False:
    path = "/home/spark159/../../media/spark159/sw/slide_seq_data(2021.07.14)/"
    fname = 'insertionlib_NMF_size_shl_meanweight1'
    note = "insertionlib_NMF"

    with open(path+fname + ".pickle", "rb") as f:
        size_shl_weight = pickle.load(f)

    size = 1
    basis = 2
    offset = 0
    chain_resi_value = {}
    base_pairs = []
    values = []
    for shl in size_shl_weight[size]:
        value = size_shl_weight[size][shl][basis]

        shl += offset # shift nucleosome location by offset
        
        if 'I' not in chain_resi_value:
            chain_resi_value['I'] = {}
        chain_resi_value['I'][shl] = value
        if 'J' not in chain_resi_value:
            chain_resi_value['J'] = {}
        chain_resi_value['J'][-shl] = value

        base_pair = [('I', shl), ('J', -shl)]
        base_pairs.append(base_pair)
        values.append(value)

    vmin = min(values)
    vmax = max(values)



### before sliding analysis
# set parameters
code = "6wz5"
#NA_color_list = ['red', 'green', 'blue']
NA_color_list = ['white', 'white']

# load pdb files
chain_resi_resn, chain_resi_index_atom, chain_seq, chain_type = read_pdb(code+".pdb")

#sys.exit(1)

# set RGB color by values
protein_chains = []
NA_chains = []
for chain, type in chain_type.items():
    if type == 'protein':
        protein_chains.append(chain)
    elif type in ['DNA', 'RNA']:
        NA_chains.append(chain)

#chain_resi_value = {}
#for chain in NA_chains:
#    for resi in chain_resi_resn[chain]:
#        value = random.random()
#        if chain not in chain_resi_value:
#            chain_resi_value[chain] = {}
#        chain_resi_value[chain][resi] = value

#chain_resi_RGB = spectrum(chain_resi_value, vmin=-0.5, vmax=1.25, cmap='jet', cbar=True)
#chain_resi_RGB = spectrum(chain_resi_value, vmin=-1.2, vmax=1.5, cmap='jet', cbar=True)
chain_resi_RGB = spectrum(chain_resi_value, vmin=vmin, vmax=vmax, cmap='jet', cbar=True)

# start write pml file
f = open("_".join([code, note, "block.pml"]), 'w')

# load the structure
print >> f, "reinitialize"
print >> f, "fetch %s" % (code)
print >> f, "hide all"
print >> f, ""

# draw proteins
#print >> f, "create protein, chain " + "+".join(protein_chains)
#print >> f, "set cartoon_color, purpleblue, protein"
##print >> f, "set cartoon_transparency, 0.6889, protein"
#print >> f, "set cartoon_transparency, 0.6, protein"
#print >> f, "show cartoon, protein"
#print >> f, ""


for chain in protein_chains:
    print >> f, "set cartoon_color, green, chain %s" % (chain)
    #print >> f, "set cartoon_color, white, chain %s" % (chain)
    #print >> f, "set cartoon_transparency, 0.6889, protein"
    #print >> f, "set cartoon_transparency, 0.6, chain %s" % (chain)
    print >> f, "set cartoon_transparency, 0.6889, chain %s" % (chain)
    #print >> f, "set cartoon_transparency, 0.3, chain %s" % (chain)
    #print >> f, "show cartoon, chain %s" % (chain)
    print >> f, ""

    

# draw nucleic acid backbones
for chain, color in zip(NA_chains, NA_color_list):
    print >> f, "create na_%s, chain %s" % (chain, chain)
    print >> f, "set cartoon_nucleic_acid_color, %s, na_%s" % (color, chain)
    #print >> f, "set cartoon_transparency, 0.2, na_%s" % (chain)
    print >> f, "show cartoon, na_%s" % (chain)
    print >> f, ""

print >> f, "select C3_prime, name C3'"
print >> f, "show sphere, C3_prime"
print >> f, "set sphere_scale, 0.2, C3_prime"
print >> f, "color gray90, C3_prime"
print >> f, ""

print >> f, "set cartoon_ladder_mode, 1"
print >> f, "set cartoon_ladder_radius, 0.1"
#print >> f, "set cartoon_ladder_radius, 0.05"
print >> f, "set cartoon_ladder_color, black"
print >> f, ""

#print >> f, "set cartoon_tube_radius, 0.05"
#print >> f, "set cartoon_tube_radius, 0.16889"
print >> f, "set cartoon_tube_radius, 0.5"
print >> f, "set cartoon_nucleic_acid_mode, 1"
print >> f, ""

# draw block nucleotides
#for chain in chain_resi_RGB:
#    for resi in sorted(chain_resi_RGB[chain].keys()):
#        if resi < 0:
#            resi_string = '\\' + str(resi)
#        else:
#            resi_string = str(resi)
#        RGB = chain_resi_RGB[chain][resi]
#        RGB_string = "[" + " ".join([str(round(comp,5)) for comp in RGB]) + "]"
#        print >> f, "dssr_block (chain %s and resi %s), block_color=N %s | edge black" % (chain,
#                                                                                          resi_string,
#                                                                                          RGB_string)  
#    print >> f, ""

for base_pair in base_pairs:
    chain1, resi1 = base_pair[0]
    chain2, resi2 = base_pair[1]
    RGB = chain_resi_RGB[chain1][resi1]
    RGB_string = "[" + " ".join([str(round(comp,5)) for comp in RGB]) + "]"
    
    selects = []
    for chain, resi in [[chain1, resi1], [chain2, resi2]]:
        if resi < 0:
            resi_string = '\\' + str(resi)
        else:
            resi_string = str(resi)
        select = "(chain %s and resi %s)" % (chain, resi_string) 
        selects.append(select)
    selects = " or ".join(selects)

    print >> f, "dssr_block %s, block_depth=1.2, block_color=N %s | edge black"  % (selects, RGB_string)
    #print >> f, "dssr_block %s, block_file=wc" % (selects) 
    #print >> f, "dssr_block %s, block_color=N %s | edge black" % (selects, RGB_string)
    #print >> f, "dssr_block %s, block_color=wc %s | edge black" % (selects, RGB_string)
    

# setup background and other details
print >> f, "set cartoon_highlight_color, grey50"
print >> f, "bg_color white"
print >> f, "remove solvent"
print >> f, "hide everything, hydro"
print >> f, ""

print >> f, "util.cbaw"
print >> f, "set sphere_quality, 4"
print >> f, "set stick_quality, 16"
print >> f, ""

print >> f, "set depth_cue, 0"
print >> f, "set ray_trace_fog, 0"
#print >> f, "set depth_cue, 1"
#print >> f, "set ray_trace_fog, 1"
print >> f, ""

print >> f, "set ray_shadow, off"
print >> f, "set orthoscopic, 1"
print >> f, ""

print >> f, "set antialias, 5"
#print >> f, "set antialias, 10"
#print >> f, "set antialias, 1"
print >> f, "set valence, 0"
print >> f, ""

print >> f, "set ambient, 0.68"
print >> f, "set reflect, 0"
print >> f, "set direct, 0.6"
print >> f, "set spec_direct, 0"
print >> f, "set light_count, 1"
print >> f, ""

# set view
#s = "set_view (\
#     0.324189186,    0.284314483,   -0.902258039,\
#     0.944719493,   -0.047908865,    0.324350089,\
#     0.048991807,   -0.957535326,   -0.284130692,\
#     0.000000000,    0.000000000, -581.525085449,\
#   153.774719238,  163.848922729,  165.090423584,\
#   459.467102051,  703.583435059,   20.000000000 )"
#s = "set_view (\
#     0.163609058,    0.126180306,   -0.978425682,\
#     0.966945112,    0.176092654,    0.184399977,\
#     0.195561796,   -0.976253808,   -0.093200102,\
#    -0.001033746,   -0.000534266, -457.852844238,\
#   156.874664307,  172.335540771,  159.350051880,\
#   335.751647949,  579.868103027,   20.000000000 )"

#s = "set_view (\
#    -0.786984861,    0.235012606,   -0.570466459,\
#     0.570210278,   -0.076067306,   -0.817970514,\
#    -0.235625342,   -0.969010115,   -0.074144654,\
#    -0.000735998,    0.000361085, -443.076202393,\
#   140.194595337,  176.541549683,  160.670333862,\
#   320.966125488,  565.082702637,   20.000000000 )"

#s = "set_view (\
#    -0.786379397,    0.178082779,   -0.591525733,\
#     0.577834487,   -0.126564890,   -0.806282520,\
#    -0.218449622,   -0.975840926,   -0.003376609,\
#    -0.000735998,    0.000361085, -386.940155029,\
#   140.194595337,  176.541549683,  160.670333862,\
#   264.830078125,  508.946716309,   20.000000000 )"

#s = "set_view (\
#    -0.856835425,    0.237675741,   -0.457547367,\
#     0.465591997,   -0.024561757,   -0.884660900,\
#    -0.221499100,   -0.971031368,   -0.089616254,\
#    -0.000735998,    0.000361085, -352.970520020,\
#   140.194595337,  176.541549683,  160.670333862,\
#   240.873413086,  464.964233398,   20.000000000 )"

#s = "set_view (\
#    -0.639455795,   -0.018048557,    0.768615782,\
#    -0.768812954,    0.009000443,   -0.639408052,\
#     0.004623585,   -0.999795020,   -0.019631024,\
#     0.000777408,    0.001110421, -328.629608154,\
#   254.589477539,  254.571029663,  248.570312500,\
#   259.106475830,  398.184173584,   20.000000000 )"

s = "set_view (\
    -0.639631331,   -0.010068387,    0.768615782,\
    -0.768640876,    0.018592561,   -0.639408052,\
    -0.007851653,   -0.999774754,   -0.019631024,\
     0.000777408,    0.001110421, -328.629608154,\
   254.589477539,  254.571029663,  248.570312500,\
   259.106475830,  398.184173584,   20.000000000 )"

print >> f, s

# save the image
print >> f, "set ray_trace_mode, 3"
print >> f, "set ray_trace_disco_factor, 1"

#print >> f, "ray 1000"
print >> f, "ray"
print >> f, "png %s_%s_block.png" % (code, note)
print >> f, ""

f.close()





"""
### after sliding analysis
# set parameters
code = "5O9G"
#NA_color_list = ['red', 'green', 'blue']
NA_color_list = ['white', 'white']

# load pdb files
chain_resi_resn, chain_resi_index_atom, chain_seq, chain_type = read_pdb(code+".pdb")

# set RGB color by values
protein_chains = []
NA_chains = []
for chain, type in chain_type.items():
    if type == 'protein':
        protein_chains.append(chain)
    elif type in ['DNA', 'RNA']:
        NA_chains.append(chain)

#chain_resi_value = {}
#for chain in NA_chains:
#    for resi in chain_resi_resn[chain]:
#        value = random.random()
#        if chain not in chain_resi_value:
#            chain_resi_value[chain] = {}
#        chain_resi_value[chain][resi] = value

#chain_resi_RGB = spectrum(chain_resi_value, vmin=-0.5, vmax=1.25, cmap='jet', cbar=True)
#chain_resi_RGB = spectrum(chain_resi_value, vmin=-1.2, vmax=1.5, cmap='jet', cbar=True)
chain_resi_RGB = spectrum(chain_resi_value, vmin=vmin, vmax=vmax, cmap='jet', cbar=True)

# start write pml file
f = open("_".join([code, note, "block.pml"]), 'w')

# load the structure
print >> f, "reinitialize"
print >> f, "fetch %s" % (code)
print >> f, "hide all"
print >> f, ""

# draw proteins
#print >> f, "create protein, chain " + "+".join(protein_chains)
#print >> f, "set cartoon_color, purpleblue, protein"
##print >> f, "set cartoon_transparency, 0.6889, protein"
#print >> f, "set cartoon_transparency, 0.6, protein"
#print >> f, "show cartoon, protein"
#print >> f, ""


#for chain in protein_chains:
#    if chain == 'W':
#        print >> f, "set cartoon_color, purpleblue, chain %s" % (chain)
#    else:
#        #print >> f, "set cartoon_color, green, chain %s" % (chain)
#        print >> f, "set cartoon_color, white, chain %s" % (chain)
#    #print >> f, "set cartoon_transparency, 0.6889, protein"
#    #print >> f, "set cartoon_transparency, 0.6, chain %s" % (chain)
#    print >> f, "set cartoon_transparency, 0.6889, chain %s" % (chain)
#    #print >> f, "set cartoon_transparency, 0.3, chain %s" % (chain)
#    print >> f, "show cartoon, chain %s" % (chain)
#    print >> f, ""

# draw Chd1 surface
for chain in protein_chains:
    if chain == 'W':
        print >> f, "set surface_color, lightpink, chain %s" % (chain)
        print >> f, "set surface_transparency, 0.6889, chain %s" % (chain)
        print >> f, "show surface, chain %s and resi 377-871" % (chain)
        #print >> f, "set cartoon_color, purple, chain %s" % (chain)
        #print >> f, "set cartoon_transparency, 0.6889, chain %s" % (chain)
        #print >> f, "show cartoon, chain %s and resi 377-871" % (chain)
    print >> f, "set surface_proximity, off"
    print >> f, "set surface_smooth_edges, on"
    print >> f, ""

    

# draw nucleic acid backbones
for chain, color in zip(NA_chains, NA_color_list):
    print >> f, "create na_%s, chain %s" % (chain, chain)
    print >> f, "set cartoon_nucleic_acid_color, %s, na_%s" % (color, chain)
    #print >> f, "set cartoon_transparency, 0.2, na_%s" % (chain)
    print >> f, "show cartoon, na_%s" % (chain)
    print >> f, ""

print >> f, "select C3_prime, name C3'"
print >> f, "show sphere, C3_prime"
print >> f, "set sphere_scale, 0.2, C3_prime"
print >> f, "color gray90, C3_prime"
print >> f, ""

print >> f, "set cartoon_ladder_mode, 1"
print >> f, "set cartoon_ladder_radius, 0.1"
#print >> f, "set cartoon_ladder_radius, 0.05"
print >> f, "set cartoon_ladder_color, black"
print >> f, ""

#print >> f, "set cartoon_tube_radius, 0.05"
#print >> f, "set cartoon_tube_radius, 0.16889"
print >> f, "set cartoon_tube_radius, 0.5"
print >> f, "set cartoon_nucleic_acid_mode, 1"
print >> f, ""

# draw block nucleotides
#for chain in chain_resi_RGB:
#    for resi in sorted(chain_resi_RGB[chain].keys()):
#        if resi < 0:
#            resi_string = '\\' + str(resi)
#        else:
#            resi_string = str(resi)
#        RGB = chain_resi_RGB[chain][resi]
#        RGB_string = "[" + " ".join([str(round(comp,5)) for comp in RGB]) + "]"
#        print >> f, "dssr_block (chain %s and resi %s), block_color=N %s | edge black" % (chain,
#                                                                                          resi_string,
#                                                                                          RGB_string)  
#    print >> f, ""

for base_pair in base_pairs:
    chain1, resi1 = base_pair[0]
    chain2, resi2 = base_pair[1]
    RGB = chain_resi_RGB[chain1][resi1]
    RGB_string = "[" + " ".join([str(round(comp,5)) for comp in RGB]) + "]"
    
    selects = []
    for chain, resi in [[chain1, resi1], [chain2, resi2]]:
        if resi < 0:
            resi_string = '\\' + str(resi)
        else:
            resi_string = str(resi)
        select = "(chain %s and resi %s)" % (chain, resi_string) 
        selects.append(select)
    selects = " or ".join(selects)

    print >> f, "dssr_block %s, block_depth=1.2, block_color=N %s | edge black"  % (selects, RGB_string)
    #print >> f, "dssr_block %s, block_file=wc" % (selects) 
    #print >> f, "dssr_block %s, block_color=N %s | edge black" % (selects, RGB_string)
    #print >> f, "dssr_block %s, block_color=wc %s | edge black" % (selects, RGB_string)
    

# setup background and other details
print >> f, "set cartoon_highlight_color, grey50"
print >> f, "bg_color white"
print >> f, "remove solvent"
print >> f, "hide everything, hydro"
print >> f, ""

print >> f, "util.cbaw"
print >> f, "set sphere_quality, 4"
print >> f, "set stick_quality, 16"
print >> f, ""

#print >> f, "set depth_cue, 0"
#print >> f, "set ray_trace_fog, 0"
print >> f, "set depth_cue, 1"
print >> f, "set ray_trace_fog, 1"
print >> f, ""

print >> f, "set ray_shadow, off"
print >> f, "set orthoscopic, 1"
print >> f, ""

print >> f, "set antialias, 5"
#print >> f, "set antialias, 10"
#print >> f, "set antialias, 1"
print >> f, "set valence, 0"
print >> f, ""

print >> f, "set ambient, 0.68"
print >> f, "set reflect, 0"
print >> f, "set direct, 0.6"
print >> f, "set spec_direct, 0"
print >> f, "set light_count, 1"
print >> f, ""

# set view
#s = "set_view (\
#     0.324189186,    0.284314483,   -0.902258039,\
#     0.944719493,   -0.047908865,    0.324350089,\
#     0.048991807,   -0.957535326,   -0.284130692,\
#     0.000000000,    0.000000000, -581.525085449,\
#   153.774719238,  163.848922729,  165.090423584,\
#   459.467102051,  703.583435059,   20.000000000 )"
#s = "set_view (\
#     0.163609058,    0.126180306,   -0.978425682,\
#     0.966945112,    0.176092654,    0.184399977,\
#     0.195561796,   -0.976253808,   -0.093200102,\
#    -0.001033746,   -0.000534266, -457.852844238,\
#   156.874664307,  172.335540771,  159.350051880,\
#   335.751647949,  579.868103027,   20.000000000 )"

#s = "set_view (\
#    -0.786984861,    0.235012606,   -0.570466459,\
#     0.570210278,   -0.076067306,   -0.817970514,\
#    -0.235625342,   -0.969010115,   -0.074144654,\
#    -0.000735998,    0.000361085, -443.076202393,\
#   140.194595337,  176.541549683,  160.670333862,\
#   320.966125488,  565.082702637,   20.000000000 )"

s = "set_view (\
    -0.786379397,    0.178082779,   -0.591525733,\
     0.577834487,   -0.126564890,   -0.806282520,\
    -0.218449622,   -0.975840926,   -0.003376609,\
    -0.000735998,    0.000361085, -386.940155029,\
   140.194595337,  176.541549683,  160.670333862,\
   264.830078125,  508.946716309,   20.000000000 )"

#s = "set_view (\
#    -0.856835425,    0.237675741,   -0.457547367,\
#     0.465591997,   -0.024561757,   -0.884660900,\
#    -0.221499100,   -0.971031368,   -0.089616254,\
#    -0.000735998,    0.000361085, -352.970520020,\
#   140.194595337,  176.541549683,  160.670333862,\
#   240.873413086,  464.964233398,   20.000000000 )"

print >> f, s

# save the image
print >> f, "set ray_trace_mode, 3"
print >> f, "set ray_trace_disco_factor, 1"

#print >> f, "ray 1000"
print >> f, "ray"
print >> f, "png %s_%s_block.png" % (code, note)
print >> f, ""

f.close()
"""