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

# mapping color according to the value based on the colormap
def spectrum(chain_resi_value,
             cmap='jet',
             vmin=None,
             vmax=None,
             cbar=False,
             note=""):

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
        fig = plt.figure(figsize=(2.5,1))
        plt.imshow([values], vmin=vmin, vmax=vmax, cmap=cmap)
        plt.gca().set_visible(False)
        cbar = plt.colorbar(ticks=[min(values), max(values)])
        cbar.ax.set_yticklabels([str(round(min(values), 2)), str(round(max(values), 2))], fontsize=5)
        plt.savefig("colorbar_%s.svg" % (note), format='svg', bbox_inches='tight')
        plt.close()
    
    return chain_resi_RGB

### parameters
# perturbation size
#size_choice = 2
#size_choice = 'mean'

# other parameters for NMF data type
#dyad_loc = None # dyad location in the DNA template

# other parameters for NMF data type
#dyad_shift = 40 # shift dyad location by offset
#basis_choice = 1 # NMF basis
#flip = False # flip the SHL coordinate

# input data (.pickle)
path = ""
pickle_fname = "Park97IDlib_NMF_size_shl_meanweight2.pickle"

# data type (NMF or energy landscape)
dtype = "NMF"

# PDB code
#code = "6wz5" # nucleosome without Chd1
code = "7tn2" # nucleosome with Chd1 (Bowman lab)

# basis name (NMF data)
#basis_name = {0:'hold', 1:'right-skew', 2:'partial-hold', 3:'left-skew'}
#basis_name = {0:'hold', 1:'sym', 2:'1-step', 3:'noisy'}
#basis_name = {0:'hold', 1:'noisy'}
basis_name = {0:'zero', 1:'-20', 2:'+20', 3:'-10', 4:'+10'}


# auto pymol run and quit
auto = True

# load pdb files
chain_resi_resn, chain_resi_index_atom, chain_seq, chain_type = read_pdb(code+".pdb")

# set SHL range
I_shl_st, I_shl_ed = min(chain_resi_resn['I'].keys()), max(chain_resi_resn['I'].keys())+1
J_shl_st, J_shl_ed = min(chain_resi_resn['J'].keys()), max(chain_resi_resn['J'].keys())+1
shl_st, shl_ed = max([I_shl_st, J_shl_st]), min([I_shl_ed, J_shl_ed])

# set figure list (size, basis, dyad, flip)
#fig_list = [('mean', 0, 0, False),
#            ('mean', 1, 0, False),
#            ('mean', 1, 40, False),
#            ('mean', 3, 0, False),
#            ('mean', 3, -20, False),
#            ('mean', 2, 0, False),
#            ('mean', 2, 30, False)]

#fig_list = [('mean', 0, 0, False),
#            ('mean', 2, 0, False),
#            ('mean', 2, 10, False),
#            ('mean', 1, 0, False),
#            ('mean', 1, -20, False),
#            ('mean', 1, 20, False),
#            ('mean', 3, 0, False),
#            ('mean', 3, -20, False),
#            ('mean', 3, 30, False)]

#fig_list = [('mean', 3, 0, False),
#            ('mean', 3, -10, False),
#            ('mean', 3, -20, False),
#            ('mean', 3, 30, False),
#            ('mean', 3, 40, False)]

#fig_list = [('mean', 0, 0, False),
#            ('mean', 1, 0, False)]

fig_list = [('mean', 0, 0, False),
            ('mean', 1, 0, False),
            ('mean', 1, -10, False),
            ('mean', 1, -20, False),
            ('mean', 2, 0, False),
            ('mean', 2, 10, False),
            ('mean', 2, 20, False),
            ('mean', 3, 0, False),
            ('mean', 3, -10, False),
            ('mean', 4, 0, False),
            ('mean', 4, 10, False)]


#### make pymol a pml file to plot the heatmap projcted on the structure
for size_choice, basis_choice, dyad_shift, flip in fig_list:

    ### set values for coloring
    # for energy plot
    if dtype == 'energy':

        with open(path + pickle_fname, "rb") as f:
            dyad_size_shl_values = pickle.load(f)

        chain_resi_value = {}
        base_pairs = []
        values = []

        shl_values = dyad_size_shl_values[dyad_choice][size_choice]

        #for shl in range(5):
        for shl in shl_values:
            value = np.mean(shl_values[shl])
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
    if dtype == 'NMF':

        with open(path + pickle_fname, "rb") as f:
            size_shl_weight = pickle.load(f)

        # average over all size and add it to the data
        shl_basis_weights = {}
        for size in size_shl_weight:
            for shl in size_shl_weight[size]:
                for basis in range(len(size_shl_weight[size][shl])):
                    if shl not in shl_basis_weights:
                        shl_basis_weights[shl] = []
                    if len(shl_basis_weights[shl]) <= basis:
                        shl_basis_weights[shl].append([])
                    weight = size_shl_weight[size][shl][basis]
                    shl_basis_weights[shl][basis].append(weight)

        size_shl_weight['mean'] = {}
        for shl in shl_basis_weights:
            for basis in range(len(shl_basis_weights[shl])):
                weights = shl_basis_weights[shl][basis]
                if shl not in size_shl_weight['mean']:
                    size_shl_weight['mean'][shl] = []
                size_shl_weight['mean'][shl].append(np.mean(weights))

        chain_resi_value = {}
        base_pairs = []
        values = []
        missing_data = []

        #if shl_sym:
        #    lim = max([abs(shl) for shl in size_shl_weight[size_choice]])
        #    shl_st, shl_ed = -lim, lim+1
        #else:
        #    shl_st, shl_ed = min(size_shl_weight[size]), max(size_shl_weight[size])+1

        #for shl in range(-lim, lim+1):
        #for shl in range(min(size_shl_weight[size]), max(size_shl_weight[size])+1):
        #for shl in size_shl_weight[size]:
        for shl in range(shl_st, shl_ed):
            oldshl = shl + dyad_shift
            if flip:
                oldshl = -oldshl
            try:
                value = size_shl_weight[size_choice][oldshl][basis_choice]
            except:
                missing_data.append(('I', shl))
                missing_data.append(('J', -shl))
                continue

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


    # note
    if dtype == 'NMF':
        try:
            bname = basis_name[basis_choice]
        except:
            bname = basis_choice
        note = 'NMF_%sbp_basis-%s_dyad-%d' % (size_choice, bname, dyad_shift)
        if flip:
            note += '_flipped'

    elif dtype == 'energy':
        note = 'Energy_%sbp_dyad-%d' % (size_choice, dyad_choice)


    ## nucleosome strcuture without Chd1
    if code == "6wz5":
        # set parameters
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
        chain_resi_RGB = spectrum(chain_resi_value, vmin=vmin, vmax=vmax, cmap='jet', cbar=True,
                                  note = '_'.join([code, note]))

        # start write pml file
        pml_fname = "_".join([code, note, "block.pml"])
        f = open(pml_fname, 'w')

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
            #print >> f, "set cartoon_transparency, 0.6889, chain %s" % (chain)
            #print >> f, "set cartoon_transparency, 0.3, chain %s" % (chain)
            #print >> f, "show cartoon, chain %s" % (chain)
            print >> f, ""

            #print >> f, "set surface_color, lightpink, chain %s" % (chain)
            #print >> f, "set surface_color, white, chain %s" % (chain)
            #print >> f, "set surface_color, aquamarine, chain %s" % (chain)
            #print >> f, "set surface_transparency, 0.6889, chain %s" % (chain)
            #print >> f, "show surface, chain %s" % (chain)
            #print >> f, "show surface, chain %s and resi 377-871" % (chain)
            #print >> f, "set cartoon_color, purple, chain %s" % (chain)
            #print >> f, "set cartoon_color, lightblue, chain %s" % (chain)
            #print >> f, "set cartoon_transparency, 0.6889, chain %s" % (chain)
            #print >> f, "show cartoon, chain %s and resi 377-871" % (chain)
            #print >> f, "set surface_proximity, off"
            #print >> f, "set surface_smooth_edges, on"
            #print >> f, ""



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


        # fill missing data with white color
        selects = []
        for chain, resi in missing_data:
            try:
                chain_resi_resn[chain][resi]
            except:
                continue
            if resi < 0:
                resi_string = '\\' + str(resi)
            else:
                resi_string = str(resi)
            select = "(chain %s and resi %s)" % (chain, resi_string) 
            selects.append(select)
        selects = " or ".join(selects)
        print >> f, "dssr_block %s, block_depth=1.3, block_color=N white | edge black"  % (selects)
        #print >> f, "dssr_block %s, block_depth=1.3, block_color=N [0.77 0.77 0.77] | edge black"  % (selects)


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
            try:
                chain_resi_resn[chain1][resi1]
                chain_resi_resn[chain2][resi2]
            except:
                continue
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
        #print >> f, "png %s_%s_block.png" % (code, note)
        print >> f, "png %s_%s_block.png, width=1.5 in, height=1.5 in, dpi=500, ray=1" % (code, note)
        if auto:
            print >> f, "quit"
        print >> f, ""

        f.close()


    ## nucleosome strcuture with Chd1 (Bowman lab)
    elif code == "7tn2":
        # set parameters
        #code = "6g0l"
        #NA_color_list = ['red', 'green', 'blue']
        NA_color_list = ['white', 'white']
        #NA_color_list = ['black', 'black']
        #NA_color_list = ['deepblue', 'deepblue']
        #NA_color_list = ['density', 'density']

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
        chain_resi_RGB = spectrum(chain_resi_value, vmin=vmin, vmax=vmax, cmap='jet', cbar=True,
                                  note = '_'.join([code, note]))

        # start write pml file
        pml_fname = "_".join([code, note, "block.pml"])
        f = open(pml_fname, 'w')

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
                #print >> f, "set surface_color, lightpink, chain %s" % (chain)
                #print >> f, "set surface_color, white, chain %s" % (chain)
                #print >> f, "set surface_color, aquamarine, chain %s" % (chain)
                #print >> f, "set surface_transparency, 0.6889, chain %s" % (chain)
                #print >> f, "show surface, chain %s and resi 377-871" % (chain)
                #print >> f, "set cartoon_color, purple, chain %s" % (chain)
                #print >> f, "set cartoon_color, white, chain %s" % (chain)
                #print >> f, "set cartoon_transparency, 0.1, chain %s" % (chain)
                #print >> f, "show cartoon, chain %s and resi 377-871" % (chain)
                continue
            print >> f, "set cartoon_color, green, chain %s" % (chain)
            print >> f, "set cartoon_transparency, 0.3, chain %s" % (chain)
            print >> f, "show cartoon, chain %s" % (chain)
            #print >> f, "set surface_proximity, off"
            #print >> f, "set surface_smooth_edges, on"
            print >> f, ""


        if True:
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

            # fill missing data with white color
            selects = []
            for chain, resi in missing_data:
                if resi < 0:
                    resi_string = '\\' + str(resi)
                else:
                    resi_string = str(resi)
                select = "(chain %s and resi %s)" % (chain, resi_string) 
                selects.append(select)
            selects = " or ".join(selects)
            print >> f, "dssr_block %s, block_depth=1.3, block_color=N white | edge black"  % (selects)
            #print >> f, "dssr_block %s, block_depth=1.3, block_color=N [0.77 0.77 0.77] | edge black"  % (selects)

            # fill valid data with defined colormap
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
                #print >> f, "dssr_block %s, block_depth=1.2, block_color=N gray | edge black"  % (selects)
                #print >> f, "dssr_block %s, block_file=wc" % (selects) 
                #print >> f, "dssr_block %s, block_color=N %s | edge black" % (selects, RGB_string)
                #print >> f, "dssr_block %s, block_color=wc %s | edge black" % (selects, RGB_string)

        if False:
            # draw nucleic acid backbones
            for chain, color in zip(NA_chains, NA_color_list):
                resi_range = []
                for resi in [min(chain_resi_value[chain].keys())-5, max(chain_resi_value[chain].keys())+5]:
                    if resi < 0:
                        resi_range.append('\\' + str(resi))
                    else:
                        resi_range.append(str(resi))
                select = "chain %s and resi %s-%s" % (chain, resi_range[0], resi_range[1])

                print >> f, "create na_%s, %s" % (chain, select)
                print >> f, "set cartoon_nucleic_acid_color, %s, na_%s" % (color, chain)
                print >> f, "show cartoon, na_%s" % (chain)
                print >> f, ""

                print >> f, "select C3_prime, name C3 and %s" % (select)
                print >> f, "show sphere, C3_prime"
                print >> f, "set sphere_scale, 0.2, C3_prime"
                print >> f, "color gray90, C3_prime"
                print >> f, ""


            print >> f, "set cartoon_ladder_mode, 1"
            print >> f, "set cartoon_ladder_radius, 0.1"
            print >> f, "set cartoon_ladder_color, black"
            print >> f, ""

            print >> f, "set cartoon_tube_radius, 0.5"
            print >> f, "set cartoon_nucleic_acid_mode, 1"
            print >> f, ""

            # fill missing data with white color
            selects = []
            for chain, resi in missing_data:
                if resi < 0:
                    resi_string = '\\' + str(resi)
                else:
                    resi_string = str(resi)
                select = "(chain %s and resi %s)" % (chain, resi_string) 
                selects.append(select)
            selects = " or ".join(selects)
            print >> f, "dssr_block %s, block_depth=1.3, block_color=N white | edge black"  % (selects)
            #print >> f, "dssr_block %s, block_depth=1.3, block_color=N [0.77 0.77 0.77] | edge black"  % (selects)

            # fill valid data with defined colormap
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

                print >> f, "dssr_block %s, block_depth=1.3, block_color=N %s | edge black"  % (selects, RGB_string)
                #print >> f, "dssr_block %s, block_depth=1.3, block_color=N [0.77 0.77 0.77] | edge black"  % (selects)
                #print >> f, "dssr_block %s, block_depth=1.3, block_color=N white | edge black"  % (selects)



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
        s = "set_view (\
             0.671171546,    0.063488662,   -0.738575995,\
             0.255896538,   -0.954918563,    0.150458679,\
            -0.695731640,   -0.289983302,   -0.657160223,\
            -0.000061989,    0.002155930, -348.055267334,\
           155.531921387,  160.811904907,  162.158294678,\
           253.623901367,  443.240203857,   20.000000000 )"

        print >> f, s

        # save the image
        print >> f, "set ray_trace_mode, 3"
        print >> f, "set ray_trace_disco_factor, 1"

        #print >> f, "ray 1000"
        #print >> f, "ray 3000"
        print >> f, "ray"
        print >> f, "png %s_%s_block.png" % (code, note)
        #print >> f, "png %s_%s_block.png, width=1.5 in, height=1.5 in, dpi=500, ray=1" % (code, note)
        if auto:
            print >> f, "quit"
        print >> f, ""

        f.close()


    if auto:
        cmd = ['pymol', '-Q', pml_fname]
        subprocess.call(cmd)



sys.exit(1)

### after sliding analysis (one-bound Chd1)
# set parameters
code = "5O9G"
#code = "6g0l"
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
#print >> f, "load Chd1apo_temp.pdb"
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
        #print >> f, "set surface_color, lightpink, chain %s" % (chain)
        print >> f, "set surface_color, white, chain %s" % (chain)
        #print >> f, "set surface_color, aquamarine, chain %s" % (chain)
        #print >> f, "set surface_transparency, 0.6889, chain %s" % (chain)
        #print >> f, "show surface, chain %s and resi 377-871" % (chain)
        #print >> f, "set cartoon_color, purple, chain %s" % (chain)
        #print >> f, "set cartoon_color, lightblue, chain %s" % (chain)
        #print >> f, "set cartoon_transparency, 0.6889, chain %s" % (chain)
        #print >> f, "show cartoon, chain %s and resi 377-871" % (chain)
    print >> f, "set surface_proximity, off"
    print >> f, "set surface_smooth_edges, on"
    print >> f, ""

    
if False:
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

if True:
    # draw nucleic acid backbones
    for chain, color in zip(NA_chains, NA_color_list):
        resi_range = []
        for resi in [min(chain_resi_value[chain].keys())-5, max(chain_resi_value[chain].keys())+5]:
            if resi < 0:
                resi_range.append('\\' + str(resi))
            else:
                resi_range.append(str(resi))
        select = "chain %s and resi %s-%s" % (chain, resi_range[0], resi_range[1])

        print >> f, "create na_%s, %s" % (chain, select)
        print >> f, "set cartoon_nucleic_acid_color, %s, na_%s" % (color, chain)
        print >> f, "show cartoon, na_%s" % (chain)
        print >> f, ""

        print >> f, "select C3_prime, name C3 and %s" % (select)
        print >> f, "show sphere, C3_prime"
        print >> f, "set sphere_scale, 0.2, C3_prime"
        print >> f, "color gray90, C3_prime"
        print >> f, ""


    print >> f, "set cartoon_ladder_mode, 1"
    print >> f, "set cartoon_ladder_radius, 0.1"
    print >> f, "set cartoon_ladder_color, black"
    print >> f, ""

    print >> f, "set cartoon_tube_radius, 0.5"
    print >> f, "set cartoon_nucleic_acid_mode, 1"
    print >> f, ""


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

#here
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
#print >> f, "ray 3000"
print >> f, "ray"
print >> f, "png %s_%s_block.png" % (code, note)
#print >> f, "png %s_%s_block.png, width=1.5 in, height=1.5 in, dpi=500, ray=1" % (code, note)
print >> f, ""

f.close()


### after sliding analysis (two-bound chd1)
# set parameters
#code = "5O9G"
code = "6g0l"
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

chain_resi_RGB = spectrum(chain_resi_value, vmin=vmin, vmax=vmax, cmap='jet', cbar=True)

# start write pml file
f = open("_".join([code, note, "block.pml"]), 'w')

# load the structure
print >> f, "reinitialize"
print >> f, "fetch %s" % (code)
print >> f, "hide all"
print >> f, ""


## draw Chd1 surface
#for chain in protein_chains:
#    if chain == 'W':
#        #print >> f, "set surface_color, lightpink, chain %s" % (chain)
#        #print >> f, "set surface_color, aquamarine, chain %s" % (chain)
#        #print >> f, "set transparency, 0.6889, chain %s" % (chain)
#        #print >> f, "show surface, chain %s and resi 377-871" % (chain)
#        #print >> f, "show surface, chain %s and resi 388-860" % (chain)
#
#        print >> f, "set cartoon_color, lightpink, chain %s" % (chain)
#        #print >> f, "set cartoon_color, aquamarine, chain %s" % (chain)
#        print >> f, "set transparency, 0.6889, chain %s" % (chain)
#        print >> f, "show cartoon, chain %s and resi 388-860" % (chain)
#
#    if chain == 'M':   
#        #print >> f, "set surface_color, lightpink, chain %s" % (chain)
#        #print >> f, "set surface_color, aquamarine, chain %s" % (chain)
#        #print >> f, "set transparency, 0.6889, chain %s" % (chain)
#        #print >> f, "show surface, chain %s and resi 377-871" % (chain)
#        #print >> f, "show surface, chain %s and resi 388-860" % (chain)
#
#        #print >> f, "set cartoon_color, lightpink, chain %s" % (chain)
#        print >> f, "set cartoon_color, aquamarine, chain %s" % (chain)
#        print >> f, "set transparency, 0.6889, chain %s" % (chain)
#        print >> f, "show cartoon, chain %s and resi 388-860" % (chain)
#
#    print >> f, "set surface_proximity, off"
#    print >> f, "set surface_smooth_edges, on"
#    print >> f, ""


# draw nucleic acid backbones
for chain, color in zip(NA_chains, NA_color_list):
    resi_range = []
    for resi in [min(chain_resi_value[chain].keys())-10, max(chain_resi_value[chain].keys())+10]:
        if resi < 0:
            resi_range.append('\\' + str(resi))
        else:
            resi_range.append(str(resi))
    select = "chain %s and resi %s-%s" % (chain, resi_range[0], resi_range[1])
    
    print >> f, "create na_%s, %s" % (chain, select)
    print >> f, "set cartoon_nucleic_acid_color, %s, na_%s" % (color, chain)
    print >> f, "show cartoon, na_%s" % (chain)
    print >> f, ""

    print >> f, "select C3_prime, name C3 and %s" % (select)
    print >> f, "show sphere, C3_prime"
    print >> f, "set sphere_scale, 0.2, C3_prime"
    print >> f, "color gray90, C3_prime"
    print >> f, ""

print >> f, "set cartoon_ladder_mode, 1"
print >> f, "set cartoon_ladder_radius, 0.1"
print >> f, "set cartoon_ladder_color, black"
print >> f, ""

print >> f, "set cartoon_tube_radius, 0.5"
print >> f, "set cartoon_nucleic_acid_mode, 1"
print >> f, ""


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

print >> f, "set depth_cue, 1"
print >> f, "set ray_trace_fog, 1"
print >> f, ""

print >> f, "set ray_shadow, off"
print >> f, "set orthoscopic, 1"
print >> f, ""

print >> f, "set antialias, 5"
print >> f, "set valence, 0"
print >> f, ""

print >> f, "set ambient, 0.68"
print >> f, "set reflect, 0"
print >> f, "set direct, 0.6"
print >> f, "set spec_direct, 0"
print >> f, "set light_count, 1"
print >> f, ""

# set view
s = "set_view (\
     0.738110662,   -0.309598595,    0.599450290,\
    -0.046806283,    0.862858117,    0.503274381,\
    -0.673053622,   -0.399529785,    0.622393906,\
     0.000000000,    0.000000000, -561.930480957,\
   183.504074097,  190.452285767,  185.815856934,\
   421.386383057,  702.474548340,   20.000000000 )"

#set_view (\
#     0.731664240,   -0.324979514,    0.599212110,\
#    -0.047133148,    0.852819264,    0.520074189,\
#    -0.680033147,   -0.408761948,    0.608660758,\
#     0.000000000,    0.000000000, -505.023406982,\
#   183.504074097,  190.452285767,  185.815856934,\
#   364.479309082,  645.567504883,   20.000000000 )

print >> f, s

# save the image
print >> f, "set ray_trace_mode, 3"
#print >> f, "set ray_trace_mode, 0"
print >> f, "set ray_trace_disco_factor, 1"

#print >> f, "ray 1000"
print >> f, "ray"
print >> f, "png %s_%s_block.png" % (code, note)
print >> f, ""

f.close()

