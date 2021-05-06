import sys
import math
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

fname = "IDlib_bubble_5_1rep_energy_wrt_601"
with open(fname + ".pickle", "rb") as f:
    size_dyad_shl_values = pickle.load(f)

#dyad_changes = {}
#for size in size_dyad_shl_values:
#    for dyad in size_dyad_shl_values[size]:
#        for shl in size_dyad_shl_values[size][dyad]:
#            values = size_dyad_shl_values[size][dyad][shl]
#            if dyad not in dyad_changes:
#                dyad_changes[dyad] = []
#            dyad_changes[dyad] += values

#var_dyad = sorted([(np.var(changes), dyad) for dyad, changes in dyad_changes.items()], reverse=True)


#dyad_list = [ dyad for var, dyad in var_dyad[:4] ]
#dyad_list = [112]
dyad_list = [92, 101, 112]

NCPlen = 147
    
shl_list = range(-(NCPlen/2), NCPlen/2 + 1)
    
f = open(fname + '.txt', 'w')
print >> f, '>Position with respect to dyad'
print >> f, ','.join([str(shl) for shl in shl_list])
for size in sorted(size_dyad_shl_values.keys()):
    dyad_shl_values = size_dyad_shl_values[size]
    for dyad in dyad_list:
        shl_values = dyad_shl_values[dyad]
        profile = []
        for shl in shl_list:
            try:
                value = np.mean(shl_values[shl])
            except:
                value = np.NaN
            profile.append(value)
        fig = plt.figure()
        plt.title(dyad)
        plt.plot(profile, 'o-')
        plt.show()
        plt.close()
        print >> f, '>perturbation-size:%s, dyad-location:%s' % (size, dyad)
        print >> f, ','.join([str(value) for value in profile])

f.close()
        
    
