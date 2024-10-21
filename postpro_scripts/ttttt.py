import math
import matplotlib.pyplot as plt

# periodicity of DNA in nucleosome
L = 10
# stiffness coefficient matrix
K = {'TA':12, 'TG':6, 'TT':68, 'TC':60, 'CA':6, 'CG':12, 'CT':35, 'CC':34, 'AA':68, 'AG':35, 'AT':48, 'AC':44, 'GA':60, 'GG':34, 'GT':44, 'GC':58}
# Nucleosomal DNA length
NCPlen = 147

def omega (index, phi=0):
    #phi = math.pi*float(phi)/180.0
    angle = (math.pi / L)*index - 0.5*(phi - 1.5)
    return 9*(math.sin(angle))**2

X = range(1,21)
Y = [omega(x, phi=290) for x in X]

fig = plt.figure()
plt.plot(X, Y)
plt.show()
