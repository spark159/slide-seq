import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
from scipy.special import betainc
from scipy.special import beta as betaf

def prior (x):
    assert x > 0
    output = 0.0
    j = 0
    while j <= beta -1:
        temp = ((-1)**j)*comb(beta-1,j)*(x**(2*alpha+j-1))
        if x <=1:
            temp *= betaf(((2*alpha+j+1)/float(alpha))+1, beta)
        else:
            temp *= betainc(((2*alpha+j+1)/float(alpha))+1, beta, (1.0/x)**alpha)*betaf(((2*alpha+j+1)/float(alpha))+1, beta)
        output += temp
        j +=1
    output *= alpha*beta*beta
    return output

alpha = 2
beta = 3

X = np.linspace(0.0001,0.9999, 1000)
Y = []
for x in X:
    Y.append(alpha*beta*(x**(alpha-1))*(1-x**alpha)**(beta-1))

fig = plt.figure()
plt.plot(X,Y)
plt.show()

X = np.linspace(0.0001,10, 1000)
Y = []
for x in X:
    Y.append(prior(x))

fig = plt.figure()
plt.plot(X, Y)
plt.show()
