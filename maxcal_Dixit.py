import sys
import numpy as np
import copy
import random
from scipy.optimize import fsolve
from scipy.optimize import newton_krylov
from scipy.optimize import anderson
from scipy.optimize import broyden1
from scipy.optimize import broyden2
from scipy.optimize import excitingmixing
import matplotlib
import matplotlib.pyplot as plt
#import load


def steady_maxcal (steady_dist, edge_matrix=[[None]], obs_matrix=[[None]], obs_average=None, iter_num=np.inf):
    def get_W_matrix (gamma):
        N = len(steady_dist)
        W_matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                W_matrix[i][j] = edge_matrix[i][j] * np.exp(-gamma*obs_matrix[i][j])
        return W_matrix
    def get_Beta_Lambda (W_matrix):
        def distance(vec1, vec2):
            return np.sqrt(np.sum(((vec1-vec2)/vec1)**2))
        def operator (vec):
            steady_dist_vec = np.asarray(steady_dist).reshape(vec.shape) 
            return steady_dist_vec / vec
        N = len(steady_dist)
        Lambda = np.ones((N,1))
        Beta = np.ones((N,1))
        i = 0
        while i < iter_num:
            Beta = operator(np.matmul(W_matrix, Lambda))
            new_Lambda = operator(np.matmul(W_matrix.transpose(), Beta))
            if distance(Lambda, new_Lambda) < 10**-10:
                return np.squeeze(Beta), np.squeeze(new_Lambda)
            Lambda = copy.deepcopy(new_Lambda)
            i +=1
        return np.squeeze(Beta), np.squeeze(new_Lambda)
    def get_Beta_Lambda_numeric (W_matrix): 
        def equations (Lambda):
            eqns = []
            for i in range(N):
                factor1 = 0.0
                for j in range(N):
                    factor2 = 0.0
                    for k in range(N):
                        factor2 += Lambda[k]*W_matrix[j][k]
                    factor1 += (steady_dist[j]/factor2) * W_matrix[j][i]
                eqns.append(steady_dist[i]/factor1 - Lambda[i])
            return eqns
        N = len(steady_dist)
        #Lambda = fsolve(equations, [1.0]*N)
        Lambda = newton_krylov(equations, [1.0]*N)
        Beta = []
        for i in range(N):
            factor = 0.0
            for j in range(N):
                factor += Lambda[j]* W_matrix[i][j]
            Beta.append(steady_dist[i]/factor)
        return Beta, Lambda
    def get_average(Beta, Lambda, W_matrix):
        N = len(steady_dist)
        average = 0.0
        for i in range(N):
            for j in range(N):
                average += Beta[i]*Lambda[j]*W_matrix[i][j]*obs_matrix[i][j]
        return average
    def get_steady_matrix (Beta, Lambda, W_matrix):
        N = len(steady_dist)
        steady_matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                steady_matrix[i][j] = Beta[i]*Lambda[j]*W_matrix[i][j]/steady_dist[i]
        return steady_matrix

    N = len(steady_dist)
    if np.any(np.isin(edge_matrix, None)):
        edge_matrix = np.ones((N,N))
        
    if np.any(np.isin(obs_matrix, None)) or obs_average == None:
        gamma = 0.0
        W_matrix = edge_matrix
        Beta, Lambda = get_Beta_Lambda(W_matrix)
        return get_steady_matrix (Beta, Lambda, W_matrix)
    
    gamma_guess = np.linspace(0.0, 100.0, num=1000)
    #gamma_guess = [100]
    diff = sys.maxint
    gamma_best, Beta_best, Lambda_best, W_best = None, None, None, None

    for k in range(len(gamma_guess)):
        gamma = gamma_guess[k]
        W_matrix = get_W_matrix(gamma)
        Beta, Lambda = get_Beta_Lambda(W_matrix)
        average = get_average(Beta, Lambda, W_matrix)
        #print "gamma:", gamma
        #print "W matrix:"
        #print W_matrix
        #print "Beta, Lambda:", Beta, Lambda
        #print "average:", average
        #print "steady matrix:"
        #print get_steady_matrix(Beta, Lambda, W_matrix)
        #print
        
        if abs(obs_average - average) < diff:
            diff = abs(obs_average - average)
            gamma_best = gamma
            Beta_best = Beta
            Lambda_best = Lambda
            W_best = W_matrix

    print "gamma best:", gamma_best
    print "real average:", obs_average
    print "best average:", get_average(Beta_best, Lambda_best, W_best)
    steady_matrix = get_steady_matrix(Beta_best, Lambda_best, W_best)
    return steady_matrix


def eqm_maxcal (eqm_dist, edge_matrix=[[None]], obs_matrix=[[None]], obs_average=None):
    def get_W_matrix (gamma):
        N = len(eqm_dist)
        W_matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                W_matrix[i][j] = edge_matrix[i][j] * np.exp(-0.5*gamma*(obs_matrix[i][j]+obs_matrix[j][i]))
        return W_matrix
    def get_Rho (W_matrix): # not use anymore
        def distance(vec1, vec2):
            return np.sqrt(np.sum(((vec1-vec2)/vec1)**2))
        def operator (vec):
            eqm_dist_vec = np.asarray(eqm_dist).reshape(vec.shape) 
            return eqm_dist_vec / vec
        N = len(eqm_dist)
        Rho = np.ones((N,1))
        i = 0
        while i < iter_num:
            new_Rho = operator(np.matmul(W_matrix, Rho))
            if distance(Rho, new_Rho) < 10**-10:
                return np.squeeze(new_Rho)
            Rho = copy.deepcopy(new_Rho)
            i +=1
        return np.squeeze(new_Rho)
    def get_Rho_numeric (W_matrix):
        def equations (Rho):
            eqns = []
            for i in range(N):
                factor = 0.0
                for j in range(N):
                    factor += W_matrix[i][j] * Rho[j]
                eqns.append(eqm_dist[i]/factor - Rho[i])
            return eqns
        N = len(eqm_dist)
        Rho = fsolve(equations, [1.0]*N)
        return Rho
    def get_average(Rho, W_matrix):
        N = len(eqm_dist)
        average = 0.0
        for i in range(N):
            for j in range(N):
                average += Rho[i]*Rho[j]*W_matrix[i][j]*obs_matrix[i][j]
        return average
    def get_eqm_matrix (Rho, W_matrix):
        N = len(eqm_dist)
        eqm_matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                eqm_matrix[i][j] = Rho[i]*Rho[j]*W_matrix[i][j]/eqm_dist[i]
        return eqm_matrix

    N = len(eqm_dist)
    if np.any(np.isin(edge_matrix, None)):
        edge_matrix = np.ones((N,N))

    if np.any(np.isin(obs_matrix, None)) or obs_average == None:
        gamma = 0.0
        W_matrix = edge_matrix
        Rho = get_Rho_numeric(W_matrix)
        return get_eqm_matrix (Rho, W_matrix)

    gamma_guess = np.linspace(0.0, 100.0, num=1000)
    #gamma_guess = [100]
    diff = sys.maxint
    gamma_best, Rho_best, W_best = None, None, None

    for k in range(len(gamma_guess)):
        gamma = gamma_guess[k]
        W_matrix = get_W_matrix(gamma)
        Rho = get_Rho_numeric(W_matrix)
        
        if obs_matrix == None or obs_average == None:
            return get_eqm_matrix(Rho, W_matrix)

        average = get_average(Rho, W_matrix)
        #print "gamma:", gamma
        #print "W matrix:"
        #print W_matrix
        #print "Rho:", Rho
        #print "average:", average
        #print "eqm matrix:"
        #print get_eqm_matrix(Rho, W_matrix)
        #print
        
        if abs(obs_average - average) < diff:
            diff = abs(obs_average - average)
            gamma_best = gamma
            Rho_best = Rho
            W_best = W_matrix

    print "gamma best:", gamma_best
    print "real average:", obs_average
    print "best average:", get_average(Rho_best, W_best)
    eqm_matrix = get_eqm_matrix(Rho_best, W_best)
    return eqm_matrix

#dist = [0.546]
#dist.append(1-dist[0])
#obs_matrix = [[1.0,1.0],[1.0,1.0]]
#obs_average = 1.0
#print eqm_maxcal(eqm_dist, obs_matrix, obs_average)
#print steady_maxcal(dist)
#print eqm_maxcal(dist)

def eqm_to_steady (eqm_matrix, steady_dist, edge_matrix=[[None]], obs_matrix=[[None]], obs_average=None, iter_num=np.inf):
    def get_W_matrix (gamma):
        N = len(steady_dist)
        W_matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                W_matrix[i][j] = edge_matrix[i][j] * eqm_matrix[i][j] *  np.exp(gamma*obs_matrix[i][j])
        return W_matrix
    def get_Tau_Lambda (W_matrix):
        def distance(vec1, vec2):
            return np.sqrt(np.sum(((vec1-vec2)/vec1)**2))
        def operator (vec):
            steady_dist_vec = np.asarray(steady_dist).reshape(vec.shape) 
            return steady_dist_vec / vec
        N = len(steady_dist)
        Lambda = np.ones((N,1))
        Tau = np.ones((N,1))
        i = 0
        while i < iter_num:
            Tau = operator(np.matmul(W_matrix, Lambda))
            new_Lambda = operator(np.matmul(W_matrix.transpose(), Tau))
            if distance(Lambda, new_Lambda) < 10**-10:
                return np.squeeze(Tau), np.squeeze(new_Lambda)
            Lambda = copy.deepcopy(new_Lambda)
            i +=1
        return np.squeeze(Tau), np.squeeze(new_Lambda)
    def get_Tau_Lambda_numeric (W_matrix): 
        def equations (Lambda):
            eqns = []
            for i in range(N):
                factor1 = 0.0
                for j in range(N):
                    factor2 = 0.0
                    for k in range(N):
                        factor2 += Lambda[k]*W_matrix[j][k]
                    factor1 += (steady_dist[j]/factor2) * W_matrix[j][i]
                eqns.append(steady_dist[i]/factor1 - Lambda[i])
            return eqns
        N = len(steady_dist)
        Lambda = fsolve(equations, [1.0]*N)
        #Lambda = newton_krylov(equations, [1.0]*N)
        Tau = []
        for i in range(N):
            factor = 0.0
            for j in range(N):
                factor += Lambda[j]* W_matrix[i][j]
            Tau.append(steady_dist[i]/factor)
        return Tau, Lambda
    def get_average(Tau, Lambda, W_matrix):
        N = len(steady_dist)
        average = 0.0
        for i in range(N):
            for j in range(N):
                average += Tau[i]*Lambda[j]*W_matrix[i][j]*obs_matrix[i][j]
        return average
    def get_steady_matrix (Tau, Lambda, W_matrix):
        N = len(steady_dist)
        steady_matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                steady_matrix[i][j] = Tau[i]*Lambda[j]*W_matrix[i][j]/steady_dist[i]
        return steady_matrix

    N = len(steady_dist)
    if np.any(np.isin(edge_matrix, None)):
        edge_matrix = np.ones((N,N))
        
    if np.any(np.isin(obs_matrix, None)) or obs_average == None:
        gamma = 0.0
        N = len(steady_dist)
        W_matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                W_matrix[i][j] = edge_matrix[i][j] * eqm_matrix[i][j]
        Tau, Lambda = get_Tau_Lambda(W_matrix)
        return get_steady_matrix(Tau, Lambda, W_matrix)
    
    gamma_guess = np.linspace(0.0, 100.0, num=1000)
    #gamma_guess = [100]
    diff = sys.maxint
    gamma_best, Tau_best, Lambda_best, W_best = None, None, None, None

    for k in range(len(gamma_guess)):
        gamma = gamma_guess[k]
        W_matrix = get_W_matrix(gamma)
        Tau, Lambda = get_Tau_Lambda(W_matrix)
        average = get_average(Tau, Lambda, W_matrix)
        #print "gamma:", gamma
        #print "W matrix:"
        #print W_matrix
        #print "Tau, Lambda:", Tau, Lambda
        #print "average:", average
        #print "steady matrix:"
        #print get_steady_matrix(Tau, Lambda, W_matrix)
        #print
        
        if abs(obs_average - average) < diff:
            diff = abs(obs_average - average)
            gamma_best = gamma
            Tau_best = Tau
            Lambda_best = Lambda
            W_best = W_matrix

    print "gamma best:", gamma_best
    print "real average:", obs_average
    print "best average:", get_average(Tau_best, Lambda_best, W_best)
    steady_matrix = get_steady_matrix(Tau_best, Lambda_best, W_best)
    return steady_matrix


def eqm_to_neqm (eqm_matrix, eqm_dist, neqm_dist, edge_matrix=[[None]], obs_matrix=[[None]], obs_average=None, iter_num=np.inf):
    def get_W_matrix (gamma):
        N = len(eqm_dist)
        W_matrix = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                W_matrix[i][j] = edge_matrix[i][j] * eqm_matrix[i][j]*np.exp(0.5*gamma*(obs_matrix[i][j]+obs_matrix[j][i]))/eqm_dist[j]
        return W_matrix
    def get_Mu (W_matrix):  # not use anymore
        def distance(vec1, vec2):
            return np.sqrt(np.sum(((vec1-vec2)/vec1)**2))
        def operator (vec):
            neqm_dist_vec = np.asarray(neqm_dist).reshape(vec.shape) 
            return neqm_dist_vec / vec
        N = len(eqm_dist)
        Mu = np.ones((N,1))
        i = 0
        while i < iter_num:
            new_Mu = operator(np.matmul(W_matrix, Mu))
            if distance(Mu, new_Mu) < 10**-10:
                return np.squeeze(new_Mu)
            Mu = copy.deepcopy(new_Mu)
            i +=1
        return np.squeeze(new_Mu)
    def get_Mu_numeric (W_matrix):
        def equations (Mu):
            eqns = []
            for i in range(N):
                factor = 0.0
                for j in range(N):
                    factor += W_matrix[i][j] * Mu[j]
                eqns.append(neqm_dist[i]/factor - Mu[i])
            return eqns
        N = len(eqm_dist)
        Mu = fsolve(equations, [1.0]*N)
        #Mu = newton_krylov(equations, [1.0]*N)
        #Mu = broyden1(equations, [1.0]*N)
        #Mu = broyden2(equations, [1.0]*N)
        #Mu = anderson(equations, [1.0]*N)
        #Mu = diagbroyden(equations, [1.0]*N)
        #Mu = excitingmixing(equations, [1.0]*N)
        return Mu
    def get_average(Mu, W_matrix):
        N = len(eqm_dist)
        average = 0.0
        for i in range(N):
            for j in range(N):
                average += Mu[i]*Mu[j]*W_matrix[i][j]*obs_matrix[i][j]
        return average
    def get_neqm_matrix (Mu, W_matrix):
        N = len(eqm_dist)
        eqm_nmatrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                eqm_nmatrix[i][j] = Mu[i]*Mu[j]*W_matrix[i][j]/neqm_dist[i]
        return eqm_nmatrix

    N = len(eqm_dist)

    if np.any(np.isin(edge_matrix, None)):
        edge_matrix = np.ones((N,N))

    if np.any(np.isin(obs_matrix, None)) or obs_average == None:
        gamma = 0.0
        W_matrix = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                W_matrix[i][j] = edge_matrix[i][j] * eqm_matrix[i][j] /eqm_dist[j]
        Mu = get_Mu_numeric(W_matrix)
        return get_neqm_matrix (Mu, W_matrix)

    gamma_guess = np.linspace(0.0, 100.0, num=1000)
    #gamma_guess = [100]
    diff = sys.maxint
    gamma_best, Mu_best, W_best = None, None, None

    for k in range(len(gamma_guess)):
        gamma = gamma_guess[k]
        W_matrix = get_W_matrix(gamma)
        Mu = get_Mu_numeric(W_matrix)
        average = get_average(Mu, W_matrix)
        #print "gamma:", gamma
        #print "W matrix:"
        #print W_matrix
        #print "Mu:", Mu
        #print "average:", average
        #print "neqm matrix:"
        #print get_neqm_matrix(Mu, W_matrix)
        #print
        
        if abs(obs_average - average) < diff:
            diff = abs(obs_average - average)
            gamma_best = gamma
            Mu_best = Mu
            W_best = W_matrix

    print "gamma best:", gamma_best
    print "real average:", obs_average
    print "best average:", get_average(Mu_best, W_best)
    neqm_matrix = get_neqm_matrix(Mu_best, W_best)
    return neqm_matrix


#eqm_dist = [0.2, 0.8]
#eqm_matrix = [[0.2,0.8],[0.2,0.8]]
#neqm_dist = [3.0/11, 8.0/11]
#obs_matrix = [[1.0,1.0],[1.0,1.0]]
#obs_average = 5.0
#print eqm_to_neqm(eqm_matrix, eqm_dist, neqm_dist)
#print eqm_to_steady(eqm_matrix, neqm_dist)

"""
ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
NCP_len = 147
filenames1 = ["../../Illumina/plusoneHS/data/Plslib-HS_S1_L001_R.sort"]
filenames2 = ["../../Illumina/plusoneHS/data/Plslib-HS-30min_S2_L001_R.sort"]
key_slider1 = load.load_files(filenames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)
key_slider2 = load.load_files(filenames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)
#net_flow = key_slider1['8'].eqm_flux()
#net_flow2 = key_slider2['8'].eqm_flux()


dyadmap = key_slider1['8'].KDE(band_width=0.1)[NCP_len/2:ref_length-NCP_len/2]
force = key_slider1['8'].force_profile(kT=0.5)[NCP_len/2:ref_length-NCP_len/2]
dyadmap2 = key_slider2['8'].KDE(band_width=0.1)[NCP_len/2:ref_length-NCP_len/2]
#print dyadmap

N = len(dyadmap)
total = sum(dyadmap)
eqm_dist = [float(value) / total for value in dyadmap]
total = sum(dyadmap2)
neqm_dist = [float(value) / total for value in dyadmap2]

edge_matrix = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        if abs(i-j) > 1:
            continue
        edge_matrix[i][j] = 1.0
#print edge_matrix


eqm_matrix = eqm_maxcal(eqm_dist, edge_matrix=edge_matrix)
print eqm_matrix
#neqm_matrix = eqm_to_neqm(eqm_matrix, eqm_dist, neqm_dist, edge_matrix=edge_matrix)
#eqm_matrix = eqm_maxcal(neqm_dist, edge_matrix=edge_matrix)
#neqm_matrix = eqm_maxcal(neqm_dist, edge_matrix=edge_matrix)
neqm_matrix = eqm_to_steady(eqm_matrix, neqm_dist, edge_matrix = edge_matrix)
#neqm_matrix = steady_maxcal(neqm_dist, edge_matrix = edge_matrix)
print neqm_matrix
#eqm_matrix = eqm_maxcal(neqm_dist, edge_matrix=edge_matrix)
#print eqm_matrix

net_flow = []
net_flow2 = []
for i in range(N):
    if i == 0:
        net_flow.append(eqm_matrix[i][i+1])
        net_flow2.append(neqm_matrix[i][i+1])
    elif i == N-1:
        net_flow.append(-eqm_matrix[i][i-1])
        net_flow2.append(-neqm_matrix[i][i-1])
    else:
        net_flow.append(eqm_matrix[i][i+1] - eqm_matrix[i][i-1])
        net_flow2.append(neqm_matrix[i][i+1] - neqm_matrix[i][i-1])

#change_flow = [ (net_flow2[i] - net_flow[i]) / abs(net_flow[i]) for i in range(len(net_flow))]
change_flow = [ (net_flow2[i] - net_flow[i])  for i in range(len(net_flow))]        


fig = plt.figure()
#plt.plot(range(len(eqm_dist)), eqm_dist)
#plt.plot(range(len(neqm_dist)), neqm_dist)
plt.plot(range(len(net_flow)), net_flow, label="MaxCal")
plt.plot(range(len(net_flow)), net_flow2, label="MaxCal2")
plt.plot(range(len(net_flow)), change_flow, label="Diff")
#plt.plot(range(len(net_flow)), force, label="Force")
#net_flow = [int(round(value*1000)) for value in net_flow]
#plt.scatter(range(len(eqm_dist)), eqm_dist, s=3, c=net_flow, cmap='bwr', vmin=-500, vmax=500)
plt.legend()
plt.savefig("maxcal.png")
plt.show()
plt.close()

"""
