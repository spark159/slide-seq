import sys
import numpy as np
import copy
import random
from scipy.optimize import fsolve

def steady_maxcal (steady_dist, obs_matrix, obs_average, iter_num=1000):
    def distance(vec1, vec2):
        temp = 0.0
        for i in range(len(vec1)):
            temp += (vec1[i] - vec2[i])**2
        return np.sqrt(temp)
    def update (Lambda, gammma, steady_dist, obs_matrix):
        N = len(Lambda)
        new_Lambda = []
        for i in range(N):
            factor1 = 0.0
            for j in range(N):
                factor2 = 0.0
                for k in range(N):
                    factor2 += Lambda[k]*np.exp(-gamma*obs_matrix[j][k])
                factor1 += steady_dist[j]*np.exp(-gamma*obs_matrix[j][i])/factor2
            new_Lambda.append(steady_dist[i]/factor1)
        return new_Lambda
    def get_Beta (Lambda, gamma, steady_dist, obs_matrix):
        N = len(Lambda)
        Beta = []
        for i in range(N):
            factor = 0.0
            for j in range(N):
                factor += Lambda[j] * np.exp(-gamma*obs_matrix[i][j])
            Beta.append(steady_dist[i] / factor)
        return Beta
    def get_matrix (Lambda, gamma, steady_dist, obs_matrix):
        N = len(Lambda)
        Beta = get_Beta(Lambda, gamma, steady_dist, obs_matrix)
        matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                matrix[i][j] = Beta[i]*Lambda[j]*np.exp(-gamma*obs_matrix[i][j])/steady_dist[i]
        return matrix
    def get_average (Lambda, gamma, steady_dist, obs_matrix):
        Beta = get_Beta(Lambda, gamma, steady_dist, obs_matrix)
        average = 0.0
        N = len(Lambda)
        for i in range(N):
            for j in range(N):
                average += Beta[i] * Lambda[j] * np.exp(-gamma*obs_matrix[i][j])*obs_matrix[i][j]
        return average
    
    N = len(steady_dist)
    gamma_guess = np.linspace(0.0, 10.0, num=10)
    diff = sys.maxint
    gamma_best, Lambda_best = None, None
    
    for k in range(len(gamma_guess)):
        gamma = gamma_guess[k]
        Lambda = [1.0] * N
        i = 0
        while i < iter_num:
            new_Lambda = update(Lambda, gamma, steady_dist, obs_matrix)
            #print new_Lambda
            print Lambda
            print new_Lambda
            print distance(new_Lambda, Lambda)
            print 
            if distance(new_Lambda, Lambda) < 10**-15:
                break
            Lambda = copy.deepcopy(new_Lambda)        
            i += 1
        average = get_average(Lambda, gamma, steady_dist, obs_matrix)
        #print abs(obs_average - average)
        #print get_matrix(new_Lambda, gamma, steady_dist, obs_matrix)
        if abs(obs_average - average) < diff:
            diff = abs(obs_average - average)
            gamma_best = gamma
            Lambda_best = new_Lambda
        
    steady_matrix = get_matrix(Lambda_best, gamma_best, steady_dist, obs_matrix)
    return steady_matrix
    
def eqm_maxcal (eqm_dist, obs_matrix, obs_average, iter_num=1000):
    def distance(vec1, vec2):
        temp = 0.0
        for i in range(len(vec1)):
            temp += (vec1[i] - vec2[i])**2
        return np.sqrt(temp)
    def update (Rho, gammma, eqm_dist, obs_matrix):
        N = len(Rho)
        new_Rho = []
        for i in range(N):
            factor = 0.0
            for j in range(N):
                factor += Rho[j] * np.exp(-0.5*gamma*(obs_matrix[i][j]+obs_matrix[j][i]))
            new_Rho.append(eqm_dist[i]/factor)
        return new_Rho
    def get_matrix (Rho, gamma, eqm_dist, obs_matrix):
        N = len(Rho)
        matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                matrix[i][j] = Rho[i]*Rho[j]*np.exp(-0.5*gamma*(obs_matrix[i][j]+obs_matrix[j][i]))/eqm_dist[i]
        return matrix
    def get_average (Rho, gamma, eqm_dist, obs_matrix):
        average = 0.0
        N = len(Rho)
        for i in range(N):
            for j in range(N):
                average += Rho[i]*Rho[j]*np.exp(-0.5*gamma*(obs_matrix[i][j]+obs_matrix[j][i]))*obs_matrix[i][j]
        return average
    
    N = len(eqm_dist)
    #gamma_guess = np.linspace(0.0, 10.0, num=1000)
    gamma_guess = [1.381]
    diff = sys.maxint
    gamma_best, Rho_best = None, None
    
    for k in range(len(gamma_guess)):
        gamma = gamma_guess[k]
        Rho = [1.0] * N
        i = 0
        while i < iter_num:
            new_Rho = update(Rho, gamma, eqm_dist, obs_matrix)
            print Rho
            print new_Rho
            print distance(new_Rho, Rho)
            print 
            if distance(new_Rho, Rho) < 10**-15:
                break
            Rho = copy.deepcopy(new_Rho)        
            i += 1
        average = get_average(Rho, gamma, eqm_dist, obs_matrix)
        #print abs(obs_average - average)
        #print get_matrix(new_Rho, gamma, eqm_dist, obs_matrix)
        if abs(obs_average - average) < diff:
            diff = abs(obs_average - average)
            gamma_best = gamma
            Rho_best = new_Rho
    print gamma_best
    eqm_matrix = get_matrix(Rho_best, gamma_best, eqm_dist, obs_matrix)
    return eqm_matrix
    
eqm_dist = [0.5, 0.5]
obs_matrix = [[1,1],[1,1]]
obs_average = 1
print eqm_maxcal(eqm_dist, obs_matrix, obs_average)
#print steady_maxcal(eqm_dist, obs_matrix, obs_average)

def eqm_to_neqm (eqm_matrix, neqm_dist, iter_num=1000):
    def distance(vec1, vec2):
        temp = 0.0
        for i in range(len(vec1)):
            temp += (vec1[i] - vec2[i])**2
        return np.sqrt(temp)
    def equations (Omega):
        N = len(Omega)
        eqns = []
        for i in range(N):
            factor1 = 0.0
            for j in range(N):
                factor2 = 0.0
                for k in range(N):
                    factor2 += neqm_dist[k]*eqm_matrix[k][j]*Omega[k]
                factor1 += eqm_matrix[i][j]*neqm_dist[j]/factor2
            eqns.append(Omega[i] - 1.0/factor1)
        return eqns
    def get_matrix (Omega):
        N = len(Omega)
        matrix = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                factor = 0.0
                for k in range(N):
                    factor += neqm_dist[k]*eqm_matrix[k][j]*Omega[k]
                matrix[i][j] = eqm_matrix[i][j]*neqm_dist[j]*Omega[i]/factor
        return matrix

    N = len(neqm_dist)
    Omega = fsolve(equations, [1.0]*N)
    neqm_matrix = get_matrix(Omega)
    return neqm_matrix

def eqm_to_steady (eqm_matrix, steady_dist, iter_num=1000):
    def distance(vec1, vec2):
        temp = 0.0
        for i in range(len(vec1)):
            temp += (vec1[i] - vec2[i])**2
        return np.sqrt(temp)
    def update(Omega, eqm_matrix, steady_dist):
        N = len(Omega)
        new_Omega = []
        for i in range(N):
            factor1 = 0.0
            for j in range(N):
                factor2 = 0.0
                for k in range(N):
                    factor2 += steady_dist[k]*eqm_matrix[k][j]*Omega[k]
                factor1 += eqm_matrix[i][j]*steady_dist[j]/factor2
            new_Omega.append(1.0/factor1)
        return new_Omega
    def get_matrix(Omega, eqm_matrix, steady_dist):
        N = len(Omega)
        temp = []
        for i in range(N):
            factor = 0.0
            for j in range(N):
                factor += steady_dist[j] * eqm_matrix[j][i] * Omega[j]
            temp.append(factor/steady_dist[i])
        matrix = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                matrix[i][j] = eqm_matrix[i][j] * Omega[i] * (1.0/temp[j])
        return matrix

    N = len(steady_dist)
    Omega = [1.0] * N
    i = 0
    while i < iter_num:
        new_Omega = update(Omega, eqm_matrix, steady_dist)
        #print new_Omega
        if distance(new_Omega, Omega) < 10**-15:
            return get_matrix(Omega, eqm_matrix, steady_dist)
        Omega = copy.deepcopy(new_Omega)        
        i += 1
    return get_matrix(Omega, eqm_matrix, steady_dist)
