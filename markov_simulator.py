import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import maxcal_Dixit as maxcal

def Markov_chain (trans_matrix, pop_size, time_len):
    def update (Traces):
        for i in range(len(Traces)):
            cur = Traces[i][-1]
            prob_list = trans_matrix[cur]
            choice = np.random.choice(len(prob_list), 1, p=prob_list)[0]
            Traces[i].append(choice)
        return Traces
    Traces = []
    for i in range(pop_size):
        Traces.append([0])
    for i in range(time_len):
        update(Traces)
    return Traces

def pop_average (Traces, state_obs):
    average = 0.0 
    for i in range(len(Traces)):
        state = Traces[i][-1]
        obs = state_obs[state]
        average += obs
    return average / len(Traces)

def time_average (Traces, state_obs):
    average = 0.0
    for i in range(len(Traces[0])):
        state = Traces[0][i]
        obs = state_obs[state]
        average += obs
    return average / len(Traces[0])

def get_dist (Traces, trans_matrix, option):
    state_num = len(trans_matrix)
    state_obs = []
    for i in range(state_num):
        obs = np.zeros(state_num)
        obs[i] = 1
        state_obs.append(obs)
    if option == "pop":
        return pop_average(Traces, state_obs)
    elif option == "time":
        return time_average(Traces, state_obs)        

def get_random_chain (state_num=None, edge_matrix=None):
    if state_num == None and edge_matrix != None:
        state_num = len(edge_matrix)
    elif state_num != None and edge_matrix == None:
        edge_matrix = np.ones((state_num, state_num))
    else:
        assert state_num == len(edge_matrix)
    trans_matrix = []
    for i in range(state_num):
        row = []
        total = 0.0
        for j in range(state_num):
            if edge_matrix[i][j] <= 0:
                row.append(0.0)
                continue
            choice = random.random()
            total += choice
            row.append(choice)
        row = [ value/total for value in row ]
        trans_matrix.append(row)
    return trans_matrix

def display_traces (Traces):
    fig = plt.figure()
    for i in range(len(Traces)):
        trace = Traces[i]
        x = range(len(trace))
        plt.plot(x, trace)
    plt.xlabel("Time steps")
    plt.ylabel("States")
    plt.show()
    plt.close()
    return None

trans_matrix = get_random_chain(state_num=1000)
#print trans_matrix
#trans_matrix = [[0.5,0.5],[0.5,0.5]]
Traces = Markov_chain(trans_matrix, 10, 30)
display_traces (Traces)
#print get_dist(Traces, trans_matrix, "pop")
