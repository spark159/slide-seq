import random
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import load

def kinetic_MC (trans_rate_matrix, pop_size, time_len, init_state=0):
    def update (cur_state, cur_time):
        state_list, rate_list = [], []
        for state, rate in trans_rate_matrix[cur_state].items():
            state_list.append(state)
            rate_list.append(rate)
        prob_list = [ float(rate)/sum(rate_list) for rate in rate_list ]
        new_state = state_list[np.random.choice(len(prob_list), 1, p=prob_list)[0]]
        new_time = cur_time + np.log(1.0/(1.0 - random.random()))/sum(rate_list)
        return new_state, new_time
    State_traces, Time_traces = [], []
    for i in range(pop_size):
        state_trace, time_trace = [init_state], [0]
        while True:
            cur_state, cur_time = state_trace[-1], time_trace[-1]
            new_state, new_time = update(cur_state, cur_time)
            state_trace.append(new_state)
            if new_time >= time_len:
                time_trace.append(time_len)
                break
            time_trace.append(new_time)
        State_traces.append(state_trace)
        Time_traces.append(time_trace)
    return State_traces, Time_traces

def pop_average (State_traces, Time_traces, state_obs, time):
    average = 0.0
    for i in range(len(State_traces)):
        idx = np.searchsorted(Time_traces[i], time, side='right') - 1
        state = State_traces[i][idx]
        obs = state_obs[state]
        average += obs
    return average / len(State_traces)

def time_average (State_traces, Time_traces, state_obs, replica):
    state_trace, time_trace = State_traces[replica], Time_traces[replica]
    average, total_time = 0.0, 0.0
    for i in range(len(time_trace)-1):
        time_step = time_trace[i+1] - time_trace[i]
        state = state_trace[i]
        obs = state_obs[state]
        average += obs*time_step
        total_time += time_step
    assert total_time == time_trace[-1]
    return average / float(total_time)

def get_dist (State_traces, Time_traces, state_num, replica=None, time=None):
    state_obs = []
    for i in range(state_num):
        obs = np.zeros(state_num)
        obs[i] = 1.0
        state_obs.append(obs)
    if replica != None and time == None:
        return time_average (State_traces, Time_traces, state_obs, replica)
    elif replica == None and time != None:
        return pop_average (State_traces, Time_traces, state_obs, time)
    else:
        return None

def display_traces (State_traces, Time_traces):
    fig = plt.figure()
    for i in range(len(State_traces)):
        state_trace = State_traces[i]
        time_trace = Time_traces[i]
        plt.plot(time_trace, state_trace)
    plt.xlabel("Time")
    plt.ylabel("States")
    plt.show()
    plt.close()
    return None

def uniform_trans_rate_matrix (state_num, rate, boundary='block'):
    trans_rate_matrix = {}
    for i in range(state_num):
        if i not in trans_rate_matrix:
            trans_rate_matrix[i] = {}
        if boundary == 'periodic' and i == 0:
            trans_rate_matrix[i][state_num-1] = rate
        else:
            trans_rate_matrix[i][i-1] = rate
        if boundary == 'periodic' and i == state_num - 1:
            trans_rate_matrix[i][0] = rate
        else:
            trans_rate_matrix[i][i+1] = rate
    return trans_rate_matrix

def guess_trans_rate_matrix (energy_profile, diff_const):
    state_num = len(energy_profile)
    trans_rate_matrix = {}
    for i in range(state_num):
        if i not in trans_rate_matrix:
            trans_rate_matrix[i] = {}
        if i < state_num - 1:
            trans_rate_matrix[i][i+1] = diff_const*np.exp(0.5*(energy_profile[i]-energy_profile[i+1]))
            #trans_rate_matrix[i][i+1] = diff_const
        if i > 0:
            trans_rate_matrix[i][i-1] = diff_const*np.exp(0.5*(energy_profile[i]-energy_profile[i-1]))
            #trans_rate_matrix[i][i-1] = diff_const
    return trans_rate_matrix

def modify_trans_rate_matrix (trans_rate_matrix, perted_sites, offsets = [19, 20, 21], scale=0.1):
    new_matrix = copy.deepcopy(trans_rate_matrix)
    for perted_site in perted_sites:
        for offset in offsets:
            try:
                new_matrix[perted_site+offset][perted_site+offset-1] = scale*new_matrix[perted_site+offset][perted_site+offset-1]
            except:
                pass
            try:
                new_matrix[perted_site-offset][perted_site-offset+1] = scale*new_matrix[perted_site-offset][perted_site-offset+1]
            except:
                pass
    return new_matrix

ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
NCP_len = 147

fnames1 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_before_.combined.sort"]
fnames2 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_after_.combined.sort"]
Control1 = load.load_files(fnames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)['601']
Control2 = load.load_files(fnames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)['601'] 

#trans_rate_matrix = guess_trans_rate_matrix (Control2.energy_profile(), 10)

#state_num = ref_length
state_num = 100
trans_rate_matrix = uniform_trans_rate_matrix (state_num, rate=10, boundary='periodic')

forward_list, backward_list = [] ,[]
for i in range(state_num):
    if i < state_num-1:
        forward = trans_rate_matrix[i][i+1]
    else:
        forward = 0.0
    forward_list.append(forward)
    if i > 0:
        backward = trans_rate_matrix[i][i-1]
    else:
        backward = 0.0
    backward_list.append(backward)

fig = plt.figure()
plt.plot(forward_list, '.')
plt.plot(backward_list, '.')
#plt.show()
plt.close()

st, ed = state_num/2, state_num/2 + 2
#scale_list = [1.0 - 0.1*i for i in range(4)]
#offset_list = [20 for i in range(4)]

scale_list = [0.8 for i in range(4)]
offset_list = [10*(i+1) for i in range(4)]

matrix_list = []
for i in range(4):
    scale = scale_list[i]
    offset = offset_list[i]
    new_matrix = modify_trans_rate_matrix (trans_rate_matrix, perted_sites=range(st, ed), offsets = range(offset-1, offset+2), scale=scale)
    matrix_list.append(new_matrix)

dist_list = []
for matrix in matrix_list:
    State_traces, Time_traces = kinetic_MC (matrix, 1, 20000, init_state=state_num/2)
    dist = get_dist (State_traces, Time_traces, state_num, replica=0, time=None)
    dist_list.append(dist)

fig = plt.figure()
uniform = 1.0/state_num
for i in range(len(dist_list)):
    dist = dist_list[i]
    scale = scale_list[i]
    offset = offset_list[i]
    plt.plot(dist, label='scale:'+str(scale)+', offset:'+str(offset))
    plt.axhline(y=uniform, linestyle = '--', color='k', alpha=0.25)
#plt.yscale("log")
plt.title("Kinetic Monte Carlo simulation")
plt.xlabel("Position")
plt.ylabel("Probability")
plt.ylim([0.0, uniform*3.0])
plt.axvspan(st, ed-1, alpha=0.5, color='red')
plt.legend()
plt.show()
plt.close()
