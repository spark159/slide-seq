import random
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

def uniform_trans_rate_matrix (state_num, rate):
    trans_rate_matrix = {}
    for i in range(state_num):
        if i not in trans_rate_matrix:
            trans_rate_matrix[i] = {}
        if i < state_num - 1:
            trans_rate_matrix[i][i+1] = rate
        if i > 0:
            trans_rate_matrix[i][i-1] = rate
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
    for perted_site in perted_sites:
        for offset in offsets:
            try:
                trans_rate_matrix[perted_site+offset][perted_site+offset-1] = scale*trans_rate_matrix[perted_site+offset][perted_site+offset-1]
            except:
                pass
            try:
                trans_rate_matrix[perted_site-offset][perted_site-offset+1] = scale*trans_rate_matrix[perted_site-offset][perted_site-offset+1]
            except:
                pass
    return

ref_length = 225
dyad_axis = ref_length/2
dyad_offset = 52
NCP_len = 147

fnames1 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_before_.combined.sort"]
fnames2 = ["/home/spark159/../../media/spark159/sw/AscanlibFinal/601_after_.combined.sort"]
Control1 = load.load_files(fnames1, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)['601']
Control2 = load.load_files(fnames2, ref_length, dyad_axis, dyad_offset, filter_num = 10, fill=None)['601'] 

#trans_rate_matrix = guess_trans_rate_matrix (Control2.energy_profile(), 10)

state_num = ref_length
trans_rate_matrix = uniform_trans_rate_matrix (state_num, rate=10)

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
    
State_traces, Time_traces = kinetic_MC (trans_rate_matrix, 1, 100000, init_state=state_num/2)
#display_traces (State_traces, Time_traces)
dist1 = get_dist (State_traces, Time_traces, state_num, replica=0, time=None)

st, ed = state_num/2, state_num/2 +1
modify_trans_rate_matrix (trans_rate_matrix, perted_sites=range(st,ed), offsets = [19, 20, 21], scale=0.5)
State_traces, Time_traces = kinetic_MC (trans_rate_matrix, 1, 100000, init_state=state_num/2)
#display_traces (State_traces, Time_traces)
dist2 = get_dist (State_traces, Time_traces, state_num, replica=0, time=None)



fig = plt.figure()
plt.plot(dist1)
plt.plot(dist2)
plt.yscale("log")
plt.axvspan(st, ed-1, alpha=0.5, color='red')
plt.show()
plt.close()
