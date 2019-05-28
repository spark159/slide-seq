#!/usr/bin/env python
import sys
import math

class MarkovModel:
    def __init__(self, order):
        self.order = order

    # Train Markove Model from genome and a list of DNA sequences
    def train(self, genome, seq_list):
        # Return if no input
        if len(seq_list) <= 0:
            return

        # Compute background signal
        self.bg_prob = [{} for i in range(self.order+1)]
        for seq in genome.values():
            for i in range(len(seq)):
                for j in range(self.order + 1):
                    if i < j:
                        continue
                    prev = seq[i-j:i]
                    cur = seq[i]
                    assert cur in "ACGT"
                    prob_j = self.bg_prob[j]
                    if prev not in prob_j:
                        prob_j[prev] = {}
                    prob_prev = prob_j[prev]
                    if cur not in prob_prev:
                        prob_prev[cur] = 1
                    else:
                        prob_prev[cur] += 1

        # Normalize and take the log values of them
        for prob_order in self.bg_prob:
            for prev, prob_prev in prob_order.items():
                total_count = sum(prob_prev.values())
                total_count = float(total_count)
                for cur, count in prob_prev.items():
                    prob_prev[cur] = math.log(count / total_count)

        # Make sure that sequences are of the same length
        self.seq_len = len(seq_list[0])
        for seq in seq_list[1:]:
            assert self.seq_len == len(seq)

        # self.prob[order][position][previous_bases + current_base]
        self.prob = [[{} for j in range(self.seq_len)] for i in range(self.order+1)]
        for seq in seq_list:
            for i in range(self.seq_len):
                for j in range(self.order + 1):
                    if i < j:
                        continue
                    prev = seq[i-j:i]
                    cur = seq[i]
                    assert cur in "ACGT"
                    prob_i = self.prob[j][i]
                    if prev not in prob_i:
                        prob_i[prev] = {}
                    prob_prev = prob_i[prev]
                    if cur not in prob_prev:
                        prob_prev[cur] = 1
                    else:
                        prob_prev[cur] += 1

        # Normalize and take the log values of them
        for prob_order in self.prob:
            for prob_pos in prob_order:
                for prev, prob_prev in prob_pos.items():
                    total_count = sum(prob_prev.values())
                    total_count = float(total_count)
                    for cur, count in prob_prev.items():
                        prob_prev[cur] = math.log(count / total_count)

    # Calculate log-probability of a given sequence
    def predict(self, seq, background = True):
        if self.seq_len != len(seq):
            return

        score = 0.0
        for i in range(self.seq_len):
            cur = seq[i]
            if i >= self.order:
                prev = seq[i-self.order:i]
                order = self.order
            else:
                prev = seq[:i]
                order = i

            if prev not in self.prob[order][i] or cur not in self.prob[order][i][prev]:
                score = (-sys.float_info.max / 2)
            else:
                score += self.prob[order][i][prev][cur]

            if background:
                if prev not in self.bg_prob[order] or cur not in self.bg_prob[order][prev]:
                    score = (-sys.float_info.max / 2)
                else:
                    score -= self.bg_prob[order][prev][cur]
                
        return score

    # Calculate the "single" NCP positioning probability profile of given sequence
    def single_profile (self, seq):
        profile = []
        total = 0.0
        for i in range(len(seq)):
            if i < (self.seq_len / 2) or i > len(seq) - 1 - (self.seq_len / 2):
                score = 0.0
            else:
                score = self.predict(seq[i - (self.seq_len / 2) : i + (self.seq_len / 2) + 1])
            profile.append(score)
            total += score
        profile = [s / total for s in profile]
        return profile

    # Calculate the sum of weight-factor of given sub sequence (forward)
    def get_forward(self, seq):
        N = len(seq)
        F = [0.0] * (N + 1)
        bg_prob = self.bg_prob[0]['']
        for i in range(N):
            bp = seq[i]
            bg = bg_prob[bp]
            a = F[i] + bg
            b = 0.0
            if i >= self.seq_len:
                b = F[i+1-self.seq_len] + self.predict(seq[i-self.seq_len+1:i+1], False)
                b = math.log(1 + math.exp(b - a))
            F[i+1] = a + b
                  
        return F

    # Calculate the sum of weight-factor of given sub sequence (reverse)
    def get_reverse(self, seq):
        N = len(seq)
        R = [0.0] * (N + 1)
        bg_prob = self.bg_prob[0]['']
        for i in reversed(range(N)):
            bp = seq[i]
            bg = bg_prob[bp]
            a = R[i+1] + bg
            b = 0.0
            if i + self.seq_len <= N:
                b = R[i + self.seq_len] + self.predict(seq[i:i+self.seq_len], False)
                b = math.log(1 + math.exp(b - a))
            R[i] = a + b
        return R

    # Calculate the probabilty of NCP would be on ith position of given sequence
    def logprob_i(self, seq, i, F, R):
        N = len(seq)
        i -= (self.seq_len / 2)
        if i < 0 or i + self.seq_len > N:
            return 0.0

        return F[i] + self.predict(seq[i:i+self.seq_len], False) + R[i + self.seq_len] - R[0]

    # Calculate the positioning probability profile of given sequence
    def logprob_profile(self, seq):
        profile = []
        F = self.get_forward(seq)
        R = self.get_reverse(seq)
        for i in range(len(seq)):
            profile.append(self.logprob_i(seq, i, F, R))
        return profile
        
    # Display some detail about the model
    def help(self):
        print >> sys.stderr, "%d-order Markov Model" % (self.order)
        for order in range(self.order + 1):
            for pos in range(len(self.prob[order])):
                if pos % 20 != 0:
                    continue
                if len(self.prob[order][pos]) == 0:
                    continue
                print >> sys.stderr, "%d-order at pos %d" % (order, pos)
                print >> sys.stderr, "\t", self.prob[order][pos]


# To be implemented
class InterpolatedMarkovModel:
    def __init__(self):
        None
        

# To be implemented
class HiddenMarkovModel:
    def __init__(self):
        None
