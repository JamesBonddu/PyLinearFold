import numpy as np
import heapq
import time
from collections import defaultdict
from itertools import product
from utility import *
from utility_v import * 

# Constants
NOTON = 5
NOTOND = 25
NOTONT = 125
SINGLE_MAX_LEN = 30
SINGLE_MIN_LEN = 0
EXPLICIT_MAX_LEN = 4
HAIRPIN_MAX_LEN = 30
INTERNAL_MAX_LEN = 30
SYMMETRIC_MAX_LEN = 15
ASYMMETRY_MAX_LEN = 28
BULGE_MAX_LEN = 30
MIN_CUBE_PRUNING_SIZE = 20
VALUE_MIN = -float('inf')

# Manner enum
MANNER_NONE = 0
MANNER_H = 1
MANNER_HAIRPIN = 2
MANNER_SINGLE = 3
MANNER_HELIX = 4
MANNER_MULTI = 5
MANNER_MULTI_eq_MULTI_plus_U = 6
MANNER_P_eq_MULTI = 7
MANNER_M2_eq_M_plus_P = 8
MANNER_M_eq_M2 = 9
MANNER_M_eq_M_plus_U = 10
MANNER_M_eq_P = 11
MANNER_C_eq_C_plus_U = 12
MANNER_C_eq_C_plus_P = 13

# BestTypes enum
TYPE_C = 0
TYPE_H = 1
TYPE_P = 2
TYPE_M = 3
TYPE_Multi = 4
TYPE_M2 = 5

class State:
    def __init__(self, score=VALUE_MIN, manner=MANNER_NONE, trace=None):
        self.score = score
        self.manner = manner
        self.trace = trace if trace is not None else TraceInfo()

class TraceInfo:
    def __init__(self, split=0, paddings=None):
        self.split = split
        self.paddings = paddings if paddings is not None else Paddings()

class Paddings:
    def __init__(self, l1=0, l2=0):
        self.l1 = l1
        self.l2 = l2

class BeamCKYParser:
    def __init__(self, beam_size=100, no_sharp_turn=True, is_verbose=False, use_constraints=False, zuker_subopt=False, energy_delta=5.0, shape_file_path="", is_fasta=False, dangle_model=2):
        self.beam = beam_size
        self.no_sharp_turn = no_sharp_turn
        self.is_verbose = is_verbose
        self.use_constraints = use_constraints
        self.zuker = zuker_subopt
        self.zuker_energy_delta = energy_delta
        self.use_shape = shape_file_path != ""
        self.m = 1.8
        self.b = -0.6
        self.is_fasta = is_fasta
        self.dangle_model = dangle_model

        self.initialize()

    def initialize(self):
        self.bestH = defaultdict(dict)
        self.bestP = defaultdict(dict)
        self.bestM2 = defaultdict(dict)
        self.bestMulti = defaultdict(dict)
        self.bestM = defaultdict(dict)
        self.bestC = [State(VALUE_MIN, MANNER_NONE, TraceInfo(0, Paddings(0, 0))) for _ in range(SINGLE_MAX_LEN + 1)]
        self.nucs = []
        self.scores = []

    def parse(self, seq, cons=None):
        self.prepare(len(seq))
        self.nucs = [self.GET_ACGU_NUM(nuc) for nuc in seq]

        for j in range(len(seq)):
            self.beam_step(j)

        result = self.get_result()
        return result

    def prepare(self, length):
        self.seq_length = length
        self.bestH.clear()
        self.bestP.clear()
        self.bestM2.clear()
        self.bestMulti.clear()
        self.bestM.clear()
        self.bestC = [State(VALUE_MIN, MANNER_NONE, TraceInfo(0, Paddings(0, 0))) for _ in range(length)]
        self.nucs = [0] * length
        self.scores = []

    def beam_step(self, j):
        pass  # Placeholder for the actual beam search logic

    def get_result(self):
        result = ''.join(['.' for _ in range(self.seq_length)])
        return result

    def GET_ACGU_NUM(self, nuc):
        return 'ACGU'.index(nuc)

    def update_if_better(self, state, newscore, manner, split=0, paddings=None):
        if paddings is None:
            paddings = Paddings(0, 0)
        if state.score < newscore:
            state.score = newscore
            state.manner = manner
            state.trace = TraceInfo(split, paddings)

    def beam_prune(self, beamstep):
        scores = [(item[1].score, item[0]) for item in beamstep.items()]
        if len(scores) <= self.beam:
            return VALUE_MIN
        threshold = self.quickselect(scores, 0, len(scores) - 1, len(scores) - self.beam)
        for score, key in scores:
            if score < threshold:
                del beamstep[key]
        return threshold

    def quickselect(self, scores, lower, upper, k):
        if lower == upper:
            return scores[lower][0]
        pivot_index = self.partition(scores, lower, upper)
        if k == pivot_index:
            return scores[k][0]
        elif k < pivot_index:
            return self.quickselect(scores, lower, pivot_index - 1, k)
        else:
            return self.quickselect(scores, pivot_index + 1, upper, k)

    def partition(self, scores, lower, upper):
        pivot = scores[upper][0]
        i = lower
        for j in range(lower, upper):
            if scores[j][0] < pivot:
                scores[i], scores[j] = scores[j], scores[i]
                i += 1
        scores[i], scores[upper] = scores[upper], scores[i]
        return i

    def get_parentheses(self, result, seq):
        result = ['.' for _ in range(self.seq_length)]
        stack = [(0, self.seq_length - 1, self.bestC[self.seq_length - 1])]

        while stack:
            i, j, state = stack.pop()

            if state.manner == MANNER_H:
                continue
            elif state.manner == MANNER_HAIRPIN:
                result[i] = '('
                result[j] = ')'
                if self.is_verbose:
                    tetra_hex_tri = -1
                    if j - i - 1 == 4:
                        tetra_hex_tri = self.if_tetraloops[i]
                    elif j - i - 1 == 6:
                        tetra_hex_tri = self.if_hexaloops[i]
                    elif j - i - 1 == 3:
                        tetra_hex_tri = self.if_triloops[i]
                    newscore = v_score_hairpin(i, j, self.nucs[i], self.nucs[i + 1], self.nucs[j - 1], self.nucs[j], tetra_hex_tri)
                    print(f"Hairpin loop ( {i + 1}, {j + 1}) {seq[i]}{seq[j]} : {newscore / -100.0}")
            elif state.manner == MANNER_SINGLE:
                result[i] = '('
                result[j] = ')'
                p = i + state.trace.paddings.l1
                q = j - state.trace.paddings.l2
                stack.append((p, q, self.bestP[q][p]))
                if self.is_verbose:
                    newscore = v_score_single(i, j, p, q, self.nucs[i], self.nucs[i + 1], self.nucs[j - 1], self.nucs[j], self.nucs[p - 1], self.nucs[p], self.nucs[q], self.nucs[q + 1])
                    print(f"Interior loop ( {i + 1}, {j + 1}) {seq[i]}{seq[j]}; ( {p + 1}, {q + 1}) {seq[p]}{seq[q]} : {newscore / -100.0}")
            elif state.manner == MANNER_HELIX:
                result[i] = '('
                result[j] = ')'
                stack.append((i + 1, j - 1, self.bestP[j - 1][i + 1]))
                if self.is_verbose:
                    p = i + 1
                    q = j - 1
                    newscore = v_score_single(i, j, p, q, self.nucs[i], self.nucs[i + 1], self.nucs[j - 1], self.nucs[j], self.nucs[p - 1], self.nucs[p], self.nucs[q], self.nucs[q + 1])
                    print(f"Interior loop ( {i + 1}, {j + 1}) {seq[i]}{seq[j]}; ( {p + 1}, {q + 1}) {seq[p]}{seq[q]} : {newscore / -100.0}")
            elif state.manner == MANNER_MULTI:
                p = i + state.trace.paddings.l1
                q = j - state.trace.paddings.l2
                stack.append((p, q, self.bestM2[q][p]))
            elif state.manner == MANNER_MULTI_eq_MULTI_plus_U:
                p = i + state.trace.paddings.l1
                q = j - state.trace.paddings.l2
                stack.append((p, q, self.bestM2[q][p]))
            elif state.manner == MANNER_P_eq_MULTI:
                result[i] = '('
                result[j] = ')'
                stack.append((i, j, self.bestMulti[j][i]))
            elif state.manner == MANNER_M2_eq_M_plus_P:
                k = state.trace.split
                stack.append((i, k, self.bestM[k][i]))
                stack.append((k + 1, j, self.bestP[j][k + 1]))
            elif state.manner == MANNER_M_eq_M2:
                stack.append((i, j, self.bestM2[j][i]))
            elif state.manner == MANNER_M_eq_M_plus_U:
                stack.append((i, j - 1, self.bestM[j - 1][i]))
            elif state.manner == MANNER_M_eq_P:
                stack.append((i, j, self.bestP[j][i]))
            elif state.manner == MANNER_C_eq_C_plus_U:
                k = j - 1
                if k != -1:
                    stack.append((0, k, self.bestC[k]))
            elif state.manner == MANNER_C_eq_C_plus_P:
                k = state.trace.split
                if k != -1:
                    stack.append((0, k, self.bestC[k]))
                    stack.append((k + 1, j, self.bestP[j][k + 1]))
                else:
                    stack.append((i, j, self.bestP[j][i]))

        return ''.join(result)

    def get_parentheses_inside_real_backtrace(self, i, j, state, global_visited_inside, window_visited):
        manner = state.manner
        if manner == MANNER_H:
            return ""
        elif manner == MANNER_HAIRPIN:
            if (TYPE_P, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_P, i, j)]
            inner = '.' * (j - i - 1)
            global_visited_inside[(TYPE_P, i, j)] = '(' + inner + ')'
            self.window_fill(window_visited, i, j, self.seq_length, self.window_size)
            return global_visited_inside[(TYPE_P, i, j)]
        elif manner == MANNER_SINGLE:
            if (TYPE_P, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_P, i, j)]
            p = i + state.trace.paddings.l1
            q = j - state.trace.paddings.l2
            left = '.' * (state.trace.paddings.l1 - 1)
            right = '.' * (state.trace.paddings.l2 - 1)
            if (TYPE_P, p, q) not in global_visited_inside:
                global_visited_inside[(TYPE_P, p, q)] = self.get_parentheses_inside_real_backtrace(p, q, self.bestP[q][p], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_P, i, j)] = '(' + left + global_visited_inside[(TYPE_P, p, q)] + right + ')'
            self.window_fill(window_visited, i, j, self.seq_length, self.window_size)
            return global_visited_inside[(TYPE_P, i, j)]
        elif manner == MANNER_HELIX:
            if (TYPE_P, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_P, i, j)]
            if (TYPE_P, i + 1, j - 1) not in global_visited_inside:
                global_visited_inside[(TYPE_P, i + 1, j - 1)] = self.get_parentheses_inside_real_backtrace(i + 1, j - 1, self.bestP[j - 1][i + 1], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_P, i, j)] = '(' + global_visited_inside[(TYPE_P, i + 1, j - 1)] + ')'
            self.window_fill(window_visited, i, j, self.seq_length, self.window_size)
            return global_visited_inside[(TYPE_P, i, j)]
        elif manner == MANNER_MULTI:
            if (TYPE_Multi, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_Multi, i, j)]
            p = i + state.trace.paddings.l1
            q = j - state.trace.paddings.l2
            left = '.' * (state.trace.paddings.l1 - 1)
            right = '.' * (state.trace.paddings.l2 - 1)
            if (TYPE_M2, q, p) not in global_visited_inside:
                global_visited_inside[(TYPE_M2, q, p)] = self.get_parentheses_inside_real_backtrace(p, q, self.bestM2[q][p], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_Multi, i, j)] = left + global_visited_inside[(TYPE_M2, q, p)] + right
            return global_visited_inside[(TYPE_Multi, i, j)]
        elif manner == MANNER_MULTI_eq_MULTI_plus_U:
            if (TYPE_Multi, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_Multi, i, j)]
            p = i + state.trace.paddings.l1
            q = j - state.trace.paddings.l2
            left = '.' * (state.trace.paddings.l1 - 1)
            right = '.' * (state.trace.paddings.l2 - 1)
            if (TYPE_M2, p, q) not in global_visited_inside:
                global_visited_inside[(TYPE_M2, p, q)] = self.get_parentheses_inside_real_backtrace(p, q, self.bestM2[q][p], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_Multi, i, j)] = left + global_visited_inside[(TYPE_M2, p, q)] + right
            return global_visited_inside[(TYPE_Multi, i, j)]
        elif manner == MANNER_P_eq_MULTI:
            if (TYPE_P, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_P, i, j)]
            if (TYPE_Multi, i, j) not in global_visited_inside:
                global_visited_inside[(TYPE_Multi, i, j)] = self.get_parentheses_inside_real_backtrace(i, j, self.bestMulti[j][i], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_P, i, j)] = '(' + global_visited_inside[(TYPE_Multi, i, j)] + ')'
            self.window_fill(window_visited, i, j, self.seq_length, self.window_size)
            return global_visited_inside[(TYPE_P, i, j)]
        elif manner == MANNER_M2_eq_M_plus_P:
            if (TYPE_M2, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_M2, i, j)]
            k = state.trace.split
            if (TYPE_M, i, k) not in global_visited_inside:
                global_visited_inside[(TYPE_M, i, k)] = self.get_parentheses_inside_real_backtrace(i, k, self.bestM[k][i], global_visited_inside, window_visited)
            if (TYPE_P, k + 1, j) not in global_visited_inside:
                global_visited_inside[(TYPE_P, k + 1, j)] = self.get_parentheses_inside_real_backtrace(k + 1, j, self.bestP[j][k + 1], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_M2, i, j)] = global_visited_inside[(TYPE_M, i, k)] + global_visited_inside[(TYPE_P, k + 1, j)]
            return global_visited_inside[(TYPE_M2, i, j)]
        elif manner == MANNER_M_eq_M2:
            if (TYPE_M, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_M, i, j)]
            if (TYPE_M2, i, j) not in global_visited_inside:
                global_visited_inside[(TYPE_M2, i, j)] = self.get_parentheses_inside_real_backtrace(i, j, self.bestM2[j][i], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_M, i, j)] = global_visited_inside[(TYPE_M2, i, j)]
            return global_visited_inside[(TYPE_M, i, j)]
        elif manner == MANNER_M_eq_M_plus_U:
            if (TYPE_M, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_M, i, j)]
            if (TYPE_M, i, j - 1) not in global_visited_inside:
                global_visited_inside[(TYPE_M, i, j - 1)] = self.get_parentheses_inside_real_backtrace(i, j - 1, self.bestM[j - 1][i], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_M, i, j)] = global_visited_inside[(TYPE_M, i, j - 1)] + '.'
            return global_visited_inside[(TYPE_M, i, j)]
        elif manner == MANNER_M_eq_P:
            if (TYPE_M, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_M, i, j)]
            if (TYPE_P, i, j) not in global_visited_inside:
                global_visited_inside[(TYPE_P, i, j)] = self.get_parentheses_inside_real_backtrace(i, j, self.bestP[j][i], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_M, i, j)] = global_visited_inside[(TYPE_P, i, j)]
            return global_visited_inside[(TYPE_M, i, j)]
        elif manner == MANNER_C_eq_C_plus_U:
            assert i == 0
            if (TYPE_C, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_C, i, j)]
            k = j - 1
            if (TYPE_C, i, k) not in global_visited_inside:
                if k != -1:
                    global_visited_inside[(TYPE_C, i, k)] = self.get_parentheses_inside_real_backtrace(0, k, self.bestC[k], global_visited_inside, window_visited)
                else:
                    global_visited_inside[(TYPE_C, i, k)] = ""
            global_visited_inside[(TYPE_C, i, j)] = global_visited_inside[(TYPE_C, i, k)] + '.'
            return global_visited_inside[(TYPE_C, i, j)]
        elif manner == MANNER_C_eq_C_plus_P:
            assert i == 0
            if (TYPE_C, i, j) in global_visited_inside:
                return global_visited_inside[(TYPE_C, i, j)]
            k = state.trace.split
            if (TYPE_C, i, k) not in global_visited_inside:
                if k != -1:
                    global_visited_inside[(TYPE_C, i, k)] = self.get_parentheses_inside_real_backtrace(0, k, self.bestC[k], global_visited_inside, window_visited)
                else:
                    global_visited_inside[(TYPE_C, i, k)] = ""
            if (TYPE_P, k + 1, j) not in global_visited_inside:
                global_visited_inside[(TYPE_P, k + 1, j)] = self.get_parentheses_inside_real_backtrace(k + 1, j, self.bestP[j][k + 1], global_visited_inside, window_visited)
            global_visited_inside[(TYPE_C, i, j)] = global_visited_inside[(TYPE_C, i, k)] + global_visited_inside[(TYPE_P, k + 1, j)]
            return global_visited_inside[(TYPE_C, i, j)]
        else:
            print(f"wrong manner inside at {i}, {j}: manner {state.manner}", file=sys.stderr)
            assert False

    def window_fill(self, window_visited, i, j, seq_length, window_size):
        for ii in range(max(0, i - window_size), min(seq_length - 1, i + window_size) + 1):
            for jj in range(max(0, j - window_size), min(seq_length - 1, j + window_size) + 1):
                if ii < jj:
                    window_visited.add((ii, jj))
                    

parser = BeamCKYParser()
# seq = "ACGUACGUACGU"
seq = "GUUUUUAUCUUACACACGCUUGUGUAAGAUAGUUA"
result = parser.parse(seq)
print(result)