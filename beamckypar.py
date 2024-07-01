import sys
import time
from collections import defaultdict
import heapq

DOUBLE_MIN = float('-inf')
MIN_CUBE_PRUNING_SIZE = 20

class Manner:
    NONE = 0  # 0: empty
    H = 1     # 1: hairpin candidate
    HAIRPIN = 2  # 2: hairpin
    SINGLE = 3  # 3: single
    HELIX = 4  # 4: helix
    MULTI = 5  # 5: multi = ..M2. [30 restriction on the left and jump on the right]
    MULTI_eq_MULTI_plus_U = 6  # 6: multi = multi + U
    P_eq_MULTI = 7  # 7: P = (multi)
    M2_eq_M_plus_P = 8  # 8: M2 = M + P
    M_eq_M2 = 9  # 9: M = M2
    M_eq_M_plus_U = 10  # 10: M = M + U
    M_eq_P = 11  # 11: M = P
    C_eq_C_plus_U = 12  # 12: C = C + U
    C_eq_C_plus_P = 13  # 13: C = C + P

class State:
    def __init__(self, score=DOUBLE_MIN, manner=Manner.NONE, trace=None):
        self.score = score
        self.manner = manner
        self.trace = trace if trace is not None else {}

    def set(self, score, manner, split=None, l1=None, l2=None):
        self.score = score
        self.manner = manner
        if split is not None:
            self.trace['split'] = split
        if l1 is not None and l2 is not None:
            self.trace['paddings'] = (l1, l2)

class BeamCKYParser:
    def __init__(self, beam_size=0, vienna=False, candidate_list=True, nosharpturn=True, cube_pruning=True, verbose=False):
        self.beam = beam_size
        self.use_vienna = vienna
        self.is_candidate_list = candidate_list
        self.no_sharp_turn = nosharpturn
        self.is_cube_pruning = cube_pruning
        self.is_verbose = verbose
        self.seq_length = 0
        self.bestH = []
        self.bestP = []
        self.bestM2 = []
        self.bestM = []
        self.bestC = []
        self.bestMulti = []
        self.sorted_bestM = []
        self.bestC = []
        self.nucs = []
        self.scores = []

    def parse(self, seq):
        start_time = time.time()
        self.prepare(len(seq))
        self.nucs = [self.GET_ACGU_NUM(ch) for ch in seq]
        next_pair = self.compute_next_pair(self.nucs)

        if self.use_vienna:
            if self.seq_length > 0:
                self.bestC[0].set(-self.v_score_external_unpaired(0, 0), Manner.C_eq_C_plus_U)
            if self.seq_length > 1:
                self.bestC[1].set(-self.v_score_external_unpaired(0, 1), Manner.C_eq_C_plus_U)
        else:
            if self.seq_length > 0:
                self.bestC[0].set(self.score_external_unpaired(0, 0), Manner.C_eq_C_plus_U)
            if self.seq_length > 1:
                self.bestC[1].set(self.score_external_unpaired(0, 1), Manner.C_eq_C_plus_U)

        for j in range(self.seq_length):
            self.process_beam_H(j, next_pair)
            self.process_beam_Multi(j, next_pair)
            self.process_beam_P(j, next_pair)
            self.process_beam_M2(j, next_pair)
            self.process_beam_M(j, next_pair)
            self.process_beam_C(j, next_pair)

        viterbi = self.bestC[self.seq_length - 1]
        result = self.get_parentheses(seq)
        elapsed_time = time.time() - start_time

        return {
            "structure": result,
            "score": viterbi.score,
            "num_states": sum(len(beam) for beam in self.bestH + self.bestP + self.bestM2 + self.bestM + self.bestC + self.bestMulti),
            "time": elapsed_time
        }

    def prepare(self, length):
        self.seq_length = length
        self.bestH = [defaultdict(State) for _ in range(length)]
        self.bestP = [defaultdict(State) for _ in range(length)]
        self.bestM2 = [defaultdict(State) for _ in range(length)]
        self.bestM = [defaultdict(State) for _ in range(length)]
        self.bestC = [State() for _ in range(length)]
        self.bestMulti = [defaultdict(State) for _ in range(length)]
        self.sorted_bestM = [[] for _ in range(length)]
        self.nucs = [0] * length
        self.scores = []

    def compute_next_pair(self, nucs):
        next_pair = [[] for _ in range(len(nucs))]
        for nuci in range(len(nucs)):
            next_pair[nuci] = [-1] * len(nucs)
            next = -1
            for j in range(len(nucs) - 1, -1, -1):
                next_pair[nuci][j] = next
                if self._allowed_pairs[nuci][nucs[j]]:
                    next = j
        return next_pair
    
    def beam_prune(self, beamstep):
        scores = []
        for i, state in beamstep.items():
            k = i - 1
            newscore = (self.bestC[k].score if k >= 0 else 0) + state.score
            scores.append((newscore, i))
        if len(scores) <= self.beam:
            return DOUBLE_MIN
        threshold = self.quickselect(scores, 0, len(scores) - 1, len(scores) - self.beam)
        for p in scores:
            if p[0] < threshold:
                beamstep.pop(p[1], None)
        return threshold
    
    
    def quickselect(self, scores, lower, upper, k):
        if lower == upper:
            return scores[lower][0]
        split = self.quickselect_partition(scores, lower, upper)
        length = split - lower + 1
        if length == k:
            return scores[split][0]
        elif k < length:
            return self.quickselect(scores, lower, split - 1, k)
        else:
            return self.quickselect(scores, split + 1, upper, k - length)

    def quickselect_partition(self, scores, lower, upper):
        pivot = scores[upper][0]
        while lower < upper:
            while scores[lower][0] < pivot:
                lower += 1
            while scores[upper][0] > pivot:
                upper -= 1
            if scores[lower][0] == scores[upper][0]:
                lower += 1
            elif lower < upper:
                scores[lower], scores[upper] = scores[upper], scores[lower]
        return upper

    def update_if_better(self, state, newscore, manner, split=None, l1=None, l2=None):
        if newscore > state.score:
            state.set(newscore, manner, split, l1, l2)
            
    def sortM(self, threshold, beamstepM, sorted_bestM):
        sorted_items = []
        for i, state in beamstepM.items():
            if state.score > threshold:
                sorted_items.append((state.score, i))
        sorted_items.sort(reverse=True)
        sorted_bestM.clear()
        for score, i in sorted_items:
            sorted_bestM.append((score, i))

    def process_beam_H(self, j, next_pair):
        beamstepH = self.bestH[j]
        if len(beamstepH) > self.beam:
            self.beam_prune(beamstepH)

        nucj = self.nucs[j]
        jnext = next_pair[nucj][j]
        if jnext != -1 and (not self.no_sharp_turn or jnext - j >= 4):
            nucjnext = self.nucs[jnext]
            nucjnext_1 = self.nucs[jnext - 1] if jnext - 1 > -1 else -1
            newscore = self.v_score_hairpin(j, jnext, nucj, self.nucs[j + 1] if j + 1 < self.seq_length else -1, nucjnext_1, nucjnext) if self.use_vienna else self.score_hairpin(j, jnext, nucj, self.nucs[j + 1] if j + 1 < self.seq_length else -1, nucjnext_1, nucjnext)
            self.update_if_better(self.bestH[jnext][j], newscore, Manner.H)

        for i, state in beamstepH.items():
            nuci = self.nucs[i]
            jnext = next_pair[nuci][j]
            if jnext != -1:
                nuci1 = self.nucs[i + 1] if i + 1 < self.seq_length else -1
                nucjnext = self.nucs[jnext]
                nucjnext_1 = self.nucs[jnext - 1] if jnext - 1 > -1 else -1
                newscore = self.v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext) if self.use_vienna else self.score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext)
                self.update_if_better(self.bestH[jnext][i], newscore, Manner.H)
            self.update_if_better(self.bestP[j][i], state.score, Manner.HAIRPIN)

    def process_beam_Multi(self, j, next_pair):
        beamstepMulti = self.bestMulti[j]
        if len(beamstepMulti) > self.beam:
            self.beam_prune(beamstepMulti)

        for i, state in beamstepMulti.items():
            nuci = self.nucs[i]
            nuci1 = self.nucs[i + 1] if i + 1 < self.seq_length else -1
            jnext = next_pair[nuci][j]
            if jnext != -1:
                new_l1 = state.trace['paddings'][0]
                new_l2 = state.trace['paddings'][1] + jnext - j
                if new_l1 + new_l2 <= 30:
                    newscore = state.score - self.v_score_multi_unpaired(j, jnext - 1) if self.use_vienna else state.score + self.score_multi_unpaired(j, jnext - 1)
                    self.update_if_better(self.bestMulti[jnext][i], newscore, Manner.MULTI_eq_MULTI_plus_U, new_l1, new_l2)
            newscore = state.score - self.v_score_multi(i, j, nuci, nuci1, self.nucs[j - 1], nucj, self.seq_length) if self.use_vienna else state.score + self.score_multi(i, j, nuci, nuci1, self.nucs[j - 1], nucj, self.seq_length)
            self.update_if_better(self.bestP[j][i], newscore, Manner.P_eq_MULTI)

    def process_beam_P(self, j, next_pair):
        beamstepP = self.bestP[j]
        if len(beamstepP) > self.beam:
            self.beam_prune(beamstepP)

        use_cube_pruning = self.is_cube_pruning and self.beam > MIN_CUBE_PRUNING_SIZE and len(beamstepP) > MIN_CUBE_PRUNING_SIZE

        for i, state in beamstepP.items():
            nuci = self.nucs[i]
            nuci_1 = self.nucs[i - 1] if i - 1 > -1 else -1

            if i > 0 and j < self.seq_length - 1:
                precomputed = self.score_junction_B(j, i, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, nuci_1, nuci) if not self.use_vienna else 0
                for p in range(i - 1, max(i - 30, -1), -1):
                    nucp = self.nucs[p]
                    q = next_pair[nucp][j]
                    while q != -1 and (i - p) + (q - j) <= 30:
                        if p == i - 1 and q == j + 1:
                            newscore = -self.v_score_single(p, q, i, j, nucp, self.nucs[p + 1] if p + 1 < self.seq_length else -1, self.nucs[q - 1] if q - 1 > -1 else -1, self.nucs[q], nuci_1, nuci, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1) + state.score if self.use_vienna else self.score_helix(nucp, self.nucs[p + 1] if p + 1 < self.seq_length else -1, self.nucs[q - 1] if q - 1 > -1 else -1, self.nucs[q]) + state.score
                            self.update_if_better(self.bestP[q][p], newscore, Manner.HELIX)
                        else:
                            newscore = -self.v_score_single(p, q, i, j, nucp, self.nucs[p + 1] if p + 1 < self.seq_length else -1, self.nucs[q - 1] if q - 1 > -1 else -1, self.nucs[q], nuci_1, nuci, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1) + state.score if self.use_vienna else self.score_junction_B(p, q, nucp, self.nucs[p + 1] if p + 1 < self.seq_length else -1, self.nucs[q - 1] if q - 1 > -1 else -1, self.nucs[q]) + precomputed + self.score_single_without_junctionB(p, q, i, j, nuci_1, nuci, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1) + state.score
                            self.update_if_better(self.bestP[q][p], newscore, Manner.SINGLE, i - p, q - j)
                        q = next_pair[nucp][q]

            if i > 0 and j < self.seq_length - 1:
                newscore = -self.v_score_M1(i, j, j, nuci_1, nuci, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + state.score if self.use_vienna else self.score_M1(i, j, j, nuci_1, nuci, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + state.score
                self.update_if_better(self.bestM[i], newscore, Manner.M_eq_P)

            if not use_cube_pruning:
                k = i - 1
                if k > 0 and not self.bestM[k].empty():
                    M1_score = -self.v_score_M1(i, j, j, nuci_1, nuci, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + state.score if self.use_vienna else self.score_M1(i, j, j, nuci_1, nuci, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + state.score
                    if not self.is_candidate_list or self.bestM2[i].manner == Manner.NONE or M1_score > self.bestM2[i].score:
                        for m in self.bestM[k]:
                            newi = m.first
                            newscore = M1_score + m.second.score
                            self.update_if_better(self.bestM2[newi], newscore, Manner.M2_eq_M_plus_P, k)
            else:
                valid_Ps = []
                M1_scores = []
                for item in beamstepP:
                    i = item.first
                    state = item.second
                    nuci = self.nucs[i]
                    nuci_1 = self.nucs[i - 1] if i - 1 > -1 else -1
                    k = i - 1
                    if k > 0 and not self.bestM[k].empty():
                        M1_score = -self.v_score_M1(i, j, j, nuci_1, nuci, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + state.score if self.use_vienna else self.score_M1(i, j, j, nuci_1, nuci, self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + state.score
                        if not self.is_candidate_list or self.bestM2[i].manner == Manner.NONE or M1_score > self.bestM2[i].score:
                            valid_Ps.append(i)
                            M1_scores.append(M1_score)

                heap = []
                for p in range(len(valid_Ps)):
                    i = valid_Ps[p]
                    k = i - 1
                    heap.append((M1_scores[p] + self.sorted_bestM[k][0].first, (p, 0)))
                    heapq.heapify(heap)

                filled = 0
                prev_score = DOUBLE_MIN
                current_score = DOUBLE_MIN
                while filled < self.beam or current_score == prev_score:
                    top = heapq.heappop(heap)
                    prev_score = current_score
                    current_score = top[0]
                    index_P = top[1][0]
                    index_M = top[1][1]
                    i = valid_Ps[index_P]
                    k = i - 1
                    newi = self.sorted_bestM[k][index_M].second
                    newscore = M1_scores[index_P] + self.bestM[k][newi].score
                    if self.bestM2[newi].manner == Manner.NONE:
                        self.update_if_better(self.bestM2[newi], newscore, Manner.M2_eq_M_plus_P, k)
                        filled += 1
                    index_M += 1
                    while index_M < len(self.sorted_bestM[k]):
                        candidate_score = M1_scores[index_P] + self.sorted_bestM[k][index_M].first
                        candidate_newi = self.sorted_bestM[k][index_M].second
                        if self.bestM2[candidate_newi].manner == Manner.NONE:
                            heapq.heappush(heap, (candidate_score, (index_P, index_M)))
                            break
                        index_M += 1

            k = i - 1
            if k >= 0:
                prefix_C = self.bestC[k]
                if prefix_C.manner != Manner.NONE:
                    newscore = -self.v_score_external_paired(k + 1, j, self.nucs[k] if k > -1 else -1, self.nucs[k + 1], self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + prefix_C.score + state.score if self.use_vienna else self.score_external_paired(k + 1, j, self.nucs[k] if k > -1 else -1, self.nucs[k + 1], self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + prefix_C.score + state.score
                    self.update_if_better(self.beamstepC, newscore, Manner.C_eq_C_plus_P, k)
            else:
                newscore = -self.v_score_external_paired(0, j, -1, self.nucs[0], self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + state.score if self.use_vienna else self.score_external_paired(0, j, -1, self.nucs[0], self.nucs[j], self.nucs[j + 1] if j + 1 < self.seq_length else -1, self.seq_length) + state.score
                self.update_if_better(self.beamstepC, newscore, Manner.C_eq_C_plus_P, -1)

    def process_beam_M2(self, j, next_pair):
        beamstepM2 = self.bestM2[j]
        if len(beamstepM2) > self.beam:
            self.beam_prune(beamstepM2)

        for i, state in beamstepM2.items():
            for p in range(i - 1, max(i - 30, -1), -1):
                nucp = self.nucs[p]
                q = next_pair[nucp][j]
                while q != -1 and (i - p) + (q - j) <= 30:
                    newscore = state.score - self.v_score_multi_unpaired(p + 1, i - 1) - self.v_score_multi_unpaired(j + 1, q - 1) + state.score if self.use_vienna else state.score + self.score_multi_unpaired(p + 1, i - 1) + self.score_multi_unpaired(j + 1, q - 1) + state.score
                    self.update_if_better(self.bestMulti[q][p], newscore, Manner.MULTI, i - p, q - j)
                    q = next_pair[nucp][q]

            self.update_if_better(self.bestM[i], state.score, Manner.M_eq_M2)

    def process_beam_M(self, j, next_pair):
        beamstepM = self.bestM[j]
        threshold = DOUBLE_MIN
        if len(beamstepM) > self.beam:
            threshold = self.beam_prune(beamstepM)

        if self.is_cube_pruning:
            self.sortM(threshold, beamstepM, self.sorted_bestM[j])

        for i, state in beamstepM.items():
            if j < self.seq_length - 1:
                newscore = state.score - self.v_score_multi_unpaired(j + 1, j + 1) + state.score if self.use_vienna else state.score + self.score_multi_unpaired(j + 1, j + 1) + state.score
                self.update_if_better(self.bestM[j + 1][i], newscore, Manner.M_eq_M_plus_U)


    def process_beam_C(self, j, next_pair):
        if j < self.seq_length - 1:
            newscore = -self.v_score_external_unpaired(j + 1, j + 1) + self.bestC[j].score if self.use_vienna else self.score_external_unpaired(j + 1, j + 1) + self.bestC[j].score
            self.update_if_better(self.bestC[j + 1], newscore, Manner.C_eq_C_plus_U)


    def get_parentheses(self, seq):
        result = ['.' for _ in range(self.seq_length)]
        stack = [(0, self.seq_length - 1, self.bestC[self.seq_length - 1])]

        while stack:
            i, j, state = stack.pop()
            if state.manner == Manner.HAIRPIN:
                result[i] = '('
                result[j] = ')'
            elif state.manner == Manner.SINGLE:
                result[i] = '('
                result[j] = ')'
                p, q = i + state.trace['paddings'][0], j - state.trace['paddings'][1]
                stack.append((p, q, self.bestP[q][p]))
            elif state.manner == Manner.HELIX:
                result[i] = '('
                result[j] = ')'
                stack.append((i + 1, j - 1, self.bestP[j - 1][i + 1]))
            elif state.manner == Manner.MULTI:
                p, q = i + state.trace['paddings'][0], j - state.trace['paddings'][1]
                stack.append((p, q, self.bestM2[q][p]))
            elif state.manner == Manner.MULTI_eq_MULTI_plus_U:
                p, q = i + state.trace['paddings'][0], j - state.trace['paddings'][1]
                stack.append((p, q, self.bestM2[q][p]))
            elif state.manner == Manner.P_eq_MULTI:
                result[i] = '('
                result[j] = ')'
                stack.append((i, j, self.bestMulti[j][i]))
            elif state.manner == Manner.M2_eq_M_plus_P:
                k = state.trace['split']
                stack.append((i, k, self.bestM[k][i]))
                stack.append((k + 1, j, self.bestP[j][k + 1]))
            elif state.manner == Manner.M_eq_M2:
                stack.append((i, j, self.bestM2[j][i]))
            elif state.manner == Manner.M_eq_M_plus_U:
                stack.append((i, j - 1, self.bestM[j - 1][i]))
            elif state.manner == Manner.M_eq_P:
                stack.append((i, j, self.bestP[j][i]))
            elif state.manner == Manner.C_eq_C_plus_U:
                k = j - 1
                if k != -1:
                    stack.append((0, k, self.bestC[k]))
            elif state.manner == Manner.C_eq_C_plus_P:
                k = state.trace['split']
                if k != -1:
                    stack.append((0, k, self.bestC[k]))
                    stack.append((k + 1, j, self.bestP[j][k + 1]))
                else:
                    stack.append((i, j, self.bestP[j][i]))

        return ''.join(result)

    def GET_ACGU_NUM(self, ch):
        return 'ACGU'.find(ch)

    def v_score_external_unpaired(self, i, j):
        return 0.0

    def score_external_unpaired(self, i, j):
        return 0.0

    def v_score_hairpin(self, i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri):
        return 0.0

    def score_hairpin(self, i, j, nuci, nuci1, nucj_1, nucj):
        return 0.0

    def v_score_single(self, p, q, i, j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1):
        return 0.0

    def score_single(self, p, q, i, j, nucp, nucp1, nucq_1, nucq, nuci_1, nuci, nucj, nucj1):
        return 0.0

    def v_score_multi(self, i, j, nuci, nuci1, nucj_1, nucj, seq_length):
        return 0.0

    def score_multi(self, i, j, nuci, nuci1, nucj_1, nucj, seq_length):
        return 0.0

    def v_score_multi_unpaired(self, i, j):
        return 0.0

    def score_multi_unpaired(self, i, j):
        return 0.0

    def v_score_M1(self, i, j, j, nuci_1, nuci, nucj, nucj1, seq_length):
        return 0.0

    def score_M1(self, i, j, j, nuci_1, nuci, nucj, nucj1, seq_length):
        return 0.0

    def v_score_external_paired(self, k, j, nuck, nuck1, nucj, nucj1, seq_length):
        return 0.0

    def score_external_paired(self, k, j, nuck, nuck1, nucj, nucj1, seq_length):
        return 0.0

    def _allowed_pairs(self, nuci, nucj):
        return (nuci, nucj) in [(0, 3), (1, 2), (2, 1), (3, 0)]

def main():
    beamsize = 100
    use_vienna = False
    is_cube_pruning = True
    is_candidate_list = True
    sharpturn = False
    is_verbose = False
    regularprint = False

    for line in sys.stdin:
        if line.strip() == '' or line[0] in [';', '>']:
            print(line.strip())
            continue

        if not line[0].isalpha():
            print(f"Unrecognized sequence: {line.strip()}")
            continue

        if regularprint:
            print(f"seq:\n{line.strip()}")
        else:
            print(line.strip())

        parser = BeamCKYParser(beamsize, use_vienna, is_candidate_list, not sharpturn, is_cube_pruning, is_verbose)
        result = parser.parse(line.strip())

        if regularprint:
            if use_vienna:
                print(f"Energy(kcal/mol): {result['score'] / -100.0:.2f}")
            else:
                print(f"Viterbi score: {result['score']}")
            print(f"Time: {result['time']:.2f} len: {len(line.strip())} score {result['score']}")
            print(f">structure\n{result['structure']}\n")
        else:
            print(f"{result['structure']} ({result['score'] / -100.0:.2f})")

if __name__ == "__main__":
    main()