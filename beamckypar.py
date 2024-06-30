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

    def process_beam_H(self, j, next_pair):
        pass

    def process_beam_Multi(self, j, next_pair):
        pass

    def process_beam_P(self, j, next_pair):
        pass

    def process_beam_M2(self, j, next_pair):
        pass

    def process_beam_M(self, j, next_pair):
        pass

    def process_beam_C(self, j, next_pair):
        pass

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