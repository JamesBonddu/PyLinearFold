from enum import Enum
from typing import List, Dict, Union, Tuple
import sys
import time

DOUBLE_MIN = -sys.float_info.max
MIN_CUBE_PRUNING_SIZE = 20

class Manner(Enum):
    MANNER_NONE = 0  # 0: empty
    MANNER_H = 1     # 1: hairpin candidate
    MANNER_HAIRPIN = 2  # 2: hairpin
    MANNER_SINGLE = 3  # 3: single
    MANNER_HELIX = 4  # 4: helix
    MANNER_MULTI = 5  # 5: multi = ..M2. [30 restriction on the left and jump on the right]
    MANNER_MULTI_eq_MULTI_plus_U = 6  # 6: multi = multi + U
    MANNER_P_eq_MULTI = 7  # 7: P = (multi)
    MANNER_M2_eq_M_plus_P = 8  # 8: M2 = M + P
    MANNER_M_eq_M2 = 9  # 9: M = M2
    MANNER_M_eq_M_plus_U = 10  # 10: M = M + U
    MANNER_M_eq_P = 11  # 11: M = P
    MANNER_C_eq_C_plus_U = 12  # 12: C = C + U
    MANNER_C_eq_C_plus_P = 13  # 13: C = C + P


class State:
    def __init__(self, score=DOUBLE_MIN, manner=Manner.MANNER_NONE):
        self.score = score
        self.manner = manner
        self.trace = {"split": None, "paddings": None}

    def set(self, score, manner):
        self.score = score
        self.manner = manner

    def set_with_split(self, score, manner, split):
        self.score = score
        self.manner = manner
        self.trace["split"] = split
        self.trace["paddings"] = None

    def set_with_paddings(self, score, manner, l1, l2):
        self.score = score
        self.manner = manner
        self.trace["split"] = None
        self.trace["paddings"] = (l1, l2)

    def get_trace(self):
        if self.trace["split"] is not None:
            return self.trace["split"]
        elif self.trace["paddings"] is not None:
            return self.trace["paddings"]
        return None


class BeamCKYParser:
    def __init__(self, beam_size=0, 
                use_vienna=False, candidate_list=True, no_sharp_turn=True, 
                cube_pruning=True, verbose=False):
        self.beam = beam_size
        self.use_vienna = use_vienna
        self.is_candidate_list = candidate_list
        self.no_sharp_turn = no_sharp_turn
        self.is_cube_pruning = cube_pruning
        self.is_verbose = verbose
        self.initialize()
        if not use_vienna:
            self.initialize_cachesingle()

        self.seq_length = 0
        self.bestH = []
        self.bestP = []
        self.bestM2 = []
        self.bestMulti = []
        self.bestM = []
        self.sorted_bestM = []
        self.bestC = []
        self.nucs = []
        self.if_tetraloops = []
        self.if_hexaloops = []
        self.if_triloops = []
        self.scores = []

    def parse(self, seq: str) -> Dict:
        self.prepare(len(seq))
        # Implement the parsing logic here
        # This is a placeholder for the actual parsing logic
        structure = ""
        score = 0.0
        num_states = 0
        time_taken = 0.0
        return {
            "structure": structure,
            "score": score,
            "num_states": num_states,
            "time": time_taken
        }

    def prepare(self, length: int):
        self.seq_length = length
        self.bestH = [{} for _ in range(length)]
        self.bestP = [{} for _ in range(length)]
        self.bestM2 = [{} for _ in range(length)]
        self.bestMulti = [{} for _ in range(length)]
        self.bestM = [{} for _ in range(length)]
        self.sorted_bestM = [[] for _ in range(length)]
        self.bestC = [State() for _ in range(length + 1)]
        self.nucs = [ord(nuc) for nuc in seq]
        self.if_tetraloops = [0] * length
        self.if_hexaloops = [0] * length
        self.if_triloops = [0] * length
        self.scores = []

    def update_if_better(self, state: State, newscore: float, manner: Manner):
        if state.score < newscore or state.manner == Manner.MANNER_NONE:
            state.set(newscore, manner)

    def update_if_better_with_split(self, state: State, newscore: float, manner: Manner, split: int):
        if state.score < newscore or state.manner == Manner.MANNER_NONE:
            state.set_with_split(newscore, manner, split)

    def update_if_better_with_paddings(self, state: State, newscore: float, manner: Manner, l1: int, l2: int):
        if state.score < newscore or state.manner == Manner.MANNER_NONE:
            state.set_with_paddings(newscore, manner, l1, l2)

    def beam_prune(self, beamstep: Dict[int, State]) -> float:
        scores = [(state.score, key) for key, state in beamstep.items()]
        scores.sort(reverse=True)
        pruning_score = scores[self.beam][0] if self