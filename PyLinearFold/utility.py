"""
utility.py
provides feature functions.

author: Kai Zhao, Dezhong Deng
edited by: 02/2018
"""
from fastcky_w import *

import numpy as np

INF = 1000000007

NOTON = 5  # NUM_OF_TYPE_OF_NUCS
NOTOND = 25
NOTONT = 125

EXPLICIT_MAX_LEN = 4
SINGLE_MIN_LEN = 0
SINGLE_MAX_LEN = 30  # NOTE: *must* <= sizeof(char), otherwise modify State::TraceInfo accordingly

HAIRPIN_MAX_LEN = 30
BULGE_MAX_LEN = SINGLE_MAX_LEN
INTERNAL_MAX_LEN = SINGLE_MAX_LEN
SYMMETRIC_MAX_LEN = 15
ASYMMETRY_MAX_LEN = 28

def GET_ACGU_NUM(x):
    return {'A': 0, 'C': 1, 'G': 2, 'U': 3}.get(x, 4)

_allowed_pairs = np.zeros((NOTON, NOTON), dtype=bool)
_helix_stacking = np.zeros((NOTON, NOTON, NOTON, NOTON), dtype=bool)
cache_single = np.zeros((SINGLE_MAX_LEN + 1, SINGLE_MAX_LEN + 1))

def initialize_cachesingle():
    global cache_single
    cache_single = np.zeros((SINGLE_MAX_LEN + 1, SINGLE_MAX_LEN + 1))
    for l1 in range(SINGLE_MIN_LEN, SINGLE_MAX_LEN + 1):
        for l2 in range(SINGLE_MIN_LEN, SINGLE_MAX_LEN + 1):
            if l1 == 0 and l2 == 0:
                continue

            # bulge
            if l1 == 0:
                cache_single[l1][l2] += bulge_length[l2]
            elif l2 == 0:
                cache_single[l1][l2] += bulge_length[l1]
            else:
                # internal
                cache_single[l1][l2] += internal_length[min(l1 + l2, INTERNAL_MAX_LEN)]

                # internal explicit
                if l1 <= EXPLICIT_MAX_LEN and l2 <= EXPLICIT_MAX_LEN:
                    cache_single[l1][l2] += internal_explicit[l1 <= l2 and l1 * EXPLICIT_MAX_LEN + l2 or l2 * EXPLICIT_MAX_LEN + l1]

                # internal symmetry
                if l1 == l2:
                    cache_single[l1][l2] += internal_symmetric_length[min(l1, SYMMETRIC_MAX_LEN)]

                else:  # internal asymmetry
                    diff = abs(l1 - l2)
                    cache_single[l1][l2] += internal_asymmetry[min(diff, ASYMMETRY_MAX_LEN)]
    return

def initialize():
    global _allowed_pairs, _helix_stacking
    _allowed_pairs[GET_ACGU_NUM('A')][GET_ACGU_NUM('U')] = True
    _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('A')] = True
    _allowed_pairs[GET_ACGU_NUM('C')][GET_ACGU_NUM('G')] = True
    _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('C')] = True
    _allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('U')] = True
    _allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('G')] = True

    HELIX_STACKING_OLD('A', 'U', 'A', 'U', True)
    HELIX_STACKING_OLD('A', 'U', 'C', 'G', True)
    HELIX_STACKING_OLD('A', 'U', 'G', 'C', True)
    HELIX_STACKING_OLD('A', 'U', 'G', 'U', True)
    HELIX_STACKING_OLD('A', 'U', 'U', 'A', True)
    HELIX_STACKING_OLD('A', 'U', 'U', 'G', True)
    HELIX_STACKING_OLD('C', 'G', 'A', 'U', True)
    HELIX_STACKING_OLD('C', 'G', 'C', 'G', True)
    HELIX_STACKING_OLD('C', 'G', 'G', 'C', True)
    HELIX_STACKING_OLD('C', 'G', 'G', 'U', True)
    HELIX_STACKING_OLD('C', 'G', 'U', 'G', True)
    HELIX_STACKING_OLD('G', 'C', 'A', 'U', True)
    HELIX_STACKING_OLD('G', 'C', 'C', 'G', True)
    HELIX_STACKING_OLD('G', 'C', 'G', 'U', True)
    HELIX_STACKING_OLD('G', 'C', 'U', 'G', True)
    HELIX_STACKING_OLD('G', 'U', 'A', 'U', True)
    HELIX_STACKING_OLD('G', 'U', 'G', 'U', True)
    HELIX_STACKING_OLD('G', 'U', 'U', 'G', True)
    HELIX_STACKING_OLD('U', 'A', 'A', 'U', True)
    HELIX_STACKING_OLD('U', 'A', 'G', 'U', True)
    HELIX_STACKING_OLD('U', 'G', 'G', 'U', True)

def HELIX_STACKING_OLD(x, y, z, w, value):
    global _helix_stacking
    _helix_stacking[GET_ACGU_NUM(x)][GET_ACGU_NUM(y)][GET_ACGU_NUM(z)][GET_ACGU_NUM(w)] = value

# ------------- nucs based scores -------------

def base_pair_score(nuci, nucj):
    return base_pair[nucj * NOTON + nuci]

def helix_stacking_score(nuci, nuci1, nucj_1, nucj):
    return helix_stacking[nuci * NOTONT + nucj * NOTOND + nuci1 * NOTON + nucj_1]

def helix_closing_score(nuci, nucj):
    return helix_closing[nuci * NOTON + nucj]

def terminal_mismatch_score(nuci, nuci1, nucj_1, nucj):
    return terminal_mismatch[nuci * NOTONT + nucj * NOTOND + nuci1 * NOTON + nucj_1]

def bulge_nuc_score(nuci):
    return bulge_0x1_nucleotides[nuci]

def internal_nuc_score(nuci, nucj):
    return internal_1x1_nucleotides[nuci * NOTON + nucj]

def dangle_left_score(nuci, nuci1, nucj):
    return dangle_left[nuci * NOTOND + nucj * NOTON + nuci1]

def dangle_right_score(nuci, nucj_1, nucj):
    return dangle_right[nuci * NOTOND + nucj * NOTON + nucj_1]

# ------------- length based scores -------------

def hairpin_score(i, j):
    return hairpin_length[min(j - i - 1, HAIRPIN_MAX_LEN)]

def internal_length_score(l):
    return internal_length[min(l, INTERNAL_MAX_LEN)]

def internal_explicit_score(l1, l2):
    l1_ = min(l1, EXPLICIT_MAX_LEN)
    l2_ = min(l2, EXPLICIT_MAX_LEN)
    return internal_explicit[l1_ <= l2_ and l1_ * NOTON + l2_ or l2_ * NOTON + l1_]

def internal_sym_score(l):
    return internal_symmetric_length[min(l, SYMMETRIC_MAX_LEN)]

def internal_asym_score(l1, l2):
    diff = abs(l1 - l2)
    return internal_asymmetry[min(diff, ASYMMETRY_MAX_LEN)]

def bulge_length_score(l):
    return bulge_length[min(l, BULGE_MAX_LEN)]

def hairpin_at_least_score(l):
    return hairpin_length_at_least[min(l, HAIRPIN_MAX_LEN)]

def buldge_length_at_least_score(l):
    return bulge_length_at_least[min(l, BULGE_MAX_LEN)]

def internal_length_at_least_score(l):
    return internal_length_at_least[min(l, INTERNAL_MAX_LEN)]

#-----------------------------------------------------
def score_junction_A(i, j, nuci, nuci1, nucj_1, nucj, len):
    return helix_closing_score(nuci, nucj) + \
           (0 if i >= len - 1 else dangle_left_score(nuci, nuci1, nucj)) + \
           (0 if j <= 0 else dangle_right_score(nuci, nucj_1, nucj))

def score_junction_B(i, j, nuci, nuci1, nucj_1, nucj):
    return helix_closing_score(nuci, nucj) + terminal_mismatch_score(nuci, nuci1, nucj_1, nucj)

def score_hairpin_length(len):
    return hairpin_length[min(len, HAIRPIN_MAX_LEN)]

def score_hairpin(i, j, nuci, nuci1, nucj_1, nucj):
    return hairpin_length[min(j - i - 1, HAIRPIN_MAX_LEN)] + \
           score_junction_B(i, j, nuci, nuci1, nucj_1, nucj)

def score_helix(nuci, nuci1, nucj_1, nucj):
    return helix_stacking_score(nuci, nuci1, nucj_1, nucj) + base_pair_score(nuci1, nucj_1)

def score_single_nuc(i, j, p, q, nucp_1, nucq1):
    l1 = p - i - 1
    l2 = j - q - 1
    if l1 == 0 and l2 == 1:
        return bulge_nuc_score(nucq1)
    if l1 == 1 and l2 == 0:
        return bulge_nuc_score(nucp_1)
    if l1 == 1 and l2 == 1:
        return internal_nuc_score(nucp_1, nucq1)
    return 0

def score_single(i, j, p, q, len, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1):
    l1 = p - i - 1
    l2 = j - q - 1
    return cache_single[l1][l2] + \
           base_pair_score(nucp, nucq) + \
           score_junction_B(i, j, nuci, nuci1, nucj_1, nucj) + \
           score_junction_B(q, p, nucq, nucq1, nucp_1, nucp) + \
           score_single_nuc(i, j, p, q, nucp_1, nucq1)

# score_single without socre_junction_B
def score_single_without_junctionB(i, j, p, q, nucp_1, nucp, nucq, nucq1):
    l1 = p - i - 1
    l2 = j - q - 1
    return cache_single[l1][l2] + \
           base_pair_score(nucp, nucq) + \
           score_single_nuc(i, j, p, q, nucp_1, nucq1)

def score_multi(i, j, nuci, nuci1, nucj_1, nucj, len):
    return score_junction_A(i, j, nuci, nuci1, nucj_1, nucj, len) + \
           multi_paired + multi_base

def score_multi_unpaired(i, j):
    return (j - i + 1) * multi_unpaired

def score_M1(i, j, k, nuci_1, nuci, nuck, nuck1, len):
    return score_junction_A(k, i, nuck, nuck1, nuci_1, nuci, len) + \
           score_multi_unpaired(k + 1, j) + base_pair_score(nuci, nuck) + multi_paired

def score_external_paired(i, j, nuci_1, nuci, nucj, nucj1, len):
    return score_junction_A(j, i, nucj, nucj1, nuci_1, nuci, len) + \
           external_paired + base_pair_score(nuci, nucj)

def score_external_unpaired(i, j):
    return (j - i + 1) * external_unpaired