"""
utility_v.py
provides feature functions for vienna model.

author: Kai Zhao, Dezhong Deng
edited by: 02/2018
"""

import math
from energy_parameter import *
from shared import *
from intl11 import *
from intl21 import *
from intl22 import *

# pairs: 0:NP 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA 7:NN
# nucleotides: CONTRAfold: 0:A 1:C 2:G 3:U 4:N ; Vienna: 0:N 1:A 2:C 3:G 4:U

def nuc_to_pair(x, y):
    if x == 1:
        return 5 if y == 4 else 0
    elif x == 2:
        return 1 if y == 3 else 0
    elif x == 3:
        return 2 if y == 2 else (3 if y == 4 else 0)
    elif x == 4:
        return 4 if y == 3 else (6 if y == 1 else 0)
    else:
        return 0

# Constants and parameters
MAXLOOP = 30

def min2(a, b):
    return a if a <= b else b

def max2(a, b):
    return a if a >= b else b

def v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops):
    # TetraLoops
    if_tetraloops.resize(max(0, seq_length - 5), -1)
    for i in range(seq_length - 5):
        if seq[i] == 'C' and seq[i + 5] == 'G':
            ts = Tetraloops.find(seq[i:i + 6])
            if ts != -1:
                if_tetraloops[i] = ts // 7

    # Triloops
    if_triloops.resize(max(0, seq_length - 4), -1)
    for i in range(seq_length - 4):
        if (seq[i] == 'C' and seq[i + 4] == 'G') or (seq[i] == 'G' and seq[i + 4] == 'C'):
            ts = Triloops.find(seq[i:i + 5])
            if ts != -1:
                if_triloops[i] = ts // 6

    # Hexaloops
    if_hexaloops.resize(max(0, seq_length - 7), -1)
    for i in range(seq_length - 7):
        if seq[i] == 'A' and seq[i + 7] == 'U':
            ts = Hexaloops.find(seq[i:i + 8])
            if ts != -1:
                if_hexaloops[i] = ts // 9

def v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri_index=-1):
    size = j - i - 1
    type = nuc_to_pair(nuci, nucj)
    energy = hairpin37[min(size, 30)] if size <= 30 else hairpin37[30] + int(lxc37 * math.log(size / 30.0))

    if size < 3:
        return energy

    if size == 4 and tetra_hex_tri_index > -1:
        return Tetraloop37[tetra_hex_tri_index]
    elif size == 6 and tetra_hex_tri_index > -1:
        return Hexaloop37[tetra_hex_tri_index]
    elif size == 3:
        if tetra_hex_tri_index > -1:
            return Triloop37[tetra_hex_tri_index]
        return energy + (TerminalAU37 if type > 2 else 0)

    energy += mismatchH37[type][nuci1][nucj_1]
    return energy

def v_score_single(i, j, p, q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1):
    type = nuc_to_pair(nuci, nucj)
    type_2 = nuc_to_pair(nucq, nucp)
    n1 = p - i - 1
    n2 = j - q - 1
    nl, ns = (n1, n2) if n1 > n2 else (n2, n1)
    energy = 0

    if nl == 0:
        return stack37[type][type_2]

    if ns == 0:
        energy = bulge37[min(nl, 30)] if nl <= 30 else bulge37[30] + int(lxc37 * math.log(nl / 30.0))
        if nl == 1:
            energy += stack37[type][type_2]
        else:
            if type > 2:
                energy += TerminalAU37
            if type_2 > 2:
                energy += TerminalAU37
        return energy

    if ns == 1:
        if nl == 1:
            return int11_37[type][type_2][nuci1][nucj_1]
        if nl == 2:
            if n1 == 1:
                energy = int21_37[type][type_2][nuci1][nucq1][nucj_1]
            else:
                energy = int21_37[type_2][type][nucq1][nuci1][nucp_1]
            return energy
        energy = internal_loop37[min(nl + 1, 30)] if nl + 1 <= 30 else internal_loop37[30] + int(lxc37 * math.log((nl + 1) / 30.0))
        energy += min(MAX_NINIO, (nl - ns) * ninio37)
        energy += mismatch1nI37[type][nuci1][nucj_1] + mismatch1nI37[type_2][nucq1][nucp_1]
        return energy

    if ns == 2:
        if nl == 2:
            return int22_37[type][type_2][nuci1][nucp_1][nucq1][nucj_1]
        if nl == 3:
            energy = internal_loop37[5] + ninio37
            energy += mismatch23I37[type][nuci1][nucj_1] + mismatch23I37[type_2][nucq1][nucp_1]
            return energy

    u = nl + ns
    energy = internal_loop37[min(u, 30)] if u <= 30 else internal_loop37[30] + int(lxc37 * math.log(u / 30.0))
    energy += min(MAX_NINIO, (nl - ns) * ninio37)
    energy += mismatchI37[type][nuci1][nucj_1] + mismatchI37[type_2][nucq1][nucp_1]
    return energy

def E_MLstem(type, si1, sj1, dangle_model):
    energy = 0
    if dangle_model != 0:
        if si1 >= 0 and sj1 >= 0:
            energy += mismatchM37[type][si1][sj1]

    if type > 2:
        energy += TerminalAU37

    energy += ML_intern37
    return energy

def v_score_M1(i, j, k, nuci_1, nuci, nuck, nuck1, len, dangle_model):
    p = i
    q = k
    tt = nuc_to_pair(nuci, nuck)
    return E_MLstem(tt, nuci_1, nuck1, dangle_model)

def v_score_multi_unpaired(i, j):
    return 0

def v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, len, dangle_model):
    tt = nuc_to_pair(nucj, nuci)
    return E_MLstem(tt, nucj_1, nuci1, dangle_model) + ML_closing37

def v_score_external_paired(i, j, nuci_1, nuci, nucj, nucj1, len, dangle_model):
    type = nuc_to_pair(nuci, nucj)
    energy = 0
    if dangle_model != 0:
        if nuci_1 >= 0 and nucj1 >= 0:
            energy += mismatchExt37[type][nuci_1][nucj1]
        elif nuci_1 >= 0:
            energy += dangle5_37[type][nuci_1]
        elif nucj1 >= 0:
            energy += dangle3_37[type][nucj1]

    if type > 2:
        energy += TerminalAU37
    return energy

def v_score_external_unpaired(i, j):
    return 0