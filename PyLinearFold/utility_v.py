"""
 *utility_v.h*
 provides feature functions for vienna model.

 author: Kai Zhao, Dezhong Deng
 edited by: 02/2018
"""
from energy_par import (
    bulge37, hairpin37, stack37, lxc37,  
    Tetraloops, Triloops, Hexaloops,
    Tetraloop37, Triloop37, Hexaloop37, 
    TerminalAU37,
    internal_loop37, 
    mismatchH37, mismatch1nI37, mismatch23I37,
    mismatchI37,
    mismatchM37,
    MAX_NINIO, ninio37,
    dangle5_37,
    dangle3_37,
    ML_intern37,
    ML_closing37,
    mismatchExt37
)
from intl11 import *
from intl21 import *
from intl22 import *

from math import log

MAXLOOP = 30

# 定义NUM_TO_NUC宏的Python函数
def NUM_TO_NUC(x):
    return -1 if x == -1 else (0 if x == 4 else x + 1)

# 定义NUM_TO_PAIR宏的Python函数
def NUM_TO_PAIR(x, y):
    if x == 0:
        return 5 if y == 3 else 0
    if x == 1:
        return 1 if y == 2 else 0
    if x == 2:
        return 2 if y == 1 else 3 if y == 3 else 0
    if x == 3:
        return 4 if y == 2 else 6 if y == 0 else 0
    return 0  # 默认情况，如果x不是0, 1, 2, 或3

# 定义NUC_TO_PAIR宏的Python函数
def NUC_TO_PAIR(x, y):
    if x == 1:
        return 5 if y == 4 else 0
    if x == 2:
        return 1 if y == 3 else 0
    if x == 3:
        return 2 if y == 2 else 3 if y == 4 else 0
    if x == 4:
        return 4 if y == 3 else 6 if y == 1 else 0
    return 0  # 默认情况，如果x不是1, 2, 3, 或4


# 定义MIN2和MAX2的内联函数
def MIN2(a, b):
    return a if a < b else b

def MAX2(a, b):
    return a if a > b else b


def v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops):
    # TetraLoops
    if_tetraloops = [-1] * max(0, seq_length - 5)
    for i in range(seq_length - 5):
        if seq[i] != 'C' or seq[i + 5] != 'G':
            continue
        tl = seq[i:i + 6]
        ts = Tetraloops.find(tl)
        if ts != -1:
            if_tetraloops[i] = ts // 7

    # Triloops
    if_triloops = [-1] * max(0, seq_length - 4)
    for i in range(seq_length - 4):
        if not ((seq[i] == 'C' and seq[i + 4] == 'G') or (seq[i] == 'G' and seq[i + 4] == 'C')):
            continue
        tl = seq[i:i + 5]
        ts = Triloops.find(tl)
        if ts != -1:
            if_triloops[i] = ts // 6

    # Hexaloops
    if_hexaloops = [-1] * max(0, seq_length - 7)
    for i in range(seq_length - 7):
        if seq[i] != 'A' or seq[i + 7] != 'U':
            continue
        tl = seq[i:i + 8]
        ts = Hexaloops.find(tl)
        if ts != -1:
            if_hexaloops[i] = ts // 9

    return if_tetraloops, if_hexaloops, if_triloops



def v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri_index=-1):
    size = j - i - 1
    type = NUM_TO_PAIR(nuci, nucj)
    si1 = NUM_TO_NUC(nuci1)
    sj1 = NUM_TO_NUC(nucj_1)

    if size <= 30:
        energy = hairpin37[size]
    else:
        energy = hairpin37[30] + int(lxc37 * log(size / 30.0))

    if size < 3:
        return energy  # should only be the case when folding alignments

    # if(special_hp):
    if size == 4 and tetra_hex_tri_index > -1:
        return Tetraloop37[tetra_hex_tri_index]
    elif size == 6 and tetra_hex_tri_index > -1:
        return Hexaloop37[tetra_hex_tri_index]
    elif size == 3:
        if tetra_hex_tri_index > -1:
            return Triloop37[tetra_hex_tri_index]
        return energy + (0 if type <= 2 else TerminalAU37)

    energy += mismatchH37[type][si1][sj1]

    return energy

def v_score_single(i, j, p, q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1):
    si1 = NUM_TO_NUC(nuci1)
    sj1 = NUM_TO_NUC(nucj_1)
    sp1 = NUM_TO_NUC(nucp_1)
    sq1 = NUM_TO_NUC(nucq1)
    type = NUM_TO_PAIR(nuci, nucj)
    type_2 = NUM_TO_PAIR(nucq, nucp)
    n1 = p - i - 1
    n2 = j - q - 1
    nl, ns, u, energy = 0, 0, 0, 0

    if n1 > n2:
        nl = n1
        ns = n2
    else:
        nl = n2
        ns = n1

    if nl == 0:
        return stack37[type][type_2]  # stack

    if ns == 0:  # bulge
        energy = bulge37[nl] if nl <= MAXLOOP else bulge37[30] + int(lxc37 * log(nl / 30.0))
        if nl == 1:
            energy += stack37[type][type_2]
        else:
            if type > 2:
                energy += TerminalAU37
            if type_2 > 2:
                energy += TerminalAU37
        return energy
    else:  # interior loop
        if ns == 1:
            if nl == 1:  # 1x1 loop
                return int11_37[type][type_2][si1][sj1]
            if nl == 2:  # 2x1 loop
                if n1 == 1:
                    energy = int21_37[type][type_2][si1][sq1][sj1]
                else:
                    energy = int21_37[type_2][type][sq1][si1][sp1]
                return energy
            else:  # 1xn loop
                energy = internal_loop37[nl + 1] if nl + 1 <= MAXLOOP else internal_loop37[30] + int(lxc37 * log((nl + 1) / 30.0))
                energy += min(MAX_NINIO, (nl - ns) * ninio37)
                energy += mismatch1nI37[type][si1][sj1] + mismatch1nI37[type_2][sq1][sp1]
                return energy
        elif ns == 2:
            if nl == 2:  # 2x2 loop
                return int22_37[type][type_2][si1][sp1][sq1][sj1]
            elif nl == 3:  # 2x3 loop
                energy = internal_loop37[5] + ninio37
                energy += mismatch23I37[type][si1][sj1] + mismatch23I37[type_2][sq1][sp1]
                return energy
        else:  # generic interior loop (no else here!)
            u = nl + ns
            energy = internal_loop37[u] if u <= MAXLOOP else internal_loop37[30] + int(lxc37 * log(u / 30.0))
            energy += min(MAX_NINIO, (nl - ns) * ninio37)
            energy += mismatchI37[type][si1][sj1] + mismatchI37[type_2][sq1][sp1]
            return energy


def E_MLstem(type, si1, sj1):
    energy = 0

    if si1 >= 0 and sj1 >= 0:
        energy += mismatchM37[type][si1][sj1]
    elif si1 >= 0:
        energy += dangle5_37[type][si1]
    elif sj1 >= 0:
        energy += dangle3_37[type][sj1]

    if type > 2:
        energy += TerminalAU37

    energy += ML_intern37

    return energy

def v_score_M1(i, j, k, nuci_1, nuci, nuck, nuck1, len):
    p = i
    q = k
    tt = NUM_TO_PAIR(nuci, nuck)
    sp1 = NUM_TO_NUC(nuci_1)
    sq1 = NUM_TO_NUC(nuck1)

    return E_MLstem(tt, sp1, sq1)

def v_score_multi_unpaired(i, j):
    return 0

def v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, len):
    tt = NUM_TO_PAIR(nucj, nuci)
    si1 = NUM_TO_NUC(nuci1)
    sj1 = NUM_TO_NUC(nucj_1)

    return E_MLstem(tt, sj1, si1) + ML_closing37

def v_score_external_paired(i, j, nuci_1, nuci, nucj, nucj1, len):
    type = NUM_TO_PAIR(nuci, nucj)
    si1 = NUM_TO_NUC(nuci_1)
    sj1 = NUM_TO_NUC(nucj1)
    energy = 0

    if si1 >= 0 and sj1 >= 0:
        energy += mismatchExt37[type][si1][sj1]
    elif si1 >= 0:
        energy += dangle5_37[type][si1]
    elif sj1 >= 0:
        energy += dangle3_37[type][sj1]

    if type > 2:
        energy += TerminalAU37

    return energy

def v_score_external_unpaired(i, j):
    return 0