from typing import List
from dataclasses import dataclass
import math
import numpy as np
from enum import Enum
import random
from Crypto.Hash import keccak
from project.hash import *
from project.polynomial import *
from project.parameters import *

# parameters
#ELLBITS = 256
#ELLSIZE = int(ELLBITS / 8)
#RBITS = 40973
#DV = 137
#T1 = 264
#NBITS = RBITS * 2
#RSIZE = math.ceil(RBITS / 8)
#NSIZE = math.ceil(NBITS / 8)
#tau = 3
#nb_iter = 5
#R_DQWORDS = math.ceil(RSIZE / 16)
#SHAKE256_BLOCK_SIZE = 136
#SHAKE256_STATE_SIZE = 200


# datatypes
@dataclass
class PublicKey:
    value: np.ndarray = np.zeros(RSIZE, dtype=np.uint8)
    raw: np.ndarray = np.zeros(RSIZE, dtype=np.uint8)


@dataclass
class SecretKey:
    val0: np.ndarray = np.zeros(RSIZE, dtype=np.uint8)
    val1: np.ndarray = np.zeros(RSIZE, dtype=np.uint8)
    raw: np.ndarray = np.zeros(2 * RSIZE, dtype=np.uint8)
    sigma: np.ndarray = np.zeros(ELLSIZE, dtype=np.uint8)


@dataclass
class CipherText:
    val0: np.ndarray = np.zeros(RSIZE, dtype=np.uint8)
    val1: np.ndarray = np.zeros(ELLSIZE, dtype=np.uint8)
    raw: np.ndarray = np.zeros(RSIZE + ELLSIZE, dtype=np.uint8)


@dataclass
class Syndrome:
    raw: np.ndarray = np.zeros(RBITS, dtype=np.uint8)


@dataclass
class SharedSecret:
    raw: np.ndarray = np.zeros(RBITS, dtype=np.uint8)


@dataclass
class Seed:
    raw: np.ndarray = np.zeros(32, dtype=np.uint8)
    qwords: np.ndarray = np.zeros(4, dtype=np.uint8)


@dataclass
class DoubleSeed:
    s1: Seed = Seed()
    s2: Seed = Seed()
    raw: np.ndarray = np.zeros(64, dtype=np.uint8)


@dataclass
class SHAstate:
    buffer: np.ndarray = np.zeros(SHAKE256_STATE_SIZE, dtype=np.uint8)
    pos: int = 0


class SeedID(Enum):
    G_SEED = 0
    H_SEED = 1
    M_SEED = 2
    E_SEED = 3


class SeedPurpose(Enum):
    KEYGEN_SEEDS = 0
    ENCAPS_SEEDS = 1
    DECAPS_SEEDS = 2


# utility functions
def var_th_fct(x):
    return max(17.8785 + 0.00402312 * x, 69)


def bit(len: int):
    return 1 << len


def mask(len: int):
    return bit(len) - 1


def bit_scan_reverse(val: int):
    index = 0

    while (val != 0):
        val >>= 1
        index += 1

    return index


def check_bit(array, position):
    idx = int(position // 8)
    pos = position % 8
    return (array[idx] >> pos) & 0x01


def set_bit(array, position):
    idx = int(position // 8)
    pos = position % 8
    array[idx] = array[idx] | 1 << pos


def get_seeds(seeds_type: SeedPurpose) -> DoubleSeed:
    d_seeds = DoubleSeed()
    for i in range(0, 32):
        d_seeds.s1.raw[i] = np.uint8(random.randint(0, 256))
        d_seeds.s2.raw[i] = np.uint8(random.randint(0, 256))
    return d_seeds


def convert2compact(a: np.ndarray) -> np.ndarray:
    out = np.zeros(DV, dtype=np.uint8)
    idx = 0
    for i in range(RSIZE):
        for j in range(8):
            if i * 8 + j == RBITS:
                break
            if a[i] >> j & 1:
                out[idx] = i * 8 + j
                idx += 1
    return out


def convertbin2byte(in_array: np.ndarray, length) -> np.ndarray:
    padding_len = length % 8
    num_bytes = length // 8
    if padding_len == 0:
        num_bytes = length // 8
    else:
        num_bytes = 1 + length // 8
    out = np.zeros(num_bytes, dtype=np.uint8)

    for i in range(num_bytes):
        for j in range(8):
            if (i * 8 + j) == length:
                break
            if in_array[i * 8 + j]:
                out[i] = out[i] | (1 << j)
    return out


def convertbyte2bin(in_array: np.ndarray, length: int) -> np.ndarray:
    paddingLen = length % 8
    num_bytes = len(in_array)
    out = np.zeros(length, dtype=np.uint8)
    if paddingLen == 0:
        num_bytes = length // 8
    else:
        num_bytes = 1 + length // 8

    for i in range(num_bytes):
        for j in range(8):
            if i * 8 + j == length:
                break
            if (in_array[i] >> j) & 1:
                out[i * 8 + j] = 1

    return out


def transpose(in_array: np.ndarray) -> np.ndarray:
    out = np.zeros(RBITS, dtype=np.uint8)
    out[0] = in_array[0]
    for i in range(1, RBITS):
        out[i] = in_array[RBITS - i]
    return out


def get_hamming_weight(in_array: np.ndarray, len) -> int:
    count = 0
    for i in range(0, len):
        count += in_array[i]
    return count


def ctr(h_compact_col: np.ndarray, pos: int, s: np.ndarray) -> int:
    count = 0
    for i in range(0, DV):
        if s[(h_compact_col[i] + pos) % RBITS]:
            count += 1
    return count


# random generator
def shake256_init(in_char, in_len) -> SHAstate:
    s = SHAstate()
    r = 1088
    c = 256
    sfx = 0x1F
    R = int(r / 8)
    i, b = 0, 0
    while (in_len > 0):
        if in_len < R:
            b = in_len
        else:
            b = R
        for i in range(0, b):
            s.buffer[i] ^= in_char[i]
        in_char += b
        in_len -= b
        if b == R:
            s.buffer = KeccakF1600(bytearray(s.buffer))
            b = 0
    s.buffer[b] ^= sfx
    if ((sfx & 0x80) and (b == (R - 1))):
        s.buffer = KeccakF1600(bytearray(s.buffer))
    s.buffer[R - 1] ^= 0x80
    s.buffer = KeccakF1600(bytearray(s.buffer))
    s.pos = 0
    return s


def shake256_prng(state: SHAstate, len):
    if len + state.pos <= SHAKE256_BLOCK_SIZE:
        res = state.buffer[state.pos:len]
        state.pos += len
        return res
    else:
        idx = SHAKE256_BLOCK_SIZE - state.pos
        a = state.buffer[state.pos:idx]
        state.pos = 0
        state.buffer = KeccakF1600(bytearray(state.buffer))
        state.pos = len - idx
        a = a + state.buffer[:state.pos]
        return a


def get_rand_mod_len_keccak(len, state: SHAstate):
    res = shake256_prng(state, 32 // 8)
    r = int.from_bytes(res, 'big')
    res = (r * len)
    return res >> 32


def generate_sparse_rep_keccak(weight, len, shake256_state):
    rand_pos = 0
    r = np.zeros(math.ceil(len / 8), dtype=np.uint8)
    for i in range(weight - 1, -1, -1):
        rand_pos = get_rand_mod_len_keccak(len - i, shake256_state)
        rand_pos += 1
        if check_bit(r, rand_pos):
            rand_pos = i
        set_bit(r, rand_pos)
    return r


# hash functions
def function_h(m: np.ndarray) -> np.ndarray:
    seed_for_hash = Seed()
    seed_for_hash.raw = m[:ELLSIZE]
    sha_state = shake256_init(seed_for_hash.raw, ELLSIZE)
    e = generate_sparse_rep_keccak(T1, NBITS, sha_state)
    return e


def function_l(e: np.ndarray) -> np.ndarray:
    e0, e1 = split_polynomial(e)
    hash_value = SHA3_384(bytearray(np.concatenate((e0, e1))))
    hv = hash_value[:ELLSIZE]
    return np.array([x for x in hv])


def function_k(m: np.ndarray, c0: np.ndarray, c1: np.ndarray) -> np.ndarray:
    tmp = np.concatenate((m, c0, c1))
    hash = SHA3_384(tmp)
    h = hash[:ELLSIZE]
    return np.array([x for x in h])


# syndrome
def compute_syndrome(ct: CipherText, sk: SecretKey) -> Syndrome:
    syndrome = Syndrome()
    s0 = mod_mul(sk.val0, ct.val0)
    bytes = convertbyte2bin(s0, RBITS)
    syndrome.raw = transpose(bytes)
    return syndrome


def recompute_syndrome(syndrome: np.ndarray, pos: int, h0_compact: np.ndarray, h1_compact: np.ndarray) -> np.ndarray:
    s = Syndrome()
    if pos < RBITS:
        for j in range(0, DV):
            if h0_compact[j] <= pos:
                s.raw[pos - h0_compact[j]] = syndrome[pos - h0_compact[j]] ^ 1
            else:
                s.raw[RBITS - h0_compact[j]] = syndrome[RBITS - h0_compact[j]] ^ 1
    else:
        for j in range(0, DV):
            if h1_compact[j] <= pos - RBITS:
                s.raw[(pos - RBITS) - h1_compact[j]] = syndrome[(pos - RBITS) - h1_compact[j]] ^ 1
            else:
                s.raw[RBITS - h1_compact[j] + (pos - RBITS)] = syndrome[RBITS - h1_compact[j] + (pos - RBITS)] ^ 1
    return s.raw


# threshold
def threshold(s, len):
    x = get_hamming_weight(s, len)
    return np.floor(var_th_fct(x))


# decoder
def getCol(row: np.array):
    col = np.zeros(DV, dtype=np.uint8)
    if row[0] == 0:
        col[0] = 0
        for i in range(1, DV):
            col[i] = RBITS - row[DV - i]
    else:
        for i in range(0, DV):
            col[i] = RBITS - row[DV - 1 - i]
    return col


def flip_adjusted_error_position(e, pos):
    adjusted_position = pos
    if pos != 0 and pos != RBITS:
        if pos > RBITS:
            adjusted_position = (NBITS - pos) + RBITS
        else:
            adjusted_position = RBITS - pos
    e[adjusted_position] = e[adjusted_position] ^ 1
    return e


def BFIter(e: np.ndarray, black: np.ndarray, gray: np.ndarray, s: np.ndarray, T, h0_compact, h1_compact, col_h0_compact,
           col_h1_compact):
    pos = np.zeros(2 * RBITS, dtype=np.uint8)
    for j in range(RBITS):
        counter = ctr(col_h0_compact, j, s)
        if counter >= T:
            e = flip_adjusted_error_position(e, j)
            pos[j] = 1
            black[j] = 1
        elif counter >= T - tau:
            gray[j] = 1
    for j in range(RBITS):
        counter = ctr(col_h1_compact, j, s)
        if counter >= T:
            e = flip_adjusted_error_position(e, RBITS + j)
            pos[RBITS + j] = 1
            black[RBITS + j] = 1
        elif counter >= T - tau:
            gray[RBITS + j] = 1
    for j in range(2 * RBITS):
        if pos[j] == 1:
            s = recompute_syndrome(s, j, h0_compact, h1_compact)
    return e, black, gray, s


def BFMaskedIter(e, s, mask, T, h0_compact, h1_compact, col_h0_compact, col_h1_compact):
    pos = np.zeros(2 * RBITS)
    for j in range(RBITS):
        counter = ctr(col_h0_compact, j, s)
        if counter >= T and mask[j]:
            e = flip_adjusted_error_position(e, j)
            pos[j] = 1
    for j in range(RBITS):
        counter = ctr(col_h1_compact, j, s)
        if counter >= T and mask[RBITS + j]:
            e = flip_adjusted_error_position(e, RBITS + j)
            pos[RBITS + j] = 1
    for j in range(RBITS * 2):
        if pos[j] == 1:
            s = recompute_syndrome(s, j, h0_compact, h1_compact)
    return e, s


def BGF_decoder(syndrome: Syndrome, h0_compact: np.ndarray, h1_compact: np.ndarray):
    e = np.zeros(2 * RBITS, dtype=np.uint8)
    compact_col_h0 = getCol(h0_compact)
    compact_col_h1 = getCol(h1_compact)

    black = np.zeros(2 * RBITS, dtype=np.uint8)
    gray = np.zeros(2 * RBITS, dtype=np.uint8)

    for i in range(5):
        black = np.zeros(2 * RBITS, dtype=np.uint8)
        gray = np.zeros(2 * RBITS, dtype=np.uint8)

        T = math.floor(var_th_fct(get_hamming_weight(syndrome.raw, RBITS)))

        e, black, gray, syndrome.raw = BFIter(e, black, gray, syndrome.raw, T, h0_compact, h1_compact, compact_col_h0,
                                              compact_col_h1)
        if i == 1:
            e, syndrome.raw = BFMaskedIter(e, syndrome.raw, black, (DV + 1) / 2 + 1, h0_compact, h1_compact,
                                           compact_col_h0,
                                           compact_col_h1)
            e, syndrome.raw = BFMaskedIter(e, syndrome.raw, gray, (DV + 1) / 2 + 1, h0_compact, h1_compact,
                                           compact_col_h0,
                                           compact_col_h1)

    res = 0
    if get_hamming_weight(syndrome.raw, RBITS) == 0:
        res = 0
    else:
        res = 1
    return res, e


# kem mechanism
def generate_keys() -> (PublicKey, SecretKey):
    cond = True
    while cond:
        double_seeds = get_seeds(SeedPurpose.KEYGEN_SEEDS)
        state = shake256_init(double_seeds.s1.raw, ELLSIZE)
        h0 = generate_sparse_rep_keccak(DV, RBITS, state)
        h1 = generate_sparse_rep_keccak(DV, RBITS, state)
        d, inv_h0 = mod_inv(h0)
        if d == 1:
            cond = False

    h = mod_mul(h1, inv_h0)
    sigma = double_seeds.s2.raw

    return PublicKey(value=np.array(h), raw=np.array(h)), SecretKey(val0=np.array(h0), val1=np.array(h1),
                                                                    raw=np.concatenate((h0, h1)),
                                                                    sigma=sigma)


def encapsulate(pk: PublicKey) -> (CipherText, SharedSecret):
    double_seeds = get_seeds(SeedPurpose.ENCAPS_SEEDS)

    m = double_seeds.s1.raw

    e = function_h(m)
    e0, e1 = split_polynomial(e)

    c0 = poly_add(e0, mod_mul(e1, pk.value))
    c1 = np.zeros(ELLSIZE, dtype=np.uint8)
    tl = function_l(e)
    for i in range(ELLSIZE):
        c1[i] = tl[i] ^ m[i]

    k = function_k(m, c0, c1)

    return CipherText(val0=c0, val1=c1, raw=np.concatenate((c0, c1))), SharedSecret(raw=k)


def decapsulate(ct: CipherText, sk: SecretKey) -> SharedSecret:
    h0_compact = convert2compact(sk.val0)
    h1_compact = convert2compact(sk.val1)

    syndrome = compute_syndrome(ct, sk)
    res_e, e_tmp = BGF_decoder(syndrome, h0_compact, h1_compact)

    e_prime = convertbin2byte(e_tmp, 2 * RBITS)
    e_l = function_l(e_prime)

    ss = SharedSecret()

    m_prime = np.zeros(ELLSIZE, dtype=np.uint8)
    for i in range(ELLSIZE):
        m_prime[i] = ct.val1[i] ^ e_l[i]

    e = function_h(m_prime)

    eq = True
    for i in range(len(e)):
        if e[i] != e_prime[i]:
            eq = False

    if eq:
        k = function_k(m_prime, ct.val0, ct.val1)
        ss.raw = k
    else:
        k = function_k(sk.sigma, ct.val0, ct.val1)
        ss.raw = k
    return ss
