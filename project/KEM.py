from typing import List
from dataclasses import dataclass
import math
import numpy as np
from enum import Enum
import random

# parameters
ELLBITS = 256
ELLSIZE = int(ELLBITS / 8)
RBITS = 40973
DV = 137
T1 = 264
NBITS = RBITS * 2
RSIZE = math.ceil(RBITS / 8)
NSIZE = math.ceil(NBITS / 8)
tau = 3
nb_iter = 5
R_DQWORDS = math.ceil(RSIZE / 16)
SHAKE256_BLOCK_SIZE = 136
SHAKE256_STATE_SIZE = 200


# datatypes
@dataclass
class PublicKey:
    value: np.ndarray = np.zeros(RSIZE)
    raw: np.ndarray = np.zeros(RSIZE)


@dataclass
class SecretKey:
    val0: np.ndarray = np.zeros(RSIZE)
    val1: np.ndarray = np.zeros(RSIZE)
    raw: np.ndarray = np.zeros(2 * RSIZE)
    sigma: np.ndarray = np.zeros(ELLSIZE)


@dataclass
class CipherText:
    val0: np.ndarray = np.zeros(RSIZE)
    val1: np.ndarray = np.zeros(ELLSIZE)
    raw: np.ndarray = np.zeros(RSIZE + ELLSIZE)


@dataclass
class Syndrome:
    raw: np.ndarray = np.zeros(RBITS)


@dataclass
class SharedSecret:
    raw: np.ndarray = np.zeros(RBITS)


@dataclass
class Seed:
    raw: np.ndarray = np.zeros(32)
    qwords: np.ndarray = np.zeros(4)


@dataclass
class DoubleSeed:
    s1: Seed = Seed()
    s2: Seed = Seed()
    raw: np.ndarray = np.zeros(64)


@dataclass
class SHAstate:
    buffer: np.ndarray = np.zeros(SHAKE256_STATE_SIZE)
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


def get_seeds(seeds_type: SeedPurpose) -> DoubleSeed:
    d_seeds = DoubleSeed()
    for i in range(0, 32):
        d_seeds.s1.raw = random.randint(0, 256)
        d_seeds.s2.raw = random.randint(0, 256)
    return d_seeds


def convert2compact(a: np.ndarray) -> np.ndarray:
    out = np.zeros(DV)
    idx = 0
    for i in range(RSIZE):
        for j in range(8):
            if i*8 + j == RBITS:
                break
            if a[i] >> j & 1:
                out[idx+1] = i*8+j
    return out


def convertbin2byte(e: np.ndarray, len) -> np.ndarray:
    out = np.zeros(NSIZE)
    padding_len = len % 8
    num_bytes = len / 8
    if padding_len != 0:
        num_bytes += 1

    for i in range(num_bytes):
        for j in range(8):
            if (i*8 + j) == len:
                break
            if e[i*8 +j]:
                out[i] = out[i] | int(1 << j)
    return out


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
    return s


def generate_sparse_rep_keccak(weight, len, shake256_state):
    raise NotImplementedError


# modulo 2 field operations

def mod_mul(b1, b2):
    raise NotImplementedError


def mod_inv(b):
    raise NotImplementedError


def mod_add(b1, b2):
    raise NotImplementedError


def split_polynomial(e):
    raise NotImplementedError


# hash functions
def function_h(m: np.ndarray) -> np.ndarray:
    seed_for_hash = Seed()

    raise NotImplementedError


def function_l(e: np.ndarray) -> np.ndarray:
    raise NotImplementedError


def function_k(m: np.ndarray, c0: np.ndarray, c1: np.ndarray) -> np.ndarray:
    raise NotImplementedError


# syndrome
def compute_syndrome(ct: CipherText, sk: SecretKey) -> Syndrome:
    raise NotImplementedError


# decoder
def BGF_decoder(syndrome: Syndrome, h0_compact: np.ndarray, h1_compact: np.ndarray):
    raise NotImplementedError


# kem mechanism
def generate_keys() -> (PublicKey, SecretKey):
    double_seeds = get_seeds(SeedPurpose.KEYGEN_SEEDS)
    state = shake256_init(double_seeds.s1.raw, ELLSIZE)
    h0 = generate_sparse_rep_keccak(DV, RBITS, state)
    h1 = generate_sparse_rep_keccak(DV, RBITS, state)

    inv_h0 = mod_inv(h0)
    h = mod_mul(h1, inv_h0)
    sigma = double_seeds.s2.raw

    return PublicKey(value=h, raw=h), SecretKey(val0=h0, val1=h1, raw=h0 + h1, sigma=sigma)


def encapsulate(pk: PublicKey) -> (CipherText, SharedSecret):
    double_seeds = get_seeds(SeedPurpose.ENCAPS_SEEDS)

    m = double_seeds.s1.raw

    e = function_h(m)
    e0, e1 = split_polynomial(e)

    c0 = mod_add(e0, mod_mul(e1, pk.value))
    c1 = m ^ function_l(e)

    k = function_k(m, c0, c1)

    return CipherText(val0=c0, val1=c1, raw=c0 + c1), SharedSecret(raw=k)


def decapsulate(ct: CipherText, sk: SecretKey) -> SharedSecret:
    h0_compact = convert2compact(sk.val0)
    h1_compact = convert2compact(sk.val1)

    syndrome = compute_syndrome(ct, sk)
    e_tmp = BGF_decoder(syndrome, h0_compact, h1_compact)

    e_prime = convertbin2byte(e_tmp, 2*RBITS)
    e_l = function_l(e_prime)

    ss = SharedSecret()
    m_prime = ct.val1 ^ e_l
    e = function_h(m_prime)
    if e_prime == e:
        k = function_k(m_prime, ct.val0, ct.val1)
        ss.raw = k
    else:
        k = function_k(sk.sigma, ct.val0, ct.val1)
        ss.raw = k
    return ss
