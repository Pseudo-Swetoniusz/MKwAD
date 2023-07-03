"""
Following functions has been copied and adapted from:
https://gist.github.com/mildsunrise/e21ae2b1649532813f2594932f9e9371
"""

# polynomials over modulo 2 field operations
import numpy as np
import math
from project.parameters import *


def mod_mul(b1, b2):
    b1_int = int.from_bytes(bytearray(b1), 'big')
    b2_int = int.from_bytes(bytearray(b2), 'big')
    modulus = np.zeros(RSIZE, dtype=np.uint8)
    modulus[0] = 128
    modulus[RSIZE - 1] = 1
    mod_int = int.from_bytes(bytearray(modulus), 'big')
    result = 0
    a = b1_int
    b = b2_int
    deg = RBITS
    result = p_mod(poly_mul(b1_int, b2_int), mod_int)
    res_byte = bytearray(result.to_bytes(RSIZE, 'big'))
    res = [x for x in res_byte]
    return res


def poly_mul(a: int, b: int) -> int:
    result = 0
    while a and b:
        if a & 1:
            result ^= b
        a >>= 1
        b <<= 1
    return result


def p_mod(a: int, b: int) -> int:
    bl = b.bit_length()
    while True:
        shift = a.bit_length() - bl
        if shift < 0:
            return a
        a ^= b << shift


def p_divmod(a: int, b: int) -> tuple[int, int]:
    q = 0
    bl = b.bit_length()
    while True:
        shift = a.bit_length() - bl
        if shift < 0:
            return q, a
        q ^= 1 << shift
        a ^= b << shift


def p_egcd(b1: int, b2: int) -> tuple[int, int, int]:
    a = (b1, 1, 0)
    b = (b2, 0, 1)
    while True:
        q, r = p_divmod(a[0], b[0])
        if not r:
            return b
        a, b = b, (r, a[1] ^ poly_mul(q, b[1]), a[2] ^ poly_mul(q, b[2]))


def mod_inv(b):
    b1_int = int.from_bytes(bytearray(b), 'big')
    modulus = np.zeros(RSIZE, dtype=np.uint8)
    modulus[0] = 128
    modulus[RSIZE - 1] = 1
    mod_int = int.from_bytes(bytearray(modulus), 'big')
    d, x, y = p_egcd(b1_int, mod_int)
    result = bytearray(x.to_bytes(RSIZE, 'big'))
    return d, [i for i in result]


def poly_add(b1, b2):
    l1 = len(b1)
    l2 = len(b2)
    x = []

    if l1 == l2:
        for i in range(l1):
            x.append(b1[i] ^ b2[i])
        return np.array(x)

    elif l1 < l2:
        for i in range(l1):
            x.append(b1[i] ^ b2[i])
        for j in range(l1, l2):
            x.append(b2[j])
        return np.array(x)

    else:
        for i in range(l2):
            x.append(b1[i] ^ b2[i])
        for j in range(l2, l1):
            x.append(b1[j])
        return np.array(x)


def split_polynomial(e):
    idx = len(e) - RSIZE
    e0 = e[idx:]
    e1 = e[:idx]
    return e1, e0
