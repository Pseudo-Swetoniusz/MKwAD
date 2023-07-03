import math

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
