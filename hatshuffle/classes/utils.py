import elgamal
import random


system_random = random.SystemRandom()


def secure_shuffle(lst):
    random.shuffle(lst, random=system_random.random)


def random_permutation(n):
    s = list(range(n))
    secure_shuffle(s)
    return s


def mk_t_randoms(n, q):
    return [q.random() for i in range(n)]


def encrypt_messages(order, pk, messages):
    return [encrypt(order, pk, message) for message in messages]


def decrypt_messages(secret, table, ciphertexts):
    return [decrypt(cs, secret, table) for cs in ciphertexts]


def encrypt(q, pk, m):
    s = q.random()
    return elgamal.enc(pk, s, m)


def make_table(pk, n):
    return elgamal.make_table(pk[0], n)


def decrypt(ctext, secret, table):
    return elgamal.dec(secret, ctext, table)


def inverse_perm(s):
    r = [None] * len(s)
    for index, value in enumerate(s):
        r[value] = index
    return r
