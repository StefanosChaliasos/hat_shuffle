from pprint import pprint
import sys
import random
import datetime

from hatshuffle.classes.crs import CRS
from hatshuffle.classes import encdec
from hatshuffle.classes.prover import Prover
from hatshuffle.classes.verifier import Verifier

n = int(sys.argv[1])


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
    return [encdec.encrypt(order, pk, message) for message in messages]


def decrypt_messages(secret, table, ciphertexts):
    return [encdec.decrypt(cs, secret, table) for cs in ciphertexts]


def demo(n, messages, batch=False):
    crs = CRS(n)
    secret = crs.sk
    pk = crs.pk
    ciphertexts = encrypt_messages(crs.order, pk, messages)

    start = datetime.datetime.now()
    sigma = random_permutation(crs.n)
    print("SIGMA = %s" % sigma)
    t_randoms = mk_t_randoms(crs.n, crs.order)
    prover = Prover(crs)
    pi_sh = prover.prove(ciphertexts, sigma, t_randoms)
    shuffled_ciphertexts = pi_sh['M_primes']
    verifier = Verifier(crs, prover)
    verify = verifier.verify_batched if batch else verifier.verify_non_batched
    perm_ok, sm_ok, cons_ok = verify(ciphertexts, pi_sh)
    print("VERIFY: %s %s %s" % (perm_ok, sm_ok, cons_ok))
    end = datetime.datetime.now()

    TABLE = encdec.make_table(pk, crs.n)
    shuffled_ms = decrypt_messages(secret, TABLE, shuffled_ciphertexts)
    #  print(shuffled_ms)
    #  print("elapsed: %s" % (end - start))


#  DEFAULT_N = 10
#  parser = argparse.ArgumentParser(description='Hat shuffler')
#  parser.add_argument('n', metavar='N', type=int, nargs='?', default=DEFAULT_N,
                    #  help='number of messages')
#  parser.add_argument('--nb', action='store_true',
                    #  default=False,
                    #  help="Don't batch verification (default is batching)")


if __name__ == '__main__':
    #  args = parser.parse_args()
    messages = list(range(n))
    demo(len(messages), messages)
    #  demo(len(messages), messages, batch=not args.nb)
