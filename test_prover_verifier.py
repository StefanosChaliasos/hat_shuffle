import unittest

from hatshuffle.classes.crs import CRS
from hatshuffle.classes.prover import Prover
from hatshuffle.classes.verifier import Verifier
from hatshuffle.classes.utils import random_permutation, encrypt_messages, \
        mk_t_randoms


class ProverVerifierTest(unittest.TestCase):
    def setUp(self):
        n = 10

        self.crs = CRS(n)
        self.ciphertexts = encrypt_messages(self.crs.order, self.crs.pk,
                                            list(range(n)))
        self.sigma = random_permutation(self.crs.n)
        self.t_randoms = mk_t_randoms(self.crs.n, self.crs.order)
        self.prover = Prover(self.crs)
        self.pi_sh = self.prover.prove(self.ciphertexts, self.sigma, self.t_randoms)


    def test_verifier(self):
        verifier = Verifier(self.crs, self.prover)
        verify = verifier.verify_batched # if batch else veriefier.verify_non_batched

        self.assertEqual(verify(self.ciphertexts, self.pi_sh), (True, True, True))

        self.crs.g2alpha = self.crs.g2alpha * 2
        self.prover = Prover(self.crs)
        self.pi_sh = self.prover.prove(self.ciphertexts, self.sigma, self.t_randoms)
        verifier = Verifier(self.crs, self.prover)
        verify = verifier.verify_batched # if batch else veriefier.verify_non_batched

        self.assertEqual(verify(self.ciphertexts, self.pi_sh), (False, False, False))
