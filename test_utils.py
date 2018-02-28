import unittest

from libffpy import LibffPy, BigNum
from hatshuffle.classes.utils import make_table, encrypt, decrypt


class UtilTest(unittest.TestCase):
    def setUp(self):
        self.MAX = 1000
        self.g = LibffPy(self.MAX)


    def key_gen(self, gen2, order):
        priv = order.random()
        pub = (gen2, priv * gen2)
        return (pub, priv)


    def test_make_table(self):
        order = self.g.order()
        gen2 = self.g.gen2()
        pk, sk = self.key_gen(gen2, order)
        table = make_table(pk, self.MAX)

        self.assertEqual(table[500 * pk[0]], 500)
        self.assertNotEqual(table[2 * pk[0]], 500)


    def test_encdec(self):
        order = self.g.order()
        gen2 = self.g.gen2()
        pk, sk = self.key_gen(gen2, order)

        c = encrypt(order, pk, 500)

        table = make_table(pk, self.MAX)

        self.assertEqual(decrypt(c, sk, table), 500)

        import random
        ps = random.sample(range(self.MAX), 100)
        for i in range(100):
            c = encrypt(order, pk, ps[i])
            self.assertEqual(decrypt(c, sk, table), ps[i])
