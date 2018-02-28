from libffpy import G1Py, G2Py
from utils import inverse_perm
from elgamal import enc


class Prover:

    def __init__(self, crs):
        self.crs = crs


    def get_infs(self):
        inf1 = G1Py.inf()
        inf2 = G2Py.inf()
        return inf1, inf2


    def tuple_map(self, func, tpl):
        return tuple(map(func, tpl))


    def tuple_add(self, tpl1, tpl2):
        zipped = zip(tpl1, tpl2)
        return tuple(z[0] + z[1] for z in zipped)


    def step1a(self):
        randoms = [self.crs.order.random() for i in range(self.crs.n - 1)]
        g1_randoms = list(map(lambda random: random * self.crs.g1_rho, randoms))
        return randoms, g1_randoms


    def step1b(self, sigma, randoms, g1_randoms):
        A = []
        B = []
        A_hat = []
        inverted_sigma = inverse_perm(sigma)
        for perm_i, g1_random, random in zip(inverted_sigma, g1_randoms, randoms):
            A.append(self.crs.g1_polys[perm_i] + g1_random)
            B.append(self.crs.g2_polys[perm_i] + random * self.crs.g2_rho)
            A_hat.append(self.crs.g1_poly_hats[perm_i] +
                         random * self.crs.g1_rhohat)
        return A, B, A_hat


    def step2(self, A, B):
        inf1, inf2 = self.get_infs()
        A.append(self.crs.g1_sum - sum(A, inf1))
        B.append(self.crs.g2_sum - sum(B, inf2))
        return A, B


    def step3(self, A_hat):
        inf1, inf2 = self.get_infs()
        A_hat.append(self.crs.g1_hat_sum - sum(A_hat, inf1))
        return A_hat


    def step4(self, randoms, g1_randoms):
        #  rand_n = - sum(randoms) % self.crs.order
        rand_n = - sum(randoms)
        randoms.append(rand_n)
        g1_randoms.append(rand_n * self.crs.g1_rho)
        return randoms, g1_randoms


    def step5a(self, randoms, A, g1_randoms, sigma):
        C = []
        inverted_sigma = inverse_perm(sigma)
        for random, a, g1_random, perm_i in zip(
                randoms, A, g1_randoms, inverted_sigma):
            C.append(random * (2 * (a + self.crs.g1_poly_zero) - g1_random) +
                     self.crs.g1_poly_squares[perm_i])
        return C


    def step5b(self, sigma, randoms):
        D = []
        inverted_sigma = inverse_perm(sigma)
        for perm_i, random in zip(inverted_sigma, randoms):
            D.append(self.crs.g1_beta_polys[perm_i] +
                     random * self.crs.g1_beta_rhos)
        return D


    def step6(self, t_randoms):
        rt = self.crs.order.random()
        inf1, _ = self.get_infs()
        g1_t = rt * self.crs.g1_rhohat
        for t_random, g1_poly_hat in zip(t_randoms, self.crs.g1_poly_hats):
            g1_t += t_random * g1_poly_hat
        return rt, g1_t


    def step7(self, sigma, t_randoms, ciphertexts):
        _, inf2 = self.get_infs()
        M_primes = []
        for perm_i, t_random in zip(sigma, t_randoms):
            M_primes.append(self.tuple_add(ciphertexts[perm_i],
                                      enc(self.crs.pk, t_random, 0)))
        return M_primes


    def step8(self, rt, randoms, ciphertexts):
        N = self.tuple_map(lambda elem: rt * elem, self.crs.pk)
        for random, ciphertext in zip(randoms, ciphertexts):
            N = self.tuple_add(N, self.tuple_map(lambda elem: random * elem, ciphertext))
        return N


    def step9(self, B, C):
        return {'B': B[:-1], 'C': C}


    def step10(self, D):
        return {'D': D}


    def step11(self, A_hat, g1_t, N):
        return {'A_hat': A_hat[:-1], 'g1_t': g1_t, 'N': N}


    def step12(self, M_primes, A, pi_1sp, pi_sm, pi_con):
        return {'M_primes': M_primes, 'A': A[:-1], 'pi_1sp': pi_1sp,
                'pi_sm': pi_sm, 'pi_con': pi_con}


    def prove(self, ciphertexts, sigma, t_randoms):
        randoms, g1_randoms = self.step1a()
        A, B, A_hat = self.step1b(sigma, randoms, g1_randoms)
        A, B = self.step2(A, B)
        A_hat = self.step3(A_hat)
        randoms, g1_randoms = self.step4(randoms, g1_randoms)
        C = self.step5a(randoms, A, g1_randoms, sigma)
        D = self.step5b(sigma, randoms)
        rt, g1_t = self.step6(t_randoms)
        M_primes = self.step7(sigma, t_randoms, ciphertexts)
        N = self.step8(rt, randoms, ciphertexts)
        pi_1sp = self.step9(B, C)
        pi_sm = self.step10(D)
        pi_con = self.step11(A_hat, g1_t, N)
        pi_sh = self.step12(M_primes, A, pi_1sp, pi_sm, pi_con)
        return pi_sh
