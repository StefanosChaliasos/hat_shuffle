from libffpy import GTPy
from prover import Prover
from functools import reduce


class Verifier:

    def __init__(self, crs):
        self.crs


    def get_infT(self):
        return GTPy.one()


    def mexp(self, scals, elems, inf):  # find optimized version in lib
        return reduce(lambda x, y: x + y,
                      map(lambda tup: tup[0] * tup[1], zip(scals, elems)),
                      inf)


    def step1(self, pi_sh):
        M_prime, g1_a, pi_1sp, pi_sm, pi_con = pi_sh
        g2_b, g1_c = pi_1sp
        g1_d = pi_sm
        g1_hat_a, g1_t, N = pi_con
        return M_prime, g1_a, pi_1sp, pi_sm, pi_con, \
            g2_b, g1_c, g1_d, g1_hat_a, g1_t, N


    def step2(self, g1_a, g2_b, g1_hat_a):
        g1_a, g2_b = prover.step2(self.crs.gk, g1_a, g2_b, self.crs.g1_sum, \
                                  self.crs.g2_sum)
        g1_hat_a = prover.step3(self.crs.gk, g1_hat_a, self.crs.g1_hat_sum)
        return g1_a, g2_b, g1_hat_a


    def step3(self, g1_a, g2_b, g1_c):
        for a, b, c in zip(g1_a, g2_b, g1_c):
            left = self.crs.g.pair(a + self.crs.g1_alpha_poly_zero, b +
                                 self.crs.g2alpha)
            right = self.crs.gk.e(c, self.crs.g2_rho) * self.crs.pair_alpha
            if left != right:
                return False
        return True


    def step3_batched(self, g1_a, g2_b, g1_c):
        y1 = [self.crs.order.random() for x in range(self.crs.n - 1)] + [1]
        infT = get_infT()
        inf1, _ = Prover.get_infs()
        right = (self.crs.g.pair(mexp(y1, g1_c, inf1), self.crs.g2_rho) *
                 (self.crs.pair_alpha ** (sum(y1) % self.crs.order)))
        left = infT
        for yi, a, b in zip(y1, g1_a, g2_b):
            left *= self.crs.g.pair(yi * (a + self.crs.g1_alpha_poly_zero), b +
                                  self.crs.g2alpha)
        return left == right


    def step4(self, g1_d, g1_a, g1_hat_a):
        g2_beta = crs_sm.g2_beta
        g2_beta_hat = crs_sm.g2_betahat

        for d, a, ah in zip(g1_d, g1_a, g1_hat_a):
            left = self.crs.g.pair(d, self.crs.g2)
            right = self.crs.g.pair(a, self.crs.g2_beta) * \
                    self.crs.g.pair(ah, self.crs.g2_beta_hat)
            if left != right:
                return False
        return True


    def step4_batched(self, g1_d, g1_a, g1_hat_a):
        y2 = [self.crs.order.random() for x in range(self.crs.n - 1)] + [1]
        inf1, _ = Prover.get_infs()
        left = self.crs.g.pair(mexp(y2, g1_d, inf1), self.crs.g2)
        right = (self.crs.g.pair(mexp(y2, g1_a, inf1), self.crs.g2_beta) *
                 self.crs.g.pair(mexp(y2, g1_hat_a, inf1), self.crs.g2_beta_hat))
        return left == right


    def step5(self, g1_t, g1_a_hat, N, M, MP):
        infT = get_infT()
        right = [self.crs.g.pair(g1_t, self.crs.pk[i]) *
                 self.crs.g.pair(self.crs.g1_hat_rho, N[i]).inv()
                 for i in range(2)]
        left = [infT, infT]
        for i in range(2):
            for ph, mp, a, m in zip(self.crs.g1_poly_hats, MP, g1_a_hat, M):
                left[i] *= self.crs.g.pair(ph, mp[i]) * \
                           self.crs.g.pair(a, m[i]).inv()
        return left == right


    def step5_batched(self, g1_t, g1_a_hat, N, M, MP):
        infT = get_infT()
        _, inf2 = Prover.get_infs()
        y3 = [self.crs.order.random(), 1]

        q = self.crs.g.pair(g1_t, mexp(y3, self.crs.pk, inf2))
        right = q * self.crs.g.pair(self.crs.g1_hat_rho, mexp(y3, N, inf2)).inv()

        left = infT
        for ph, mp, a, m in zip(self.crs.g1_poly_hats, MP, g1_a_hat, M):
            left *= self.crs.g.pair(ph, mexp(y3, mp, inf2)) * \
                    self.crs.g.pair(a, mexp(y3, m, inf2)).inv()
        return left == right


    def verify_batched(self, M, pi_sh):
        print("Using Batching")
        M_prime, g1_a, pi_1sp, pi_sm, \
            pi_con, g2_b, g1_c, g1_d, g1_hat_a, g1_t, N = step1(pi_sh)

        g1_a, g2_b, g1_hat_a = step2(g1_a, g2_b, g1_hat_a)

        perm_ok = step3_batched(g1_a, g2_b, g1_c)
        sm_ok = step4_batched(g1_d, g1_a, g1_hat_a)
        cons_ok = step5_batched(g1_t, g1_hat_a, N, M, M_prime)

        return perm_ok, sm_ok, cons_ok


    def verify_non_batched(self, M, pi_sh):
        print("Not Using Batching")

        M_prime, g1_a, pi_1sp, pi_sm, \
            pi_con, g2_b, g1_c, g1_d, g1_hat_a, g1_t, N = step1(pi_sh)

        g1_a, g2_b, g1_hat_a = step2(g1_a, g2_b, g1_hat_a)

        perm_ok = step3(g1_a, g2_b, g1_c)
        sm_ok = step4(g1_d, g1_a, g1_hat_a)
        cons_ok = step5(g1_t, g1_hat_a, N, M, M_prime)

        return perm_ok, sm_ok, cons_ok
