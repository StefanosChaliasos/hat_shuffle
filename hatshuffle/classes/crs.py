from libffpy import LibffPy, BigNum


class CRS:
    def __init__(self, n):
        self.n = n

        #  mk_gk
        self.g = LibffPy(n)  # Generator
        self.order = self.g.order()
        self.g1 = self.g.gen1()
        self.g2 = self.g.gen2()
        self.gt = self.g.pair(self.g1, self.g2)
        self.pair = self.g.pair

        # mk_Chi
        self.chi = self.order.random()
        self.alpha = self.order.random()
        self.beta = self.order.random()
        self.betahat = self.order.random()
        self.rho = 1 + (self.order - 1).random()
        self.rhohat = 1 + (self.order - 1).random()
        self.sk = self.order.random()

        # mk_crs
        self.polys_all = self.generate_pis()
        self.poly_zero = self.polys_all[0]
        self.polys = self.polys_all[1:]
        self.poly_hats = self.generate_pi_hats()

        self.pk = (self.g2, self.chi * self.g2)
        self.g1_polys = list(map(lambda poly: poly * self.g1, self.polys))
        self.g1_rho = self.rho * self.g1
        self.g2_polys = list(map(lambda poly: poly * self.g2, self.polys))
        self.g2_rho = self.rho * self.g2

        # crs_sm
        self.g1_beta_polys = []
        for poly, poly_hat in zip(self.polys, self.poly_hats):
            self.g1_beta_polys.append(
                (self.beta * poly + self.betahat * poly_hat) * self.g1)
        self.g1_beta_rhos = (self.beta * self.rho + self.betahat * self.rhohat)\
                 * self.g1
        self.g2_beta = self.beta * self.g2
        self.g2_betahat = self.betahat * self.g2

        # crs_1sp
        self.poly_sum = sum(self.polys)
        self.poly_hat_sum = sum(self.poly_hats)
        #  G1
        self.g1_alpha_poly_zero = (self.alpha + self.poly_zero) * self.g1
        self.g1_poly_zero = self.poly_zero * self.g1
        self.inv_rho = self.rho.mod_inverse()
        self.g1_poly_squares = []
        for poly in self.polys:
            numer = (poly + self.poly_zero) ** 2 - 1
            self.g1_poly_squares.append((numer * self.inv_rho) * self.g1)
        self.g1_sum = self.poly_sum * self.g1
        self.g1_hat_sum = self.poly_hat_sum * self.g1
        #  G2
        self.g2alpha = (-self.alpha + self.poly_zero) * self.g2
        self.g2_sum = self.poly_sum * self.g2
        #  GT
        #  self.pair_alpha = self.gt ** ((1 - self.alpha ** 2) % self.order)
        self.pair_alpha = self.gt ** (1 - self.alpha ** 2)

        #  crs_con
        self.g1_poly_hats = list(map(lambda poly_hat: poly_hat * self.g1,
                                     self.poly_hats))
        self.g1_rhohat = self.rhohat * self.g1

        self.trapdoor = (self.chi, self.rhohat)


    def generate_pis(self):
        """Computes vector of elements P_i(chi) for i = 0, ..., n.

        P_0 is special, defined as ln_plus1(chi) - 1, where ln_plus1 is the
        (n+1)th Lagrange basis polynomial.
        The rest P_i are defined as 2*ln_i(ch) + ln_plus1(chi), where ln_i is the
        ith Lagrange basis polynomial.
        It returns a list with the values of the n+1 polynomials P_i(chi).
        Uses Lagrange basis polynomials with distinct points 1, ..., n + 1.
        """
        if self.chi <= self.n + 1:
            raise ValueError("chi should be greater than n + 1, chi=%s n+1=%s"
                             % (self.chi, self.n+1))
        pis = []

        prod = BigNum(1)
        # prod = (x - w_1) (x - w_2) ... (x - w_{n+1})
        for j in range(1, self.n + 2):
            prod *= (self.chi - j)

        # denoms[0] = 1 / (w_1 - w_2) (w_1 - w_3) ... (w_1 - w_{n + 1})
        # denoms[1] = 1 / (w_2 - w_1) (w_2 - w_3) ... (w_2 - w_{n + 1})
        # denoms[n] = 1 / (w_{n+1}- w_1) (w_{n+1} - w_2) ... (w_{n+1} - w_n)
        denoms = self.compute_denominators()

        missing_factor = self.chi - (self.n + 1)

        ln_plus1 = prod * missing_factor.mod_inverse()
        ln_plus1 *= denoms[self.n].mod_inverse()

        # P_0 is special
        pis.append(ln_plus1 - BigNum(1))

        two = BigNum(2)
        for i in range(1, self.n + 1):
            missing_factor = self.chi - i
            l_i = prod * missing_factor.mod_inverse()
            l_i *= denoms[i - 1].mod_inverse()
            pis.append(two * l_i + ln_plus1)
        return pis


    def generate_pi_hats(self):
        """Computes vector of elements \hat{P}_i(chi) for i = 1, ..., n.

        These are simply the polynomials x^{(i+1)(n+1)}.
        It returns a list with the values of the n polynomials \hat{P}_i(chi).
        """
        hpis = []
        #  hpi = self.chi**(self.n+1) % self.order
        #  hpis = [ hpi**(i+1) % self.order for i in range(1, self.n+1)]
        hpi = self.chi**(self.n+1)
        hpis = [ hpi**(i+1) for i in range(1, self.n+1)]

        return hpis


    def compute_denominators(self):
        """Computes denominators for Lagrange basis polynomials.

        Uses distinct points 1, ...,k
        k -- number of basis polynomials
        """
        k = self.n + 1
        denominators = []
        temp = BigNum(1)
        for i in range(1, k+1):
            if i == 1:
                for j in range(2, k+1):
                    elem = i - j;
                    temp *= elem
            elif i == k:
                elem = 1 - k;
                temp *= elem
            else:
                inverse = BigNum(i - 1 - k)
                inverse = inverse.mod_inverse()
                elem = i - 1
                temp *= elem
                temp *= inverse
            denominators.append(temp)
        return denominators
