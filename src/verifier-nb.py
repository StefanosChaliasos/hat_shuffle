import prover
from bplib.bp import GTElem

def step1(pi_sh):
    M_prime,g1_a,pi_1sp,pi_sm,pi_con=pi_sh
    g2_b,g1_c=pi_1sp
    g1_d=pi_sm
    g1_hat_a,g1_t,N=pi_con
    return M_prime,g1_a,pi_1sp,pi_sm,pi_con,g2_b,g1_c,g1_d,g1_hat_a,g1_t,N


def step2(g1_a, g2_b, g1_hat_a,g1_sum,g2_sum,g1_hat_sum):
    
    g1_an=g1_sum-sum(g1_a)
    g1_a.append(g1_an)
    
    g2_bn=g2_sum-sum(g2_b)
    g2_a.append(g2_bn)
    
    g1_hat_a_an=g1_hat_sum-sum(g1_hata_a)
    g1_hat_a.append(g1_hat_an)
    
    return g1_a, g2_b, g1_hat_a


def get_infT(gk):
    return GTElem.one(gk.G)


def step3(gk,g1_a, g2_b, g1_c,crs_1sp,g2_rho):
# CRS_1SP_T(g1_alpha_poly_zero, g1_poly_zero, g1_poly_squares,
#            g1_sum, g1_hat_sum, g2alpha, g2_sum, pair_alpha)  


    g1_alpha_poly_zero=crs_1sp.g1_alpha_poly_zero
    g2_alpha=crs_1sp.g2_alpha
    pair_alpha=crs_1sp.pair_alpha
    
    for a,b,c in zip(g1_a,g2_b,g1_c):
        left=gk.e(a+g1_alpha_poly_zero,b+g2_alpha)
        right=gk.e(c,g2_rho)*pair_alpha
        if not (left==right) return false
    return true


def step4(gk,g1_d,g1_a,g1_hat_a,crs_sm):
# CRS_SM_T(g1_beta_polys, g1_beta_rhos,g2_beta, g2_betahat)


    g2_beta=crs_sm.g2_beta
    g2_beta_hat=crs_sm.g2_betahat
    
    for d,a,ah in zip (g1_d,g1_a, g1_hat_a):
        left=gk.e(d,gk.g2)
        right=gk.e(a,g2_beta) * gk.e(ah,g2_beta_hat)
        if not (left==right) return false
    return true


def step5(gk,g1_t,crs_con,g1_a_hat,N,M,MP):
# CRS_CON_T(g1_poly_hats, g1rhohat)


    infT=get_infT(gk)
    g1_hat_rho=crs_con.g1rhohat
    g1_poly_hats=crs_con.g1_poly_hats
    
    right=[e(g1_t,pk[i]) * e(g1_hat_rho,N[i]).inv() for i in range(2)]
    left=[infT,infT]
    for i in range(2):
        for ph,mp,a,m in zip(g1_poly_hats,zip(*MP)[i],g1_a_hat,zip(*M)[i]):
            left[i]*=e(ph,mp)
            left[i]*=e(a,m).inv()
    if not (left==right) return false
    return true

     
def verify(crs,pi_sh):
# gk_T(G, q, g1, g2, gt, G.pair)    
# CRS_T(gk, pk, g1_polys, g1rho, g2_polys, g2rho, crs_sm, crs_1sp, crs_con)
    
    
    M_prime,g1_a,pi_1sp,pi_sm,pi_con,g2_b,g1_c,g1_d,g1_hat_a,g1_t,N = step1(pi_sh)
    
    g1_a, g2_b, g1_hat_a = step2(g1_a, g2_b, g1_hat_a,crs.crs_1sp.g1_sum,  \
    crs.crs_1sp.g2_sum,crs.crs_1sp.g1_hat_sum)
    
    perm_ok = step3(gk,g1_a, g2_b, g1_c,crs.crs_1sp,crs.g2rho)
    
    sm_ok   = step4(gk,g1_d,g1_a, g1_hat_a,crs.crs_sm)
    cons_ok = step5(gk,g1_t,crs.crs_con,g1_hat_a,N,M,M_prime)
    return perm_ok,sm_ok,cons_ok