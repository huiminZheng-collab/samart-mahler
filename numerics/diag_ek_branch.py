# Diagnostic: why does EK(tau') != n2(s)?
# Answer (numerical): Samart's Prop 2.1(i) is proven for Im tau >= 1/2 only.
# tau' = (1+s7)/8 has Im = sqrt7/8 < 1/2, so EK(tau') is an analytic-
# continuation value on a different branch of mtilde.  The second preimage
# of s under s2 on X_0(4) (w4 partner tau_w = -1/(4 tau') = (-1+sqrt(-7))/4,
# Im = sqrt7/4 > 1/2) satisfies s2(tau_w) = s and EK(tau_w) = n2(s) exactly.
# Moreover the branch gap has a closed form in {M7, d7}:
#   EK(tau') = (8/7)(44 M7 - d7)   (pslq below, 40 digits).

from mpmath import mp, mpf, mpc, pi, exp, sqrt, sin, cos, nstr, pslq, log, catalan

mp.dps = 50

def eta(tau):
    q = exp(2 * pi * 1j * tau)
    p = mpc(1); qn = q
    for n in range(1, 500):
        p *= (1 - qn); qn *= q
        if abs(qn) < mpf(10) ** (-60):
            break
    return exp(pi * 1j * tau / 12) * p

def s2(tau):
    return -(eta(tau + mpf(1) / 2)) ** 24 / (eta(2 * tau + 1)) ** 24

def U(j, tau):
    tot = mpc(0)
    for m in range(1, 3000):
        z = j * pi * m * tau
        term = cos(z) / sin(z) ** 3 / m
        tot += term
        if abs(term) < mpf(10) ** (-60):
            break
    return 2 * pi ** 3 * tot

def EK(tau):
    return (2 * pi * tau + (2 / pi ** 3) * (U(1, tau) - 4 * U(4, tau))).imag

s7 = sqrt(mpc(-7))
s_t = (47 + 45 * s7) / 2
tau_p = (1 + s7) / 8
tau_w = -1 / (4 * tau_p)      # = (-1+sqrt(-7))/4
print("tau_w =", nstr(tau_w, 20), " Im =", nstr(tau_w.imag, 10))
print("|s2(tau_w) - s| =", nstr(abs(s2(tau_w) - s_t), 5))
EKw, EKp = EK(tau_w), EK(tau_p)
print("EK(tau_w) =", nstr(EKw, 45))
print("EK(tau')  =", nstr(EKp, 45))

M7 = mpf("0.10267160777890201121045659489829291399889482708922")
d7 = mpf(7) * sqrt(mpf(7)) / (4 * pi) * mpf("1.1519254705444910471016923973205499647978214")
rhs = mpf(4) / 7 * (54 * M7 + d7)
print("rhs       =", nstr(rhs, 45))
print("EK(tau_w) - rhs =", nstr(EKw - rhs, 5))

# pslq on the branch gap against a pool of constants
gap = rhs - EKp   # = n2(s) - EK(tau')
cands = [mpf(1), M7, d7, pi, log(2), catalan,
         mpf("1.1519254705444910471016923973205499647978214"), pi ** 2]
rel = pslq([gap] + cands, tol=mpf(10) ** (-40), maxcoeff=10 ** 6, maxsteps=1000)
print("pslq(gap, 1, M7, d7, pi, log2, Catalan, L(chi,2), pi^2) =", rel)
print("gap - (12 d7 - 136 M7)/7 =", nstr(gap - (12 * d7 - 136 * M7) / 7, 5))
print("EK(tau') - (8/7)(44 M7 - d7) =", nstr(EKp - mpf(8) / 7 * (44 * M7 - d7), 5))
