# n4_m144_true_mahler.py -- true-Mahler side of n4(-144) at high precision.
#
# c = c(tau1) = 2 sqrt(3) e^{-i pi/4} lies OUTSIDE the closed astroid
# Ast = P(T^3) (astroid functional 3.634 > 4^{2/3} = 2.520, cert2), so
# P + c has no zeros on T^3: no root ever touches |z| = 1, the Jensen
# integrand is real-analytic and periodic in (t1, t2), and the plain
# trapezoid rule converges spectrally (no cusp subtraction needed,
# unlike the c = 3 machine n4_81_final30.py).
#
#   n4(-144) = 4 m(P + c),
#   m(P+c) = (2 pi)^{-2} int int inner(t1,t2) dt1 dt2,
#   inner = sum_{roots z of z^4 + A z + B, |z|>1} log|z|   (Jensen),
#   A = c x y,  B = x^4 + y^4 + 1,  x = e^{i t1}, y = e^{i t2}.
#
# Roots: float64 numpy seeds + 3 Newton steps in mpc (50 dps).
# Result is compared against the EK4 side (10/3)(4 M12 + d3) with
# M12, d3 taken from verify_P1_n4_m144.py (60-dps Mellin values).

import numpy as np
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log

mp.dps = 50

c = 2*sqrt(3)*exp(-1j*pi/4)

# EK4-side reference, 60 digits (verify_P1_n4_m144.py [L3], [D2]):
M12 = mpf("0.3016149874129407464690529311477683998854")
d3  = mpf("0.3230659472194505140936365107238063940722")
EK4_ref = mpf(10)/3*(4*M12 + d3)

def inner(t1, t2):
    x = exp(1j*t1)
    y = exp(1j*t2)
    A = c*x*y
    B = x**4 + y**4 + 1
    seeds = np.roots([1, 0, 0, complex(A), complex(B)])
    tot = mpf(0)
    for s0 in seeds:
        z = mpc(s0.real, s0.imag)
        for _ in range(3):           # Newton, quadratic from 1e-16 seed
            z2 = z*z
            z = z - (z2*z2 + A*z + B)/(4*z2*z + A)
        a = abs(z)
        if a > 1:
            tot += log(a)
    return tot

def mahler_trap(N):
    """trapezoid on [0, 2 pi)^2, N points per axis (spectral for
    analytic periodic integrands)."""
    h = 2*pi/N
    tot = mpf(0)
    for i in range(N):
        for j in range(N):
            tot += inner(i*h, j*h)
    return tot*(h/(2*pi))**2

if __name__ == "__main__":
    import time
    for N in (64, 96, 128):
        t0 = time.time()
        m = mahler_trap(N)
        n4 = 4*m
        print("N = %3d   n4 = %s" % (N, mp.nstr(n4, 40)))
        print("          |n4 - EK4_ref| = %.2e   (%.1fs)"
              % (abs(n4 - EK4_ref), time.time()-t0))
    print()
    print("EK4 side (10/3)(4 M12 + d3) =", mp.nstr(EK4_ref, 40))
