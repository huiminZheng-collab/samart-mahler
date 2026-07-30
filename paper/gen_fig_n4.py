# gen_fig_n4.py
#
# Data generator for paper/formal2.tex, Figure 2 (the wrong-sheet
# departure of Samart's EK4 formula below Im tau = 1/sqrt(2)).
#
# Writes n4_departure.dat with three whitespace-separated columns
#     y   ek4 = EK4(i y)   fourm = 4*m(P + c(i y))
# for y = 0.30, 0.32, ..., 1.50 (61 ordinates), where
#   P = (x^4+y^4+z^4+1)/(xyz),
#   c(tau) = (eta(2t)/eta(t))^6 (16 A^4 + A^-4),
#   A = eta(t) eta(4t)^2 / eta(2t)^3   (natural holomorphic branch of
#   s4^{1/4}; c(tau)^4 = s4(tau), checked at every sample),
#   EK4(tau) = Im[ 2 pi tau + (10/pi^3) (U1 - 2 U2) ],
#   U_j = 2 pi^3 sum_{m>=1} cos(j pi m tau) / (m sin^3(j pi m tau))
#   (same implementation as diag_n4_astroid.py / verify_P1_n4_81.py).
#
# On the imaginary axis c(i y) is real and >= 4 (s4(i y) >= 256 always,
# Samart's Lemma 2.2), so the z-Jensen integrand is smooth and the direct
# integral is reliable everywhere on the grid.
#
# Precision plan: eta / c(tau) / s4 / EK4 at mpmath 30 dps.  The direct
# integral uses the z-Jensen reduction (for fixed x,y the z-polynomial is
# z^4 + c x y z + (x^4+y^4+1); m = torus mean of sum_{|root|>1} log|root|)
# evaluated by the midpoint trapezoid rule with batched numpy float64
# roots (companion eigvals), which converges spectrally for the smooth
# periodic integrand.  Validated against the mpmath nested tanh-sinh
# quadrature mahler4() of diag_n4_astroid.py: agreement 1.5e-11 at
# c = 4.9421... and 2.7e-7 at the near-critical c = 4.0283... ; the grid
# is refined (N up to 1024) near the critical value c = 4 (y ~ 1/sqrt(2)),
# where the error stays <~1e-4, invisible at plot scale.
#
# Anchors (run first): EK4(i) = 4 m(P + 648^{1/4}) = 6.4332830658...
# (Samart's proved region, s4(i) = 648); and EK4(0.3 i) - 4m = +11.369...
# reproducing the diag_n4_astroid.py discrepancy at tau = 0.3 i.
#
# New file; does not modify anything else.

import time
import numpy as np
from mpmath import mp, mpf, mpc, pi, exp, cos, sin

mp.dps = 30

# ============================ eta, s4, c(tau) ============================
def eta(tau, nterms=400):
    q = exp(2 * pi * 1j * tau)
    p = mpc(1)
    qn = q
    for n in range(1, nterms + 1):
        p *= (1 - qn)
        qn *= q
        if abs(qn) < mpf(10) ** (-40):
            break
    return exp(pi * 1j * tau / 12) * p

def s4(tau):
    """Delta-quotient definition (verify_P1_n4_81.py [A0])."""
    e1, e2, e4 = eta(tau), eta(2 * tau), eta(4 * tau)
    r = e1 * e4 ** 2 / e2 ** 3
    return (e2 / e1) ** 24 * (16 * r ** 4 + r ** -4) ** 4

def c_of_tau(tau):
    """Natural holomorphic branch of s4^{1/4}."""
    e1, e2, e4 = eta(tau), eta(2 * tau), eta(4 * tau)
    A = e1 * e4 ** 2 / e2 ** 3
    return (e2 / e1) ** 6 * (16 * A ** 4 + A ** -4)

# ============================ EK4 U-series ============================
def U(j, tau):
    tot = mpc(0)
    for m in range(1, 3000):
        z = j * pi * m * tau
        term = cos(z) / sin(z) ** 3 / m
        tot += term
        if abs(term) < mpf(10) ** (-40):
            break
    return 2 * pi ** 3 * tot

def EK4(tau):
    return (2 * pi * tau + (10 / pi ** 3) * (U(1, tau) - 2 * U(2, tau))).imag

# ============================ direct integral (numpy, float64) ============================
def four_mahler_trap(c, N):
    """4*m(x^4+y^4+z^4+1+c x y z), real c > 4, midpoint rule, N x N."""
    t = np.pi * (np.arange(N) + 0.5) / N
    T1, T2 = np.meshgrid(t, t, indexing="ij")
    X = np.exp(1j * T1)
    Y = np.exp(1j * T2)
    A = c * X * Y
    B = X ** 4 + Y ** 4 + 1
    # companion matrices of z^4 + A z + B (batched)
    M = np.zeros(T1.shape + (4, 4), dtype=complex)
    M[..., 0, 1] = 1
    M[..., 1, 2] = 1
    M[..., 2, 3] = 1
    M[..., 3, 0] = -B
    M[..., 3, 1] = -A
    R = np.linalg.eigvals(M)
    return 4.0 * np.log(np.maximum(np.abs(R), 1.0)).sum(-1).mean()

def four_mahler(c):
    """Grid refinement near the critical value c = 4."""
    dc = float(c) - 4.0
    N = 256 if dc > 0.05 else (512 if dc > 0.005 else 1024)
    return four_mahler_trap(float(c), N)

# ============================ main ============================
def main():
    t_start = time.time()
    # ---- control 1: tau = i (proved region, s4(i) = 648) ----
    ek_i = EK4(mpc(0, 1))
    c_i = c_of_tau(mpc(0, 1)).real
    m4_i = four_mahler(c_i)
    print("control tau=i: c(i) = %.12f (648^(1/4) = %.12f)" %
          (mpf(c_i), mpf(648) ** mpf("0.25")))
    print("control tau=i: EK4 = %.15f  4m = %.15f  |diff| = %.2e" %
          (ek_i, m4_i, abs(mpf(ek_i) - m4_i)))
    assert abs(mpf(ek_i) - m4_i) < mpf(10) ** (-8), "control at tau=i failed"
    # ---- control 2: tau = 0.3 i (wrong sheet; diag_n4_astroid.py) ----
    ek3 = EK4(mpc(0, mpf("0.3")))
    c3 = c_of_tau(mpc(0, mpf("0.3"))).real
    m43 = four_mahler(c3)
    print("control tau=0.3i: c = %.10f  EK4 - 4m = %.6f (diag: +11.4)" %
          (mpf(c3), ek3 - m43))

    # ---- grid ----
    Y0, Y1, NY = mpf("0.30"), mpf("1.50"), 61
    rows = []
    print("%-6s %-16s %-18s %-18s %-10s" % ("y", "c(iy)", "EK4(iy)", "4m(P+c)", "EK4-4m"))
    for k in range(NY):
        y = Y0 + (Y1 - Y0) * k / (NY - 1)
        tau = mpc(0, y)
        c = c_of_tau(tau)
        assert abs(c.imag) < mpf(10) ** (-20), "c(iy) should be real"
        assert abs(c ** 4 - s4(tau)) < mpf(10) ** (-20), "c^4 != s4"
        ek = EK4(tau)
        mv = four_mahler(c.real)
        rows.append((y, mpf(ek), mv))
        if k % 10 == 0 or k == NY - 1:
            print("%-6.3f %-16.10f %-18.12f %-18.12f %-10.4f"
                  % (y, c.real, ek, mv, ek - mv))

    with open("n4_departure.dat", "w") as f:
        f.write("# n4_departure.dat -- Figure 2 data for formal2.tex\n")
        f.write("# generated by gen_fig_n4.py: 61 ordinates y in [0.30,1.50];\n")
        f.write("# eta/c/EK4 at mpmath 30 dps; direct integral by z-Jensen +\n")
        f.write("# midpoint trapezoid (numpy float64 batched roots, validated\n")
        f.write("# vs mpmath quad to 1.5e-11 at c=4.94, 2.7e-7 at c=4.028).\n")
        f.write("y ek4 fourm\n")
        for y, ek, mv in rows:
            f.write("%s %s %s\n" % (mp.nstr(y, 6), mp.nstr(ek, 25),
                                    mp.nstr(mv, 18)))

    print("wrote n4_departure.dat (%d rows), total %.1f s"
          % (len(rows), time.time() - t_start))

if __name__ == "__main__":
    main()
