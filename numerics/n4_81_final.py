# n4_81_final.py -- FINAL high-precision direct computation of
#     n4(81) = 4 m(C_81),  C_81 = x^4+y^4+z^4+1+3xyz   (c = 81^{1/4} = 3).
#
# Structure (see NOTES.md for the full derivation):
#   * Jensen in z: m = E_{t1,t2}[ g ],  g = sum_j log^+ |z_j|,
#     z^4 + 3xy z + (x^4+y^4+1) = 0,  x = e^{i t1}, y = e^{i t2}.
#   * g is even in t1, t2, pi/2-PERIODIC in t1 (z -> i z maps roots of the
#     polynomial at x -> i x onto those at x, preserving |z|), and symmetric
#     about t1 = pi/4.  Hence  m = (4/pi^2) int_0^{pi/4} inner(t1) dt1.
#   * inner(t1) has |x|-type cusps exactly where the kink count in t2
#     changes.  Kinks (some |z_j| = 1) are found EXACTLY via the resultant
#     S(W; X) of z^4+Az+B and z^4 conj(P)(1/z), W = y^4, X = x^4
#     (kink_resultant_W.py, degree 8 in W): no grid misses.
#   * Critical t1 in (0, pi/4) by exact-count bisection (float64; the
#     split-point error ~1e-10 contributes ~slope*delta^2 ~ 1e-20):
#         t1a = 0.36136712(3)  (count 6 -> 2; bisection band width 2.3e-10)
#         t1b = pi/6 exactly   (count 2 -> 0)
#   * inner: tanh-sinh maxdegree=7 on each kink-free t2 piece
#     (deg7 vs deg8 differ by 1.3e-25 at the kink-free point t1 = 1.0,
#     bench_inner.py; pieces are analytic up to their endpoints).
#   * g eval: numpy seeds + 2 Newton steps at 42 dps, err 2e-41 (measured).
#
# Expected error budget: inner deg7 ~1e-25/eval; kink positions ~1e-12
# (float64 resultant roots) -> ~1e-24 residual per piece; critical split
# points ~1e-10 -> ~1e-20; outer tanh-sinh truncation measured by the
# deg5..deg8 sequence below (endpoint cusps only -> fast convergence).

import numpy as np
from mpmath import mp, mpf, mpc, pi, exp, log, quad, polyroots
import time

mp.dps = 42
c = mpf(3)

from kink_resultant_W import COEFFS

# ---------------- g evaluation ----------------

def make_g(t1):
    x = exp(1j * t1)
    x4 = x ** 4
    cx = c * x
    def g(t2):
        y = exp(1j * t2)
        A = cx * y
        B = x4 + y ** 4 + 1
        seeds = np.roots([1, 0, 0, complex(A), complex(B)])
        tot = mpf(0)
        for s0 in seeds:
            z = mpc(s0.real, s0.imag)
            for _ in range(2):
                z2 = z * z
                z = z - (z2 * z2 + A * z + B) / (4 * z2 * z + A)
            a = abs(z)
            if a > 1:
                tot += log(a)
        return tot
    return g

# ---------------- exact kink detection via resultant in W = y^4 ----------------

def S_coeffs_np(t1f):
    X = np.exp(4j * t1f)
    out = []
    for entry in reversed(COEFFS):   # descending order for np.roots
        if entry is None:
            out.append(0j)
            continue
        lo, cl = entry
        v = 0j
        for ci in reversed(cl):
            v = v * X + ci
        out.append(v * X ** lo)
    return np.array(out, dtype=complex)

def kinks_t2(t1, g):
    """All t2 in (0, pi) where some |z_j(t1, t2)| = 1.  Exact detection."""
    t1f = float(t1)
    r = np.roots(S_coeffs_np(t1f))
    out = []
    for wj in r:
        if abs(abs(wj) - 1) < 1e-6:
            th = np.angle(wj) % (2 * np.pi)
            for k in range(2):       # k=0,1 give t2 in (0, pi)
                t2 = (th + 2 * np.pi * k) / 4
                if 1e-10 < t2 < np.pi - 1e-10:
                    out.append(t2)
    out = sorted(set(out))
    # verify at mp precision: a z-root must sit on the unit circle
    ver = []
    for t2 in out:
        y = exp(1j * mpf(t2))
        x = exp(1j * t1)
        ok = False
        for z in polyroots([1, 0, 0, c * x * y, x ** 4 + y ** 4 + 1]):
            if abs(abs(z) - 1) < mpf('1e-6'):
                ok = True
                break
        if ok:
            ver.append(mpf(t2))
    return ver

# ---------------- inner / outer quadrature ----------------

T1A = mpf('0.361367123929')   # count 6 -> 2 (bisection band +-1.2e-10)
T1B = pi / 6                  # count 2 -> 0 (exact special value)

def inner(t1, deg=7):
    g = make_g(t1)
    ks = kinks_t2(t1, g)
    pts = [mpf(0)] + ks + [pi]
    return quad(g, pts, maxdegree=deg)

def mahler(deg_out, deg_in=7):
    f = lambda t1: inner(t1, deg_in)
    v = quad(f, [mpf(0), T1A, T1B, pi / 4], maxdegree=deg_out)
    return 4 * v / pi ** 2

if __name__ == "__main__":
    t_start = time.time()
    # symmetry sanity: inner(t) = inner(pi/2 - t) at low degree
    for t in ('0.7', '1.1'):
        a = inner(mpf(t), 5)
        b = inner(pi / 2 - mpf(t), 5)
        print("symmetry t=%s: |inner(t)-inner(pi/2-t)| = %.3e"
              % (t, float(abs(a - b))), flush=True)
    prev = None
    for deg_out in (5, 6, 7, 8):
        t0 = time.time()
        m = mahler(deg_out)
        print("deg_out=%d  m  = %s  [%.0fs]"
              % (deg_out, m, time.time() - t0), flush=True)
        print("deg_out=%d n4(81) = 4m = %s" % (deg_out, 4 * m), flush=True)
        if prev is not None:
            print("   change from deg %d: %.3e"
                  % (deg_out - 1, float(abs(m - prev))), flush=True)
        prev = m
    M7 = mpf("0.10267160777890201121045659489829291399889482708922")
    print("40 M7 =", 40 * M7, flush=True)
    print("4m(deg8) - 40 M7 =", 4 * prev - 40 * M7, flush=True)
    print("total %.0fs" % (time.time() - t_start), flush=True)
