# diag_n4_astroid.py
#
# Go/no-go diagnostic: does  4 m(P + c(tau)) = EK4(tau)  hold for the n4
# family  P = (x^4+y^4+z^4+1)/(xyz)  inside the astroid disc
#     D = { |Re c|^{2/3} + |Im c|^{2/3} <= 4^{2/3} }   (= P(T^3)),
# away from the cross  CROSS = [-4,4] u i[-4,4] ?
#
#   c(tau) = (eta(2t)/eta(t))^6 (16 A^4 + A^-4),  A = eta(t) eta(4t)^2/eta(2t)^3
# is the natural holomorphic branch of s4(tau)^{1/4}  (c(tau)^4 = s4(tau),
# checked numerically against the Delta-quotient definition to 40 digits).
#
# EK4(tau): agent-7's U-series implementation (verify_P1_n4_81.py),
#     EK4(tau) = Im[ 2 pi tau + (10/pi^3)(U1 - 2 U2) ],
#     U_j = 2 pi^3 sum_{m>=1} cos(j pi m tau)/(m sin^3(j pi m tau)),
# convergent for every Im tau > 0.
#
# Direct side: m(P+c) = m(x^4+y^4+z^4+1+c x y z)  (m(xyz)=0);  for fixed
# x,y on the torus the z-polynomial is  z^4 + (c x y) z + (x^4+y^4+1),
# Jensen via numerical polyroots (generalizes mahler_m.py; same scheme as
# mahler4() in verify_P1_n4_81.py, extended to complex c).
#
# New file; does not modify anything else.  mpmath 50 dps throughout,
# except the 20-point path scan where the direct integral runs at 35 dps
# / tanh-sinh degree 4 for speed (noted in the output).

import time
from mpmath import (mp, mpf, mpc, pi, exp, cos, sin, log, quad, polyroots)

mp.dps = 50

# ============================ eta, s4, c(tau) ============================
def eta(tau, nterms=400):
    q = exp(2 * pi * 1j * tau)
    p = mpc(1)
    qn = q
    for n in range(1, nterms + 1):
        p *= (1 - qn)
        qn *= q
        if abs(qn) < mpf(10) ** (-60):
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
        if abs(term) < mpf(10) ** (-60):
            break
    return 2 * pi ** 3 * tot

def EK4(tau):
    return (2 * pi * tau + (10 / pi ** 3) * (U(1, tau) - 2 * U(2, tau))).imag

# ============================ astroid / cross ============================
AST = mpf(4) ** (mpf(2) / 3)

def astroid_v(c):
    return abs(c.real) ** (mpf(2) / 3) + abs(c.imag) ** (mpf(2) / 3)

def on_cross(c, tol=mpf(10) ** (-25)):
    return (abs(c.imag) < tol and abs(c.real) <= 4) or \
           (abs(c.real) < tol and abs(c.imag) <= 4)

# ============================ direct Jensen integral ============================
def mahler4(c, maxdegree):
    """m(x^4+y^4+z^4+1+c x y z) by z-Jensen + tanh-sinh over [0,pi]^2
    (integrand is pi-periodic: t -> t+pi sends r -> -r)."""
    def g(t1, t2):
        x = exp(1j * t1)
        yy = exp(1j * t2)
        roots = polyroots([1, 0, 0, c * x * yy, x ** 4 + yy ** 4 + 1])
        s = mpf(0)
        for r in roots:
            aa = abs(r)
            if aa > 1:
                s += log(aa)
        return s
    val = quad(lambda t1: quad(lambda t2: g(t1, t2), [0, pi], maxdegree=maxdegree),
               [0, pi], maxdegree=maxdegree)
    return val / pi ** 2

def fmtc(z, n=8):
    f = "%+." + str(n) + "f"
    return "(" + (f % z.real) + (f % z.imag) + "j)"

# ============================ Part 1: branch check ============================
print("=" * 78)
print("PART 1: c(tau)^4 = s4(tau) (Delta-quotient), target 40 digits")
print("=" * 78)
for tau in [mpc(mpf(1) / 8, mpf("0.4")), mpc(-mpf(1) / 4, mpf("0.55")),
            mpc(0, mpf("0.7")), mpc(mpf(3) / 8, mpf("0.33")),
            mpc(-mpf(3) / 8, mpf("0.6"))]:
    d = abs(c_of_tau(tau) ** 4 - s4(tau))
    print("  tau = %s   |c^4 - s4| = %.2e   %s" %
          (fmtc(tau, 4), mpf(d), "OK(>40d)" if d < mpf(10) ** (-40) else "FAIL"))

# ============================ Part 2: grid scan, pick points ============================
print()
print("=" * 78)
print("PART 2: grid scan of c(tau); classification (astroid / cross)")
print("=" * 78)
RES = [0, mpf(1) / 8, -mpf(1) / 8, mpf(1) / 4, -mpf(1) / 4, mpf(3) / 8, -mpf(3) / 8]
IMS = [mpf("0.3"), mpf("0.4"), mpf("0.5"), mpf("0.6")]
interior_pts, exterior_pts, cross_pts = [], [], []
for re_ in RES:
    for im_ in IMS:
        tau = mpc(re_, im_)
        c = c_of_tau(tau)
        v = astroid_v(c)
        cr = on_cross(c)
        tag = ("IN " if v < AST else "out") + (" CROSS" if cr else "")
        print("  tau=(%6s,%4s)  c=%s  ast=%8.5f  %s" %
              (re_, im_, fmtc(c, 5), v, tag))
        if v < AST and not cr:
            interior_pts.append(tau)
        elif v >= AST and not cr:
            exterior_pts.append(tau)
        else:
            cross_pts.append(tau)

# selection: >=3 interior non-cross (take all, capped at 5);
# >=2 exterior non-cross, mixed: one near-real c plus complex ones.
sel_in = interior_pts[:5]
ext_real = [t for t in exterior_pts if abs(c_of_tau(t).imag) < mpf("0.01")]
ext_cplx = [t for t in exterior_pts if abs(c_of_tau(t).imag) >= mpf("0.01")]
sel_out = ext_real[:1] + ext_cplx[:2]
if len(sel_out) < 2:
    sel_out = exterior_pts[:3]
print()
print("selected interior non-cross points:", [str(t) for t in sel_in])
print("selected exterior non-cross points:", [str(t) for t in sel_out])
assert len(sel_in) >= 3, "need >=3 interior non-cross sample points"
assert len(sel_out) >= 2, "need >=2 exterior non-cross sample points"

# ============================ Part 3: EK4 vs 4 m(P+c) ============================
# Precision plan (calibrated by _probe3_time.py): EK4/eta at 50 dps; the
# direct Jensen integral at 35 dps / tanh-sinh degree 4 (~60-90 s per
# point).  Interior points have a log-singular integrand (zeros on T^3),
# so the integral itself is good to ~1e-3 there; exterior smooth points to
# ~1e-15.  A machinery control at tau = i (s4 = 648, proved region) is
# included first.
print()
print("=" * 78)
print("PART 3: EK4(tau) vs 4 m(P + c(tau))")
print("        eta/EK4 at 50 dps; direct integral at 35 dps / deg 4")
print("=" * 78)
t0 = time.time()
c_i = c_of_tau(mpc(0, 1))
ek_i = EK4(mpc(0, 1))
with mp.workdps(35):
    mv_i = 4 * mahler4(mpc(c_i), 4)
print("control tau=i (s4=648, PROVED region): EK4 = %.18f  4m = %.18f  |diff| = %.2e  (%.0fs)"
      % (ek_i, mv_i, abs(ek_i - mv_i), time.time() - t0))
print()
print("%-22s %-9s %-26s %-26s %-10s" % ("c(tau)", "class", "EK4(tau)", "4m(P+c)", "|diff|"))
results = []
for tau in sel_in + sel_out:
    c = c_of_tau(tau)
    v = astroid_v(c)
    cls = ("IN " if v < AST else "out") + ("/CROSS" if on_cross(c) else "")
    t0 = time.time()
    ek = EK4(tau)
    with mp.workdps(35):
        mv = 4 * mahler4(mpc(c), 4)
    dt = time.time() - t0
    d = abs(ek - mv)
    results.append((tau, c, cls, ek, mv, d))
    print("%-22s %-9s %-26s %-26s %-10.2e   (%.0fs)" %
          (fmtc(c, 6), cls, ek, mv, mpf(d), dt))

# ============================ Part 4: vertical path scan ============================
print()
print("=" * 78)
print("PART 4: path scan  Re tau = 1/4, Im tau = 0.70 -> 0.20 (20 points)")
print("        EK4 at 50 dps; direct integral at 30 dps / deg 3 (speed)")
print("=" * 78)
print("%-6s %-24s %-9s %-24s %-24s %-10s" %
      ("Im", "c(tau)", "class", "EK4(tau)", "4m(P+c)", "|diff|"))
path_rows = []
NP = 20
for k in range(NP):
    im_ = mpf("0.70") - (mpf("0.70") - mpf("0.20")) * k / (NP - 1)
    tau = mpc(mpf(1) / 4, im_)
    c = c_of_tau(tau)
    v = astroid_v(c)
    cls = ("IN " if v < AST else "out") + ("/CROSS" if on_cross(c, mpf(10) ** (-8)) else "")
    ek = EK4(tau)
    with mp.workdps(30):
        cl = mpc(c)
        mv = 4 * mahler4(cl, 3)
    d = abs(ek - mv)
    path_rows.append((im_, c, cls, ek, mv, d))
    print("%-6.3f %-24s %-9s %-24s %-24s %-10.2e" %
          (im_, fmtc(c, 5), cls, ek, mv, mpf(d)))

# trajectory summary: crossings
print()
print("trajectory summary:")
prev = None
for im_, c, cls, ek, mv, d in path_rows:
    state = "IN" if astroid_v(c) < AST else "out"
    near_cross = min(abs(c.real), abs(c.imag)) < mpf("0.05")
    if prev is not None and (prev[0] != state or (near_cross and not prev[1])):
        print("  near Im tau = %.3f: state %s, c = %s%s" %
              (im_, state, fmtc(c, 5), "  (close to cross)" if near_cross else ""))
    prev = (state, near_cross)
worst = max(path_rows, key=lambda r: r[5])
print("  largest |EK4 - 4m| on path: %.2e at Im tau = %.3f, c = %s (%s)" %
      (mpf(worst[5]), worst[0], fmtc(worst[1], 5), worst[2]))

# ============================ Part 5: verdict ============================
print()
print("=" * 78)
print("PART 5: VERDICT")
print("=" * 78)
in_diffs = [d for (t, c, cl, e, m, d) in results if cl.startswith("IN")]
out_diffs = [d for (t, c, cl, e, m, d) in results if cl.startswith("out")]
print("interior non-cross sample points: worst |diff| = %.2e" % max(in_diffs))
print("exterior non-cross sample points: worst |diff| = %.2e" % max(out_diffs))
print("path scan: worst |diff| = %.2e" % max(r[5] for r in path_rows))
