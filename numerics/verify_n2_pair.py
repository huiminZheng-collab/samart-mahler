# Verify Samart's Table-4 conjecture (still unproved):
#   n2(s) := 2 m( (x+1/x)(y+1/y)(z+1/z) + sqrt(s) )
#        =? (4/7) ( 54 M7 + d7 ),
#   s = (47 + 45 sqrt(-7))/2,   CM point tau' = (1 + sqrt(-7))/8,
#   conjugate pair tau'' = (-1 + sqrt(-7))/8  <->  s-bar.
#   M7 = L'(g7,0),  d7 = L'(chi_{-7}, -1).
#
# All numerics mpmath at 50 dps; report >= 40 digits.

from mpmath import (mp, mpf, mpc, pi, exp, log, sqrt, sin, cos, quad,
                    zeta, bernoulli, fac, nstr)
from fractions import Fraction

mp.dps = 50
DG = 45  # digits to print

# ----------------------------------------------------------------------
# Part 1: s2(tau) = -Delta(tau+1/2)/Delta(2 tau+1) at tau', tau''
# ----------------------------------------------------------------------
def eta(tau, nterms=500):
    q = exp(2 * pi * 1j * tau)
    p = mpc(1)
    qn = q
    for n in range(1, nterms + 1):
        p *= (1 - qn)
        qn *= q
        if abs(qn) < mpf(10) ** (-70):
            break
    return exp(pi * 1j * tau / 12) * p

def s2(tau):
    return -(eta(tau + mpf(1) / 2)) ** 24 / (eta(2 * tau + 1)) ** 24

s7 = sqrt(mpc(-7))
tau_p = (1 + s7) / 8
tau_m = (-1 + s7) / 8
s_target = (47 + 45 * s7) / 2

s2_p = s2(tau_p)
s2_m = s2(tau_m)
print("=== Part 1: s2 at the CM points ===")
print("tau'  =", nstr(tau_p, 20))
print("s2(tau')          =", nstr(s2_p, DG))
print("(47+45sqrt(-7))/2 =", nstr(s_target, DG))
print("diff =", nstr(abs(s2_p - s_target), 5))
print("s2(tau'')          =", nstr(s2_m, DG))
print("(47-45sqrt(-7))/2 =", nstr((47 - 45 * s7) / 2, DG))
print("diff =", nstr(abs(s2_m - (47 - 45 * s7) / 2), 5))
print("|s| =", nstr(abs(s_target), 10))

# ----------------------------------------------------------------------
# Part 2: direct Mahler integral with complex c = sqrt(s).
# Jensen over t3: with A = 8 cos t1 cos t2 (real) and complex c,
#   (1/2pi) int log|A cos t + c| dt = m_z( (A/2) z^2 + c z + A/2 )
#   = log(max(|c + sqrt(c^2 - A^2)|, |c - sqrt(c^2 - A^2)|) / 2).
# (Product of the two roots z+ z- = 1, so |z_max| >= 1 and the formula
#  reduces to log(|A|/2) + log max(1, |z_max|) = the above.  For real
#  c > |A| it gives log((c + sqrt(c^2 - A^2))/2), matching mahler_m.py.)
# Reciprocity m(f+c) = m(f-c) swaps the two branches, so any branch of
# sqrt(s) (and of the inner sqrt) gives the same answer.
# Since Im c^2 = 45 sqrt(7)/2 != 0 and A^2 is real, c^2 - A^2 never
# vanishes: the integrand is analytic, no cusp splitting needed.
# ----------------------------------------------------------------------
c = sqrt(s_target)          # either branch OK by reciprocity
print("\n=== Part 2: direct Mahler integral, c = sqrt(s) ===")
print("c =", nstr(c, 20), "  |c| =", nstr(abs(c), 10))

def Jc(A):
    w = sqrt(c * c - A * A)
    return log(max(abs(c + w), abs(c - w)) / 2)

def inner(t1):
    cc = cos(t1)
    return quad(lambda t2: Jc(8 * cc * cos(t2)), [0, pi / 2])

mval = 4 / pi ** 2 * quad(inner, [0, pi / 2])
n2_direct = 2 * mval
print("n2(s) = 2 m(f + sqrt s) =", nstr(n2_direct, DG))

# ----------------------------------------------------------------------
# Part 3: EK formula at tau' (Samart Prop 2.1(i), Poisson trick)
# ----------------------------------------------------------------------
print("\n=== Part 3: EK(tau') ===")
def U(j, tau):
    tot = mpc(0)
    for m in range(1, 3000):
        z = j * pi * m * tau
        term = cos(z) / sin(z) ** 3 / m
        tot += term
        if abs(term) < mpf(10) ** (-65):
            break
    return 2 * pi ** 3 * tot

def EK(tau):
    return (2 * pi * tau + (2 / pi ** 3) * (U(1, tau) - 4 * U(4, tau))).imag

EK_p = EK(tau_p)
EK_m = EK(tau_m)
print("EK(tau')  =", nstr(EK_p, DG))
print("EK(tau'') =", nstr(EK_m, DG), " (conjugate check: same, since n2(s-bar)=n2(s))")
print("EK(tau') - n2(s) =", nstr(EK_p - n2_direct, 5))

# --- Part 3b: the VALID CM point for the EK identity ------------------------
# Samart's Prop 2.1(i) is proven only for Im tau >= 1/2, and tau' has
# Im = sqrt7/8 < 1/2: EK(tau') is a continuation value on another branch
# (numerically EK(tau') = (8/7)(44 M7 - d7) != n2(s), see diag_ek_branch.py).
# The second preimage of s under s2 on X_0(4) (s2 has degree 2; partner under
# w4: tau -> -1/(4 tau), since lam(-1/(2t)) = 1 - lam(2t) leaves lam(1-lam)
# invariant) is tau_w = -1/(4 tau') = (-1+sqrt(-7))/4 with Im = sqrt7/4 > 1/2,
# INSIDE the proven region, and s2(tau_w) = s exactly.  Check EK there.
print("\n=== Part 3b: EK at tau_w = -1/(4 tau') = (-1+sqrt(-7))/4 (Im > 1/2) ===")
tau_w = -1 / (4 * tau_p)
print("tau_w =", nstr(tau_w, 20), " Im =", nstr(tau_w.imag, 10))
print("|s2(tau_w) - s| =", nstr(abs(s2(tau_w) - s_target), 5))
EK_w = EK(tau_w)
print("EK(tau_w)        =", nstr(EK_w, DG))
print("EK(tau_w) - n2(s) =", nstr(EK_w - n2_direct, 5))

# ----------------------------------------------------------------------
# Part 4: d7 = L'(chi_{-7}, -1), independently.
# chi_{-7} odd, conductor 7. Completed L:  Lam(s) = (7/pi)^{(s+1)/2}
# Gamma((s+1)/2) L(s),  Lam(s) = Lam(1-s).
# L(chi,-1) = -B_{2,chi}/2 = 0 (parity), so near s = -1:
#   Lam(-1) = 2 L'(chi,-1)   [residue of Gamma((s+1)/2) is 2]
#   Lam(2)  = (7/pi)^{3/2} Gamma(3/2) L(chi,2) = (7^{3/2} sqrt(pi)/4) L(2)
# => d7 = L'(chi,-1) = 7 sqrt(7)/(4 pi) * L(chi_{-7}, 2).
# ----------------------------------------------------------------------
print("\n=== Part 4: d7 = L'(chi_{-7}, -1) ===")
def chi(n):
    n %= 7
    if n == 0:
        return 0
    return 1 if pow(n, 3, 7) == 1 else -1   # (n/7) = n^3 mod 7
assert [chi(n) for n in range(7)] == [0, 1, 1, -1, 1, -1, -1]

# --- sanity: L(chi,-1) = -B_{2,chi}/2 must be 0 (parity) ---
B2chi = sum(Fraction(chi(a)) * (Fraction(a * a, 49) - Fraction(a, 7)
                                 + Fraction(1, 6)) for a in range(1, 8)) * 7
print("B_{2,chi} =", B2chi, " -> L(chi,-1) =", -B2chi / 2, "(expect 0)")

# --- L(chi_{-7}, 2), two independent ways ---
# (a) Hurwitz zeta decomposition:  L = 7^{-2} sum_a chi(a) zeta(2, a/7)
L2_hur = mpf(7) ** (-2) * sum(chi(a) * zeta(2, mpf(a) / 7) for a in range(1, 8))

# (b) truncation at 7K + Euler-Maclaurin tail per residue class
K, J = 150, 15
direct = sum(mpf(chi(n)) / mpf(n) ** 2 for n in range(1, 7 * K + 1))
tail = mpf(0)
for a in range(1, 8):
    if chi(a) == 0:
        continue
    u = 7 * K + a
    ta = mpf(1) / (7 * u) + mpf(1) / (2 * mpf(u) ** 2)   # mpf: avoid float64 division
    for j in range(1, J + 1):
        ta += mpf(bernoulli(2 * j)) * 7 ** (2 * j - 1) * mpf(u) ** (-(2 * j + 1))
    tail += chi(a) * ta
L2_em = direct + tail
print("L(chi,2) Hurwitz  =", nstr(L2_hur, DG))
print("L(chi,2) E-M      =", nstr(L2_em, DG))
print("diff =", nstr(abs(L2_hur - L2_em), 5))

# --- extra FE sanity at a generic point: Lam(0.3) ?= Lam(0.7) ---
def Lchi(s):
    return mpf(7) ** (-s) * sum(chi(a) * zeta(s, mpf(a) / 7) for a in range(1, 8))
def Lam(s):
    from mpmath import gamma
    return (7 / pi) ** ((s + 1) / 2) * gamma((s + 1) / 2) * Lchi(s)
print("FE check |Lam(0.3)-Lam(0.7)| =", nstr(abs(Lam(mpf("0.3")) - Lam(mpf("0.7"))), 5))
# and L(chi,0) = -B_{1,chi} = 1 (= 2h(-7)/w):
print("L(chi,0) =", nstr(Lchi(0), 10), " (expect 1)")

d7 = 7 * sqrt(mpf(7)) / (4 * pi) * L2_hur
print("d7 = L'(chi_{-7},-1) = 7 sqrt7/(4 pi) L(chi,2) =", nstr(d7, DG))

# ----------------------------------------------------------------------
# Part 5: final comparison
# ----------------------------------------------------------------------
print("\n=== Part 5: final comparison ===")
M7 = mpf("0.10267160777890201121045659489829291399889482708922")
rhs = mpf(4) / 7 * (54 * M7 + d7)
print("n2(s)              =", nstr(n2_direct, DG))
print("EK(tau_w)          =", nstr(EK_w, DG))
print("EK(tau')           =", nstr(EK_p, DG), " (continuation branch, != n2(s))")
print("(4/7)(54 M7 + d7)  =", nstr(rhs, DG))
print("n2(s) - rhs        =", nstr(n2_direct - rhs, 5))
print("EK(tau_w) - rhs    =", nstr(EK_w - rhs, 5))
print("n2(s) - EK(tau_w)  =", nstr(n2_direct - EK_w, 5))
print("EK(tau') - rhs     =", nstr(EK_p - rhs, 5))
print("n2(s) - EK(tau')   =", nstr(n2_direct - EK_p, 5))
# branch gap has a closed form: pslq gives 7*(n2-EK(tau')) = 12 d7 - 136 M7,
# i.e. EK(tau') = (8/7)(44 M7 - d7) -- itself in the span of {M7, d7}:
print("(n2-EK(tau')) - (12 d7 - 136 M7)/7 =",
      nstr((n2_direct - EK_p) - (12 * d7 - 136 * M7) / 7, 5))
