#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cert_fricke.py -- exact root numbers w = +1 for all 15 newforms of
Theorems C, D, E, G, H, I, J, K, via ONE-ORDINATE INTERVAL Fricke locks.

Mathematical basis (paper, "Levels and root numbers: exact determination"):
  * Each form F is a newform of S_3(Gamma_0(N), chi) with quadratic chi
    and REAL coefficients (exact constructions below).  By Atkin--Lehner,
    F|_3 W_N = w F with w in {+1,-1}, and w is the root number.
    Hence the Fricke ratio
        R_N(y) = F(i/(sqrt(N) y)) / (y^3 F(i y/sqrt(N)))
    equals w for EVERY y > 0 (exactly, not approximately).
  * A single rigorous interval enclosure of R_N(y0) of width < 1
    therefore determines the integer w exactly.
  * The LEVELS are exact independently (Hecke's theorem: theta series of
    a grossencharacter of conductor f over K has level |d_K| N(f);
    conductors (1), (2), (4) are determined by exact mod arithmetic --
    see the paper's conductor table).  The two-ordinate numeric Fricke
    scan of earlier versions is superseded; here only the SIGN is locked.

Tracks:
  [C*]  exact coefficient constructions (pure integer / Fraction
        arithmetic), including the integrity pins of the verify_P1_*
        scripts (theta identities to q^60, integrality, a(1)=1, sign
        conventions, Hecke recursion at small primes).
  [R*]  interval Fricke locks: R_N(y) computed with mpmath.iv (dps 70)
        with SELF-CONTAINED tail bounds:
          eta products:  |a_n| <= (n+5)^5 2^{n/2}
          theta (disc -24):  |a(n)| <= 6 n^3
          theta (disc -15, -84): |a(n)| <= 30 n^3
        (same bounds as the [V*] tracks of verify_P1_*.py; quoted there
        with proofs).  Check: 1 in enclosure, half-width < 1e-20.
  [X]   mp 60-dps cross-check of every ratio (NOT part of the proof).

Output: one line per (form, ordinate) + ALL CHECKS PASSED.
"""

from mpmath import mp, mpf, sqrt
from mpmath import iv
from fractions import Fraction as Fr

mp.dps = 60
iv.dps = 70

FAILS = []
def check(name, cond):
    print(("  [ok] " if cond else "  [FAIL] ") + name)
    if not cond:
        FAILS.append(name)

def check_lock(name, ivl, target=1):
    hw = (ivl.b - ivl.a) / 2
    ok = (ivl.a <= target <= ivl.b) and hw < mpf(10)**(-20)
    print(("  [ok] " if ok else "  [FAIL] ")
          + "%s : contains %d, half-width = %.3e" % (name, target, mpf(hw)))
    if not ok:
        FAILS.append(name)

NMAX = 800          # coefficients and series truncation

# =====================================================================
# PART 0: exact coefficient constructions (integers / Fractions)
# =====================================================================

# --- 0a. eta products: g8 = eta^2 eta(2) eta(4) eta(8)^2, g16 = eta(4)^6
def Pser(d, N):
    res = [0]*(N+1); res[0] = 1
    k = 1
    while True:
        e1, e2 = k*(3*k-1)//2*d, k*(3*k+1)//2*d
        if e1 > N and e2 > N: break
        if e1 <= N: res[e1] += (-1)**k
        if e2 <= N: res[e2] += (-1)**k
        k += 1
    return res

def s_mul(A, B, N):
    C = [0]*(N+1)
    for i in range(N+1):
        if A[i] == 0: continue
        for j in range(N+1-i):
            if B[j]: C[i+j] += A[i]*B[j]
    return C

def s_pow(A, e, N):
    R = [0]*(N+1); R[0] = 1
    for _ in range(e): R = s_mul(R, A, N)
    return R

g8ser = s_mul(s_mul(s_pow(Pser(1, NMAX), 2, NMAX), Pser(2, NMAX), NMAX),
              s_mul(Pser(4, NMAX), s_pow(Pser(8, NMAX), 2, NMAX), NMAX), NMAX)
a8 = [0] + g8ser[:-1]
g16ser = s_pow(Pser(4, NMAX), 6, NMAX)
a16 = [0] + g16ser[:-1]
chi8t = [0, 1, 0, -1, 0, -1, 0, 1]        # chi_8 = (2/.), conductor 8
a8tw = [chi8t[n % 8]*a8[n] for n in range(NMAX+1)]
a16tw = [chi8t[n % 8]*a16[n] for n in range(NMAX+1)]

# g12 = eta(2t)^3 eta(6t)^3  (Theorem B; K = Q(sqrt-3), conductor (2))
a12ser = s_mul(s_pow(Pser(2, NMAX), 3, NMAX),
               s_pow(Pser(6, NMAX), 3, NMAX), NMAX)
a12 = [0] + a12ser[:-1]

# theta pins (exact, q^60 > Sturm bounds 3, 6)
NT = 60
a_th8 = [0]*(NT+1)
for a in range(-9, 10):
    for b in range(-9, 10):
        Nm = a*a + 2*b*b
        if 1 <= Nm <= NT: a_th8[Nm] += a*a - 2*b*b
check("[C1] 2 g8 = theta over Q(sqrt-2) to q^60 (exact)",
      all(a_th8[n] == 2*a8[n] for n in range(1, NT+1)))
a_th16 = [0]*(NT+1)
for a in range(-9, 10):
    for b in range(-9, 10):
        Nm = a*a + b*b
        if 1 <= Nm <= NT and a % 2 == 1 and b % 2 == 0:
            a_th16[Nm] += a*a - b*b
check("[C2] 2 g16 = theta over Q(i), a==1(2), to q^60 (exact)",
      all(a_th16[n] == 2*a16[n] for n in range(1, NT+1)))
a_th12 = [0]*(NT+1)
for a in range(-9, 10):
    for b in range(-9, 10):
        # alpha = a + b w over Q(sqrt-3), N = a^2 - ab + b^2;
        # alpha == 1 mod 2: a odd, b even; Re(alpha^2) = a^2 - b^2
        Nm = a*a - a*b + b*b
        if 1 <= Nm <= NT and (a % 2, b % 2) == (1, 0):
            a_th12[Nm] += a*a - b*b
check("[C2b] 2 g12 = theta over Q(sqrt-3), a==1(2), to q^60 (exact)",
      all(a_th12[n] == 2*a12[n] for n in range(1, NT+1)))

# --- 0b. disc -24: g1 = P - Q, g2 = P + Q
RB = int(NMAX**0.5) + 2
P = [0]*(NMAX+1); Q = [0]*(NMAX+1)
for a in range(-RB, RB+1):
    for b in range(-RB, RB+1):
        n = a*a + 6*b*b
        if 1 <= n <= NMAX: P[n] += a*a - 6*b*b
for x in range(-RB, RB+1):
    for y in range(-RB, RB+1):
        n = 2*x*x + 3*y*y
        if 1 <= n <= NMAX: Q[n] += 4*x*x - 6*y*y
ok = all(P[n] % 2 == 0 and Q[n] % 4 == 0 for n in range(1, NMAX+1))
P = [p//2 for p in P]; Q = [q//4 for q in Q]
a1_24 = [P[n] - Q[n] for n in range(NMAX+1)]
a2_24 = [P[n] + Q[n] for n in range(NMAX+1)]
check("[C3] disc -24: P,Q integrality, a1(1)=a2(1)=1 (exact)",
      ok and a1_24[1] == a2_24[1] == 1)
chi8m = lambda n: 0 if n % 2 == 0 else (1 if n % 8 in (1, 3) else -1)
a1_24tw = [chi8m(n)*a1_24[n] for n in range(NMAX+1)]
a2_24tw = [chi8m(n)*a2_24[n] for n in range(NMAX+1)]

# --- 0c. disc -15: a1 = P + (ombar/8) R, a2 = P - (ombar/8) R
RB2 = int((8*NMAX)**0.5) + 2
Pc = [Fr(0)]*(NMAX+1); Rr = [Fr(0)]*(NMAX+1); Ri = [Fr(0)]*(NMAX+1)
for c in range(-RB2, RB2+1):
    for d in range(-RB2, RB2+1):
        nn = c*c + 15*d*d
        if nn % 4 == 0 and 1 <= nn//4 <= NMAX and (c-d) % 2 == 0:
            Pc[nn//4] += Fr(c*c - 15*d*d, 8)
        if nn % 8 == 0 and 1 <= nn//8 <= NMAX and (c-d) % 4 == 0:
            Rr[nn//8] += Fr(c*c - 15*d*d, 4)
            Ri[nn//8] += Fr(2*c*d, 4)
a1re = [Pc[n] + (Rr[n] + 15*Ri[n])/16 for n in range(NMAX+1)]
a2re = [Pc[n] - (Rr[n] + 15*Ri[n])/16 for n in range(NMAX+1)]
a1im = [(Ri[n] - Rr[n])/16 for n in range(NMAX+1)]
a2im = [-(Ri[n] - Rr[n])/16 for n in range(NMAX+1)]
ok = (all(v == 0 for v in a1im[1:]) and all(v == 0 for v in a2im[1:])
      and all(v.denominator == 1 for v in a1re[1:])
      and all(v.denominator == 1 for v in a2re[1:])
      and a1re[1] == a2re[1] == 1)
a1_15 = [int(v) for v in a1re]
a2_15 = [int(v) for v in a2re]
check("[C4] disc -15: a1,a2 rational integral, a(1)=1 (exact)", ok)
check("[C5] disc -15: a1(2)=+1, a2(2)=-1, a1(5)=+5, a2(5)=-5 (exact)",
      a1_15[2] == 1 and a2_15[2] == -1 and a1_15[5] == 5 and a2_15[5] == -5)

# --- 0d. disc -84: four forms g(e2,e3) = P0 + e2 P2 + e3 P3 + e2e3 P6
RB3 = int((2*NMAX)**0.5) + 3
P0 = [Fr(0)]*(NMAX+1); P2 = [Fr(0)]*(NMAX+1)
P3 = [Fr(0)]*(NMAX+1); P6 = [Fr(0)]*(NMAX+1)
for a in range(-RB3, RB3+1):
    for b in range(-RB3, RB3+1):
        n = a*a + 21*b*b
        if 1 <= n <= NMAX: P0[n] += Fr(a*a - 21*b*b, 2)
        n = 3*a*a + 7*b*b
        if 1 <= n <= NMAX: P3[n] += Fr(3*a*a - 7*b*b, 2)
for x in range(-RB3, RB3+1):
    for y in range(-RB3, RB3+1):
        n = 2*x*x + 2*x*y + 11*y*y
        if 1 <= n <= NMAX:
            v = (2*x+y)**2 - 21*y*y
            assert v % 4 == 0; P2[n] += Fr(v, 4)
        n = 6*x*x + 6*x*y + 5*y*y
        if 1 <= n <= NMAX:
            v = (6*x+3*y)**2 - 21*y*y
            assert v % 12 == 0; P6[n] += Fr(v, 12)
ok = (all(v.denominator == 1 for PP in (P0, P2, P3, P6) for v in PP[1:])
      and P0[1] == 1)
a84 = {}
for e2 in (1, -1):
    for e3 in (1, -1):
        a84[(e2, e3)] = [int(P0[n] + e2*P2[n] + e3*P3[n] + e2*e3*P6[n])
                         for n in range(NMAX+1)]
check("[C6] disc -84: P0,P2,P3,P6 integral, a(1)=1 all four (exact)", ok)
check("[C7] disc -84: sign pins a_(1,1)(2)=+2, a_(1,1)(3)=+3 (exact)",
      a84[(1, 1)][2] == 2 and a84[(1, 1)][3] == 3
      and a84[(-1, 1)][2] == -2 and a84[(1, -1)][3] == -3)

# Hecke recursion pins at small primes (exact; guards against transcription)
def chi_3m(n): return 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)
def chi_8m(n): return 0 if n % 2 == 0 else (1 if n % 8 in (1, 3) else -1)
def chi_4m(n): return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def chi_24m(n): return 0 if n % 2 == 0 or n % 3 == 0 else \
    (1 if n % 24 in (1, 5, 7, 11) else -1)
def chi_15(n): return 0 if n % 3 == 0 or n % 5 == 0 else \
    (1 if n % 15 in (1, 2, 4, 8) else -1)
def chi_84(n):
    if n % 2 == 0 or n % 3 == 0 or n % 7 == 0: return 0
    return 1 if n % 84 in (1, 5, 11, 17, 19, 23, 25, 31, 37, 41, 55, 71) else -1
PR = []
for n in range(2, NMAX+1):
    if all(n % p for p in PR if p*p <= n): PR.append(n)
ok = True
for aa, ch, Nlv in ((a12, chi_3m, 12), (a8, chi_8m, 8), (a16, chi_4m, 16),
                    (a1_24, chi_24m, 24), (a2_24, chi_24m, 24),
                    (a1_15, chi_15, 15), (a2_15, chi_15, 15),
                    (a84[(1, 1)], chi_84, 84), (a84[(-1, -1)], chi_84, 84)):
    for p in PR:
        if p*p > NMAX: break
        if Nlv % p == 0:
            continue            # U_p at p | N: different recursion
        if aa[p*p] != aa[p]*aa[p] - ch(p)*p*p: ok = False
check("[C8] Hecke recursion a(p^2)=a(p)^2-chi(p)p^2, p^2<=800, p coprime "
      "to level (exact)", ok)

# =====================================================================
# PART 1: interval Fricke locks (mpmath.iv, dps = 70)
# =====================================================================
rt2_iv = iv.sqrt(2)

def qser_iv(a, t, kind, n0):
    """Enclosure of sum a_n e^{-t n}; t > 0 an iv.mpf.
    kind: 'eta' -> |a_n| <= (n+5)^5 2^{n/2};
          't6'  -> |a(n)| <= 6 n^3;  't30' -> |a(n)| <= 30 n^3."""
    tot = iv.mpf(0)
    for n in range(1, min(n0, len(a)-1)+1):
        if a[n]:
            tot += a[n] * iv.exp(-t*n)
    if kind == 'eta':
        r = rt2_iv * iv.exp(-t)
        rho = (r * ((iv.mpf(n0+7))/(iv.mpf(n0+6)))**5).b
        assert rho < 1
        T = ((iv.mpf(n0+6))**5 * r**(n0+1) / (1 - rho)).b
    else:
        C = 6 if kind == 't6' else 30
        r = iv.exp(-t)
        rho = (r * ((iv.mpf(n0+2))/(iv.mpf(n0+1)))**3).b
        assert rho < 1 and n0*t > 3      # n^3 e^{-tn} decreasing
        T = (C * (iv.mpf(n0+1))**3 * r**(n0+1) / (1 - rho)).b
    return tot + iv.mpf([-T, T])

def fricke_iv(a, N, yv, kind, n0=NMAX):
    xN = iv.sqrt(iv.mpf(N))
    y = iv.mpf(yv)
    num = qser_iv(a, 2*iv.pi/(xN*y), kind, n0)
    den = y**3 * qser_iv(a, 2*iv.pi*y/xN, kind, n0)
    return num/den

def fricke_mp(a, N, yv):
    xN = sqrt(mpf(N))
    def qs(t):
        return mpf(0) + sum(a[n]*mp.exp(-t*n)
                            for n in range(1, NMAX+1) if a[n])
    return qs(2*mp.pi/(xN*yv))/(yv**3*qs(2*mp.pi*yv/xN))

FORMS = [
    ("g8",            a8,       8,  'eta'),
    ("g8 x chi8",     a8tw,     32, 'eta'),
    ("g16",           a16,      16, 'eta'),
    ("g16 x chi8",    a16tw,    64, 'eta'),
    ("g12",           a12,      12, 'eta'),
    ("g1 (disc -24)", a1_24,    24, 't6'),
    ("g2 (disc -24)", a2_24,    24, 't6'),
    ("g1tw (N=96)",   a1_24tw,  96, 't6'),
    ("g2tw (N=96)",   a2_24tw,  96, 't6'),
    ("g1 (disc -15)", a1_15,    15, 't30'),
    ("g2 (disc -15)", a2_15,    15, 't30'),
    ("g84(+,+)",      a84[(1, 1)],   84, 't30'),
    ("g84(+,-)",      a84[(1, -1)],  84, 't30'),
    ("g84(-,+)",      a84[(-1, 1)],  84, 't30'),
    ("g84(-,-)",      a84[(-1, -1)], 84, 't30'),
]

print()
print("=== [R] interval Fricke locks (iv.dps = 70, n0 = %d) ===" % NMAX)
for label, a, N, kind in FORMS:
    for ys in ("0.6", "1.1"):
        R = fricke_iv(a, N, mpf(ys), kind)
        check_lock("[R:%s, N=%d, y=%s]" % (label, N, ys), R, 1)

print()
print("=== [X] mp 60-dps cross-check (NOT part of the proof) ===")
worst = mpf(0)
for label, a, N, kind in FORMS:
    r = fricke_mp(a, N, mpf("0.6"))
    dev = abs(r - 1)
    worst = max(worst, dev)
    print("  %s : |R - 1| = %.3e" % (label, dev))
check("[X] all mp Fricke ratios within 1e-40 of +1", worst < mpf(10)**(-40))

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
else:
    print("ALL FRICKE LOCKS PASSED: w = +1 for all 15 forms "
          "(interval-certified, half-widths < 1e-20)")
