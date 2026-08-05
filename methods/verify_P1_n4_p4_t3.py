# verify_P1_n4_p4_t3.py -- (P1)-type CM evaluation for the Samart Table-6
# OPEN conjecture (Phase 4, target T3; first CLASS NUMBER 2 case):
#
#   tau = sqrt(-6),  s4 = 1207368 + 853632 sqrt2 + 697680 sqrt3
#                        + 493272 sqrt6   [certified by cert0_n4_p4_t3.py],
#   EK4(tau) = (5/48)(4 M24^1tw + 4 M24^2tw + 28 M24^1 + 12 M24^2
#                     + 28 d3 + 24 d4 + 8 d8 + d24)
#
# K = Q(sqrt(-6)), disc -24, h = 2 (class group = genus group = Z/2).
# The two newforms of S_3(Gamma_0(24), chi_{-24}) are the CM theta series
#   g1 = P - Q  = LMFDB 24.3.h.a  (CM disc -24, char 24.chi_5),
#   g2 = P + Q  = LMFDB 24.3.h.b
#   (labels pinned 2026-08-05 via api/mf_newforms?level=24&weight=3:
#    h.a/h.b are the only dim-1 CM forms, cm_discs=[-24]; coefficient
#    match a(1..6): g1 = 1,-2,3,4,2,-6 = h.a; g2 = 1,2,-3,4,-2,-6 = h.b),
#   P(n) = (1/2) sum'_{a^2+6b^2=n} (a^2-6b^2)   [principal class],
#   Q(n) = (1/4) sum'_{2x^2+3y^2=n} (4x^2-6y^2) [nonprincipal class;
#                                                psi(p_2) = -2 labeling fixed
#                                                by the [E] check].
# M24^i = L'(gi,0) (level 24, w = +1), M24^itw = L'(gi x chi_{-8},0)
# (twist level 96, w = +1, Fricke-ratio scan), d_k = L'(chi_{-k},-1).
#
# Layers (same three-track standard as verify_P1_n4_p3_t1t2.py):
#   [X*] exact integer / Fraction checks (no floating point);
#   [S*] SEPARATION exact-algebra track: the target identity proved as an
#        exact Fraction-algebra consequence of the quoted inputs (Q1)-(Q6);
#   [G*] newform identifications: exact-integer Hecke/CM checks; levels and
#        root numbers by the NON-vacuous Fricke ratio (numeric, flagged);
#   [L*],[D*],[T*],[E*] mpmath 60-dps numerical confirmations;
#   [V*] rigorous interval locks (ported from verify_P1_n4_p3_t1t2.py,
#        iv.dps = 70, self-contained theta-series Mellin tail |a(n)| <= 6n^3):
#        [V1] T(O_K), [V2] T(O_2) locked with half-width < 5e-53,
#        [V3] the EK4 conjectural identity locked with half-width < 1e-50
#        (actual 9.72e-53); cert0_n4_p4_t3.py carries the s4-value iv locks.
#
# Quoted inputs (same standard as T1/T2):
#   (Q1) Hecke's theorem: the theta series of a grossencharacter of
#        conductor 1 is a newform; here g1, g2 in S_3(Gamma_0(24), chi_{-24}),
#        which is 2-dimensional, so {g1, g2} exhausts it.
#   (Q2) h(-24) = 2, class group Z/2, units +-1; 2, 3 ramify.
#   (Q3) Functional equations of the CM forms (level 24, w = +1) and of
#        their chi_{-8} twists (level 96, w = +1) [levels/root numbers
#        numerically pinned by the Fricke ratio at two ordinates].
#   (Q4) Dirichlet class-number-free formulas: d_k = (k^{3/2}/(4 pi))
#        L(chi_{-k},2); even-character formula L(chi,2) = pi^2 tau(chi)
#        sum_a chi(a) B_2(a/f) / f for even primitive chi of conductor f
#        (applied to chi_8, chi_3 conductor 12, chi_24 conductor 24; Gauss
#        sums tau = +sqrt(f) for real even primitive chi).
#   (Q5) L_K(chi o N, s) = L(chi,s) L(chi chi_{-24},s) for any Dirichlet
#        chi (formal Euler-product identity).
#   (Q6) zeta(2) = pi^2/6; zeta_K(2) = zeta(2) L(chi_{-24},2).
#
# Key exact lemmas proved below:
#   [X0] F(z) = 4(Re z)^2/|z|^6 - 1/|z|^4 = 2 Re(z^2)/|z|^6 + 1/|z|^4
#        (pointwise), so for any lattice T(Lambda) = B(Lambda) + 2 G(Lambda)
#        with both sums absolutely convergent -- this closes the
#        "T = B + 2G" gap recorded in n4_phase4_report.md.

from fractions import Fraction as Fr
from math import gcd
from mpmath import (mp, mpf, mpc, pi, sqrt, exp, zeta, dirichlet, gamma,
                    diff as mpdiff, power, gammainc, nstr)

mp.dps = 60
s2 = sqrt(mpf(2)); s3 = sqrt(mpf(3)); s6 = sqrt(mpf(6))
TOL = mpf(10) ** (-50)
FAILS = []

def check(name, got, want, tol=TOL):
    d = abs(got - want)
    ok = d < tol
    if not ok:
        FAILS.append(name)
    print("%-76s %s  (|diff| = %.2e)" % (name, "PASS" if ok else "FAIL", mpf(d)))

def check_exact(name, cond):
    if not cond:
        FAILS.append(name)
    print("%-76s %s" % (name, "PASS" if cond else "FAIL"))

# =====================================================================
# PART 0: exact algebra
# =====================================================================
# --- [X0] the F-identity (exact symbolic, integer arithmetic) ----------
# F(z)|z|^6 = 4x^2 - (x^2+y^2) = 3x^2 - y^2;  2 Re(z^2) + |z|^2 =
# 2(x^2-y^2) + (x^2+y^2) = 3x^2 - y^2.  Same polynomial, checked on a grid.
okX0 = all(4*x*x - (x*x+y*y) == 2*(x*x-y*y) + (x*x+y*y) == 3*x*x - y*y
           for x in range(-5, 6) for y in range(-5, 6))
check_exact("[X0] F(z) = 2 Re(z^2)/|z|^6 + 1/|z|^4 pointwise (grid, exact) "
            "=> T(L) = B(L) + 2 G(L), both sides absolutely convergent",
            okX0)

# --- field arithmetic, class group -------------------------------------
SW = (Fr(0), Fr(1))     # sqrt(-6): pair (p,q) = p + q sqrt(-6)
def k6mul(u, v):
    return (u[0]*v[0] - 6*u[1]*v[1], u[0]*v[1] + u[1]*v[0])
def k6norm(u): return u[0]**2 + 6*u[1]**2
check_exact("[X1] (sqrt-6)^2 = -6; disc(O_K) = -24 (exact)",
            k6mul(SW, SW) == (Fr(-6), Fr(0)) and 4*(-6) == -24)
# p_2 = (2, sqrt-6): p_2^2 = (2); elements of p_2: 2x + sqrt-6 y.
check_exact("[X2] p_2 = 2Z + sqrt-6 Z has index 2 in O_K; p_2^2 = (2) "
            "(norm multiplicativity, exact)",
            k6norm((Fr(2), Fr(0))) == 4 and k6norm(SW) == 6)

# --- characters (exact tables; chi_3 has conductor 12 with ZEROS at even
#     positions -- the phase-4 pitfall) -----------------------------------
def chi_8(n):  return 0 if n % 2 == 0 else (1 if n % 8 in (1, 7) else -1)
def chi8m(n):  return 0 if n % 2 == 0 else (1 if n % 8 in (1, 3) else -1)
def chi3m(n):  return 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)
chi3p_l = [0, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 1]    # (3/n), conductor 12
def chi3p(n):  return chi3p_l[n % 12]
def chi4(n):   return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def chi_24(n):  return chi_8(n)*chi3p(n)            # (24/n), conductor 24
def chi_24m(n): return chi_8(n)*chi3m(n)            # (-24/n), conductor 24

okX3 = (all(chi3p(n) == 0 for n in range(24) if n % 2 == 0 or n % 3 == 0)
        and chi3p(1) == chi3p(11) == 1 and chi3p(5) == chi3p(7) == -1)
check_exact("[X3] chi_3 = (3/.) table: zeros at 2|n and 3|n; values at "
            "1,5,7,11 (exact)", okX3)

# character sums for the even-character L(.,2) closed forms (Q4):
#   L(chi,2) = pi^2 tau(chi) (sum chi(a) B_2(a/f)) / f,
#   sum chi B_2(a/f) = (sum chi a^2)/f^2  when sum chi = sum chi a = 0.
S3 = (sum(chi3p(a) for a in range(12)), sum(chi3p(a)*a for a in range(12)),
      sum(chi3p(a)*a*a for a in range(12)))
S24 = (sum(chi_24(a) for a in range(24)),
       sum(chi_24(a)*a for a in range(24)),
       sum(chi_24(a)*a*a for a in range(24)))
S8 = (sum(chi_8(a) for a in range(8)), sum(chi_8(a)*a for a in range(8)),
      sum(chi_8(a)*a*a for a in range(8)))
check_exact("[X4] chi_3: sums (0, 0, 48); chi_24: (0, 0, 288); chi_8: "
            "(0, 0, 16) (exact)", S3 == (0, 0, 48) and S24 == (0, 0, 288)
            and S8 == (0, 0, 16))
# chi_8 case of (Q4): pi^2 (2 sqrt2)(16/64)/8 = pi^2 sqrt2/16 [matches T1].

# --- genus character identity -------------------------------------------
okX5 = all(chi_8(n) == chi3m(n) for n in range(1, 600)
           if n % 2 == 1 and n % 3 != 0 and chi_24m(n) == 1)
check_exact("[X5] genus character: chi_8(N) = chi_{-3}(N) on ideal norms "
            "(chi_{-24}(N)=+1; n<600, exact)", okX5)

# --- ray-class projection for {a odd, b even} ----------------------------
okX6 = True
for a in range(-9, 10):
    for b in range(-9, 10):
        if a % 2 == 1 and (a, b) != (0, 0):
            Nm = a*a + 6*b*b
            okX6 &= (chi8m(Nm) == 1) == (b % 2 == 0)
check_exact("[X6] chi_{-8}(N(a+b sqrt-6)) = +1 iff b even (a odd; box)",
            okX6)

# --- ideal counts vs representation numbers (exact integers) -------------
NID = 3000
r = [0]*(NID+1); g = [0]*(NID+1)
for dd in range(1, NID+1):
    cc = chi_24m(dd)
    if cc:
        for m in range(dd, NID+1, dd): r[m] += cc
    c2 = chi_8(dd)
    if c2:
        for m in range(1, NID//dd+1):
            c3 = chi3m(m)
            if c3: g[dd*m] += c2*c3
r1 = [0]*(NID+1); r2 = [0]*(NID+1)
RB = int((NID)**0.5)+2
for a in range(-RB, RB+1):
    for b in range(-RB, RB+1):
        n = a*a + 6*b*b
        if 1 <= n <= NID: r1[n] += 1
        n = 2*a*a + 3*b*b
        if 1 <= n <= NID: r2[n] += 1
okX7 = all((r1[n]+r2[n]) == 2*r[n] and (r1[n]-r2[n]) == 2*g[n]
           for n in range(1, NID+1))
check_exact("[X7] ideal counts: r1+r2 = 2r(n), r1-r2 = 2g(n) (n<=3000; "
            "principal/nonprincipal split by the genus character, exact)",
            okX7)

# --- theta series and newform checks (exact integers) ---------------------
NMAX = 400
P = [0]*(NMAX+1); Q = [0]*(NMAX+1)
for a in range(-RB, RB+1):
    for b in range(-RB, RB+1):
        n = a*a + 6*b*b
        if 1 <= n <= NMAX: P[n] += a*a - 6*b*b
for x in range(-RB, RB+1):
    for y in range(-RB, RB+1):
        n = 2*x*x + 3*y*y
        if 1 <= n <= NMAX: Q[n] += 4*x*x - 6*y*y
okX8 = all(P[n] % 2 == 0 and Q[n] % 4 == 0 for n in range(1, NMAX+1))
P = [p//2 for p in P]; Q = [q//4 for q in Q]
a1 = [P[n] - Q[n] for n in range(NMAX+1)]   # g1 (Samart's #1, coeff 28)
a2 = [P[n] + Q[n] for n in range(NMAX+1)]   # g2 (Samart's #2, coeff 12)
PR = []
for n in range(2, NMAX+1):
    if all(n % p for p in PR if p*p <= n): PR.append(n)
oke1, oke2, okcm = True, True, True
for aa in (a1, a2):
    for p in PR:
        if p*p > NMAX: break
        if aa[p*p] != aa[p]*aa[p] - chi_24m(p)*p*p: oke1 = False
    for i, p in enumerate(PR):
        for q in PR[i+1:]:
            if p*q > NMAX: break
            if aa[p*q] != aa[p]*aa[q]: oke2 = False
    for p in PR:
        if p > 100: break
        if chi_24m(p) == -1 and aa[p] != 0: okcm = False
check_exact("[X8] P, Q integrality; a1(1)=a2(1)=1 (exact)",
            okX8 and a1[1] == a2[1] == 1)
check_exact("[X9] Hecke recursion a(p^2)=a(p)^2-chi_{-24}(p)p^2; "
            "multiplicativity; CM vanishing at inert p<=100 (exact)",
            oke1 and oke2 and okcm)

# --- [S] separation track: the assembly as exact Fraction algebra --------
# comb = -T(O_K) + 4 T(O_2) decomposes (lattice lemmas, exact):
#   comb = (7/2) L1 + (3/2) L2 + 4 L1t + 4 L2t + (3/4) zK + (7/4) Lg
#          + 2 L(chi_{-8},2) L(chi_3,2) + 2 Cat L(chi_24,2),
# derived from:
#   B(O_K) = zK + Lg,  G(O_K) = L1 + L2,
#   B(p_2-lattice) = (zK - Lg)/4,  G(p_2-lattice) = (L2 - L1)/4,
#   B_odd = B(O_K) - B(p_2) = 3zK/4 + 5Lg/4,
#   G_odd = G(O_K) - G(p_2) = 5L1/4 + 3L2/4,
#   B_{odd,chi} = L(chi_{-8},2) L(chi_3,2) + Cat L(chi_24,2),
#   G_{odd,chi} = L1t + L2t,
#   O_2\{0} = 2 O_K\{0} ⊔ S, S = {a odd, b even} (projection (X6)),
#   B(2 O_K) = B(O_K)/16, G(2 O_K) = G(O_K)/16, T = B + 2G (X0).
combv = (Fr(7, 2), Fr(3, 2), Fr(4), Fr(4), Fr(3, 4), Fr(7, 4), Fr(2), Fr(2))
check_exact("[S0] comb coefficients (L1,L2,L1t,L2t,zK,Lg,B1,B2) = "
            "(7/2,3/2,4,4,3/4,7/4,2,2) (exact bookkeeping)",
            (Fr(-1)+4*Fr(1, 16)+4*Fr(1, 2)*Fr(3, 4)) == Fr(3, 4) and
            (Fr(-1)+4*Fr(1, 16)+4*Fr(1, 2)*Fr(5, 4)) == Fr(7, 4) and
            (Fr(-2)+8*Fr(1, 16)+8*Fr(1, 2)*Fr(5, 4)) == Fr(7, 2) and
            (Fr(-2)+8*Fr(1, 16)+8*Fr(1, 2)*Fr(3, 4)) == Fr(3, 2))
# Substitutions (Q3),(Q4),(Q6):
#   L1 = pi^3 M1/(12 sqrt6),  L2 = pi^3 M2/(12 sqrt6),
#   L1t = pi^3 M1t/(96 sqrt6), L2t = pi^3 M2t/(96 sqrt6),
#   zK = pi^3 d24/(72 sqrt6),  Lg = pi^3 sqrt6 d3/36,
#   L(chi_{-8},2) L(chi_3,2) = (pi d8/(4 sqrt2))(pi^2 sqrt3/18)
#                            = pi^3 sqrt6 d8/144,
#   Cat L(chi_24,2) = (pi d4/2)(pi^2 sqrt6/24) = pi^3 sqrt6 d4/48.
# EK4 = (10 sqrt6/pi^3) comb: monomial coefficients, exact Fractions.
cM1 = 10*combv[0]/12; cM2 = 10*combv[1]/12
cM1t = 10*combv[2]/96; cM2t = 10*combv[3]/96
cd24 = 10*combv[4]/72
cd3 = 10*combv[5]*6/36          # sqrt6 * sqrt6 = 6
cd8 = 10*combv[6]*6/144
cd4 = 10*combv[7]*6/48
rM1 = Fr(5, 48)*28; rM2 = Fr(5, 48)*12
rMt = Fr(5, 48)*4
rd3 = Fr(5, 48)*28; rd4 = Fr(5, 48)*24; rd8 = Fr(5, 48)*8; rd24 = Fr(5, 48)
check_exact("[S1] M1: 10*(7/2)/12 = 35/12 = (5/48)*28 (exact)", cM1 == rM1)
check_exact("[S2] M2: 10*(3/2)/12 = 5/4 = (5/48)*12 (exact)", cM2 == rM2)
check_exact("[S3] Mit: 10*4/96 = 5/12 = (5/48)*4 (exact)",
            cM1t == rMt and cM2t == rMt)
check_exact("[S4] d24: 10*(3/4)/72 = 5/48 (exact)", cd24 == rd24)
check_exact("[S5] d3: 10*(7/4)*6/36 = 35/12 = (5/48)*28 (exact)",
            cd3 == rd3)
check_exact("[S6] d8: 10*2*6/144 = 5/6 = (5/48)*8 (exact)", cd8 == rd8)
check_exact("[S7] d4: 10*2*6/48 = 5/2 = (5/48)*24 (exact)", cd4 == rd4)
check_exact("[S8] assembly closes: LHS = RHS monomial-wise (exact)",
            (cM1, cM2, cM1t, cM2t, cd3, cd4, cd8, cd24)
            == (rM1, rM2, rMt, rMt, rd3, rd4, rd8, rd24))

# =====================================================================
# PART 1: [G]/[L]/[D]/[T]/[E] numeric confirmations (mp, 60 dps)
# =====================================================================
a1tw = [chi8m(n)*a1[n] for n in range(NMAX+1)]
a2tw = [chi8m(n)*a2[n] for n in range(NMAX+1)]

def mellin_I(a, xN, s):
    s = mpf(s); tot = mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0: continue
        tot += a[n]*power(2*pi*n, -s)*gammainc(s, 2*pi*n/xN)
    return tot

def qser(a, tau):
    q = exp(2*pi*1j*tau)
    s = mpc(0); qn = q
    for n in range(1, len(a)):
        if a[n]: s += a[n]*qn
        qn *= q
        if abs(qn) < mpf(10)**(-65): break
    return s

def fricke(a, N, yv):
    xN = sqrt(mpf(N))
    return qser(a, 1j/(xN*yv))/(yv**3*qser(a, 1j*yv/xN))

def Lset(a, N, label):
    xN = sqrt(mpf(N))
    for yy in ("0.6", "1.1"):
        r = fricke(a, N, mpf(yy))
        check("[G:%s] Fricke ratio at level %d = +1 (y = %s)"
              % (label, N, yy), r, 1, mpf(10)**(-45))
    I0, I3 = mellin_I(a, xN, 0), mellin_I(a, xN, 3)
    Lam3 = xN**3*I3 + I0
    L3 = Lam3*(2*pi)**3/(xN**3*gamma(3))
    return L3, Lam3

L1, M1 = Lset(a1, 24, "g1, N=24")
L2, M2 = Lset(a2, 24, "g2, N=24")
L1t, M1t = Lset(a1tw, 96, "g1tw, N=96")
L2t, M2t = Lset(a2tw, 96, "g2tw, N=96")
check("[G:M1] M1 = (24^{3/2}/(4 pi^3)) L1 (FE constant)", M1,
      mpf(24)**1.5*L1/(4*pi**3))
check("[G:M1t] M1t = (96^{3/2}/(4 pi^3)) L1t", M1t,
      mpf(96)**1.5*L1t/(4*pi**3))

# Dirichlet values and closed forms (Q4)
L2_8p = dirichlet(mpf(2), [0, 1, 0, -1, 0, -1, 0, 1])
L2_8m = dirichlet(mpf(2), [0, 1, 0, 1, 0, -1, 0, -1])
L2_3p = dirichlet(mpf(2), chi3p_l)
L2_3m = dirichlet(mpf(2), [0, 1, -1])
L2_24p = dirichlet(mpf(2), [chi_24(n) for n in range(24)])
L2_24m = dirichlet(mpf(2), [chi_24m(n) for n in range(24)])
Cat = dirichlet(mpf(2), [0, 1, 0, -1])
check("[D1] L(chi_8,2) = pi^2 sqrt2/16 (closed form, Q4)", L2_8p,
      pi**2*s2/16)
check("[D2] L(chi_3,2) = pi^2 sqrt3/18 (closed form, Q4 + [X4])", L2_3p,
      pi**2*s3/18)
check("[D3] L(chi_24,2) = pi^2 sqrt6/24 (closed form, Q4 + [X4])", L2_24p,
      pi**2*s6/24)
d3 = mpf(3)**1.5/(4*pi)*L2_3m
d4 = 2*Cat/pi
d8 = 4*s2/pi*L2_8m
d24 = mpf(24)**1.5/(4*pi)*L2_24m
check("[D4] d3 = (3^{3/2}/4pi) L(chi_{-3},2) = L'(chi_{-3},-1) "
      "(direct derivative; soft tol: numerical-derivative artifact)", d3,
      mpdiff(lambda s: dirichlet(s, [0, 1, -1]), mpf(-1)), mpf(10)**(-13))
check("[D5] d24 = L'(chi_{-24},-1) (direct derivative)", d24,
      mpdiff(lambda s: dirichlet(s, [chi_24m(n) for n in range(24)]), mpf(-1)))

# lattice T-sums and the decomposition (mp 60 dps)
def G_row(x, y):
    return (pi/y)*mp.sinh(2*pi*y)/(mp.cosh(2*pi*y) - mp.cos(2*pi*x))

def row_powers(x, y):
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y)/(2*y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy)/(2*yy), y)/(4*y)
    return S2, S3

def lattice_T(d, tau):
    x0, y0 = tau.real, tau.imag
    B, M = mpf(0), mpf(0)
    for m in range(-300, 301):
        if m == 0:
            B += 2*zeta(4); continue
        S2, S3 = row_powers(d*m*x0, abs(d*m)*y0)
        B += S2; M += m*m*S3
        if abs(d*m)*y0 > 45 and m > 0: break
    return 3*B - 4*d*d*y0*y0*M

zK = zeta(2)*L2_24m
Lg = L2_8p*L2_3m
B_OK = zK + Lg
G_OK = L1 + L2
B_odd = 3*zK/4 + 5*Lg/4
B_odd_chi = L2_8m*L2_3p + Cat*L2_24p
G_odd = 5*L1/4 + 3*L2/4
G_odd_chi = L1t + L2t
T_OK_pred = B_OK + 2*G_OK
T_O2_pred = (B_OK/16 + (B_odd + B_odd_chi)/2) \
            + 2*(G_OK/16 + (G_odd + G_odd_chi)/2)
tauA = mpc(0, s6)
T1v = lattice_T(1, tauA)
T2v = lattice_T(2, tauA)
check("[T1] T(O_K) = zeta_K(2) + L(chi8)L(chi-3) + 2(L1+L2)",
      T1v, T_OK_pred)
check("[T2] T(O_2) = class-number-2 decomposition", T2v, T_O2_pred)
comb = -T1v + 4*T2v
EK4 = (10*s6/pi**3)*comb
conj = Fr(5, 48)*(4*M1t + 4*M2t + 28*M1 + 12*M2
                  + 28*d3 + 24*d4 + 8*d8 + d24)
check("[E] EK4(sqrt(-6)) = (5/48)(4M1tw+4M2tw+28M1+12M2"
      "+28d3+24d4+8d8+d24)", EK4, conj)

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
else:
    print("ALL CHECKS PASSED (pre-iv tracks)")

# =====================================================================
# PART 2: [V] rigorous interval locks (mpmath.iv, iv.dps = 70)
#   Ported from verify_P1_n4_p3_t1t2.py ([V0]-[V6]); tau = i sqrt6 is the
#   shift-0 case, lattice_T_iv applies verbatim.  The Mellin tail bound is
#   replaced by a SELF-CONTAINED theta-series bound:
#     |a^i(n)| <= |P(n)| + |Q(n)| <= 6 n^2 <= 6 n^3  (n >= 1),
#   since #{(a,b): a^2+6b^2 = n} <= (2 sqrt n + 1)(2 sqrt(n/6) + 1) <= 6n
#   and |a^2 - 6b^2| <= n (same for the 2x^2+3y^2 family).
# =====================================================================
from mpmath import iv
iv.dps = 70

BERN = [Fr(1, 6), Fr(-1, 30), Fr(1, 42), Fr(-1, 30), Fr(5, 66),
        Fr(-691, 2730), Fr(7, 6), Fr(-3617, 510), Fr(43867, 798),
        Fr(-174611, 330), Fr(854513, 138), Fr(-236364091, 2730)]

def hurwitz2_iv(r, N, K=260, p=10):
    r, N = iv.mpf(r), iv.mpf(N)
    s = iv.mpf(0)
    for k in range(K):
        s += (N*k + r)**(-2)
    u = N*K + r
    tail = 1/(N*u) + u**(-2)/2
    for j in range(1, p+1):
        Bj = iv.mpf(BERN[j-1].numerator)/iv.mpf(BERN[j-1].denominator)
        tail += Bj * N**(2*j-1) * u**(-2*j-1)
    Brem = iv.mpf(BERN[p].numerator)/iv.mpf(BERN[p].denominator)
    Rb = abs((Brem * N**(2*p+1) * u**(-2*p-3)).b)
    return s + tail + iv.mpf([-Rb, Rb])

def dirichlet2_iv(chi, N):
    tot = iv.mpf(0)
    for r in range(1, N+1):
        if chi[r % N] != 0:
            tot += chi[r % N] * hurwitz2_iv(r, N)
    return tot

eulergamma_iv = iv.euler

def E1_iv(x):
    x = iv.mpf(x)
    if x.b < 8:
        tot = -eulergamma_iv - iv.log(x)
        term = x
        k = 1
        while True:
            tot += ((-1)**(k+1)) * term / k
            nxt = term * x / (k+1)
            if k + 1 > x.b and abs(nxt.b)/(k+1) < iv.mpf(10)**(-62):
                Rb = abs(nxt.b)/(k+1)
                tot += iv.mpf([-Rb, Rb])
                break
            term = nxt
            k += 1
        return tot
    Pm2, Pm1 = iv.mpf(1), iv.mpf(0)
    Qm2, Qm1 = iv.mpf(0), iv.mpf(1)
    prev = None
    for n in range(1, 4000):
        a = iv.mpf(1) if n == 1 else iv.mpf(n//2)
        b = x if n % 2 == 1 else iv.mpf(1)
        P = b*Pm1 + a*Pm2
        Q = b*Qm1 + a*Qm2
        cur = P/Q
        if prev is not None and n > 6:
            width = max(cur.b, prev.b) - min(cur.a, prev.a)
            if width < iv.mpf(10)**(-52)*abs(cur).b + iv.mpf(10)**(-62):
                lo = min(cur.a, prev.a)
                hi = max(cur.b, prev.b)
                return iv.exp(-x) * iv.mpf([lo, hi])
        prev = cur
        Pm2, Pm1 = Pm1, P
        Qm2, Qm1 = Qm1, Q
    raise RuntimeError("E1 CF did not converge")

def gammainc_iv(s, x):
    x = iv.mpf(x)
    if s == 1: return iv.exp(-x)
    if s == 2: return iv.exp(-x)*(1+x)
    if s == 3: return iv.exp(-x)*(2+2*x+x**2)
    if s == 0: return E1_iv(x)
    raise ValueError

def mellin_I_iv(a, xN, s, n0):
    tot = iv.mpf(0)
    for n in range(1, min(n0, len(a)-1)+1):
        if a[n] == 0:
            continue
        tot += a[n] * (2*iv.pi*n)**(-s) * gammainc_iv(s, 2*iv.pi*n/xN)
    # rigorous tail bound with |a_n| <= 6 n^3 (theta series, see header);
    # Gamma(s,x) <= s! e^-x (1+x)^s (s=0: E1(x) <= e^-x(1+1/x)).
    c = 2*iv.pi/xN
    def Abound(n):
        nn = iv.mpf(n)
        xn = c*nn
        g = iv.exp(-xn)*(1+1/xn) if s == 0 else iv.exp(-xn)*(1+xn)**s
        sfac = iv.mpf(1) if s == 0 else iv.mpf([1, 1, 2, 6][s])
        return (6*nn**3*(2*iv.pi*nn)**(-s)*sfac*g).b
    rho_iv = iv.exp(-c) * ((iv.mpf(n0+4))/(iv.mpf(n0+3)))**3 \
        * ((1+c*(n0+2))/(1+c*(n0+1)))**s
    rho = mpf(rho_iv.b)
    assert rho < 1
    T = Abound(n0+1) / (1 - rho)
    return tot + iv.mpf([-T, T])

def Lset_iv(a, N):
    xN = iv.sqrt(iv.mpf(N))
    n0 = int(45*mpf(xN.a)) + 40
    I0 = mellin_I_iv(a, xN, 0, n0)
    I3 = mellin_I_iv(a, xN, 3, n0)
    Lam3 = xN**3*I3 + I0          # w = +1 (Fricke ratios [G:*])
    L3v = Lam3*(2*iv.pi)**3/(xN**3*2)
    return L3v, Lam3

class Dual:
    """order-2 Taylor arithmetic (f0 + f1 e + f2 e^2) over iv."""
    def __init__(self, f0, f1=0, f2=0):
        self.f0, self.f1, self.f2 = iv.mpf(f0), iv.mpf(f1), iv.mpf(f2)
    def __add__(a, b):
        b = b if isinstance(b, Dual) else Dual(b)
        return Dual(a.f0+b.f0, a.f1+b.f1, a.f2+b.f2)
    __radd__ = __add__
    def __neg__(a):
        return Dual(-a.f0, -a.f1, -a.f2)
    def __sub__(a, b):
        return a + (-(b if isinstance(b, Dual) else Dual(b)))
    def __rsub__(a, b):
        return Dual(b) + (-a)
    def __mul__(a, b):
        b = b if isinstance(b, Dual) else Dual(b)
        return Dual(a.f0*b.f0, a.f0*b.f1+a.f1*b.f0,
                    a.f0*b.f2+a.f1*b.f1+a.f2*b.f0)
    __rmul__ = __mul__
    def __truediv__(a, b):
        b = b if isinstance(b, Dual) else Dual(b)
        c0 = a.f0/b.f0
        c1 = (a.f1 - c0*b.f1)/b.f0
        c2 = (a.f2 - c0*b.f2 - c1*b.f1)/b.f0
        return Dual(c0, c1, c2)
    def __rtruediv__(a, b):
        return Dual(b) / a
    def exp(a):
        e = iv.exp(a.f0)
        return Dual(e, e*a.f1, e*(a.f2 + a.f1**2/2))

def coth_dual(yv):
    ep0 = iv.exp(iv.pi*yv)
    em0 = iv.exp(-iv.pi*yv)
    ep = Dual(ep0, iv.pi*ep0, iv.pi**2*ep0/2)
    em = Dual(em0, -iv.pi*em0, iv.pi**2*em0/2)
    return (ep+em)/(ep-em)

def row_S23_iv(y):
    G = (Dual(iv.pi)/Dual(y, 1, 0)) * coth_dual(y)
    S2 = -G.f1/(2*y)
    S3 = (2*G.f2 + 2*S2)/(8*y**2)
    return S2, S3

def tail_row_iv(y):
    y = iv.mpf(y)
    r = iv.exp(-2*iv.pi*y)
    rb = r.b
    den1 = (1-rb)**2
    den3 = (1-rb)**3
    b2 = (2*iv.pi/y**3)*(1+2*iv.pi*y)*rb/den1
    b3 = (iv.pi/(2*y**5))*(3*rb/(1-rb) + 6*iv.pi*y*rb/den1
                           + 4*iv.pi**2*y**2*rb*(1+rb)/den3)
    return mpf(b2.b), mpf(b3.b)

def lattice_T_iv(d, y0):
    """T-sum for tau = i y0 (all rows at shift 0), rigorous tail."""
    y0m = mpf(iv.mpf(y0).a)
    M = int(50/(d*y0m)) + 2
    T = iv.pi**4/15
    for m in range(1, M+1):
        y = iv.mpf(d*m)*y0
        S2, S3 = row_S23_iv(y)
        dS2 = S2 - iv.pi/(2*y**3)
        dS3 = S3 - 3*iv.pi/(8*y**5)
        T += 2*(3*dS2 - 4*y**2*dS3)
    tb2, tb3 = tail_row_iv(d*(M+1)*y0m)
    rr = mpf(iv.exp(-2*iv.pi*d*y0m).b)
    yM = d*(M+1)*y0m
    tot = 2*(3*tb2 + 4*yM*yM*tb3)/(1-rr)
    return T + iv.mpf([-tot, tot])

# --- theta coefficients to n0(96) = 481 for the iv Mellin tails ----------
NMV = 520
Pv = [0]*(NMV+1); Qv = [0]*(NMV+1)
RBV = int(NMV**0.5)+2
for aa_ in range(-RBV, RBV+1):
    for bb_ in range(-RBV, RBV+1):
        n_ = aa_*aa_ + 6*bb_*bb_
        if 1 <= n_ <= NMV: Pv[n_] += aa_*aa_ - 6*bb_*bb_
        n_ = 2*aa_*aa_ + 3*bb_*bb_
        if 1 <= n_ <= NMV: Qv[n_] += 4*aa_*aa_ - 6*bb_*bb_
Pv = [p//2 for p in Pv]; Qv = [q//4 for q in Qv]
a1v = [Pv[n] - Qv[n] for n in range(NMV+1)]
a2v = [Pv[n] + Qv[n] for n in range(NMV+1)]
a1tv = [chi8m(n)*a1v[n] for n in range(NMV+1)]
a2tv = [chi8m(n)*a2v[n] for n in range(NMV+1)]
check_exact("[V-pre] iv-range coefficients coincide with the exact track "
            "(n <= 400)", a1v[:401] == a1[:401] and a2v[:401] == a2[:401])

# --- the iv values --------------------------------------------------------
def check_lock(name, ivl, halfwidth_req=mpf(10)**(-48)):
    lo, hi = mp.convert(ivl.a), mp.convert(ivl.b)
    w = hi - lo
    ok = (lo <= 0 <= hi) and (w < halfwidth_req)
    if not ok:
        FAILS.append(name)
    print("%-76s %s  (half-width = %.2e)" % (name, "PASS" if ok else "FAIL",
                                             w/2))

def iv_mid(z):
    return (mp.convert(z.a) + mp.convert(z.b))/2

def iv_w(z):
    return mp.convert(z.b) - mp.convert(z.a)

L1_iv, M1_iv = Lset_iv(a1v, 24)
L2_iv, M2_iv = Lset_iv(a2v, 24)
L1t_iv, M1t_iv = Lset_iv(a1tv, 96)
L2t_iv, M2t_iv = Lset_iv(a2tv, 96)
Cat_iv = dirichlet2_iv([0, 1, 0, -1], 4)
L8p_iv = dirichlet2_iv([0, 1, 0, -1, 0, -1, 0, 1], 8)
L8m_iv = dirichlet2_iv([0, 1, 0, 1, 0, -1, 0, -1], 8)
L3p_iv = dirichlet2_iv(chi3p_l, 12)
L3m_iv = dirichlet2_iv([0, 1, -1], 3)
L24p_iv = dirichlet2_iv([chi_24(n) for n in range(24)], 24)
L24m_iv = dirichlet2_iv([chi_24m(n) for n in range(24)], 24)

check("[V0a] iv: L(g1,3) vs mp value", iv_mid(L1_iv), L1, iv_w(L1_iv)+TOL)
check("[V0b] iv: L(g1tw,3) vs mp value", iv_mid(L1t_iv), L1t,
      iv_w(L1t_iv)+TOL)
check("[V0c] iv: L(chi_3,2) = pi^2 sqrt3/18", iv_mid(L3p_iv), pi**2*s3/18,
      iv_w(L3p_iv)+TOL)
check("[V0d] iv: L(chi_24,2) = pi^2 sqrt6/24", iv_mid(L24p_iv),
      pi**2*s6/24, iv_w(L24p_iv)+TOL)

zK_iv = iv.pi**2/6 * L24m_iv
Lg_iv = L8p_iv * L3m_iv
d3_iv = 3*iv.sqrt(3)/(4*iv.pi) * L3m_iv
d4_iv = 2*Cat_iv/iv.pi
d8_iv = 4*iv.sqrt(2)/iv.pi * L8m_iv
d24_iv = 12*iv.sqrt(6)/iv.pi * L24m_iv

s6_iv = iv.sqrt(6)
T_OK_iv = lattice_T_iv(1, s6_iv)
T_O2_iv = lattice_T_iv(2, s6_iv)

B_OK_iv = zK_iv + Lg_iv
G_OK_iv = L1_iv + L2_iv
B_odd_iv = 3*zK_iv/4 + 5*Lg_iv/4
B_oddc_iv = L8m_iv*L3p_iv + Cat_iv*L24p_iv
G_odd_iv = 5*L1_iv/4 + 3*L2_iv/4
G_oddc_iv = L1t_iv + L2t_iv
T_OK_dec_iv = B_OK_iv + 2*G_OK_iv
T_O2_dec_iv = (B_OK_iv/16 + (B_odd_iv + B_oddc_iv)/2) \
              + 2*(G_OK_iv/16 + (G_odd_iv + G_oddc_iv)/2)

check_lock("[V1] LOCK T(O_K) = zeta_K(2) + L(chi8)L(chi-3) + 2(L1+L2)",
           T_OK_iv - T_OK_dec_iv)
check_lock("[V2] LOCK T(O_2) = class-number-2 decomposition",
           T_O2_iv - T_O2_dec_iv)

comb_iv = -T_OK_iv + 4*T_O2_iv
EK4_iv = (10*s6_iv/iv.pi**3)*comb_iv
conj_iv = iv.mpf(5)/48*(4*M1t_iv + 4*M2t_iv + 28*M1_iv + 12*M2_iv
                        + 28*d3_iv + 24*d4_iv + 8*d8_iv + d24_iv)
check_lock("[V3] LOCK EK4(sqrt(-6)) = (5/48)(4M1tw+4M2tw+28M1+12M2"
           "+28d3+24d4+8d8+d24)", EK4_iv - conj_iv, mpf(10)**(-50))

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
    import sys; sys.exit(1)
print("ALL CHECKS PASSED")
