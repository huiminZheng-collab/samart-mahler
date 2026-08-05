# verify_P1_n5_e78.py -- (P1)-type CM evaluation for the Samart Table-6
# OPEN conjectures #7 and #8 (Phase 8, targets #7/#8; FIFTH field, the
# FIRST with class number 4):
#
#   #7: tau7 = (3 + sqrt(-21))/6 = 1/2 + i sqrt21/6  (y ~ 0.7638),
#       s4(tau7) = -893952 + 516096 sqrt3,
#       EK4(tau7) = (20/7)(M84^3 - M84^4 + 8 d3 - 4 d4)
#   #8: tau8 = (1 + sqrt(-21))/2 = 1/2 + i sqrt21/2  (y ~ 2.2913),
#       s4(tau8) = -893952 - 516096 sqrt3,
#       EK4(tau8) = (20/21)(M84^3 + M84^4 + 8 d3 + 4 d4)
#   [numeric locks 6e-61 / 1.2e-60 in n5_p8_e78_fit.py; s4 mp locks below,
#    s4 iv locks in cert0_n5_e78.py (parallel work)]
#
# K = Q(sqrt-21), O_K = Z[sqrt-21] (-21 == 3 (4)), disc -84, h = 4, class
# group = genus group = (Z/2)^2: classes 1, [p_2], [p_3], [p_6] with norm
# forms x^2+21y^2, 2x^2+2xy+11y^2, 3x^2+7y^2, 6x^2+6xy+5y^2; p_2^2 = (2),
# p_3^2 = (3), p_6 = p_2 p_3, p_6^2 = (6).  The four newforms of
# S_3(Gamma_0(84), chi_{-84}) are the CM theta series
#   g(e2,e3) = P0 + e2 P2 + e3 P3 + e2 e3 P6      (e2, e3 in {+-1}),
# with class theta sums P0,P2,P3,P6 (grossencharacter values REAL:
# psi(p_2) = +-2, psi(p_3) = +-3, psi(p_6) = psi(p_2)psi(p_3) = +-6;
# convention psi(p_2) = +2, psi(p_3) = +3 for g(1,1), so chi_{(e2,e3)}:
# p_2 |-> e2, p_3 |-> e3).  M84^3 = L'(g(-1,+1),0), M84^4 = L'(g(-1,-1),0)
# (Samart's #3/#4 of the level-84 block), d_k = L'(chi_{-k},-1).
#
# Anchors (fitted < 1e-40 in n5_p8_e78_fit.py; locked in [V1]-[V4] here;
# zK = zeta(2)L(chi_{-84},2), A4 = L(chi_{-4},2)L(chi_21,2),
# A3 = L(chi_{-3},2)L(chi_28,2), A7 = L(chi_{-7},2)L(chi_12,2),
# L1..L4 = L(g(1,1),3), L(g(1,-1),3), L(g(-1,1),3), L(g(-1,-1),3);
# class eigenvalues on (A4,A3,A7): p2: (-,-,+), p3: (-,+,-), p6: (+,-,-)):
#   T(O_K) = (1/2)(zK+A4+A3+A7) + (L1+L2+L3+L4)
#   T(p_2) = (1/8)(zK-A4-A3+A7) + (1/4)(L1+L2-L3-L4)
#   T(p_3) = (1/18)(zK-A4+A3-A7) + (1/9)(L1-L2+L3-L4)
#   T(p_6) = (1/72)(zK+A4-A3-A7) + (1/36)(L1-L2-L3+L4)
#
# Lattice relations (h = 4 novelty: the EK4 lattices are class lattices,
# NOT orders; no Euler-factor removal is needed):
#   #8: T(Lambda_1) = 16 T(p_2)   (tau = 1/2 + i sqrt21/2),
#       T(Lambda_2) = T(O_K)      (Z + 2 tau Z = Z[sqrt-21]);
#   #7: T(Lambda_1) = 1296 T(p_6) (tau = 1/2 + i sqrt21/6),
#       T(Lambda_2) = 81 T(p_3)   (Z + 2 tau Z = Z + Z sqrt-21/3).
#   EK4(tau) = (10 y0/pi^3)(-T(Lambda_1) + 4 T(Lambda_2)).
# Assembly (zK and A7 cancel in both combinations):
#   #7: comb = -1296 T(p_6) + 324 T(p_3) = 36(A3-A4) + 72(L3-L4),
#       EK4 = (5 sqrt21/(3 pi^3)) comb
#   #8: comb = -16 T(p_2) + 4 T(O_K)   = 4(A3+A4) + 8(L3+L4),
#       EK4 = (5 sqrt21/pi^3) comb
# Closure uses the Fricke relation (w = +1, level 84)
#   M = L'(g,0) = (84^{3/2}/(4 pi^3)) L(g,3) = (42 sqrt21/pi^3) L(g,3)
# and the even-character Gauss-sum closed forms (exact integer sums [X5])
#   L(chi_28,2) = 2 sqrt7 pi^2/49,  L(chi_21,2) = 8 sqrt21 pi^2/441,
# giving A3 = (8 sqrt21 pi^3/441) d3, A4 = (4 sqrt21 pi^3/441) d4.
#
# Layers (same three-track standard as verify_P1_n5_e6.py):
#   [X*]/[S*] exact integer / Fraction checks (no floating point);
#   [G*] newform identifications: exact-integer Hecke/CM checks; root
#        number by the non-vacuous Fricke ratio (numeric, flagged);
#   [L*],[D*],[T*],[E*],[S4*] mpmath 60-dps numerical confirmations;
#   [V*] rigorous interval locks (iv.dps = 70, self-contained theta tail
#        |a(n)| <= 30 n^3 -- see bound derivation in PART 2):
#        [V1]-[V4] the four class-lattice anchors, [V5]/[V6] the two EK4
#        identities, all locked with half-width < 1e-48.
#
# Quoted inputs (same standard as e6):
#   (Q0) n4(s4(tau)) = EK4(tau) at tau7, tau8: both lie on the certified
#        vertical leg x = 1/2, y in [0.702, sqrt21/2] of n5_line_cert.py
#        (#8's endpoint is certified interior; #7 also lies on leg B of
#        n4_p4_t4_cert.py); the propagation theorem then gives the
#        identity at both CM points.
#   (Q1) Hecke's theorem: theta series of conductor-1 grossencharacters
#        are newforms; S_3(Gamma_0(84), chi_{-84}) is 4-dimensional
#        (LMFDB: the four dim-1 CM forms 84.3.d.*, cm_discs = [-84]),
#        so {g(e2,e3)} exhausts it.
#   (Q2) h(-84) = 4, class group = genus group = (Z/2)^2, units +-1;
#        2, 3, 7 ramify.
#   (Q3) Functional equations of the CM forms (level 84, w = +1).
#   (Q4) Dirichlet class-number-free formulas d_k = (k^{3/2}/(4 pi))
#        L(chi_{-k},2) (k = 3, 4); even-character Gauss formula
#        L(chi,2) = (pi^2/d^{5/2}) sum_{a=1}^d chi(a) a^2 for even
#        primitive real chi (chi_28, chi_21; Gauss sums +sqrt28, +sqrt21).
#   (Q5) Genus theory: L_K(chi o N, s) = L(chi,s) L(chi chi_{-84},s) for
#        chi = chi_{-4}, chi_{-3}, chi_{-7} (products [X4] exact).
#   (Q6) zeta(2) = pi^2/6; zeta_K(2) = zeta(2) L(chi_{-84},2).
#   (Q7) s4 values: s4(tau7), s4(tau8) are the roots of
#        X^2 + 1787904 X + (893952^2 - 3*516096^2) [S7 exact; S4-1/2 mp;
#        iv locks in cert0_n5_e78.py]; q7 = -e^{-pi sqrt21/3},
#        q8 = -e^{-pi sqrt21} negative real [S4-3/4].

from fractions import Fraction as Fr
from mpmath import (mp, mpf, mpc, pi, sqrt, zeta, dirichlet, gamma,
                    diff as mpdiff, power, gammainc, sinh, cosh, cos, exp)

mp.dps = 60
s3 = sqrt(mpf(3)); s7 = sqrt(mpf(7)); s21 = sqrt(mpf(21))
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
okX0 = all(4*x*x - (x*x+y*y) == 2*(x*x-y*y) + (x*x+y*y) == 3*x*x - y*y
           for x in range(-5, 6) for y in range(-5, 6))
check_exact("[X0] F(z) = 2 Re(z^2)/|z|^6 + 1/|z|^4 pointwise (grid, exact) "
            "=> T(L) = B(L) + 2 G(L)", okX0)

# --- ideal arithmetic (pair form alpha = (c, d) = c + d sqrt-21) ---------
def k21norm(u): return u[0]**2 + 21*u[1]**2
S21P = (0, 1)          # sqrt-21
# p_2 = (2, 1+sqrt-21): p_2^2 = (2): generators 4, 2(1+s), (1+s)^2 =
# -20 + 2s, and 2 = 6*4 - 2(1+s) + (1+s)^2 (integer combination).
check_exact("[X1] p_2^2 = (2): (1+s)^2 = -20+2s in (2) and "
            "2 = 6*4 - 2(1+s) + (1+s)^2 (exact); N(p_2)^2 = 4 = N((2))",
            (1 - 21, 2) == (-20, 2)
            and 24 - 2 - 20 == 2 and (-2) + 2 == 0
            and 2**2 == k21norm((2, 0)))
# p_3 = (3, sqrt-21): p_3^2 = (3): generators 9, 3s, -21; 3 = -(-21) - 2*9.
check_exact("[X2] p_3^2 = (3): 3 = 21 - 2*9 (exact); N(p_3)^2 = 9 = N((3))",
            21 - 18 == 3 and 3**2 == k21norm((3, 0)))
# Class group (Z/2)^2: p_6 = p_2 p_3, p_6^2 = (6); genus eigenvalues on
# (A4,A3,A7): p2: (-,-,+), p3: (-,+,-), p6: (+,-,-) = product.
eps2, eps3, eps6 = (-1, -1, +1), (-1, +1, -1), (+1, -1, -1)
check_exact("[X3] p_6 = p_2 p_3: eigenvalue multiplicativity eps(p6) = "
            "eps(p2)*eps(p3); p_6^2 = (6): N = 36 = N((6)) (exact)",
            all(eps6[i] == eps2[i]*eps3[i] for i in range(3))
            and 6**2 == k21norm((6, 0)))

# --- characters (tables copied verbatim from n5_p8_e78_fit.py) -----------
def jacobi(a, n):
    a %= n; t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: t = -t
        a %= n
    return t if n == 1 else 0

def chi_84(n):
    if n % 2 == 0 or n % 3 == 0 or n % 7 == 0: return 0
    return (1 if n % 4 == 1 else -1) * jacobi(84, n)
chi_4m = lambda n: 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
chi_21 = lambda n: 0 if n % 3 == 0 or n % 7 == 0 else jacobi(n, 21)
chi_3m = lambda n: 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)
chi_28 = lambda n: 0 if n % 2 == 0 or n % 7 == 0 else chi_4m(n)*jacobi(n, 7)
chi_7m = lambda n: 0 if n % 7 == 0 else jacobi(n, 7)
chi_12 = lambda n: 0 if n % 2 == 0 or n % 3 == 0 else chi_4m(n)*jacobi(n, 3)

okX4 = all((chi_4m(n)*chi_84(n) == chi_21(n)
            and chi_3m(n)*chi_84(n) == chi_28(n)
            and chi_7m(n)*chi_84(n) == chi_12(n))
           for n in range(1, 600) if n % 2 and n % 3 and n % 7)
check_exact("[X4] genus characters: chi_4m*chi_84 = chi_21, chi_3m*chi_84 "
            "= chi_28, chi_7m*chi_84 = chi_12 (n coprime 84, n<600, exact)",
            okX4)

S28 = (sum(chi_28(a) for a in range(1, 29)),
       sum(chi_28(a)*a for a in range(1, 29)),
       sum(chi_28(a)*a*a for a in range(1, 29)))
S21 = (sum(chi_21(a) for a in range(1, 22)),
       sum(chi_21(a)*a for a in range(1, 22)),
       sum(chi_21(a)*a*a for a in range(1, 22)))
check_exact("[X5a] Gauss sums: sum chi_28(a) a^k = (0,0,448); sum chi_21(a) "
            "a^k = (0,0,168) (k=0,1,2; exact)", S28 == (0, 0, 448)
            and S21 == (0, 0, 168))
# L(chi_d,2) = pi^2 S_d / d^{5/2} (Q4): closed forms reduce to integers:
# L(chi_28,2) = 2 sqrt7 pi^2/49  <=>  448*49 = 2*7*(2*28^2)  [sqrt28=2 sqrt7]
# L(chi_21,2) = 8 sqrt21 pi^2/441 <=> 168*441 = 8*21^3
check_exact("[X5b] closed-form algebra: 448*49 = 4*7*28^2 and 168*441 = "
            "8*21^3 (exact)", 448*49 == 4*7*28**2 and 168*441 == 8*21**3)

# --- genus characters constant on classes; eigenvalue table --------------
FORMS = ((lambda x, y: x*x + 21*y*y, (+1, +1, +1)),
         (lambda x, y: 2*x*x + 2*x*y + 11*y*y, eps2),
         (lambda x, y: 3*x*x + 7*y*y, eps3),
         (lambda x, y: 6*x*x + 6*x*y + 5*y*y, eps6))
okX6 = True
seen = [0]*4
for i, (f, eps) in enumerate(FORMS):
    for x in range(-9, 10):
        for y in range(-9, 10):
            n = f(x, y)
            if n >= 1 and n % 2 and n % 3 and n % 7:
                seen[i] += 1
                okX6 &= (chi_4m(n), chi_3m(n), chi_7m(n)) == eps
check_exact("[X6] genus chars on classes (box): (chi_4m,chi_3m,chi_7m) = "
            "(+,+,+),(-,-,+),(-,+,-),(+,-,-) (all forms hit: %s; exact)"
            % seen, okX6 and all(v > 0 for v in seen))

# --- ideal counts vs representation numbers (exact integers) -------------
NID = 2000
r = [0]*(NID+1)
g4 = [0]*(NID+1); g3 = [0]*(NID+1); g7 = [0]*(NID+1)
for dd in range(1, NID+1):
    cc = chi_84(dd)
    if cc:
        for m in range(dd, NID+1, dd): r[m] += cc
    c4 = chi_4m(dd)
    if c4:
        for m in range(1, NID//dd+1):
            cm = chi_21(m)
            if cm: g4[dd*m] += c4*cm
    c3 = chi_3m(dd)
    if c3:
        for m in range(1, NID//dd+1):
            cm = chi_28(m)
            if cm: g3[dd*m] += c3*cm
    c7 = chi_7m(dd)
    if c7:
        for m in range(1, NID//dd+1):
            cm = chi_12(m)
            if cm: g7[dd*m] += c7*cm
r1 = [0]*(NID+1); r2 = [0]*(NID+1); r3 = [0]*(NID+1); r6 = [0]*(NID+1)
RB = int(NID**0.5)+3
for x in range(-RB, RB+1):
    for y in range(-RB, RB+1):
        for f, rr in ((FORMS[0][0], r1), (FORMS[1][0], r2),
                      (FORMS[2][0], r3), (FORMS[3][0], r6)):
            n = f(x, y)
            if 1 <= n <= NID: rr[n] += 1
okX7 = all((r1[n]+r2[n]+r3[n]+r6[n]) == 2*r[n]
           and (r1[n]-r2[n]-r3[n]+r6[n]) == 2*g4[n]
           and (r1[n]-r2[n]+r3[n]-r6[n]) == 2*g3[n]
           and (r1[n]+r2[n]-r3[n]-r6[n]) == 2*g7[n]
           for n in range(1, NID+1))
check_exact("[X7] ideal counts: r1+r2+r3+r6 = 2r(n); genus splits = "
            "2g4,2g3,2g7 (n<=2000; exact)", okX7)

# --- theta series (exact Fraction arithmetic) -----------------------------
NMAX = 400
RB2 = int((2*NMAX)**0.5)+3
P0 = [Fr(0)]*(NMAX+1); P2 = [Fr(0)]*(NMAX+1)
P3 = [Fr(0)]*(NMAX+1); P6 = [Fr(0)]*(NMAX+1)
for a in range(-RB2, RB2+1):
    for b in range(-RB2, RB2+1):
        n = a*a + 21*b*b
        if 1 <= n <= NMAX: P0[n] += Fr(a*a - 21*b*b, 2)
        n = 3*a*a + 7*b*b
        if 1 <= n <= NMAX: P3[n] += Fr(3*a*a - 7*b*b, 2)
for x in range(-RB2, RB2+1):
    for y in range(-RB2, RB2+1):
        n = 2*x*x + 2*x*y + 11*y*y
        if 1 <= n <= NMAX:
            v = (2*x+y)**2 - 21*y*y
            assert v % 4 == 0; P2[n] += Fr(v, 4)
        n = 6*x*x + 6*x*y + 5*y*y
        if 1 <= n <= NMAX:
            v = (6*x+3*y)**2 - 21*y*y
            assert v % 12 == 0; P6[n] += Fr(v, 12)
okX8 = (all(v.denominator == 1 for P in (P0, P2, P3, P6) for v in P[1:])
        and all(P[1] == 0 for P in (P2, P3, P6)) and P0[1] == 1)
check_exact("[X8] P0,P2,P3,P6 integral; a(1) = 1 for all four forms "
            "(exact Fractions)", okX8)
aE = {}
for e2 in (1, -1):
    for e3 in (1, -1):
        aE[(e2, e3)] = [int(P0[n] + e2*P2[n] + e3*P3[n] + e2*e3*P6[n])
                        for n in range(NMAX+1)]
check_exact("[X9] sign convention psi(p2)=+2, psi(p3)=+3 for g(1,1): "
            "a_(1,1)(2)=+2, a_(1,1)(3)=+3, a_(-1,1)(2)=-2, a_(1,-1)(3)=-3 "
            "(exact)", aE[(1, 1)][2] == 2 and aE[(1, 1)][3] == 3
            and aE[(-1, 1)][2] == -2 and aE[(1, -1)][3] == -3)

PR = []
for n in range(2, NMAX+1):
    if all(n % p for p in PR if p*p <= n): PR.append(n)
oke1, oke2, okcm = True, True, True
for aa in aE.values():
    for p in PR:
        if p*p > NMAX: break
        if aa[p*p] != aa[p]*aa[p] - chi_84(p)*p*p: oke1 = False
    for i, p in enumerate(PR):
        for q in PR[i+1:]:
            if p*q > NMAX: break
            if aa[p*q] != aa[p]*aa[q]: oke2 = False
    for p in PR:
        if p > 150: break
        if chi_84(p) == -1 and aa[p] != 0: okcm = False
check_exact("[X10] Hecke recursion a(p^2)=a(p)^2-chi_{-84}(p)p^2; "
            "multiplicativity; CM vanishing at inert p<=150 (exact)",
            oke1 and oke2 and okcm)

# --- [S] separation track (basis zK, A4, A3, A7, L1, L2, L3, L4) ---------
T_OK_c = (Fr(1, 2), Fr(1, 2), Fr(1, 2), Fr(1, 2), Fr(1), Fr(1), Fr(1), Fr(1))
T_p2_c = (Fr(1, 8), Fr(-1, 8), Fr(-1, 8), Fr(1, 8),
          Fr(1, 4), Fr(1, 4), Fr(-1, 4), Fr(-1, 4))
T_p3_c = (Fr(1, 18), Fr(-1, 18), Fr(1, 18), Fr(-1, 18),
          Fr(1, 9), Fr(-1, 9), Fr(1, 9), Fr(-1, 9))
T_p6_c = (Fr(1, 72), Fr(1, 72), Fr(-1, 72), Fr(-1, 72),
          Fr(1, 36), Fr(-1, 36), Fr(-1, 36), Fr(1, 36))
# internal consistency: B(C) = (1/(2 N_C^2)) sum eps(C) A_i;
# G(C) = (1/N_C^2) sum chi(C) L_i with chi(p2)=(+1,+1,-1,-1) etc.
check_exact("[S0] anchor tuples: B parts = (1/2N^2)(zK +- A4 +- A3 +- A7), "
            "G parts = (1/N^2)(+-L1+-L2+-L3+-L4) (exact bookkeeping)",
            all(T_p2_c[i] == Fr(eps2[i-1], 8) if 0 < i < 4 else True
                for i in range(8))
            and all(T_p3_c[i] == Fr(eps3[i-1], 18) if 0 < i < 4 else True
                    for i in range(8))
            and all(T_p6_c[i] == Fr(eps6[i-1], 72) if 0 < i < 4 else True
                    for i in range(8))
            and (T_p2_c[4:] == (Fr(1, 4), Fr(1, 4), Fr(-1, 4), Fr(-1, 4)))
            and (T_p3_c[4:] == (Fr(1, 9), Fr(-1, 9), Fr(1, 9), Fr(-1, 9)))
            and (T_p6_c[4:] == (Fr(1, 36), Fr(-1, 36), Fr(-1, 36), Fr(1, 36))))
comb8_c = tuple(-16*T_p2_c[i] + 4*T_OK_c[i] for i in range(8))
check_exact("[S1] #8: comb = -16 T(p2) + 4 T(O_K) = 4A4 + 4A3 + 8L3 + 8L4 "
            "(exact)", comb8_c == (Fr(0), Fr(4), Fr(4), Fr(0),
                                   Fr(0), Fr(0), Fr(8), Fr(8)))
comb7_c = tuple(-1296*T_p6_c[i] + 324*T_p3_c[i] for i in range(8))
check_exact("[S2] #7: comb = -1296 T(p6) + 324 T(p3) = -36A4 + 36A3 + "
            "72L3 - 72L4 (exact)", comb7_c == (Fr(0), Fr(-36), Fr(36),
                                                Fr(0), Fr(0), Fr(0),
                                                Fr(72), Fr(-72)))
check_exact("[S3] zK and A7 cancel in both combinations (exact)",
            comb7_c[0] == comb7_c[3] == comb8_c[0] == comb8_c[3] == 0)

# Closure substitutions (Q3),(Q4):
#   L(g,3) = (pi^3/(42 sqrt21)) M        (Fricke, w = +1, level 84)
#   A3 = (8 sqrt21 pi^3/441) d3,  A4 = (4 sqrt21 pi^3/441) d4.
pref7 = Fr(5, 3)          # EK4 = pref7 (sqrt21/pi^3) comb7
cM3_7 = pref7*72*Fr(1, 42); cM4_7 = -pref7*72*Fr(1, 42)
cd3_7 = pref7*36*Fr(8*21, 441); cd4_7 = pref7*(-36)*Fr(4*21, 441)
check_exact("[S4] #7 closes: (M3,M4,d3,d4) = (20/7,-20/7,160/7,-80/7) = "
            "(20/7)(1,-1,8,-4) (exact)",
            (cM3_7, cM4_7, cd3_7, cd4_7)
            == (Fr(20, 7), Fr(-20, 7), Fr(160, 7), Fr(-80, 7))
            == tuple(Fr(20, 7)*v for v in (1, -1, 8, -4)))
pref8 = Fr(5)             # EK4 = pref8 (sqrt21/pi^3) comb8
cM3_8 = pref8*8*Fr(1, 42); cM4_8 = pref8*8*Fr(1, 42)
cd3_8 = pref8*4*Fr(8*21, 441); cd4_8 = pref8*4*Fr(4*21, 441)
check_exact("[S5] #8 closes: (M3,M4,d3,d4) = (20/21,20/21,160/21,80/21) "
            "= (20/21)(1,1,8,4) (exact)",
            (cM3_8, cM4_8, cd3_8, cd4_8)
            == (Fr(20, 21), Fr(20, 21), Fr(160, 21), Fr(80, 21))
            == tuple(Fr(20, 21)*v for v in (1, 1, 8, 4)))
check_exact("[S6] Fricke constant: 84^{3/2}/4 = 42 sqrt21 <=> 84^3 = "
            "16*42^2*21 (exact)", 84**3 == 16*42*42*21)
CQUAD = 893952**2 - 3*516096**2
check_exact("[S7] s4 quadratic X^2 + 1787904 X + c: 2*893952 = 1787904, "
            "disc = 12*516096^2, c = 84934656 > 0 (exact)",
            2*893952 == 1787904 and CQUAD == 84934656
            and 1787904**2 - 4*CQUAD == 12*516096**2)

print()
print("Separation track summary: entries #7/#8 proved by exact algebra from:")
print("  (Q0) n4(s4(tau)) = EK4(tau): tau7, tau8 on the certified leg B of")
print("       n5_line_cert.py (x = 1/2, y in [0.702, sqrt21/2]; interior)")
print("  (Q1) s4 values -893952 +- 516096 sqrt3: mp locks [S4-1/2] here,")
print("       iv locks in cert0_n5_e78.py (parallel)")
print("  (Q2) h(-84) = 4, class group = genus group = (Z/2)^2 [X1-X3]")
print("  (Q3) four CM theta series exhaust S_3(Gamma_0(84), chi_{-84})")
print("       [Hecke X10; LMFDB]; FEs w = +1 [Fricke ratios G:*]")
print("  (Q4) genus theory L_K(chi o N) = L(chi) L(chi chi_84) [X4];")
print("       ideal-count genus splits [X7]; anchors fitted in")
print("       n5_p8_e78_fit.py, locked [V1]-[V4]")
print("  (Q5) closed forms L(chi_28,2), L(chi_21,2) [X5 exact sums, Q4]")
print("  (Q6) Dirichlet formulas d3, d4 [D3/D4]; Fricke constant [S6]")
print()

# =====================================================================
# PART 1: [G]/[D]/[T]/[E]/[S4] numeric confirmations (mp, 60 dps)
# =====================================================================
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
        rr = fricke(a, N, mpf(yy))
        check("[G:%s] Fricke ratio at level %d = +1 (y = %s)"
              % (label, N, yy), rr, 1, mpf(10)**(-45))
    I0, I3 = mellin_I(a, xN, 0), mellin_I(a, xN, 3)
    Lam3 = xN**3*I3 + I0
    L3 = Lam3*(2*pi)**3/(xN**3*gamma(3))
    return L3, Lam3

LL = {}; ML = {}
for k in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
    LL[k], ML[k] = Lset(aE[k], 84, "g(%d,%d), N=84" % k)
    check("[G:M(%d,%d)] M = (42 sqrt21/pi^3) L(g,3) (FE constant)" % k,
          ML[k], 42*s21*LL[k]/pi**3)

L84v = dirichlet(mpf(2), [chi_84(n) for n in range(84)])
L4mv = dirichlet(mpf(2), [0, 1, 0, -1])
L21v = dirichlet(mpf(2), [chi_21(n) for n in range(21)])
L3mv = dirichlet(mpf(2), [0, 1, -1])
L28v = dirichlet(mpf(2), [chi_28(n) for n in range(28)])
L7mv = dirichlet(mpf(2), [chi_7m(n) for n in range(7)])
L12v = dirichlet(mpf(2), [chi_12(n) for n in range(12)])
zK = zeta(2)*L84v
A4 = L4mv*L21v; A3 = L3mv*L28v; A7 = L7mv*L12v
check("[D1] L(chi_28,2) = 2 sqrt7 pi^2/49 (closed form, Q4 + [X5])",
      L28v, 2*s7*pi**2/49)
check("[D2] L(chi_21,2) = 8 sqrt21 pi^2/441 (closed form, Q4 + [X5])",
      L21v, 8*s21*pi**2/441)
d3 = mpf(3)**1.5/(4*pi)*L3mv
d4 = 2*L4mv/pi
check("[D3] d3 = (3^{3/2}/4pi) L(chi_{-3},2) = L'(chi_{-3},-1) "
      "(direct derivative)", d3,
      mpdiff(lambda s: dirichlet(s, [0, 1, -1]), mpf(-1)))
check("[D4] d4 = (2/pi) L(chi_{-4},2) = L'(chi_{-4},-1) "
      "(direct derivative)", d4,
      mpdiff(lambda s: dirichlet(s, [0, 1, 0, -1]), mpf(-1)))

L1, L2, L3, L4 = LL[(1, 1)], LL[(1, -1)], LL[(-1, 1)], LL[(-1, -1)]
M3, M4 = ML[(-1, 1)], ML[(-1, -1)]

def G_row(x, y):
    return (pi/y)*sinh(2*pi*y)/(cosh(2*pi*y) - cos(2*pi*x))

def row_powers(x, y):
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y)/(2*y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy)/(2*yy), y)/(4*y)
    return S2, S3

def lattice_T(d, x0, y0):
    B, M = mpf(0), mpf(0)
    for m in range(-300, 301):
        if m == 0:
            B += 2*zeta(4); continue
        S2, S3 = row_powers(d*m*x0, abs(d*m)*y0)
        B += S2; M += m*m*S3
        if abs(d*m)*y0 > 45 and m > 0: break
    return 3*B - 4*d*d*y0*y0*M

def decomp(NC, bsgn, gsgn):
    return ((zK + bsgn[0]*A4 + bsgn[1]*A3 + bsgn[2]*A7)/(2*NC*NC)
            + (gsgn[0]*L1 + gsgn[1]*L2 + gsgn[2]*L3 + gsgn[3]*L4)/(NC*NC))

T_OKv = lattice_T(1, 0, s21)               # O_K = Z[sqrt-21]
T_OKalt = lattice_T(2, mpf(1)/2, s21/2)    # same lattice, rows of Z+2 tau8 Z
T_p2b = lattice_T(1, mpf(1)/2, s21/2)      # 16 T(p_2), tau = 1/2+i sqrt21/2
T_p3b = lattice_T(1, 0, s21/3)             # 81 T(p_3)
T_p3alt = lattice_T(2, mpf(1)/2, s21/6)    # same lattice, rows of Z+2 tau7 Z
T_p6b = lattice_T(1, mpf(1)/2, s21/6)      # 1296 T(p_6), tau = 1/2+i sqrt21/6

check("[T1] T(O_K) = (1/2)(zK+A4+A3+A7) + (L1+L2+L3+L4)",
      T_OKv, decomp(1, (1, 1, 1), (1, 1, 1, 1)))
check("[T1b] T(O_K) via rows of Z + 2 tau8 Z (same lattice Z[sqrt-21])",
      T_OKv, T_OKalt)
check("[T2] T(p_2) = (1/8)(zK-A4-A3+A7) + (1/4)(L1+L2-L3-L4)",
      T_p2b/16, decomp(2, (-1, -1, 1), (1, 1, -1, -1)))
check("[T3] T(p_3) = (1/18)(zK-A4+A3-A7) + (1/9)(L1-L2+L3-L4)",
      T_p3b/81, decomp(3, (-1, 1, -1), (1, -1, 1, -1)))
check("[T3b] 81 T(p_3) via rows of Z + 2 tau7 Z (same lattice)",
      T_p3b, T_p3alt)
check("[T4] T(p_6) = (1/72)(zK+A4-A3-A7) + (1/36)(L1-L2-L3+L4)",
      T_p6b/1296, decomp(6, (1, -1, -1), (1, -1, -1, 1)))

comb7 = -T_p6b + 4*T_p3b
comb7_dec = 36*(A3 - A4) + 72*(L3 - L4)
check("[T5] #7: comb = -1296 T(p6) + 324 T(p3) = 36(A3-A4) + 72(L3-L4)",
      comb7, comb7_dec)
comb8 = -T_p2b + 4*T_OKv
comb8_dec = 4*(A3 + A4) + 8*(L3 + L4)
check("[T6] #8: comb = -16 T(p2) + 4 T(O_K) = 4(A3+A4) + 8(L3+L4)",
      comb8, comb8_dec)

EK4_7 = (5*s21/(3*pi**3))*comb7
conj7 = Fr(20, 7)*(M3 - M4 + 8*d3 - 4*d4)
check("[E1] EK4(tau7) = (20/7)(M84^3 - M84^4 + 8 d3 - 4 d4)  [#7]",
      EK4_7, conj7)
EK4_7_dir = (10*(s21/6)/pi**3)*(-T_p6b + 4*T_p3alt)
check("[E1b] #7 via direct lattice pair (Z+2 tau7 Z rows)", EK4_7, EK4_7_dir)
check("[E3] #7 assembly from decomposed comb (FE constants only)",
      (5*s21/(3*pi**3))*comb7_dec, conj7)
EK4_8 = (5*s21/pi**3)*comb8
conj8 = Fr(20, 21)*(M3 + M4 + 8*d3 + 4*d4)
check("[E2] EK4(tau8) = (20/21)(M84^3 + M84^4 + 8 d3 + 4 d4)  [#8]",
      EK4_8, conj8)
EK4_8_dir = (10*(s21/2)/pi**3)*(-T_p2b + 4*T_OKalt)
check("[E2b] #8 via direct lattice pair (Z+2 tau8 Z rows)", EK4_8, EK4_8_dir)
check("[E4] #8 assembly from decomposed comb (FE constants only)",
      (5*s21/pi**3)*comb8_dec, conj8)

# --- s4 locks (mp track; iv locks in cert0_n5_e78.py) ----------------------
def eta_hp(tau, nterms=400):
    q = exp(2*pi*1j*tau)
    p = mpc(1); qn = q
    for n in range(1, nterms+1):
        p *= (1-qn); qn *= q
        if abs(qn) < mpf(10)**(-65): break
    return exp(pi*1j*tau/12)*p

def s4_hp(tau):
    e1, e2, e4 = eta_hp(tau), eta_hp(2*tau), eta_hp(4*tau)
    W = e1*e4**2/e2**3
    return (e2/e1)**24 * (16*W**4 + W**(-4))**4

t7 = mpc(mpf(1)/2, s21/6)
t8 = mpc(mpf(1)/2, s21/2)
check("[S4-1] s4(tau7) = -893952 + 516096 sqrt3 (mp lock; iv: cert0_n5_e78)",
      s4_hp(t7), -893952 + 516096*s3, mpf(10)**(-45))
check("[S4-2] s4(tau8) = -893952 - 516096 sqrt3 (mp lock; iv: cert0_n5_e78)",
      s4_hp(t8), -893952 - 516096*s3, mpf(10)**(-45))
check("[S4-3] q7 = exp(2 pi i tau7) = -e^{-pi sqrt21/3} (negative real)",
      exp(2*pi*1j*t7), -exp(-pi*s21/3), mpf(10)**(-45))
check("[S4-4] q8 = exp(2 pi i tau8) = -e^{-pi sqrt21} (negative real)",
      exp(2*pi*1j*t8), -exp(-pi*s21), mpf(10)**(-45))

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
else:
    print("ALL CHECKS PASSED (pre-iv tracks)")

# =====================================================================
# PART 2: [V] rigorous interval locks (mpmath.iv, iv.dps = 70)
#   Machinery ported from verify_P1_n5_e6.py.  Mellin tail uses the
#   self-contained theta bound |a(n)| <= 30 n^3:
#     representation counts are boxed: r_C(n) <= (2 sqrt n+1)^2 <= 12 n;
#     summands obey |a^2-21b^2|/2 <= n/2 (P0), |(2x+y)^2-21y^2|/4 <= n/2
#     since (2x+y)^2+21y^2 = 2n (P2), |3a^2-7b^2|/2 <= n/2 (P3),
#     |(6x+3y)^2-21y^2|/12 <= n/2 since (6x+3y)^2+21y^2 = 6n (P6);
#     hence each |P_C(n)| <= 6 n^2 and |a(n)| <= 24 n^2 <= 30 n^3.
#   iv lattice sums: O_K (y0 = sqrt21) and p_3 (y0 = sqrt21/3) have all
#   rows at shift 0 (pure-imaginary lattices); p_2 (y0 = sqrt21/2) and
#   p_6 (y0 = sqrt21/6) sit on x = 1/2, so rows alternate shift 0 (m even)
#   and shift 1/2 (m odd) -- the lattice_T_iv_e6 pattern, parametrized.
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
    # rigorous tail bound with |a_n| <= 30 n^3 (theta series, see header);
    # Gamma(s,x) <= s! e^-x (1+x)^s (s=0: E1(x) <= e^-x(1+1/x)).
    c = 2*iv.pi/xN
    def Abound(n):
        nn = iv.mpf(n)
        xn = c*nn
        g = iv.exp(-xn)*(1+1/xn) if s == 0 else iv.exp(-xn)*(1+xn)**s
        sfac = iv.mpf(1) if s == 0 else iv.mpf([1, 1, 2, 6][s])
        return (30*nn**3*(2*iv.pi*nn)**(-s)*sfac*g).b
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
    def __init__(a, f0, f1=0, f2=0):
        a.f0, a.f1, a.f2 = iv.mpf(f0), iv.mpf(f1), iv.mpf(f2)
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

def cs_dual(yv, shift):
    ep0 = iv.exp(iv.pi*yv)
    em0 = iv.exp(-iv.pi*yv)
    ep = Dual(ep0, iv.pi*ep0, iv.pi**2*ep0/2)
    em = Dual(em0, -iv.pi*em0, iv.pi**2*em0/2)
    return (ep-em)/(ep+em) if shift else (ep+em)/(ep-em)

def row_S23_iv(y, shift=False):
    G = (Dual(iv.pi)/Dual(y, 1, 0)) * cs_dual(y, shift)
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

def lattice_T_iv_half(y0):
    """T-sum for tau = 1/2 + i y0 (lattice Z + Z tau): rows alternate
    shift 0 (m even, coth) and shift 1/2 (m odd, tanh); power terms and
    tail bounds are shift-independent (Poisson kernel k = 0 term)."""
    y0 = iv.mpf(y0)
    y0m = mpf(y0.a)
    M = int(50/y0m) + 2
    T = iv.pi**4/15
    for m in range(1, M+1):
        y = iv.mpf(m)*y0
        S2, S3 = row_S23_iv(y, shift=(m % 2 == 1))
        dS2 = S2 - iv.pi/(2*y**3)
        dS3 = S3 - 3*iv.pi/(8*y**5)
        T += 2*(3*dS2 - 4*y**2*dS3)
    tb2, tb3 = tail_row_iv((M+1)*y0m)
    rr = mpf(iv.exp(-2*iv.pi*y0m).b)
    yM = (M+1)*y0m
    tot = 2*(3*tb2 + 4*yM*yM*tb3)/(1-rr)
    return T + iv.mpf([-tot, tot])

# --- theta coefficients to n0(84) = 452 for the iv Mellin tails ----------
NMV = 500
P0r = [0]*(NMV+1); P2r = [0]*(NMV+1)
P3r = [0]*(NMV+1); P6r = [0]*(NMV+1)
RBV = int((2*NMV)**0.5)+3
for a_ in range(-RBV, RBV+1):
    for b_ in range(-RBV, RBV+1):
        n_ = a_*a_ + 21*b_*b_
        if 1 <= n_ <= NMV: P0r[n_] += a_*a_ - 21*b_*b_
        n_ = 3*a_*a_ + 7*b_*b_
        if 1 <= n_ <= NMV: P3r[n_] += 3*a_*a_ - 7*b_*b_
for x_ in range(-RBV, RBV+1):
    for y_ in range(-RBV, RBV+1):
        n_ = 2*x_*x_ + 2*x_*y_ + 11*y_*y_
        if 1 <= n_ <= NMV: P2r[n_] += ((2*x_+y_)**2 - 21*y_*y_)//4
        n_ = 6*x_*x_ + 6*x_*y_ + 5*y_*y_
        if 1 <= n_ <= NMV: P6r[n_] += ((6*x_+3*y_)**2 - 21*y_*y_)//12
assert all(v % 2 == 0 for v in P0r[1:]) and all(v % 2 == 0 for v in P3r[1:])
P0r = [v//2 for v in P0r]; P3r = [v//2 for v in P3r]
aV = {}
for e2 in (1, -1):
    for e3 in (1, -1):
        aV[(e2, e3)] = [P0r[n] + e2*P2r[n] + e3*P3r[n] + e2*e3*P6r[n]
                        for n in range(NMV+1)]
check_exact("[V-pre] iv-range coefficients coincide with the exact track "
            "(n <= 400)", all(aV[k][:NMAX+1] == aE[k] for k in aE))

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

L1_iv, M1_iv = Lset_iv(aV[(1, 1)], 84)
L2_iv, M2_iv = Lset_iv(aV[(1, -1)], 84)
L3_iv, M3_iv = Lset_iv(aV[(-1, 1)], 84)
L4_iv, M4_iv = Lset_iv(aV[(-1, -1)], 84)
L84_iv = dirichlet2_iv([chi_84(n) for n in range(84)], 84)
L4m_iv = dirichlet2_iv([0, 1, 0, -1], 4)
L21_iv = dirichlet2_iv([chi_21(n) for n in range(21)], 21)
L3m_iv = dirichlet2_iv([0, 1, -1], 3)
L28_iv = dirichlet2_iv([chi_28(n) for n in range(28)], 28)
L7m_iv = dirichlet2_iv([chi_7m(n) for n in range(7)], 7)
L12_iv = dirichlet2_iv([chi_12(n) for n in range(12)], 12)

check("[V0a] iv: L(g(1,1),3) vs mp value", iv_mid(L1_iv), L1, iv_w(L1_iv)+TOL)
check("[V0b] iv: L(g(1,-1),3) vs mp value", iv_mid(L2_iv), L2,
      iv_w(L2_iv)+TOL)
check("[V0c] iv: L(g(-1,1),3) vs mp value", iv_mid(L3_iv), L3,
      iv_w(L3_iv)+TOL)
check("[V0d] iv: L(g(-1,-1),3) vs mp value", iv_mid(L4_iv), L4,
      iv_w(L4_iv)+TOL)
check("[V0e] iv: L(chi_28,2) = 2 sqrt7 pi^2/49", iv_mid(L28_iv),
      2*s7*pi**2/49, iv_w(L28_iv)+TOL)
check("[V0f] iv: L(chi_21,2) = 8 sqrt21 pi^2/441", iv_mid(L21_iv),
      8*s21*pi**2/441, iv_w(L21_iv)+TOL)

zK_iv = iv.pi**2/6 * L84_iv
A4_iv = L4m_iv * L21_iv
A3_iv = L3m_iv * L28_iv
A7_iv = L7m_iv * L12_iv
d3_iv = 3*iv.sqrt(3)/(4*iv.pi) * L3m_iv
d4_iv = 2*L4m_iv/iv.pi

s21_iv = iv.sqrt(21)
T_OK_iv = lattice_T_iv(1, s21_iv)          # Lambda_2 of #8 (= O_K)
T_p2_iv = lattice_T_iv_half(s21_iv/2)      # Lambda_1 of #8 (= 16 T(p_2))
T_p3_iv = lattice_T_iv(1, s21_iv/3)        # Lambda_2 of #7 (= 81 T(p_3))
T_p6_iv = lattice_T_iv_half(s21_iv/6)      # Lambda_1 of #7 (= 1296 T(p_6))

check("[V0g] iv T(O_K) contains the 60-dps mp value",
      iv_mid(T_OK_iv), T_OKv, iv_w(T_OK_iv) + TOL)
check("[V0h] iv 16 T(p_2) contains the 60-dps mp value",
      iv_mid(T_p2_iv), T_p2b, iv_w(T_p2_iv) + TOL)
check("[V0i] iv 81 T(p_3) contains the 60-dps mp value",
      iv_mid(T_p3_iv), T_p3b, iv_w(T_p3_iv) + TOL)
check("[V0j] iv 1296 T(p_6) contains the 60-dps mp value",
      iv_mid(T_p6_iv), T_p6b, iv_w(T_p6_iv) + TOL)

T_OK_dec_iv = ((zK_iv + A4_iv + A3_iv + A7_iv)/2
               + (L1_iv + L2_iv + L3_iv + L4_iv))
T_p2_dec_iv = ((zK_iv - A4_iv - A3_iv + A7_iv)/8
               + (L1_iv + L2_iv - L3_iv - L4_iv)/4)
T_p3_dec_iv = ((zK_iv - A4_iv + A3_iv - A7_iv)/18
               + (L1_iv - L2_iv + L3_iv - L4_iv)/9)
T_p6_dec_iv = ((zK_iv + A4_iv - A3_iv - A7_iv)/72
               + (L1_iv - L2_iv - L3_iv + L4_iv)/36)

check_lock("[V1] LOCK T(O_K) = (1/2)(zK+A4+A3+A7) + (L1+L2+L3+L4)",
           T_OK_iv - T_OK_dec_iv)
check_lock("[V2] LOCK T(p_2) = (1/8)(zK-A4-A3+A7) + (1/4)(L1+L2-L3-L4)",
           T_p2_iv/16 - T_p2_dec_iv)
check_lock("[V3] LOCK T(p_3) = (1/18)(zK-A4+A3-A7) + (1/9)(L1-L2+L3-L4)",
           T_p3_iv/81 - T_p3_dec_iv)
check_lock("[V4] LOCK T(p_6) = (1/72)(zK+A4-A3-A7) + (1/36)(L1-L2-L3+L4)",
           T_p6_iv/1296 - T_p6_dec_iv)

comb7_iv = -T_p6_iv + 4*T_p3_iv
EK4_7_iv = (5*s21_iv/(3*iv.pi**3))*comb7_iv
conj7_iv = iv.mpf(20)/7*(M3_iv - M4_iv + 8*d3_iv - 4*d4_iv)
check_lock("[V5] LOCK #7: EK4(tau7) = (20/7)(M84^3 - M84^4 + 8 d3 - 4 d4)",
           EK4_7_iv - conj7_iv, mpf(10)**(-48))

comb8_iv = -T_p2_iv + 4*T_OK_iv
EK4_8_iv = (5*s21_iv/iv.pi**3)*comb8_iv
conj8_iv = iv.mpf(20)/21*(M3_iv + M4_iv + 8*d3_iv + 4*d4_iv)
check_lock("[V6] LOCK #8: EK4(tau8) = (20/21)(M84^3 + M84^4 + 8 d3 + 4 d4)",
           EK4_8_iv - conj8_iv, mpf(10)**(-48))

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
    import sys; sys.exit(1)
print("ALL CHECKS PASSED")
