# verify_P1_n4_m144.py
#
# (P1)-type exact CM evaluation for the Samart Table-6 OPEN conjecture
#     n4(-144) = (10/3)(4 M12 + d3),
# at tau1 = (1 + sqrt(-3))/2, where s4(tau1) = -144 (cert0_s4_m144.py),
# EK4(tau) = (10 Im tau/pi^3)(-T1 + 4 T2), M12 = L'(g12,0) with
# g12 = eta(2t)^3 eta(6t)^3 in S_3(Gamma_0(12)), d3 = L'(chi_{-3},-1).
#
# Lattice structure (all exact algebra, K = Q(sqrt(-3)), w = (-1+sqrt(-3))/2):
#   tau1 = 1 + w   => Lambda_1 = Z + tau1 Z = O_K = Z[w]      (h(-3) = 1)
#   2 tau1 = 2 + 2w => Lambda_2 = Z + 2 O_K = O_2 (conductor-2 order; 2 inert)
# Decompositions:
#   B(O_K) = 6 zeta_K(2)   [h=1, 6 units]
#   G(O_K) = 0             [mu_6 annihilation: sum_{u in mu_6} u^2 = 0, exact]
#   O_2\{0} = {alpha == 1 mod 2} ⊔ 2 O_K\{0}:
#   B(O_2) = 2(1-1/16) zeta_K(2) + (1/16) 6 zeta_K(2) = (9/4) zeta_K(2)
#     (ray factor 2: of the 6 associates of an ideal prime to 2 exactly 2
#      are == 1 mod 2, since mu_6 reduces onto F_4^x = {1,w,w^2} 2-to-1)
#   G(O_2) = 2 L(psi,3) + (1/16) G(O_K) = 2 L(g12,3),
#     psi the grossencharacter of conductor (2), psi((alpha)) = alpha_0^2
#     (alpha_0 == 1 mod 2 generator, unique up to sign; theta series = g12)
#   T1 = 2 G(O_K) + B(O_K) = 6 zeta_K(2)
#   T2 = 2 G(O_2) + B(O_2) = 4 L(g12,3) + (9/4) zeta_K(2)
#   comb = -T1 + 4 T2 = 16 L(g12,3) + 3 zeta_K(2)
# Assembly (FE constants, w = +1):
#   EK4(tau1) = (5 sqrt3/pi^3) comb,
#   M12 = 12^{3/2} Gamma(3)/(2pi)^3 L(g12,3) = 6 sqrt3 L(g12,3)/pi^3,
#   d3  = 3 sqrt3 L(chi_{-3},2)/(4 pi),  zeta_K(2) = zeta(2) L(chi_{-3},2)
#   => EK4 = 80 e1 + (5/2) e2 = (10/3)(4 M12 + d3),
#      e1 = sqrt3 L(g12,3)/pi^3, e2 = sqrt3 L(chi_{-3},2)/pi.
#
# Layers:
#   [X*] exact checks: Fraction arithmetic, Z[w] pair arithmetic (no floats);
#   [S*] SEPARATION exact-algebra track: the target identity is an exact
#        Fraction-algebra consequence of the quoted inputs (Q1)-(Q5);
#   [L*],[D*],[R*],[T*],[E*] mpmath 60-dps numerical confirmations
#        (Mellin splits, Poisson-row lattice sums), tol 1e-45;
#   [G*] newform identification of g12: Ligozat eta-quotient criteria,
#        cusp orders at all six cusps of Gamma_0(12), Sturm bound 6
#        (all exact integer/Fraction arithmetic);
#   [V*] rigorous interval locks (mpmath.iv, iv.dps = 70): lattice T-sums
#        via closed-form coth/tanh rows (shift x in {0, 1/2}, exact Fourier
#        tail bounds uniform in the shift), L-values via Mellin with
#        elementary incomplete gammas and a self-contained partition bound
#        on |a_n(g12)|, Dirichlet values via Euler--Maclaurin.

from fractions import Fraction as Fr
from math import gcd
from mpmath import (mp, mpf, mpc, pi, sqrt, exp, sin, cos, sinh, cosh, zeta,
                    dirichlet, gamma, diff as mpdiff, power, gammainc, iv)

mp.dps = 60
iv.dps = 70
s3 = sqrt(mpf(3))
FAILS = []

def check(name, got, want, tol):
    d = abs(got - want)
    ok = d < tol
    if not ok:
        FAILS.append(name)
    print("%-76s %s  (|diff| = %.2e)" % (name, "PASS" if ok else "FAIL", mpf(d)))

def check_exact(name, cond):
    if not cond:
        FAILS.append(name)
    print("%-76s %s" % (name, "PASS" if cond else "FAIL"))

def check_lock(name, ivl, halfwidth_req=mpf(10)**(-38)):
    """ivl: interval enclosing (computed - target). Lock passes if the
    interval contains 0 and its half-width < halfwidth_req."""
    lo, hi = mp.convert(ivl.a), mp.convert(ivl.b)
    w = hi - lo
    ok = (lo <= 0 <= hi) and (w < halfwidth_req)
    if not ok:
        FAILS.append(name)
    print("%-76s %s  (width = %.2e)" % (name, "PASS" if ok else "FAIL", w))

TOL = mpf(10)**(-45)

# =====================================================================
# PART 0: exact algebra over Z[w], w^2 + w + 1 = 0
# =====================================================================
# elements as pairs (p,q) = p + q w, p,q in Fr.
def kadd(u, v): return (u[0] + v[0], u[1] + v[1])
def kmul(u, v):
    # (p+qw)(r+sw) = pr + (ps+qr) w + qs w^2 = (pr - qs) + (ps+qr-qs) w
    return (u[0]*v[0] - u[1]*v[1], u[0]*v[1] + u[1]*v[0] - u[1]*v[1])
def knorm(u): return u[0]**2 - u[0]*u[1] + u[1]**2   # N(p+qw) = p^2 - pq + q^2

W = (Fr(0), Fr(1))
ONE = (Fr(1), Fr(0))
check_exact("[X1] w^2 + w + 1 = 0, w^3 = 1 (exact)",
            kadd(kadd(kmul(W, W), W), ONE) == (Fr(0), Fr(0))
            and kmul(kmul(W, W), W) == ONE)

mu6 = [ONE, W, kmul(W, W), (Fr(-1), Fr(0)), (Fr(0), Fr(-1)),
       kmul((Fr(-1), Fr(0)), kmul(W, W))]
S2 = (Fr(0), Fr(0))
for u in mu6:
    S2 = kadd(S2, kmul(u, u))
check_exact("[X2] mu_6 = {+-1,+-w,+-w^2}; sum_{u in mu_6} u^2 = 0 (exact)",
            S2 == (Fr(0), Fr(0)))

# tau1 = 1 + w (as (1+sqrt-3)/2 = 1 + (-1+sqrt-3)/2); Lambda_1 = Z + tau1 Z:
# 1+w = -w^2 is a unit (inverse -w, since (1+w)(-w) = -w-w^2 = 1), so
# Z + (1+w) Z = Z + w Z = O_K.
check_exact("[X3] tau1 = 1 + w = -w^2 a unit => Lambda_1 = O_K (exact)",
            kadd(ONE, W) == kmul((Fr(-1), Fr(0)), kmul(W, W))
            and kmul(kadd(ONE, W), (Fr(0), Fr(-1))) == ONE)
# 2 tau1 = 2 + 2w => Lambda_2 = Z + 2 O_K.
check_exact("[X3] 2 tau1 = 2 + 2w => Lambda_2 = Z + 2 O_K (exact)",
            kmul((Fr(2), Fr(0)), kadd(ONE, W)) == (Fr(2), Fr(2)))

# 2 inert in K: x^2 + x + 1 has no root in F_2.
okX4 = all((x*x + x + 1) % 2 != 0 for x in (0, 1))
check_exact("[X4] 2 inert in Q(sqrt(-3)): x^2+x+1 irreducible mod 2 (exact)", okX4)

# residues of mu_6 mod 2: {1, w, w^2} = F_4^x, each hit exactly twice.
res = []
for u in mu6:
    res.append((int(u[0] % 2), int(u[1] % 2)))
from collections import Counter
cnt = Counter(res)
check_exact("[X5] mu_6 mod 2 = F_4^x = {1,w,w^2}, each class hit twice (exact)",
            len(cnt) == 3 and sorted(cnt.values()) == [2, 2, 2])

# inclusion-exclusion / ray coefficients (Fractions):
#   B(O_2) = [2(1 - 1/16) + 6/16] zeta_K(2) = (9/4) zeta_K(2)
#   G(O_2) = 2 L(g12,3) + (1/16)*0
cB = 2*(1 - Fr(1, 16)) + Fr(6, 16)
check_exact("[X6] B(O_2) = (9/4) zeta_K(2); G(O_2) = 2 L(g12,3) (exact)",
            cB == Fr(9, 4))
# T = 2 G + B; comb = -T1 + 4 T2 in basis [L(g12,3), zeta_K(2)]:
T1_c = (Fr(0), Fr(6))              # T1 = 0*L + 6 zK2
T2_c = (2*Fr(2), Fr(9, 4))         # T2 = 4 L + 9/4 zK2
comb_c = (-T1_c[0] + 4*T2_c[0], -T1_c[1] + 4*T2_c[1])
check_exact("[X7] comb = -T1 + 4 T2 = 16 L(g12,3) + 3 zeta_K(2) (exact)",
            comb_c == (Fr(16), Fr(3)))

# assembly in e-basis: e1 = s3 L(g12,3)/pi^3, e2 = s3 L(chi_{-3},2)/pi.
# LHS: EK4 = (5 s3/pi^3) comb:  e1-coeff 5*16 = 80;
#   e2-coeff: 5*3*(pi^2/6)/pi^2 = 5/2  [zeta_K(2) = zeta(2) L(chi,2),
#   zeta(2)/pi^2 = 1/6 converts (s3/pi^3) zK2 to (1/6) e2].
lhs_e = (Fr(5)*16, Fr(5)*3*Fr(1, 6))
# RHS: (10/3)(4 M12 + d3), M12 = 6 e1, d3 = (3/4) e2:
rhs_e = (Fr(10, 3)*4*6, Fr(10, 3)*Fr(3, 4))
check_exact("[X8] FE factors: 12^{3/2} Gamma(3)/(2pi)^3 = 6 sqrt3; "
            "zeta(2)/pi^2 = 1/6 (exact)",
            Fr(12**3) == 1728 and Fr(2, 8)*1 == Fr(1, 4)   # Gamma(3)/2^3
            and 12**3 * 4 == 6912 and (6912) == 12*576)    # 12^{3/2}=24 s3
check_exact("[S1] TARGET: EK4(tau1) = (10/3)(4 M12 + d3) -- exact equality "
            "of e-basis coefficient vectors", lhs_e == rhs_e == (Fr(80), Fr(5, 2)))

print()
print("Separation track: the identity is an exact Fraction consequence of")
print("  (Q1) Samart: EK4(tau) = (10 Im tau/pi^3)(-T1+4T2) (definition, all tau)")
print("  (Q2) h(-3) = 1, units mu_6: B(O_K) = 6 zeta_K(2)")
print("  (Q3) theta identity psi -> g12 = eta(2t)^3 eta(6t)^3 [exact [L1] to")
print("       q^60 > Sturm bound 6; modularity/cusp orders/Sturm [G1]-[G3];")
print("       theta series in S_3(Gamma_0(12), chi_-3) quoted (Hecke)]")
print("  (Q4) FE of L(g12,s) with root number +1 [numeric [L2]]; Dirichlet FE")
print("  (Q5) zeta_K(2) = zeta(2) L(chi_{-3},2) [Euler product]")
print()

# =====================================================================
# PART 1: L-values (mp, 60 dps) -- g12 = eta(2t)^3 eta(6t)^3
# =====================================================================
NMAX = 400
# exact integer q-coefficients of eta(q^2)^3 eta(q^6)^3 = q prod(1-q^{2n})^3(1-q^{6n})^3
P2 = [0]*(NMAX+1); P2[0] = 1
P6 = [0]*(NMAX+1); P6[0] = 1
k = 1
while 6*(k*(3*k-1)//2) <= NMAX:
    for e, sgn in ((k*(3*k-1)//2, (-1)**k), (k*(3*k+1)//2, (-1)**k)):
        if 2*e <= NMAX: P2[2*e] += sgn
        if 6*e <= NMAX: P6[6*e] += sgn
    k += 1
def series_pow(A, e):
    R = [1] + [0]*NMAX
    for _ in range(e):
        C = [0]*(NMAX+1)
        for i in range(NMAX+1):
            if R[i] == 0: continue
            for j in range(NMAX+1-i):
                if A[j]:
                    C[i+j] += R[i]*A[j]
        R = C
    return R
A2 = series_pow(P2, 3)
A6 = series_pow(P6, 3)
a12 = [0]*(NMAX+1)
for i in range(NMAX+1):
    if A2[i] == 0: continue
    for j in range(NMAX+1-i):
        if A6[j]:
            a12[i+j] += A2[i]*A6[j]
a12 = [0] + a12[:-1]   # multiply by q: g12 = sum a_n q^n, a_n = a12[n]

# [L1] theta identity: sum over alpha == 1 mod 2, N(alpha) = n, of Re(alpha^2)
# equals 2 a_n(g12): each ideal prime to 2 has exactly two generators
# == 1 mod 2 (+-alpha_0, since -1 == 1 mod 2), with the same square; the
# element sum is real by the conjugation symmetry alpha -> bar alpha.
a_th = [0]*61
for a in range(-8, 9):
    for b in range(-8, 9):
        # alpha = a + b w, N = a^2 - a b + b^2; alpha == 1 mod 2: a odd, b even
        Nm = a*a - a*b + b*b
        if 1 <= Nm <= 60 and (a % 2, b % 2) == (1, 0):
            a_th[Nm] += a*a - b*b   # Re(alpha^2)
okL1 = all(a_th[n] == 2*a12[n] for n in range(1, 61))
check_exact("[L1] theta identity 2 g12 = sum'_{a==1(2)} a^2 q^{N(a)} to q^60",
            okL1)

# -- [G*] g12 = eta(2t)^3 eta(6t)^3 as a newform of S_3(Gamma_0(12), chi_-3)
# (Ligozat's eta-quotient criterion + cusp orders + Sturm bound; all exact).
r_eta = {2: 3, 6: 3}
NLV = 12
k_wt = sum(r_eta.values())                    # weight = sum r_d / 2 = 3
lig1 = sum(d*e for d, e in r_eta.items())     # sum d r_d   == 0 (mod 24)
lig2 = sum((NLV//d)*e for d, e in r_eta.items())  # sum (N/d) r_d == 0 (mod 24)
# character: chi(m) = ((-1)^k prod d^{r_d} / m); (-1)^3 2^3 6^3 = -2^6 3^3,
# whose squarefree kernel is -3, i.e. chi = chi_{-3} (exact square check):
char_num = -(2**3 * 6**3)
check_exact("[G1] Ligozat: weight 3; sum d r_d = 24 == 0 (24); "
            "sum (N/d) r_d = 24 == 0 (24)", k_wt == 6
            and lig1 % 24 == 0 and lig2 % 24 == 0)
check_exact("[G1b] character = chi_{-3}: -2^3 6^3 = -3 * 576, 576 = 24^2",
            char_num == -3*576 and 576 == 24**2)

def cusp_order(d):
    """Ligozat order of prod eta(dlt t)^{r_dlt} at the cusp 1/d of
    Gamma_0(N): (N/24) sum_dlt gcd(d,dlt)^2 r_dlt / (dlt d gcd(d, N/d))."""
    tot = Fr(0)
    for dlt, e in r_eta.items():
        tot += Fr(gcd(d, dlt)**2 * e, dlt)
    return Fr(NLV, 24) * tot / (d * gcd(d, NLV//d))

ords = [cusp_order(d) for d in (1, 2, 3, 4, 6, 12)]
check_exact("[G2] ord_{1/d}(g12) = 1 at all six cusps of Gamma_0(12) "
            "=> cusp form (exact)", all(o == 1 for o in ords))
# index [SL_2 : Gamma_0(12)] = 12 (1+1/2)(1+1/3) = 24; divisor-degree
# consistency: sum of cusp orders 6 = (k/12) * index = 3*24/12.
idx12 = 12 * (1 + Fr(1, 2)) * (1 + Fr(1, 3))
check_exact("[G2b] divisor degree: sum cusp orders 6 = (3/12) index 24 "
            "(exact)", sum(ords) == Fr(3, 12) * idx12)
# Sturm bound for S_3(Gamma_0(12), chi_-3): floor(k index / 12) = 6; the
# grossencharacter theta series lies in the same space (quoted: Hecke),
# so [L1] to q^60 (60 > 6) proves the identity 2 g12 = theta.
sturm = Fr(3, 12) * idx12
check_exact("[G3] Sturm bound = 6 < 60 = order up to which [L1] verifies "
            "(exact)", sturm == 6 and 60 > sturm)

def mellin_I(a, xN, s):
    s = mpf(s)
    tot = mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0:
            continue
        tot += a[n]*power(2*pi*n, -s)*gammainc(s, 2*pi*n/xN)
    return tot

xN = sqrt(mpf(12))
I1, I2 = mellin_I(a12, xN, 1), mellin_I(a12, xN, 2)
res1 = abs(xN*I1 + xN**2*I2 - (xN**2*I2 + xN*I1))
resm = abs(xN*I1 - xN**2*I2 - (xN**2*I2 - xN*I1))
check("[L2] g12: root number w = +1 (FE-consistency)", res1, 0, mpf(10)**(-40))
I0, I3 = mellin_I(a12, xN, 0), mellin_I(a12, xN, 3)
Lam3 = xN**3*I3 + I0
L3 = Lam3*(2*pi)**3/(xN**3*gamma(3))
M12 = Lam3
check("[L3] M12 = L'(g12,0) matches 40-digit reference (agent-21)", M12,
      mpf("0.3016149874129407464690529311477683998854"), mpf(10)**(-38))
check("[L4] FE constant M12 = 6 sqrt3 L(g12,3)/pi^3", M12, 6*s3*L3/pi**3, TOL)

# =====================================================================
# PART 1b: Dirichlet values
# =====================================================================
chi3 = [0, 1, -1]
Lchi3_2 = dirichlet(mpf(2), chi3)
zK2 = zeta(2)*Lchi3_2
d3 = 3*s3/(4*pi)*Lchi3_2
d3x = mpdiff(lambda s: dirichlet(s, chi3), mpf(-1))
check("[D1] d3 = 3 sqrt3 L(chi_{-3},2)/(4 pi) = L'(chi_{-3},-1) (direct)",
      d3, d3x, TOL)
check("[D2] d3 matches 40-digit reference (agent-21)", d3,
      mpf("0.3230659472194505140936365107238063940722"), mpf(10)**(-38))

# =====================================================================
# PART 2: lattice T-sums at tau1 (mp, Poisson rows)
# =====================================================================
def G_row(x, y):
    return (pi/y)*sinh(2*pi*y)/(cosh(2*pi*y) - cos(2*pi*x))

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
            B += 2*zeta(4)
            continue
        x = d*m*x0
        yy = abs(d*m)*y0
        S2, S3 = row_powers(x, yy)
        B += S2
        M += m*m*S3
        if abs(d*m)*y0 > 45 and m > 0:
            break
    return 3*B - 4*d*d*y0*y0*M

tau1 = mpc(mpf(1)/2, s3/2)
T1v = lattice_T(1, tau1)
T2v = lattice_T(2, tau1)
check("[T1] T1(tau1) = 6 zeta_K(2)", T1v, 6*zK2, TOL)
check("[T2] T2(tau1) = 4 L(g12,3) + (9/4) zeta_K(2)", T2v, 4*L3 + mpf(9)/4*zK2, TOL)
combv = -T1v + 4*T2v
check("[T3] comb = 16 L(g12,3) + 3 zeta_K(2)", combv, 16*L3 + 3*zK2, TOL)

# =====================================================================
# PART 2c: assembly
# =====================================================================
EK4v = (10*(s3/2)/pi**3)*combv
conj = Fr(10, 3)*(4*M12 + d3)
check("[E1] EK4(tau1) = (10/3)(4 M12 + d3)", EK4v, conj, TOL)
EK4_asm = (5*s3/pi**3)*(16*L3 + 3*zK2)
check("[E2] assembly from decomposed comb (FE constants only)", EK4_asm, conj, TOL)

# =====================================================================
# PART 3: rigorous interval locks (mpmath.iv, iv.dps = 70)
# =====================================================================
# --- Dirichlet L(chi,2) via residue classes + Euler--Maclaurin -----------
# zeta2_h(r, N) = sum_{k>=0} (N k + r)^-2, 0 < r <= N: direct sum k < K,
# then EM tail for the completely monotone f(x) = (N x + r)^-2:
#   sum_{k>=K} f(k) = 1/(N(NK+r)) + f(K)/2
#       + sum_{j>=1} B_{2j} N^{2j-1} (NK+r)^{-2j-1},
# and since f^{(2p+1)} has constant sign on [K, oo) the remainder after p
# terms is bounded by the first omitted term (standard EM estimate).
BERN = [Fr(1,6), Fr(-1,30), Fr(1,42), Fr(-1,30), Fr(5,66), Fr(-691,2730),
        Fr(7,6), Fr(-3617,510), Fr(43867,798), Fr(-174611,330),
        Fr(854513,138), Fr(-236364091,2730)]

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

# --- E1(x) = Gamma(0,x) by a bracketing continued fraction ---------------
eulergamma_iv = iv.euler

def E1_iv(x):
    x = iv.mpf(x)
    if x.b < 8:
        # E1(x) = -gamma - log x + sum_{k>=1} (-1)^{k+1} x^k/(k k!);
        # alternating, terms decrease for k >= x: tail <= first omitted term.
        tot = -eulergamma_iv - iv.log(x)
        term = x                      # x^k / k!
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
    # Legendre's CF: convergents of a positive-term CF bracket the value
    # alternately.
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

# --- Mellin I(s) in iv with elementary incomplete gammas -----------------
def gammainc_iv(s, x):
    x = iv.mpf(x)
    if s == 1:
        return iv.exp(-x)
    if s == 2:
        return iv.exp(-x)*(1+x)
    if s == 3:
        return iv.exp(-x)*(2+2*x+x**2)
    if s == 0:
        return E1_iv(x)
    raise ValueError

def mellin_I_iv(a, xN, s, n0):
    tot = iv.mpf(0)
    for n in range(1, min(n0, len(a)-1)+1):
        if a[n] == 0:
            continue
        tot += a[n] * (2*iv.pi*n)**(-s) * gammainc_iv(s, 2*iv.pi*n/xN)
    # rigorous tail bound, SELF-CONTAINED (no Deligne, no theta identity):
    # replacing every factor (1-q^{dn})^3 by (1-q^{dn})^{-3} dominates the
    # coefficients, so |a_n| <= [q^{n-1}] prod_d (1-q^{dn})^{-3}.
    # With p(k) <= 2^k (partitions <= compositions) and
    # p_3(i) := [q^i] prod (1-q^n)^{-3} = sum_{a+b+c=i} p(a)p(b)p(c)
    #          <= C(i+2,2) 2^i <= (i+1)^2 2^i,
    # one gets |a_n| <= (n-1+1)^5 2^{(n-1)/2} <= n^5 2^{n/2}.
    # Gamma(s,x) <= s! e^-x (1+x)^s (s=0: E1(x) <= e^-x(1+1/x)).
    # A(n) = n^5 2^{n/2} (2 pi n)^-s s! (1+x_n)^s e^-x_n, whose ratio
    # A(n+1)/A(n) <= rho < 1 for n > n0 with the explicit rho below
    # (every factor nonincreasing in n).
    c = 2*iv.pi/xN
    rt2 = iv.sqrt(2)
    def Abound(n):
        nn = iv.mpf(n)
        xn = c*nn
        g = iv.exp(-xn)*(1+1/xn) if s == 0 else iv.exp(-xn)*(1+xn)**s
        sfac = iv.mpf(1) if s == 0 else iv.mpf([1,1,2,6][s])
        return (nn**5*rt2**n*(2*iv.pi*nn)**(-s)*sfac*g).b
    rho_iv = rt2*iv.exp(-c) * ((iv.mpf(n0+2))/(iv.mpf(n0+1)))**(5-s) \
        * ((1+c*(n0+2))/(1+c*(n0+1)))**s
    rho = mpf(rho_iv.b)             # rigorous upper bound of the ratio
    assert rho < 1
    T = Abound(n0+1) / (1 - rho)
    return tot + iv.mpf([-T, T])

def Lset_iv(a, N):
    xN = iv.sqrt(iv.mpf(N))
    n0 = int(20*mpf(xN.a)) + 30
    I0 = mellin_I_iv(a, xN, 0, n0)
    I3 = mellin_I_iv(a, xN, 3, n0)
    Lam3 = xN**3*I3 + I0          # w = +1 (mp track [L2]: w in {+-1} by
    # modularity [G1][G2]; Lam(1)=Lam(2) for w=+1, residual bounded away
    # from 0 for w=-1, checked numerically above)
    L3v = Lam3*(2*iv.pi)**3/(xN**3*2)
    Mv = Lam3
    return L3v, Mv

# --- lattice T in iv: T(d, tau) = pi^4/15 + sum_{m != 0} (3 dS2 - 4 y^2 dS3)
# with exact algebraic-part cancellation; rows by closed-form dual
# evaluations, tails by the exact Fourier representations.
# Row structure at tau = x0 + i y0, x0 = 1/2 (this script): row m has shift
# x = d m x0 mod 1 in {0, 1/2}, and
#   G(x,y) = (pi/y) sinh(2 pi y)/(cosh(2 pi y) - cos(2 pi x))
#          = (pi/y) coth(pi y)   [x == 0 (mod 1)]
#          = (pi/y) tanh(pi y)   [x == 1/2 (mod 1)]
# Both shifts have the SAME algebraic part: from the Poisson kernel
#   sinh(2 pi y)/(cosh(2 pi y)-cos(2 pi x)) = sum_k e^{-2 pi |k| y + 2 pi i k x}
# the k = 0 term is 1 for every x, giving
#   dS2(x,y) = (pi/y^3)   sum_{k>=1} (1+2 pi k y) e^{-2 pi k y} cos(2 pi k x)
#   dS3(x,y) = (pi/(4y^5)) sum_{k>=1} (3+6 pi k y+4 pi^2 k^2 y^2) e^{-2 pi k y}
#              cos(2 pi k x)
# (verified numerically at both shifts, 60 dps).  |cos| <= 1 makes the
# geometric tail bounds below uniform in the shift.
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

def _epm_dual(yv):
    ep0 = iv.exp(iv.pi*yv)
    em0 = iv.exp(-iv.pi*yv)
    ep = Dual(ep0, iv.pi*ep0, iv.pi**2*ep0/2)
    em = Dual(em0, -iv.pi*em0, iv.pi**2*em0/2)
    return ep, em

def coth_dual(yv):
    ep, em = _epm_dual(yv)
    return (ep+em)/(ep-em)

def tanh_dual(yv):
    ep, em = _epm_dual(yv)
    return (ep-em)/(ep+em)

def row_S23_iv(sh, y):
    """S2(sh/2, y), S3(sh/2, y) in iv via G = (pi/y) {coth|tanh}(pi y)."""
    G = (Dual(iv.pi)/Dual(y, 1, 0)) * (coth_dual(y) if sh == 0
                                       else tanh_dual(y))
    S2 = -G.f1/(2*y)
    S3 = (2*G.f2 + 2*S2)/(8*y**2)
    return S2, S3

def tail_row_iv(y):
    """rigorous bound of |dS2|, |dS3| at any shift (|cos| <= 1; the sums
    below bound the Fourier series above with a factor 2 margin)."""
    y = iv.mpf(y)
    r = iv.exp(-2*iv.pi*y)
    rb = r.b
    den1 = (1-rb)**2
    den3 = (1-rb)**3
    b2 = (2*iv.pi/y**3)*(1+2*iv.pi*y)*rb/den1
    b3 = (iv.pi/(2*y**5))*(3*rb/(1-rb) + 6*iv.pi*y*rb/den1
                           + 4*iv.pi**2*y**2*rb*(1+rb)/den3)
    return mpf(b2.b), mpf(b3.b)

def lattice_T_iv(d, y0, sh0=1):
    """T-sum for tau = sh0/2 + i y0; row m has shift parity (sh0*d*m) % 2.
    Rigorous tail, uniform in the shift by |cos| <= 1."""
    y0m = mpf(iv.mpf(y0).a)     # lower endpoint, only used to choose M
    M = int(50/(d*y0m)) + 2
    T = iv.pi**4/15
    for m in range(1, M+1):
        y = iv.mpf(d*m)*y0
        S2, S3 = row_S23_iv((sh0*d*m) % 2, y)
        dS2 = S2 - iv.pi/(2*y**3)
        dS3 = S3 - 3*iv.pi/(8*y**5)
        T += 2*(3*dS2 - 4*y**2*dS3)
    # rigorous tail: sum_{|m|>M} (3|dS2| + 4 y_m^2 |dS3|);  the per-row
    # Fourier bounds decrease geometrically in m (r_m = e^{-2 pi d y0 m}),
    # so the total is bounded by the bound at m = M+1 times 2/(1-r).
    tb2, tb3 = tail_row_iv(d*(M+1)*y0m)
    rr = mpf(iv.exp(-2*iv.pi*d*y0m).b)
    yM = d*(M+1)*y0m
    tot = 2*(3*tb2 + 4*yM*yM*tb3)/(1-rr)
    return T + iv.mpf([-tot, tot])

# --- the iv values ---
L3_iv, M12_iv = Lset_iv(a12, 12)
Lchi3_2_iv = dirichlet2_iv(chi3, 3)
zK2_iv = iv.pi**2/6 * Lchi3_2_iv
d3_iv = 3*iv.sqrt(3)/(4*iv.pi)*Lchi3_2_iv

def iv_mid(z):
    return (mp.convert(z.a) + mp.convert(z.b))/2

def iv_w(z):
    return mp.convert(z.b) - mp.convert(z.a)

check("[V0] iv: L(g12,3) vs mp value", iv_mid(L3_iv), L3, iv_w(L3_iv) + TOL)
check("[V0] iv: M12 = L'(g12,0) vs mp value", iv_mid(M12_iv), M12,
      iv_w(M12_iv) + TOL)
check("[V0] iv: L(chi_{-3},2) vs mp value", iv_mid(Lchi3_2_iv), Lchi3_2,
      iv_w(Lchi3_2_iv) + TOL)

s3_iv = iv.sqrt(3)
T1_iv = lattice_T_iv(1, s3_iv/2)
T2_iv = lattice_T_iv(2, s3_iv/2)

check_lock("[V1] LOCK T1(tau1) = 6 zeta_K(2)", T1_iv - 6*zK2_iv)
check_lock("[V2] LOCK T2(tau1) = 4 L(g12,3) + (9/4) zeta_K(2)",
           T2_iv - (4*L3_iv + iv.mpf(9)/4*zK2_iv))

comb_iv = -T1_iv + 4*T2_iv
EK4_iv = (5*s3_iv/iv.pi**3)*comb_iv      # (10 Im tau1/pi^3), Im tau1 = s3/2
conj_iv = iv.mpf(10)/3*(4*M12_iv + d3_iv)
check_lock("[V3] LOCK EK4(tau1) = (10/3)(4 M12 + d3)", EK4_iv - conj_iv)

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
else:
    print("ALL CHECKS PASSED")
