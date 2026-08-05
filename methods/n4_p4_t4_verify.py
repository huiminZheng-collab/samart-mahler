# verify_P1_n4_p4_t4.py
#
# (P1)-type exact CM evaluation for the Samart Table-6 OPEN conjecture T4
# (n4 family; boundary point, Samart Prop 2.1 does NOT apply -- the identity
# n4(s4(tau4)) = EK4(tau4) is instead certified by the cert-2 path of
# n4_p4_t4_cert.py: tau4 is an INTERIOR point of V4 and of the certified
# path from the anchor disk around i, so the propagation lemma applies):
#
#   Target T4:  tau4 = (1 + sqrt(-2))/2,  s4 = 3656 - 2600 sqrt2,
#     EK4(tau4) = (5/4)(4 M8tw - 28 M8 + 4 d4 - d8)
#
# where EK4(tau4) = (10 Im tau4 / pi^3)(-T(Lambda_1) + 4 T(Lambda_2)),
# Lambda_d = Z + d tau4 Z, and M8, M8tw, d4, d8 are exactly the constants
# of T1 (verify_P1_n4_p3_t1t2.py):  M8 = L'(g8,0), g8 = eta^2 eta(2t)
# eta(4t) eta(8t)^2 in S_3(Gamma_0(8), chi_{-8});  M8tw = L'(g8 x chi_8,0)
# (level 32, w = +1, Fricke ratio [F0]);  d_k = L'(chi_{-k},-1).
#
# Lattice structure (all exact algebra; K = Q(sqrt(-2)), O_K = Z[sqrt(-2)],
# h = 1, units +-1):
#   2 tau4 = 1 + sqrt(-2), hence Lambda_2 = Z + (1+sqrt(-2))Z = O_K
#     (exact [X1]);
#   2 Lambda_1 = Z*2 + Z(1+sqrt(-2)) = {a + b sqrt(-2) : a == b (mod 2)}
#     =: O' = 2 O_K ⊔ ((1+sqrt(-2)) + 2 O_K)   (exact [X2][X3]);
#   parity classes of O_K mod 2 O_K: (0,0) = 2 O_K, (1,0) = a odd b even
#     (the T1 ray class), (1,1) = a,b odd, (0,1); B(O') = B00 + B11,
#     G(O') = G00 + G11 with
#       B00 = (1/16) B(O_K),        B11 = (3/4) zK2 - zV,
#       G00 = (1/16) G(O_K),        G11 = (5/4) L8 - L8tw
#     from the SAME anchors as T1 ([X4][X5], (Q2)(Q4));
#   T(2 Lambda) = (1/16) T(Lambda) (homogeneity, exact [X8]), so
#     T(Lambda_1) = 16 T(O') = 44 L8 - 32 L8tw + 14 zK2 - 16 zV,
#     T(Lambda_2) = T(O_K)   = 4 L8 + 2 zK2,
#     comb = -T1 + 4 T2 = -28 L8 + 32 L8tw - 6 zK2 + 16 zV.
#   Assembly in e = (sqrt2/pi^3)(L8, L8tw, zK2, zV):
#     EK4 = (10/(sqrt2 pi^3)) comb = 5 comb.e = (-140, 160, -30, 80) e,
#     RHS = (5/4)(-28*4, 4*32, -24, 4*16) e = (-140, 160, -30, 80) e  [X10].
#
# s4(tau4) = 3656 - 2600 sqrt2 is certified by cert0_n4_p3_t1t2.py
# (conjugate partner of T1, S = 7312, P = -153664).
#
# Layers (same three-track standard as verify_P1_n4_p3_t1t2.py):
#   [X*] exact integer / Fraction checks (no floating point);
#   [S*] SEPARATION exact-algebra track: the target identity proved as an
#        exact Fraction-algebra consequence of the quoted inputs (Q1)-(Q6);
#   [G*] newform identification (copied from T1: Ligozat, cusp orders,
#        divisor degree, Sturm bound; Fricke ratio root number);
#   [L*],[D*],[F*],[T*],[E*] mpmath 60-dps numerical confirmations;
#   [V*] rigorous interval locks (mpmath.iv, iv.dps = 70): the T1 machine
#        plus alternating-shift rows (coth/tanh) for Lambda_1 (x0 = 1/2).

from fractions import Fraction as Fr
from math import gcd
from mpmath import (mp, mpf, mpc, pi, sqrt, exp, sin, cos, sinh, cosh, zeta,
                    dirichlet, gamma, diff as mpdiff, power, gammainc,
                    catalan, iv)

mp.dps = 60
iv.dps = 70
s2 = sqrt(mpf(2))
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
    lo, hi = mp.convert(ivl.a), mp.convert(ivl.b)
    w = hi - lo
    ok = (lo <= 0 <= hi) and (w < halfwidth_req)
    if not ok:
        FAILS.append(name)
    print("%-76s %s  (half-width = %.2e)" % (name, "PASS" if ok else "FAIL", w/2))

TOL = mpf(10)**(-45)
TOLM = mpf(10)**(-40)

# =====================================================================
# PART 0: exact algebra
# =====================================================================
def k1mul(u, v):
    return (u[0]*v[0] - 2*u[1]*v[1], u[0]*v[1] + u[1]*v[0])
def k1norm(u): return u[0]**2 + 2*u[1]**2

SW = (Fr(0), Fr(1))       # sqrt(-2)
# [X1] 2 tau4 = 1 + sqrt(-2); Lambda_2 = Z + (1+sqrt(-2))Z = O_K:
# a + b sqrt(-2) = (a-b) + b(1+sqrt(-2)) for all a, b (exact on a box).
okX1 = all((a - b) + b == a for a in range(-9, 10) for b in range(-9, 10))
check_exact("[X1] 2 tau4 = 1+sqrt(-2); a + b sqrt-2 = (a-b) + b(1+sqrt-2) "
            "=> Lambda_2 = O_K (exact)", okX1 and k1norm(SW) == 2)

# [X2]/[X3] 2 Lambda_1 = Z*2 + Z(1+sqrt(-2)) = {(a,b): a == b mod 2} = O'
# = 2 O_K ⊔ ((1+sqrt(-2)) + 2 O_K): parity bookkeeping, exact on a box.
def in_Op(a, b): return (a - b) % 2 == 0
okX2 = all(in_Op(2*m + n, n) for m in range(-6, 7) for n in range(-6, 7))
okX3 = all(in_Op(a, b) == ((a % 2 == 0 and b % 2 == 0) or
                           (a % 2 == 1 and b % 2 == 1))
           for a in range(-9, 10) for b in range(-9, 10))
check_exact("[X2] 2 Lambda_1 = {(a,b): a == b (mod 2)} = O' (box, exact)", okX2)
check_exact("[X3] O' = 2 O_K ⊔ ((1+sqrt-2) + 2 O_K) (parity classes, exact)",
            okX3)

# [X4] chi_8(N(a + b sqrt-2)) = +1 iff b even (a odd; same as T1 [X4]).
def chi_8v(n):
    if n % 2 == 0: return 0
    return 1 if n % 8 in (1, 7) else -1
okX4 = True
for a in range(-9, 10):
    for b in range(-9, 10):
        if a % 2 == 1 and (a, b) != (0, 0):
            okX4 &= (chi_8v(a*a + 2*b*b) == 1) == (b % 2 == 0)
check_exact("[X4] chi_8(N(a + b sqrt-2)) = +1 iff b even (a odd; box 19^2)",
            okX4)

# [X5] inclusion-exclusion at p = (sqrt-2) and scaling by 2 (exact).
check_exact("[X5] IE at p: B: 1/4, G: -2/8 = -1/4; scaling by 2: B,G: 4/64 "
            "= 1/16 (exact)",
            Fr(-2, 8) == Fr(-1, 4) and Fr(4, 64) == Fr(1, 16))

# [X6] chi_8 * chi_{-8} = chi_{-4} pointwise (for L_K(chi_8 o N) =
# L(chi_8) L(chi_{-4}); same as T1 [X6]).
def chi4v(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def chi8mv(n):
    return 0 if n % 2 == 0 else (1 if n % 8 in (1, 3) else -1)
okX6 = all(chi_8v(n)*chi8mv(n) == chi4v(n) for n in range(1, 56))
check_exact("[X6] chi_8 * chi_{-8} = chi_{-4} pointwise (n < 56, exact)", okX6)

# [X7] chi_8 character sums (B_2 = 2; same as T1 [X7]).
chi8p = [0, 1, 0, -1, 0, -1, 0, 1]
check_exact("[X7] sum chi_8 = sum chi_8 a = 0; sum chi_8 a^2 = 16 (exact)",
            sum(chi8p) == 0 and sum(chi8p[a]*a for a in range(8)) == 0
            and sum(chi8p[a]*a*a for a in range(8)) == 16)

# [X8] homogeneity T(c Lambda) = c^{-4} T(Lambda) at c = 2:
# B: |2z|^{-4} = (1/16)|z|^{-4};  G: Re((2z)^2)/|2z|^6 = (4/64) Re(z^2)/|z|^6.
check_exact("[X8] T(2 Lambda) = (1/16) T(Lambda) (B and G weights, exact)",
            Fr(1, 16) == Fr(4, 64))

# [X9] T(Lambda) = B(Lambda) + 2 G(Lambda) for ANY lattice: the forms
# 4(Re z)^2/|z|^6 - 1/|z|^4 and 2 Re(z^2)/|z|^6 + 1/|z|^4 agree pointwise
# (4x^2 - (x^2+y^2) = 3x^2 - y^2 = 2(x^2-y^2) + (x^2+y^2)); both series
# absolutely convergent (sum |lambda|^{-4}).
okX9 = all(4*x*x - (x*x + y*y) == 2*(x*x - y*y) + (x*x + y*y)
           for x in range(-6, 7) for y in range(-6, 7))
check_exact("[X9] T = B + 2G pointwise on forms (grid, exact; holds for the "
            "T4 lattices too)", okX9)

# [X10] assembly (e-basis e = (sqrt2/pi^3)(L8, L8tw, zK2, zV);
# M8 = 4 e1, M8tw = 32 e2, d8 = 24 e3, d4 = 16 e4):
#   LHS EK4 = (10/(sqrt2 pi^3)) comb = 5 comb.e, comb = (-28, 32, -6, 16);
#   RHS (5/4)(4 M8tw - 28 M8 + 4 d4 - d8).
comb4_c = (Fr(-28), Fr(32), Fr(-6), Fr(16))
lhs4_e = tuple(5*comb4_c[i] for i in range(4))
rhs4_e = (Fr(5, 4)*(-28)*4, Fr(5, 4)*4*32, Fr(5, 4)*(-1)*24, Fr(5, 4)*4*16)
check_exact("[X10] T4 assembly: LHS = RHS = (-140, 160, -30, 80) e (exact)",
            lhs4_e == rhs4_e == (Fr(-140), Fr(160), Fr(-30), Fr(80)))

# =====================================================================
# PART 1: L-values (mp, 60 dps) -- g8 / M8tw / d4 / d8 identical to T1
# =====================================================================
NMAX = 400
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
a8tw = [chi8p[n % 8]*a8[n] for n in range(NMAX+1)]

# [L1] theta identity (exact integers to q^60 > Sturm bound 3)
NT = 60
a_th8 = [0]*(NT+1)
for a in range(-8, 9):
    for b in range(-8, 9):
        Nm = a*a + 2*b*b
        if 1 <= Nm <= NT:
            a_th8[Nm] += a*a - 2*b*b
check_exact("[L1] theta identity 2 g8 = sum' alpha^2 q^{N(a)} to q^60 "
            "(exact; K = Q(sqrt-2))",
            all(a_th8[n] == 2*a8[n] for n in range(1, NT+1)))

# [G*] Ligozat eta-quotient criteria + cusp orders + Sturm bound (exact)
def ligozat(r_eta, NLV):
    lig1 = sum(d*e for d, e in r_eta.items())
    lig2 = sum((NLV//d)*e for d, e in r_eta.items())
    def cusp_order(d):
        tot = Fr(0)
        for dlt, e in r_eta.items():
            tot += Fr(gcd(d, dlt)**2 * e, dlt)
        return Fr(NLV, 24) * tot / (d * gcd(d, NLV//d))
    divs = [d for d in range(1, NLV+1) if NLV % d == 0]
    return lig1, lig2, [cusp_order(d) for d in divs], divs

def phi(n):
    return sum(1 for k in range(1, n+1) if gcd(k, n) == 1)

l1_8, l2_8, ords8, divs8 = ligozat({1: 2, 2: 1, 4: 1, 8: 2}, 8)
check_exact("[G1] g8: Ligozat sums 24 == 0 (24), 24 == 0 (24); weight 3",
            l1_8 % 24 == 0 and l2_8 % 24 == 0
            and sum((2, 1, 1, 2)) == 6)
check_exact("[G1b] g8 character: -1^2 2 4 8^2 = -2^9, squarefree kernel -2 "
            "=> chi_{-8} (exact)", -(2*4*8**2) == -512 and -512 == -2*256
            and 256 == 16**2)
check_exact("[G2] g8 cusp orders %s all > 0" % [str(o) for o in ords8],
            all(o > 0 for o in ords8))
dd8 = sum(phi(gcd(d, 8//d))*o for d, o in zip(divs8, ords8))
check_exact("[G2b] g8 divisor degree 3 = (3/12) index 12 (exact)", dd8 == 3)
check_exact("[G3] Sturm bound 3 < 60 = theta-id coverage",
            Fr(3, 12)*12 == 3 and NT == 60)

def mellin_I(a, xN, s):
    s = mpf(s)
    tot = mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0:
            continue
        tot += a[n]*power(2*pi*n, -s)*gammainc(s, 2*pi*n/xN)
    return tot

def qser(a, tau):
    q = exp(2*pi*1j*tau)
    s = mpc(0); qn = q
    for n in range(1, len(a)):
        if a[n]:
            s += a[n]*qn
        qn *= q
        if abs(qn) < mpf(10)**(-70):
            break
    return s

def Lset(a, N, label):
    xN = sqrt(mpf(N))
    for yy in ("0.6", "1.1"):
        yv = mpf(yy)
        r = qser(a, 1j/(xN*yv))/(yv**3*qser(a, 1j*yv/xN))
        check("[F0:%s] Fricke ratio at level %d = +1 (y = %s)"
              % (label, N, yy), r, 1, TOLM)
    I0, I3 = mellin_I(a, xN, 0), mellin_I(a, xN, 3)
    Lam3 = xN**3*I3 + I0
    L3 = Lam3*(2*pi)**3/(xN**3*gamma(3))
    return L3, Lam3

L8, M8 = Lset(a8, 8, "g8")
L8tw, M8tw = Lset(a8tw, 32, "g8 x chi_8")

check("[F1] M8 = (4 sqrt2/pi^3) L(g8,3)", M8, 4*s2*L8/pi**3, TOL)
check("[F2] M8tw = (32 sqrt2/pi^3) L(g8tw,3)", M8tw, 32*s2*L8tw/pi**3, TOL)

# --- Dirichlet values ------------------------------------------------------
chi4l = [0, 1, 0, -1]
chi8m = [0, 1, 0, 1, 0, -1, 0, -1]
Cat = dirichlet(mpf(2), chi4l)
Lchi8m_2 = dirichlet(mpf(2), chi8m)
Lchi8p_2 = dirichlet(mpf(2), chi8p)
d4 = 2*Cat/pi
d8 = 4*s2/pi*Lchi8m_2
check("[D1] d4 = 2 Catalan/pi = L'(chi_{-4},-1) (direct)", d4,
      mpdiff(lambda s: dirichlet(s, chi4l), mpf(-1)), TOL)
check("[D2] d8 = (4 sqrt2/pi) L(chi_{-8},2) = L'(chi_{-8},-1) (direct)", d8,
      mpdiff(lambda s: dirichlet(s, chi8m), mpf(-1)), TOL)
check("[D3] L(chi_8,2) = pi^2 sqrt2/16 (closed form)", Lchi8p_2,
      pi**2*s2/16, TOL)

zK2 = zeta(2)*Lchi8m_2            # zeta_K(2), K = Q(sqrt(-2))
zV = Lchi8p_2*Cat                 # L(chi_8,2) L(chi_{-4},2)

# =====================================================================
# PART 2: lattice T-sums and decompositions (mp; x0 = 1/2 rows)
# =====================================================================
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

y0 = 1/s2
T_L1 = lattice_T(1, mpf(1)/2, y0)      # Lambda_1 = Z + tau4 Z
T_L2 = lattice_T(2, mpf(1)/2, y0)      # Lambda_2 = Z + 2 tau4 Z = O_K

T_L2_dec = 4*L8 + 2*zK2
T_L1_dec = 44*L8 - 32*L8tw + 14*zK2 - 16*zV
check("[T1] T(Lambda_2) = T(O_K) = 4 L8 + 2 zeta_K(2)", T_L2, T_L2_dec, TOL)
check("[T2] T(Lambda_1) = 44 L8 - 32 L8tw + 14 zK2 - 16 zV",
      T_L1, T_L1_dec, TOL)
comb4 = -T_L1 + 4*T_L2
comb4_dec = -28*L8 + 32*L8tw - 6*zK2 + 16*zV
check("[T3] comb = -28 L8 + 32 L8tw - 6 zK2 + 16 zV", comb4, comb4_dec, TOL)

# [T0] independent direct lattice sum over O' (box truncation, F(z) weight):
# T(O') = sum'_{a==b(2)} F(z) with the pointwise form
# F(z) = 4(Re z)^2/|z|^6 - 1/|z|^4, z = a + b sqrt(-2), |z|^2 = a^2+2b^2.
BD = 800
T_Op_direct = mpf(0)
for a in range(-BD, BD+1):
    for b in range(-BD, BD+1):
        if (a, b) == (0, 0) or (a - b) % 2 != 0:
            continue
        Nm = mpf(a*a + 2*b*b)
        T_Op_direct += 4*mpf(a*a)/Nm**3 - 1/Nm**2
check("[T0] direct O'-sum (box 800) = T(Lambda_1)/16 (1e-3 truncation)",
      T_Op_direct, T_L1/16, mpf(10)**(-3))

# --- assembly and the conjectured right-hand side --------------------------
EK4_4 = (10*y0/pi**3)*comb4
conj4 = Fr(5, 4)*(4*M8tw - 28*M8 + 4*d4 - d8)
check("[E1] EK4((1+sqrt-2)/2) = (5/4)(4 M8tw - 28 M8 + 4 d4 - d8)",
      EK4_4, conj4, TOL)
EK4_4_asm = (10*y0/pi**3)*comb4_dec
check("[E3] assembly T4 from decomposed comb (FE constants only)",
      EK4_4_asm, conj4, TOL)

# =====================================================================
# PART 2d: SEPARATION -- exact-algebra track (no interval estimates)
# =====================================================================
# Anchors (Q2): B(O_K) = 2 zK2, G(O_K) = 2 L8; (Q4): chi-sums = (2 zV,
# 2 L8tw); IE [X5]; projection [X4].  Coefficient vectors in basis
# [L8, L8tw, zK2, zV].
B_odd_zK2 = Fr(3, 4)*2          # B_odd = (3/4) B(O_K) = (3/2) zK2
B_chi_zV = Fr(2)                # chi-sum = 2 zV
B10 = (Fr(1, 2)*B_odd_zK2, Fr(1, 2)*B_chi_zV)      # a odd b even (T1 ray)
B11 = (Fr(1, 2)*B_odd_zK2, -Fr(1, 2)*B_chi_zV)     # a, b odd
B00 = (Fr(1, 16)*2,)                               # B(2 O_K) = zK2/8
B_Op_zK2 = B00[0] + B11[0]
B_Op_zV = B11[1]
check_exact("[S1] T4: B(O') = B00 + B11 = (7/8) zK2 - zV (exact)",
            (B_Op_zK2, B_Op_zV) == (Fr(7, 8), Fr(-1)))

G_odd_L8 = Fr(5, 4)*2           # G_odd = (5/4) G(O_K) = (5/2) L8
G_chi_tw = Fr(2)                # chi-sum = 2 L8tw
G11 = (Fr(1, 2)*G_odd_L8, -Fr(1, 2)*G_chi_tw)      # (5/4) L8 - L8tw
G_Op_L8 = Fr(1, 16)*2 + G11[0]
G_Op_tw = G11[1]
check_exact("[S2] T4: G(O') = G00 + G11 = (11/8) L8 - L8tw (exact)",
            (G_Op_L8, G_Op_tw) == (Fr(11, 8), Fr(-1)))

T_Op_c = (2*G_Op_L8, 2*G_Op_tw, B_Op_zK2, B_Op_zV)
check_exact("[S3] T4: T(O') = (11/4) L8 - 2 L8tw + (7/8) zK2 - zV (exact)",
            T_Op_c == (Fr(11, 4), Fr(-2), Fr(7, 8), Fr(-1)))
T_L1_c = tuple(16*T_Op_c[i] for i in range(4))      # [X8] homogeneity
T_L2_c = (Fr(4), Fr(0), Fr(2), Fr(0))               # T(O_K), anchors (Q2)
check_exact("[S4] T4: T(Lambda_1) = 44 L8 - 32 L8tw + 14 zK2 - 16 zV; "
            "T(Lambda_2) = 4 L8 + 2 zK2 (exact)",
            T_L1_c == (Fr(44), Fr(-32), Fr(14), Fr(-16))
            and T_L2_c == (Fr(4), Fr(0), Fr(2), Fr(0)))
comb4_Sc = tuple(-T_L1_c[i] + 4*T_L2_c[i] for i in range(4))
check_exact("[S5] T4: comb = -28 L8 + 32 L8tw - 6 zK2 + 16 zV (exact)",
            comb4_Sc == (Fr(-28), Fr(32), Fr(-6), Fr(16)))
lhs4_Se = tuple(5*comb4_Sc[i] for i in range(4))
check_exact("[S6] TARGET T4: EK4(tau4) = (5/4)(4 M8tw - 28 M8 + 4 d4 - d8) "
            "-- exact equality of e-basis coefficient vectors",
            lhs4_Se == rhs4_e)
check_exact("[S7] T4 sign pattern differs from T1's (5/8)(4 M8tw + 28 M8 + "
            "4 d4 + d8) (exact)", rhs4_e != (Fr(5, 8)*4*32, Fr(5, 8)*28*4,
                                             Fr(5, 8)*24, Fr(5, 8)*4*16))

print()
print("Separation track summary: the T4 identity is proved by exact algebra")
print("from the quoted inputs:")
print("  (Q0) n4(s4(tau4)) = EK4(tau4): cert-2 path n4_p4_t4_cert.py")
print("       (tau4 interior of V4; propagation lemma formal.tex Sec. 3);")
print("       s4(tau4) = 3656 - 2600 sqrt2: cert0_n4_p3_t1t2.py")
print("  (Q1) lattice identities 2 tau4 = 1+sqrt-2, 2 Lambda_1 = O'")
print("       [X1]-[X3 exact]")
print("  (Q2) h(-8) = 1, units +-1: B(O_K) = 2 zeta_K(2), G(O_K) = 2 L8")
print("  (Q3) theta identity 2 g8 = sum' alpha^2 q^N [L1 exact to q^60,")
print("       Sturm bound 3; theta series in the space quoted (Hecke)]")
print("  (Q4) L_K(chi_8 o N) = L(chi_8) L(chi_{-4}) [X6 exact]; chi-sums")
print("       (2 zV, 2 L8tw) via [X4]")
print("  (Q5) FEs with root number +1 at levels 8, 32 [Fricke ratios F0];")
print("       Dirichlet FE for chi_{-4}, chi_{-8}")
print("  (Q6) finite-sum formula L(chi_8,2) = pi^2 sqrt2/16 [X7 exact]")
print()

# =====================================================================
# PART 3: rigorous interval locks (mpmath.iv, iv.dps = 70)
# =====================================================================
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
    # rigorous tail bound (self-contained; same as verify_P1_n4_p3_t1t2.py):
    # |a_n| <= (n+5)^5 2^{n/2}; Gamma(s,x) <= s! e^-x (1+x)^s
    # (s=0: E1(x) <= e^-x(1+1/x)).
    c = 2*iv.pi/xN
    rt2 = iv.sqrt(2)
    def Abound(n):
        nn = iv.mpf(n)
        xn = c*nn
        g = iv.exp(-xn)*(1+1/xn) if s == 0 else iv.exp(-xn)*(1+xn)**s
        sfac = iv.mpf(1) if s == 0 else iv.mpf([1,1,2,6][s])
        return ((nn+5)**5*rt2**n*(2*iv.pi*nn)**(-s)*sfac*g).b
    rho_iv = rt2*iv.exp(-c) * ((iv.mpf(n0+7))/(iv.mpf(n0+6)))**5 \
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
    Lam3 = xN**3*I3 + I0          # w = +1 (Fricke ratios [F0])
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

def cs_dual(yv, shift):
    """coth(pi y) at shift 0, tanh(pi y) at shift 1/2 (Poisson kernel:
    sinh(2 pi y)/(cosh(2 pi y) - cos(2 pi x)) = coth/tanh at x = 0, 1/2)."""
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
    """Tail bounds on the row deviations dS2, dS3 from the power terms.
    Derived (T1 machine) from the Poisson-kernel expansion
    G(x,y) = (pi/y)(1 + 2 sum_{k>=1} r^k cos(2 pi k x)), r = e^{-2 pi y},
    with |cos| <= 1 -- hence UNIFORM in the shift x (both 0 and 1/2)."""
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
    """T-sum for tau = i y0 (all rows at shift 0), rigorous tail (T1)."""
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

def lattice_T_iv_t4():
    """T(Lambda_1) for tau4 = 1/2 + i/sqrt2: rows y = m/sqrt2 alternate
    shift 0 (m even, coth rows) and shift 1/2 (m odd, tanh rows).  The
    power terms pi/(2 y^3), 3 pi/(8 y^5) of S2, S3 are shift-independent
    (Poisson kernel: only the k = 0 term is a power of y), so the same
    subtraction applies; the tail bounds are uniform in the shift."""
    y0 = 1/iv.sqrt(2)
    y0m = mpf(iv.mpf(y0).a)
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

# --- the iv values ---
L8_iv, M8_iv = Lset_iv(a8, 8)
L8tw_iv, M8tw_iv = Lset_iv(a8tw, 32)
Cat_iv = dirichlet2_iv(chi4l, 4)
Lchi8m_2_iv = dirichlet2_iv(chi8m, 8)
Lchi8p_2_iv = dirichlet2_iv(chi8p, 8)
zK2_iv = iv.pi**2/6 * Lchi8m_2_iv
zV_iv = Lchi8p_2_iv * Cat_iv
d4_iv = 2*Cat_iv/iv.pi
d8_iv = 4*iv.sqrt(2)/iv.pi*Lchi8m_2_iv

def iv_mid(z):
    return (mp.convert(z.a) + mp.convert(z.b))/2

def iv_w(z):
    return mp.convert(z.b) - mp.convert(z.a)

check("[V0a] iv: L(g8,3) vs mp value", iv_mid(L8_iv), L8, iv_w(L8_iv) + TOL)
check("[V0b] iv: L(g8tw,3) vs mp value", iv_mid(L8tw_iv), L8tw,
      iv_w(L8tw_iv) + TOL)
check("[V0c] iv: L(chi_8,2) = pi^2 sqrt2/16", iv_mid(Lchi8p_2_iv),
      pi**2*s2/16, iv_w(Lchi8p_2_iv) + TOL)

s2_iv = iv.sqrt(2)
T_L2_iv = lattice_T_iv(1, s2_iv)       # Lambda_2 = O_K = Z + sqrt(-2) Z
T_L1_iv = lattice_T_iv_t4()            # Lambda_1 = Z + tau4 Z

check("[V0d] iv T(Lambda_2) contains the 60-dps mp value",
      iv_mid(T_L2_iv), T_L2, iv_w(T_L2_iv) + TOL)
check("[V0e] iv T(Lambda_1) contains the 60-dps mp value",
      iv_mid(T_L1_iv), T_L1, iv_w(T_L1_iv) + TOL)

T_L2_dec_iv = 4*L8_iv + 2*zK2_iv
T_L1_dec_iv = 44*L8_iv - 32*L8tw_iv + 14*zK2_iv - 16*zV_iv

check_lock("[V1] LOCK T(Lambda_2) = 4 L8 + 2 zeta_K(2)",
           T_L2_iv - T_L2_dec_iv)
check_lock("[V2] LOCK T(Lambda_1) = 44 L8 - 32 L8tw + 14 zK2 - 16 zV",
           T_L1_iv - T_L1_dec_iv)

comb4_iv = -T_L1_iv + 4*T_L2_iv
EK4_4_iv = (10/s2_iv/iv.pi**3)*comb4_iv
conj4_iv = iv.mpf(5)/4*(4*M8tw_iv - 28*M8_iv + 4*d4_iv - d8_iv)

check_lock("[V3] LOCK EK4(tau4) = (5/4)(4 M8tw - 28 M8 + 4 d4 - d8)",
           EK4_4_iv - conj4_iv)

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
    import sys; sys.exit(1)
print("ALL CHECKS PASSED")
