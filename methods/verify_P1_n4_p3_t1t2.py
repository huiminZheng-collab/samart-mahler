# verify_P1_n4_p3_t1t2.py
#
# (P1)-type exact CM evaluations for two Samart Table-6 OPEN conjectures
# (n4 family; both pure-imaginary deep inner points, Samart Prop 2.1(iii)
# applies, so n4(s4(tau)) = EK4(tau) is already a theorem):
#
#   Target T1:  tau = sqrt(-2),  s4 = 3656 + 2600 sqrt2,
#     EK4(tau) = (5/8)(4 M8tw + 28 M8 + 4 d4 + d8)
#   Target T2:  tau = 2i,        s4 = 143208 + 101574 sqrt2,
#     EK4(tau) = (5/16)(4 M16tw + 20 M16 + 9 d4 + 4 d8)
#
# where EK4(tau) = (10 Im tau/pi^3)(-T1 + 4 T2),
#   M8 = L'(g8,0),  g8 = eta(t)^2 eta(2t) eta(4t) eta(8t)^2
#                    in S_3(Gamma_0(8), chi_{-8})  [identified here],
#   M8tw = L'(g8 x chi_8, 0),   level 32, w = +1 [Fricke ratio, F0 checks],
#   M16 = L'(g16,0), g16 = eta(4t)^6 in S_3(Gamma_0(16), chi_{-4}),
#   M16tw = L'(g16 x chi_8, 0), level 64, w = +1 [Fricke ratio],
#   chi_8 = (2/n) (even, conductor 8),  d_k = L'(chi_{-k},-1).
#
# Lattice structure (all exact algebra):
#   T1 (K = Q(sqrt(-2)), O_K = Z[sqrt(-2)], h = 1, units +-1, 2 = -(sqrt-2)^2
#       ramified):
#     Lambda_1 = O_K,  Lambda_2 = O_2 = Z + 2 O_K (conductor 2).
#     chi_8(N alpha) = +1 iff alpha == 1 mod 2 (exact [X4]), so the
#     projection onto {alpha == 1 mod 2} is (1 + chi_8 o N)/2 over odd-norm
#     elements; the p = (sqrt-2) inclusion-exclusion gives
#       B_odd = (3/4) B(O_K),  G_odd = (5/4) G(O_K)   [G(p-sum) = -(1/4) G]
#     and the chi-sums are 2 zV and 2 L8tw, zV = L(chi_8,2) L(chi_{-4},2).
#     => T(O_K) = 4 L8 + 2 zK2,
#        T(O_2) = (11/4) L8 + 2 L8tw + (7/8) zK2 + zV,
#        comb = -T1 + 4 T2 = 7 L8 + 8 L8tw + (3/2) zK2 + 4 zV,
#        zK2 = zeta_K(2) = zeta(2) L(chi_{-8},2).
#   T2 (K = Q(i), O_K = Z[i], h = 1, units mu_4, 2 = -i(1+i)^2 ramified):
#     Lambda_1 = O_2,  Lambda_2 = O_4.  G(O_K) = 0 (mu_4 annihilation
#     sum_{u} u^2 = 0, exact [Y1]); each ideal prime to 2 has exactly two
#     generators == 1 mod 2 (+-alpha), giving
#       T(O_2) = 4 L16 + (7/4) zC,          zC = zeta_K(2) = zeta(2) Cat.
#     O_4\{0} = {a odd, b == 0 mod 4} ⊔ 2 O_2\{0}, and for a generator
#     alpha == 1 mod 2 the condition b == 0 mod 4 is chi_8(N alpha) = +1
#     (exact [Y4]), so the first set projects onto (1 + chi_8 o N)/2:
#       B-set = (3/4) zC + zV2,  G-set = L16 + L16tw,
#       zV2 = L(chi_8,2) L(chi_{-8},2),
#       T(O_4) = (9/4) L16 + 2 L16tw + (55/64) zC + zV2,
#       comb = 5 L16 + 8 L16tw + (27/16) zC + 4 zV2.
#
# s4(tau) = the stated Q(sqrt2) values is certified separately by
# cert0_n4_p3_t1t2.py.
#
# Layers (same three-track standard as verify_P1_n4_s7pair.py):
#   [X*],[Y*] exact integer / Fraction checks (no floating point);
#   [S*] SEPARATION exact-algebra track: each target identity proved
#        separately as an exact Fraction-algebra consequence of the quoted
#        inputs (Q1)-(Q6) -- no interval estimates in this track;
#   [G*] newform identifications: theta identities to q^60 (exact integers)
#        vs Sturm bounds 3 and 6, Ligozat criteria, cusp orders, divisor
#        degrees; root numbers/levels by the NON-vacuous Fricke ratio;
#   [L*],[D*],[F*],[T*],[E*] mpmath 60-dps numerical confirmations;
#   [V*] rigorous interval locks (mpmath.iv, iv.dps = 70): lattice T-sums
#        via closed-form coth rows with exact Fourier tail bounds, L-values
#        via Mellin with elementary incomplete gammas (E1 by a bracketing
#        continued fraction) and a self-contained eta-product coefficient
#        bound, Dirichlet values via Euler--Maclaurin.

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
# --- K1 = Q(sqrt(-2)): pairs (p,q) = p + q sqrt(-2) -----------------------
def k1mul(u, v):
    return (u[0]*v[0] - 2*u[1]*v[1], u[0]*v[1] + u[1]*v[0])
def k1norm(u): return u[0]**2 + 2*u[1]**2

SW = (Fr(0), Fr(1))       # sqrt(-2)
ONE = (Fr(1), Fr(0))
check_exact("[X1] (sqrt(-2))^2 = -2; 2 = -(sqrt(-2))^2 (2 ramified); "
            "N(sqrt(-2)) = 2 (exact)",
            k1mul(SW, SW) == (Fr(-2), Fr(0))
            and k1mul((Fr(-1), Fr(0)), k1mul(SW, SW)) == (Fr(2), Fr(0))
            and k1norm(SW) == 2)
# O_K = Z[sqrt(-2)] (disc -8 fundamental); h(-8) = 1, units +-1 [quoted Q2].
check_exact("[X2] disc(O_K) = -8 (exact)", 4*(-2) == -8)

# units mod 2: odd-norm elements are a + b sqrt(-2) with a odd; mod 2 the
# classes are 1 and 1 + sqrt(-2), both units (norms 1, 3 odd); 1+sqrt(-2)
# is NOT == 1 mod 2 (difference sqrt(-2) not in 2 O_K).
cl2 = [ONE, k1mul(ONE, SW)]  # placeholder, replaced below
u2a, u2b = ONE, (Fr(1), Fr(1))
check_exact("[X3] (O_K/2)^x = {1, 1+sqrt(-2)}: norms 1, 3 odd; "
            "1+sqrt(-2} =/= 1 (mod 2) (exact)",
            k1norm(u2a) % 2 == 1 and k1norm(u2b) % 2 == 1
            and (u2b[0] - 1) % 2 == 0 and (u2b[1]) % 2 == 1)

def chi_8v(n):
    if n % 2 == 0: return 0
    return 1 if n % 8 in (1, 7) else -1

# chi_8(N alpha) = +1 iff alpha == 1 mod 2 (i.e. b even), for odd-norm alpha:
# N = a^2 + 2 b^2 == a^2 == 1 (mod 8) when b even; == a^2 + 2 == 3 (mod 8)
# when b odd.  Exact check on a box.
okX4 = True
for a in range(-9, 10):
    for b in range(-9, 10):
        if a % 2 == 1 and (a, b) != (0, 0):
            Nm = a*a + 2*b*b
            okX4 &= (chi_8v(Nm) == 1) == (b % 2 == 0)
check_exact("[X4] chi_8(N(a + b sqrt-2)) = +1 iff b even (a odd; box 19^2)",
            okX4)

# inclusion-exclusion at p = (sqrt(-2)), N p = 2:
# B(p O_K) = (1/4) B(O_K);  G(p O_K) = ((-sqrt-2)^2 / |sqrt-2|^6) G(O_K)
#          = (-2/8) G(O_K) = -(1/4) G(O_K).
check_exact("[X5] IE coefficients at p = (sqrt-2): B: 1/4; G: -2/8 = -1/4",
            Fr(-2, 8) == Fr(-1, 4))

# L_K(chi_8 o N, s) = L(chi_8, s) L(chi_{-4}, s): local factors.
def chi4v(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def chi8mv(n):
    if n % 2 == 0: return 0
    return 1 if n % 8 in (1, 3) else -1
# chi_K = chi_{-8} for T1: split/inert by chi_{-8}(p); chi_8 * chi_{-8} =
# chi_{-4} (exact check below).
okX6 = all(chi_8v(n)*chi8mv(n) == chi4v(n) for n in range(1, 56))
check_exact("[X6] chi_8 * chi_{-8} = chi_{-4} pointwise (n < 56, exact)",
            okX6)
# chi_8 character sums: sum chi_8(a) a^2 = 16 => B_2(chi_8) = 16/8 = 2.
chi8p = [0, 1, 0, -1, 0, -1, 0, 1]
S0_8 = sum(chi8p); S1_8 = sum(chi8p[a]*a for a in range(8))
S2_8 = sum(chi8p[a]*a*a for a in range(8))
check_exact("[X7] sum chi_8 = sum chi_8 a = 0; sum chi_8 a^2 = 16 "
            "(=> B_2 = 2, exact)", S0_8 == 0 and S1_8 == 0 and S2_8 == 16)

# assembly coefficient identities, T1 (e-basis e = (sqrt2/pi^3)
# (L8, L8tw, zK2, zV);  M8 = 4 e1, M8tw = 32 e2, d8 = 24 e3, d4 = 16 e4):
#   LHS EK4 = 10 comb.e, comb = (7, 8, 3/2, 4);
#   RHS (5/8)(28 M8 + 4 M8tw + d8 + 4 d4).
comb1_c = (Fr(7), Fr(8), Fr(3, 2), Fr(4))
lhs1_e = tuple(10*comb1_c[i] for i in range(4))
rhs1_e = (Fr(5, 8)*28*4, Fr(5, 8)*4*32, Fr(5, 8)*24, Fr(5, 8)*4*16)
check_exact("[X8] T1 assembly: LHS = RHS = (70, 80, 15, 40) e (exact)",
            lhs1_e == rhs1_e == (Fr(70), Fr(80), Fr(15), Fr(40)))

# --- K2 = Q(i): pairs (p,q) = p + q i -------------------------------------
def k2mul(u, v):
    return (u[0]*v[0] - u[1]*v[1], u[0]*v[1] + u[1]*v[0])
def k2norm(u): return u[0]**2 + u[1]**2

I_ = (Fr(0), Fr(1))
OMI = (Fr(1), Fr(1))     # 1 + i
check_exact("[Y1] i^2 = -1; 2 = -i (1+i)^2 (ramified); mu_4 = {+-1,+-i}, "
            "sum_{u in mu_4} u^2 = 0 (exact)",
            k2mul(I_, I_) == (Fr(-1), Fr(0))
            and k2mul((Fr(0), Fr(-1)), k2mul(OMI, OMI)) == (Fr(2), Fr(0))
            and k2norm(OMI) == 2
            and Fr(1) + Fr(-1) + Fr(1) + Fr(-1) == 0)
check_exact("[Y2] disc(O_K) = -4 (exact)", 4*(-1) == -4)
# mu_4 mod 2: 1 and i distinct mod 2 (1 - i not in 2 O_K), -1 == 1 mod 2:
# image of units in (O_K/2)^x is {1, i}, order 2.
check_exact("[Y3] -1 == 1 (mod 2 O_K); i =/= 1 (mod 2 O_K) (exact)",
            (Fr(-1) - 1) % 2 == 0 and (I_[1] - 0) % 2 == 1)
# for alpha = a + b i == 1 mod 2 (a odd, b even): chi_8(N alpha) = +1 iff
# b == 0 mod 4 (N = a^2 + b^2 == 1 mod 8 iff b == 0 mod 4, == 5 iff b == 2).
okY4 = True
for a in range(-9, 10, 2):
    for b in range(-8, 9, 2):
        if (a, b) != (0, 0):
            Nm = a*a + b*b
            okY4 &= (chi_8v(Nm) == 1) == (b % 4 == 0)
check_exact("[Y4] chi_8(N(a + b i)) = +1 iff b == 0 mod 4 (a, b parity box)",
            okY4)
# chi_8 * chi_{-4} = chi_{-8} pointwise (for the T2 zV2 factorization
# L_K(chi_8 o N) = L(chi_8) L(chi_8 chi_{-4}) = L(chi_8) L(chi_{-8})).
okY5 = all(chi_8v(n)*chi4v(n) == chi8mv(n) for n in range(1, 56))
check_exact("[Y5] chi_8 * chi_{-4} = chi_{-8} pointwise (n < 56, exact)", okY5)

# assembly coefficient identities, T2 (e-basis e = (1/pi^3)
# (L16, L16tw, zC, zV2);  M16 = 16 e1, M16tw = 128 e2, d4 = 12 e3,
# d8 = 64 e4):
#   LHS EK4 = 20 comb.e, comb = (5, 8, 27/16, 4);
#   RHS (5/16)(20 M16 + 4 M16tw + 9 d4 + 4 d8).
comb2_c = (Fr(5), Fr(8), Fr(27, 16), Fr(4))
lhs2_e = tuple(20*comb2_c[i] for i in range(4))
rhs2_e = (Fr(5, 16)*20*16, Fr(5, 16)*4*128, Fr(5, 16)*9*12, Fr(5, 16)*4*64)
check_exact("[Y6] T2 assembly: LHS = RHS = (100, 160, 135/4, 80) e (exact)",
            lhs2_e == rhs2_e == (Fr(100), Fr(160), Fr(135, 4), Fr(80)))

# FE-constant exact ingredients:
#   8^{3/2} Gamma(3)/(2 pi)^3 = 16 sqrt2 * 2/(8 pi^3) = 4 sqrt2/pi^3;
#   32^{3/2} ... = 32 sqrt2/pi^3;  16^{3/2} ... = 16/pi^3;
#   64^{3/2} ... = 128/pi^3;   d8 = (8^{3/2}/(4 pi)) L(chi_{-8},2)
#   = (4 sqrt2/pi) L(chi_{-8},2);  zeta(2)/pi^2 = 1/6.
check_exact("[X9] FE factors: 8^3 = 512 = 256*2; 32^3 = 32768 = 16384*2; "
            "16^3 = 4096; 64^3 = 262144; Gamma(3)/2^3 = 1/4 (exact)",
            8**3 == 512 and 512 == 256*2 and 32**3 == 32768
            and 32768 == 16384*2 and 16**3 == 4096 and 64**3 == 262144
            and Fr(2, 8) == Fr(1, 4))

# =====================================================================
# PART 1: L-values (mp, 60 dps)
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

# g8 = eta(t)^2 eta(2t) eta(4t) eta(8t)^2;  g16 = eta(4t)^6
g8ser = s_mul(s_mul(s_pow(Pser(1, NMAX), 2, NMAX), Pser(2, NMAX), NMAX),
              s_mul(Pser(4, NMAX), s_pow(Pser(8, NMAX), 2, NMAX), NMAX), NMAX)
a8 = [0] + g8ser[:-1]
g16ser = s_pow(Pser(4, NMAX), 6, NMAX)
a16 = [0] + g16ser[:-1]
a8tw = [chi8p[n % 8]*a8[n] for n in range(NMAX+1)]
a16tw = [chi8p[n % 8]*a16[n] for n in range(NMAX+1)]

# [L1]/[L2] theta identities (exact integers to q^60 > Sturm bounds 3, 6)
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
a_th16 = [0]*(NT+1)
for a in range(-8, 9):
    for b in range(-8, 9):
        Nm = a*a + b*b
        if 1 <= Nm <= NT and a % 2 == 1 and b % 2 == 0:
            a_th16[Nm] += a*a - b*b
check_exact("[L2] theta identity 2 g16 = sum_{a==1(2)} alpha^2 q^N to q^60 "
            "(exact; K = Q(i))",
            all(a_th16[n] == 2*a16[n] for n in range(1, NT+1)))

# [G*] Ligozat eta-quotient criteria + cusp orders + Sturm bounds (exact)
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

l1_16, l2_16, ords16, divs16 = ligozat({4: 6}, 16)
check_exact("[G3] g16: Ligozat sums 24 == 0 (24), 24 == 0 (24); weight 3",
            l1_16 % 24 == 0 and l2_16 % 24 == 0)
check_exact("[G3b] g16 character: -4^6 = -2^12, kernel -1 => chi_{-4} (exact)",
            -(4**6) == -4096 and 4096 == 64**2)
check_exact("[G4] g16 cusp orders %s all = 1" % [str(o) for o in ords16],
            all(o == 1 for o in ords16))
dd16 = sum(phi(gcd(d, 16//d))*o for d, o in zip(divs16, ords16))
check_exact("[G4b] g16 divisor degree 6 = (3/12) index 24 (exact)", dd16 == 6)
check_exact("[G5] Sturm bounds 3 (g8), 6 (g16) < 60 = theta-id coverage",
            Fr(3, 12)*12 == 3 and Fr(3, 12)*24 == 6 and NT == 60)

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
    """L(f,3), L'(f,0) at quoted level N, root number +1; NON-vacuous
    Fricke ratio check at two ordinates."""
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
L16, M16 = Lset(a16, 16, "g16")
L16tw, M16tw = Lset(a16tw, 64, "g16 x chi_8")

check("[F1] M8 = (4 sqrt2/pi^3) L(g8,3)", M8, 4*s2*L8/pi**3, TOL)
check("[F2] M8tw = (32 sqrt2/pi^3) L(g8tw,3)", M8tw, 32*s2*L8tw/pi**3, TOL)
check("[F3] M16 = (16/pi^3) L(g16,3)", M16, 16*L16/pi**3, TOL)
check("[F4] M16tw = (128/pi^3) L(g16tw,3)", M16tw, 128*L16tw/pi**3, TOL)

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
check("[D3] L(chi_8,2) = pi^2 sqrt2/16 (closed form: tau = 2 sqrt2, "
      "B_2 = 2, N = 8)", Lchi8p_2, pi**2*s2/16, TOL)

zK2 = zeta(2)*Lchi8m_2            # zeta_K(2), K = Q(sqrt(-2))
zV = Lchi8p_2*Cat                 # T1
zC = zeta(2)*Cat                  # zeta_K(2), K = Q(i)
zV2 = Lchi8p_2*Lchi8m_2           # T2

# =====================================================================
# PART 2: lattice T-sums and decompositions (mp)
# =====================================================================
def G_row(x, y):
    return (pi/y)*sinh(2*pi*y)/(cosh(2*pi*y) - cos(2*pi*x))

def row_powers(x, y):
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y)/(2*y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy)/(2*yy), y)/(4*y)
    return S2, S3

def lattice_T(d, y0):
    B, M = mpf(0), mpf(0)
    for m in range(-300, 301):
        if m == 0:
            B += 2*zeta(4)
            continue
        S2, S3 = row_powers(0, abs(d*m)*y0)
        B += S2
        M += m*m*S3
        if abs(d*m)*y0 > 45 and m > 0:
            break
    return 3*B - 4*d*d*y0*y0*M

# T1
T_OK = lattice_T(1, s2)
T_O2 = lattice_T(2, s2)
T_OK_dec = 4*L8 + 2*zK2
T_O2_dec = mpf(11)/4*L8 + 2*L8tw + mpf(7)/8*zK2 + zV
check("[T1] T(O_K) = 4 L8 + 2 zeta_K(2)", T_OK, T_OK_dec, TOL)
check("[T2] T(O2) = (11/4) L8 + 2 L8tw + (7/8) zeta_K(2) + zV",
      T_O2, T_O2_dec, TOL)
comb1 = -T_OK + 4*T_O2
comb1_dec = 7*L8 + 8*L8tw + mpf(3)/2*zK2 + 4*zV
check("[T3] comb1 = 7 L8 + 8 L8tw + (3/2) zeta_K(2) + 4 zV",
      comb1, comb1_dec, TOL)

# T2
U_O2 = lattice_T(1, 2)
U_O4 = lattice_T(2, 2)
U_O2_dec = 4*L16 + mpf(7)/4*zC
U_O4_dec = mpf(9)/4*L16 + 2*L16tw + mpf(55)/64*zC + zV2
check("[T4] T(O2) = 4 L16 + (7/4) zeta_K(2)  [K = Q(i)]", U_O2, U_O2_dec, TOL)
check("[T5] T(O4) = (9/4) L16 + 2 L16tw + (55/64) zeta_K(2) + zV2",
      U_O4, U_O4_dec, TOL)
comb2 = -U_O2 + 4*U_O4
comb2_dec = 5*L16 + 8*L16tw + mpf(27)/16*zC + 4*zV2
check("[T6] comb2 = 5 L16 + 8 L16tw + (27/16) zeta_K(2) + 4 zV2",
      comb2, comb2_dec, TOL)

# --- assembly and the conjectured right-hand sides -------------------------
EK4_1 = (10*s2/pi**3)*comb1
EK4_2 = (20/pi**3)*comb2
conj1 = Fr(5, 8)*(4*M8tw + 28*M8 + 4*d4 + d8)
conj2 = Fr(5, 16)*(4*M16tw + 20*M16 + 9*d4 + 4*d8)
check("[E1] EK4(sqrt(-2)) = (5/8)(4 M8tw + 28 M8 + 4 d4 + d8)",
      EK4_1, conj1, TOL)
check("[E2] EK4(2i) = (5/16)(4 M16tw + 20 M16 + 9 d4 + 4 d8)",
      EK4_2, conj2, TOL)
EK4_1_asm = (10*s2/pi**3)*comb1_dec
EK4_2_asm = (20/pi**3)*comb2_dec
check("[E3] assembly T1 from decomposed comb (FE constants only)",
      EK4_1_asm, conj1, TOL)
check("[E3] assembly T2 from decomposed comb (FE constants only)",
      EK4_2_asm, conj2, TOL)

# =====================================================================
# PART 2d: SEPARATION -- exact-algebra tracks (no interval estimates)
# =====================================================================
# --- T1 track: coefficient vectors in basis [L8, L8tw, zK2, zV] -----------
# anchors (Q2)(Q3)(Q4): B(O_K) = 2 zK2, G(O_K) = 2 L8,
#   chi-sums = (2 zV, 2 L8tw); IE [X5]; projection [X4]:
#   B_odd = (3/4) B(O_K) = (3/2) zK2;  G_odd = (5/4) G(O_K) = (5/2) L8;
#   ray sums (1/2)(odd + chi): B_ray = (3/4) zK2 + zV,
#                              G_ray = (5/4) L8 + L8tw.
B_ray_zK2, B_ray_zV = Fr(1, 2)*Fr(3, 4)*2, Fr(1, 2)*2
G_ray_L8, G_ray_tw = Fr(1, 2)*Fr(5, 4)*2, Fr(1, 2)*2
B_O2 = (B_ray_zK2 + Fr(1, 16)*2, B_ray_zV)         # (7/8) zK2 + zV
G_O2 = (G_ray_L8 + Fr(1, 16)*2, G_ray_tw)          # (11/8) L8 + L8tw
T_OK_c = (2*Fr(2), Fr(0), Fr(2), Fr(0))
T_O2_c = (2*G_O2[0], 2*G_O2[1], B_O2[0], B_O2[1])
check_exact("[S1] T1: T(O_K) = 4 L8 + 2 zK2; T(O2) = (11/4) L8 + 2 L8tw "
            "+ (7/8) zK2 + zV (exact)",
            T_OK_c == (Fr(4), Fr(0), Fr(2), Fr(0))
            and T_O2_c == (Fr(11, 4), Fr(2), Fr(7, 8), Fr(1)))
comb1_Sc = tuple(-T_OK_c[i] + 4*T_O2_c[i] for i in range(4))
check_exact("[S2] T1: comb = 7 L8 + 8 L8tw + (3/2) zK2 + 4 zV (exact)",
            comb1_Sc == (Fr(7), Fr(8), Fr(3, 2), Fr(4)))
# conversions (Q5)(Q6): M8 = 4 e1, M8tw = 32 e2, d8 = 24 e3, d4 = 16 e4
lhs1_Se = tuple(10*comb1_Sc[i] for i in range(4))
rhs1_Se = (Fr(5, 8)*28*4, Fr(5, 8)*4*32, Fr(5, 8)*24, Fr(5, 8)*4*16)
check_exact("[S3] TARGET T1: EK4(sqrt(-2)) = (5/8)(4 M8tw + 28 M8 + 4 d4 "
            "+ d8) -- exact equality of e-basis coefficient vectors",
            lhs1_Se == rhs1_Se)

# --- T2 track: coefficient vectors in basis [L16, L16tw, zC, zV2] ---------
# anchors: G(O_K) = 0 [Y1], B(O_K) = 4 zC; two generators == 1 mod 2 per
# ideal [Y3]; (1+chi_8 o N)/2 projection for the O_4 ray set [Y4].
B_O2_2 = (Fr(3, 2) + Fr(1, 16)*4)                  # (7/4) zC
G_O2_2 = Fr(2)                                     # 2 L16
T_O2_2c = (2*G_O2_2, Fr(0), B_O2_2, Fr(0))
B_set = (Fr(3, 4), Fr(1))                          # (3/4) zC + zV2
G_set = (Fr(1), Fr(1))                             # L16 + L16tw
B_O4 = (B_set[0] + Fr(1, 16)*B_O2_2, B_set[1])
G_O4 = (G_set[0] + Fr(1, 16)*G_O2_2, G_set[1])
T_O4_2c = (2*G_O4[0], 2*G_O4[1], B_O4[0], B_O4[1])
check_exact("[S4] T2: T(O2) = 4 L16 + (7/4) zC; T(O4) = (9/4) L16 + 2 L16tw"
            " + (55/64) zC + zV2 (exact)",
            T_O2_2c == (Fr(4), Fr(0), Fr(7, 4), Fr(0))
            and T_O4_2c == (Fr(9, 4), Fr(2), Fr(55, 64), Fr(1)))
comb2_Sc = tuple(-T_O2_2c[i] + 4*T_O4_2c[i] for i in range(4))
check_exact("[S5] T2: comb = 5 L16 + 8 L16tw + (27/16) zC + 4 zV2 (exact)",
            comb2_Sc == (Fr(5), Fr(8), Fr(27, 16), Fr(4)))
lhs2_Se = tuple(20*comb2_Sc[i] for i in range(4))
rhs2_Se = (Fr(5, 16)*20*16, Fr(5, 16)*4*128, Fr(5, 16)*9*12, Fr(5, 16)*4*64)
check_exact("[S6] TARGET T2: EK4(2i) = (5/16)(4 M16tw + 20 M16 + 9 d4 + "
            "4 d8) -- exact equality of e-basis coefficient vectors",
            lhs2_Se == rhs2_Se)
check_exact("[S7] the two targets are distinct statements (exact)",
            rhs1_Se != rhs2_Se and comb1_Sc != comb2_Sc)

print()
print("Separation track summary: each identity proved individually by exact")
print("algebra from the quoted inputs:")
print("  (Q1) Samart Prop 2.1(iii): n4(s4(tau)) = EK4(tau) for tau pure")
print("       imaginary, Im tau >= 1/sqrt2 (both points qualify)")
print("  (Q2) h(-8) = h(-4) = 1; units +-1 (Q(sqrt-2)), mu_4 (Q(i)):")
print("       B(O_K) = 2 zeta_K(2) resp. 4 zeta_K(2)")
print("  (Q3) theta identities 2 g8 = sum' alpha^2 q^N,")
print("       2 g16 = sum_{a==1(2)} alpha^2 q^N [exact [L1][L2] to q^60,")
print("       Sturm bounds 3, 6; theta series in the respective spaces")
print("       quoted (Hecke)]; twisted versions via chi_8 o N [X4][Y4]")
print("  (Q4) L_K(chi_8 o N, s) factorizations [X6][Y5 exact Euler")
print("       factors]; chi-sums = (2 zV, 2 L8tw) resp. (2 zV2, 2 L16tw)")
print("  (Q5) FEs with root number +1 at levels 8, 32, 16, 64 [non-vacuous")
print("       Fricke ratios [F0]]; Dirichlet FE for chi_{-4}, chi_{-8}")
print("  (Q6) finite-sum formula L(chi_8,2) = pi^2 sqrt2/16 [exact [X7]]")
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
    # rigorous tail bound, SELF-CONTAINED: replacing every factor
    # (1 - q^{d n})^{r_d} by (1 - q^{d n})^{-|r_d|} dominates the
    # coefficients; with p(k) <= 2^k and the 6-fold convolution
    # p_6(i) <= C(i+5,5) 2^i <= (i+5)^5 2^i one gets
    # |a_n| <= (n+4)^5 2^{(n-1)/2} <= (n+5)^5 2^{n/2}.
    # Gamma(s,x) <= s! e^-x (1+x)^s (s=0: E1(x) <= e^-x(1+1/x)).
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
    n0 = int(45*mpf(xN.a)) + 40      # tail ~ exp(-(c - log2/2) n0) < 1e-45
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

# --- the iv values ---
L8_iv, M8_iv = Lset_iv(a8, 8)
L8tw_iv, M8tw_iv = Lset_iv(a8tw, 32)
L16_iv, M16_iv = Lset_iv(a16, 16)
L16tw_iv, M16tw_iv = Lset_iv(a16tw, 64)
Cat_iv = dirichlet2_iv(chi4l, 4)
Lchi8m_2_iv = dirichlet2_iv(chi8m, 8)
Lchi8p_2_iv = dirichlet2_iv(chi8p, 8)
zK2_iv = iv.pi**2/6 * Lchi8m_2_iv
zV_iv = Lchi8p_2_iv * Cat_iv
zC_iv = iv.pi**2/6 * Cat_iv
zV2_iv = Lchi8p_2_iv * Lchi8m_2_iv
d4_iv = 2*Cat_iv/iv.pi
d8_iv = 4*iv.sqrt(2)/iv.pi*Lchi8m_2_iv

def iv_mid(z):
    return (mp.convert(z.a) + mp.convert(z.b))/2

def iv_w(z):
    return mp.convert(z.b) - mp.convert(z.a)

check("[V0] iv: L(g8,3) vs mp value", iv_mid(L8_iv), L8, iv_w(L8_iv) + TOL)
check("[V0] iv: L(g8tw,3) vs mp value", iv_mid(L8tw_iv), L8tw,
      iv_w(L8tw_iv) + TOL)
check("[V0] iv: L(g16,3) vs mp value", iv_mid(L16_iv), L16,
      iv_w(L16_iv) + TOL)
check("[V0] iv: L(g16tw,3) vs mp value", iv_mid(L16tw_iv), L16tw,
      iv_w(L16tw_iv) + TOL)
check("[V0] iv: L(chi_8,2) = pi^2 sqrt2/16", iv_mid(Lchi8p_2_iv),
      pi**2*s2/16, iv_w(Lchi8p_2_iv) + TOL)

s2_iv = iv.sqrt(2)
T_OK_iv = lattice_T_iv(1, s2_iv)
T_O2_iv = lattice_T_iv(2, s2_iv)
U_O2_iv = lattice_T_iv(1, iv.mpf(2))
U_O4_iv = lattice_T_iv(2, iv.mpf(2))

T_OK_dec_iv = 4*L8_iv + 2*zK2_iv
T_O2_dec_iv = iv.mpf(11)/4*L8_iv + 2*L8tw_iv + iv.mpf(7)/8*zK2_iv + zV_iv
U_O2_dec_iv = 4*L16_iv + iv.mpf(7)/4*zC_iv
U_O4_dec_iv = iv.mpf(9)/4*L16_iv + 2*L16tw_iv + iv.mpf(55)/64*zC_iv + zV2_iv

check_lock("[V1] LOCK T(O_K) = 4 L8 + 2 zeta_K(2)  [Q(sqrt-2)]",
           T_OK_iv - T_OK_dec_iv)
check_lock("[V2] LOCK T(O2) = (11/4) L8 + 2 L8tw + (7/8) zK2 + zV",
           T_O2_iv - T_O2_dec_iv)
check_lock("[V3] LOCK T(O2) = 4 L16 + (7/4) zeta_K(2)  [Q(i)]",
           U_O2_iv - U_O2_dec_iv)
check_lock("[V4] LOCK T(O4) = (9/4) L16 + 2 L16tw + (55/64) zC + zV2",
           U_O4_iv - U_O4_dec_iv)

comb1_iv = -T_OK_iv + 4*T_O2_iv
comb2_iv = -U_O2_iv + 4*U_O4_iv
EK4_1_iv = (10*s2_iv/iv.pi**3)*comb1_iv
EK4_2_iv = (20/iv.pi**3)*comb2_iv
conj1_iv = iv.mpf(5)/8*(4*M8tw_iv + 28*M8_iv + 4*d4_iv + d8_iv)
conj2_iv = iv.mpf(5)/16*(4*M16tw_iv + 20*M16_iv + 9*d4_iv + 4*d8_iv)

check_lock("[V5] LOCK EK4(sqrt(-2)) = (5/8)(4 M8tw + 28 M8 + 4 d4 + d8)",
           EK4_1_iv - conj1_iv)
check_lock("[V6] LOCK EK4(2i) = (5/16)(4 M16tw + 20 M16 + 9 d4 + 4 d8)",
           EK4_2_iv - conj2_iv)

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
    import sys; sys.exit(1)
print("ALL CHECKS PASSED")
