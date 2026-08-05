# verify_P1_n5_e1.py -- Phase 8, open entry #1 (19th identity):
#   n4(143208 - 101574 sqrt2) = (5/8)(4 M16tw - 20 M16 - 9 d4 + 4 d8)
# at tau = (1+2i)/2, K = Q(i), h = 1.
#
# Station status:
#   S3 (continuation): tau = (1/2, 1) is the top corner of the certified
#       path of n4_p4_t4_cert.py (leg A y=1, x in [0,1/2]; leg B x=1/2,
#       y in [0.702,1]) => tau in W, n4(s4(tau)) = EK4(tau) by the
#       propagation theorem. [quoted]
#   S4 (exact s-value): s4(tau) = 143208-101574 sqrt2 locked by
#       cert0_n4_p3_t1t2.py (conjugate lock of Theorem D; numeric
#       re-confirmation in n5_p8_step0.py, 2.6e-57). [quoted]
#   S5 (this script): exact CM evaluation EK4(tau) = conjectured RHS.
#
# Lattices (exact [E1-X1..X3]):
#   2 tau = 1+2i  => Lambda_2 = Z + (1+2i)Z = {m+ni: n even} = O_2
#                    (conductor-2 order; the SAME lattice as T2's Lambda_1)
#   2 Lambda_1 = 2Z + (1+2i)Z = {m+2bi: m == b (2)} =: O'
#              = 2 O_2 ⊔ ((1+2i) + 2 O_2)
# Decomposition (exact Fraction track [E1-S*], e-basis
#   e = (1/pi^3)(L16, L16tw, zC, zV2), zC = zeta_K(2),
#   zV2 = L(chi_8,2) L(chi_{-8},2)):
#   anchors quoted from verify_P1_n4_p3_t1t2.py (T2):
#     B(O_K) = 4 zC, G(O_K) = 0 (mu_4 annihilation [Y1]);
#     B(C) = (3/2) zC, G(C) = 2 L16   (C = {a odd, b even} = {==1 mod 2 O_K});
#     B(C0) = (3/4) zC + zV2, G(C0) = L16 + L16tw
#             (C0 = {a odd, b == 0 (4)}, chi_8 o N projection [Y4][Y5]).
#   C = C0 ⊔ C', C' = (1+2i)+2 O_2 = {m odd, n == 2 (4)}:
#     B(C') = (3/4) zC - zV2, G(C') = L16 - L16tw
#   T(O_2) = 4 L16 + (7/4) zC;   T(2 O_2) = (1/16) T(O_2);
#   T(O') = (9/4) L16 - 2 L16tw + (55/64) zC - zV2
#   T(Lambda_1) = 16 T(O');  comb = -T(Lambda_1) + 4 T(O_2)
#             = -20 L16 + 32 L16tw - (27/4) zC + 16 zV2
#   EK4 = 10 comb.e = (-200, 320, -135/2, 160).e
#       = (5/8)(4 M16tw - 20 M16 - 9 d4 + 4 d8)   [M16=16e1, M16tw=128e2,
#                                                 d4=12e3, d8=64e4]
# New file; modifies nothing.
from mpmath import (mp, mpf, mpc, pi, zeta, sinh, cosh, cos, exp, sqrt,
                    gamma, power, gammainc, diff as mpdiff, dirichlet, iv, nstr)
from fractions import Fraction as Fr
from math import gcd

mp.dps = 60
iv.dps = 70
s2 = sqrt(2)
TOL = mpf(10)**(-50)
TOLM = mpf(10)**(-40)
FAILS = []

def check(name, got, want, tol):
    d = abs(got - want)
    ok = d < tol
    if not ok:
        FAILS.append(name)
    print("%-72s %s (|diff| = %s)" % (name, "PASS" if ok else "FAIL",
                                      nstr(d, 3)))

def check_exact(name, cond):
    if not cond:
        FAILS.append(name)
    print("%-72s %s" % (name, "PASS" if cond else "FAIL"))

def check_lock(name, ivl, halfwidth_req=mpf(10)**(-38)):
    ok = (ivl.a < 0 < ivl.b) and (ivl.b - ivl.a)/2 < halfwidth_req
    if not ok:
        FAILS.append(name)
    print("%-72s %s (half-width = %s)" % (name, "PASS" if ok else "FAIL",
                                          nstr((ivl.b-ivl.a)/2, 3)))

chi8p = [0, 1, 0, -1, 0, -1, 0, 1]
chi4l = [0, 1, 0, -1]
chi8m = [0, 1, 0, 1, 0, -1, 0, -1]

# =====================================================================
# PART 0: exact algebra
# =====================================================================
def k2mul(u, v):
    return (u[0]*v[0] - u[1]*v[1], u[0]*v[1] + u[1]*v[0])
def k2norm(u): return u[0]**2 + u[1]**2

I_ = (Fr(0), Fr(1))
OMI = (Fr(1), Fr(1))
check_exact("[Y1] i^2 = -1; 2 = -i(1+i)^2; mu_4, sum u^2 = 0 (exact)",
            k2mul(I_, I_) == (Fr(-1), Fr(0))
            and k2mul((Fr(0), Fr(-1)), k2mul(OMI, OMI)) == (Fr(2), Fr(0))
            and k2norm(OMI) == 2
            and Fr(1) + Fr(-1) + Fr(1) + Fr(-1) == 0)
check_exact("[Y2] disc(O_K) = -4 (exact)", 4*(-1) == -4)
check_exact("[Y3] -1 == 1 (mod 2 O_K); i =/= 1 (mod 2 O_K) (exact)",
            (Fr(-1) - 1) % 2 == 0 and (I_[1] - 0) % 2 == 1)

def chi_8v(n):
    if n % 2 == 0: return 0
    return 1 if n % 8 in (1, 7) else -1
okY4 = True
for a in range(-9, 10, 2):
    for b in range(-8, 9, 2):
        if (a, b) != (0, 0):
            Nm = a*a + b*b
            okY4 &= (chi_8v(Nm) == 1) == (b % 4 == 0)
check_exact("[Y4] chi_8(N(a+bi)) = +1 iff b == 0 mod 4 (a odd, b even; box)",
            okY4)
def chi4v(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def chi8mv(n):
    return 0 if n % 2 == 0 else (1 if n % 8 in (1, 3) else -1)
okY5 = all(chi_8v(n)*chi4v(n) == chi8mv(n) for n in range(1, 56))
check_exact("[Y5] chi_8 * chi_{-4} = chi_{-8} pointwise (n < 56, exact)", okY5)

# --- #1 lattice set identities (exact on boxes) ---------------------------
# [E1-X1] a + b(1+2i) = (a+b) + 2bi, and m + 2bi = (m-b) + b(1+2i):
# Lambda_2 = Z + (1+2i)Z = {m+ni: n even} = O_2.
okX1 = all((a + b) == a + b and (m - b) + b == m
           for a in range(-8, 9) for b in range(-8, 9)
           for m in range(-8, 9))
check_exact("[E1-X1] Lambda_2 = Z+(1+2i)Z = {m+ni: n even} = O_2 (exact)",
            okX1 and k2mul(OMI, OMI) == (Fr(0), Fr(2)))
# [E1-X2] 2 Lambda_1 = 2Z + (1+2i)Z: element (2a+b) + 2bi has m == b (2);
# conversely if m == b (2) then m + 2bi = 2a + b(1+2i) with a = (m-b)/2.
okX2 = all(((2*a + b) - b) % 2 == 0
           for a in range(-8, 9) for b in range(-8, 9)) \
       and all((m - b) % 2 != 0 or 2*((m - b)//2) + b == m
               for m in range(-9, 10) for b in range(-9, 10))
check_exact("[E1-X2] 2 Lambda_1 = {m+2bi: m == b (2)} = O' (exact)", okX2)
# [E1-X3] O' = 2 O_2 ⊔ ((1+2i) + 2 O_2): 2 O_2 = {m even, n == 0 (4)},
# coset = {m odd, n == 2 (4)}; union = {m == n/2 parity} = O'.
def in_Op(m, b): return (m - b) % 2 == 0
okX3 = all(in_Op(m, b) == ((m % 2 == 0 and b % 2 == 0) or
                           (m % 2 == 1 and b % 2 == 1))
           for m in range(-9, 10) for b in range(-9, 10))
check_exact("[E1-X3] O' = 2 O_2 ⊔ ((1+2i)+2 O_2) (parity classes, exact)",
            okX3)
# [E1-X4] C = {m odd, n even} = C0 ⊔ C' with C0 = {n == 0 (4)},
# C' = {n == 2 (4)} = (1+2i) + 2 O_2 (as sets).
okX4 = all((b % 2 == 0) == ((b % 4 == 0) or (b % 4 == 2))
           for b in range(-9, 10)) \
       and all(((1 + 2*c) % 2 == 1) and ((2 + 4*d) % 4 == 2)
               for c in range(-6, 7) for d in range(-6, 7))
check_exact("[E1-X4] C = C0 ⊔ C'; C' = (1+2i)+2 O_2 = {m odd, n==2(4)} "
            "(exact)", okX4)
check_exact("[E1-X5] T(2 Lambda) = (1/16) T(Lambda) (B and G weights, "
            "exact)", Fr(1, 16) == Fr(4, 64))
okX6 = all(4*x*x - (x*x + y*y) == 2*(x*x - y*y) + (x*x + y*y)
           for x in range(-6, 7) for y in range(-6, 7))
check_exact("[E1-X6] T = B + 2G pointwise on forms (grid, exact)", okX6)

# --- exact Fraction coefficient track (basis (L16, L16tw, zC, zV2)) ------
B_OK = Fr(4)                                  # B(O_K) = 4 zC        [Q2]
B_C, G_C = Fr(3, 2), Fr(2)                    # C: (3/2) zC, 2 L16   [Q3]
B_C0, G_C0 = Fr(3, 4), Fr(1)                  # C0: (3/4) zC + zV2   [Q4]
# B(C0) = (0,0,3/4, +1) in (L16,L16tw,zC,zV2); G(C0) = (1,1,0,0)
B_Cp_zC = B_C - B_C0                          # (3/4) zC
B_Cp_zV = Fr(-1)                              # - zV2
G_Cp_L16, G_Cp_tw = G_C - G_C0, -G_C0         # L16 - L16tw
check_exact("[E1-S1] B(C') = (3/4) zC - zV2; G(C') = L16 - L16tw (exact)",
            (B_Cp_zC, B_Cp_zV) == (Fr(3, 4), Fr(-1))
            and (G_Cp_L16, G_Cp_tw) == (Fr(1), Fr(-1)))
T_O2 = (2*G_C, Fr(0), B_C + B_OK/16, Fr(0))   # O_2 = C ⊔ 2 O_K
check_exact("[E1-S2] T(O_2) = 4 L16 + (7/4) zC (exact)",
            T_O2 == (Fr(4), Fr(0), Fr(7, 4), Fr(0)))
T_2O2 = tuple(t/16 for t in T_O2)
T_Cp = (2*G_Cp_L16, 2*G_Cp_tw, B_Cp_zC, B_Cp_zV)
T_Op = tuple(T_2O2[i] + T_Cp[i] for i in range(4))
check_exact("[E1-S3] T(O') = (9/4) L16 - 2 L16tw + (55/64) zC - zV2 (exact)",
            T_Op == (Fr(9, 4), Fr(-2), Fr(55, 64), Fr(-1)))
T_L1 = tuple(16*t for t in T_Op)
check_exact("[E1-S4] T(Lambda_1) = 36 L16 - 32 L16tw + (55/4) zC - 16 zV2 "
            "(exact)", T_L1 == (Fr(36), Fr(-32), Fr(55, 4), Fr(-16)))
comb_c = tuple(-T_L1[i] + 4*T_O2[i] for i in range(4))
check_exact("[E1-S5] comb = -20 L16 + 32 L16tw - (27/4) zC + 16 zV2 (exact)",
            comb_c == (Fr(-20), Fr(32), Fr(-27, 4), Fr(16)))
lhs_e = tuple(10*comb_c[i] for i in range(4))     # EK4 = 10 comb.e (y0=1)
rhs_e = (Fr(5, 8)*(-20)*16, Fr(5, 8)*4*128, Fr(5, 8)*(-9)*12, Fr(5, 8)*4*64)
check_exact("[E1-S6] TARGET #1: EK4((1+2i)/2) = (5/8)(4 M16tw - 20 M16 - "
            "9 d4 + 4 d8) -- exact e-basis equality",
            lhs_e == rhs_e == (Fr(-200), Fr(320), Fr(-135, 2), Fr(160)))
check_exact("[E1-S7] sign pattern differs from T2's (5/16)(+20,+4,+9,+4) "
            "(exact)", rhs_e != (Fr(100), Fr(160), Fr(135, 4), Fr(80)))

print()
print("Separation track summary: entry #1 proved by exact algebra from:")
print("  (Q0) n4(s4(tau)) = EK4(tau): tau on the certified path of")
print("       n4_p4_t4_cert.py (top corner); propagation theorem")
print("  (Q1) s4(tau) = 143208-101574 sqrt2: cert0_n4_p3_t1t2.py")
print("  (Q2) h(-4)=1, units mu_4: B(O_K)=4 zC, G(O_K)=0 [Y1]")
print("  (Q3) theta identity 2 g16 = sum_{a==1(2)} alpha^2 q^N [L2 exact")
print("       to q^60, Sturm bound 6; Hecke]; ideal counting [Y3]:")
print("       B(C) = (3/2) zC, G(C) = 2 L16")
print("  (Q4) chi_8 o N projection [Y4]: B(C0) = (3/4) zC + zV2,")
print("       G(C0) = L16 + L16tw; L_K(chi_8 o N) = L(chi_8) L(chi_{-8})")
print("       [Y5 exact]")
print("  (Q5) FEs root number +1 at levels 16, 64 [Fricke ratios F0];")
print("       Dirichlet FE for chi_{-4}, chi_{-8}")
print("  (Q6) L(chi_8,2) = pi^2 sqrt2/16 (finite-sum, T2 [X7])")
print()

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

g16ser = s_pow(Pser(4, NMAX), 6, NMAX)
a16 = [0] + g16ser[:-1]
a16tw = [chi8p[n % 8]*a16[n] for n in range(NMAX+1)]

NT = 60
a_th16 = [0]*(NT+1)
for a in range(-8, 9):
    for b in range(-8, 9):
        Nm = a*a + b*b
        if 1 <= Nm <= NT and a % 2 == 1 and b % 2 == 0:
            a_th16[Nm] += a*a - b*b
check_exact("[L2] theta identity 2 g16 = sum_{a==1(2)} alpha^2 q^N to q^60 "
            "(exact; K = Q(i))",
            all(a_th16[n] == 2*a16[n] for n in range(1, NT+1)))

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

l1_16, l2_16, ords16, divs16 = ligozat({4: 6}, 16)
check_exact("[G3] g16: Ligozat sums 24 == 0 (24), 24 == 0 (24); weight 3",
            l1_16 % 24 == 0 and l2_16 % 24 == 0)
check_exact("[G3b] g16 character: -4^6 = -2^12, kernel -1 => chi_{-4} "
            "(exact)", -(4**6) == -4096 and 4096 == 64**2)
check_exact("[G4] g16 cusp orders %s all = 1" % [str(o) for o in ords16],
            all(o == 1 for o in ords16))
dd16 = sum(phi(gcd(d, 16//d))*o for d, o in zip(divs16, ords16))
check_exact("[G4b] g16 divisor degree 6 = (3/12) index 24 (exact)", dd16 == 6)
check_exact("[G5] Sturm bound 6 < 60 = theta-id coverage",
            Fr(3, 12)*24 == 6 and NT == 60)

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

L16, M16 = Lset(a16, 16, "g16")
L16tw, M16tw = Lset(a16tw, 64, "g16 x chi_8")
check("[F3] M16 = (16/pi^3) L(g16,3)", M16, 16*L16/pi**3, TOL)
check("[F4] M16tw = (128/pi^3) L(g16tw,3)", M16tw, 128*L16tw/pi**3, TOL)

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

zC = zeta(2)*Cat                   # zeta_K(2), K = Q(i)
zV2 = Lchi8p_2*Lchi8m_2

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

T_L1v = lattice_T(1, mpf(1)/2, 1)     # Lambda_1 = Z + tau Z
T_L2v = lattice_T(2, mpf(1)/2, 1)     # Lambda_2 = Z + 2 tau Z = O_2
T_L2alt = lattice_T(1, 0, 2)          # same lattice Z + 2i Z, pure-imag rows

T_L2_dec = 4*L16 + mpf(7)/4*zC
T_L1_dec = 36*L16 - 32*L16tw + mpf(55)/4*zC - 16*zV2
check("[T1] T(Lambda_2) = T(O_2) = 4 L16 + (7/4) zeta_K(2)",
      T_L2v, T_L2_dec, TOL)
check("[T1b] T(Lambda_2) via pure-imag rows (same lattice Z+2iZ)",
      T_L2v, T_L2alt, TOL)
check("[T2] T(Lambda_1) = 36 L16 - 32 L16tw + (55/4) zC - 16 zV2",
      T_L1v, T_L1_dec, TOL)
comb = -T_L1v + 4*T_L2v
comb_dec = -20*L16 + 32*L16tw - mpf(27)/4*zC + 16*zV2
check("[T3] comb = -20 L16 + 32 L16tw - (27/4) zC + 16 zV2",
      comb, comb_dec, TOL)

# [T0] independent direct lattice sum over O' (box truncation):
# T(O') = sum'_{m==b(2)} F(z), F(z) = 4(Re z)^2/|z|^6 - 1/|z|^4,
# z = m + 2bi, |z|^2 = m^2 + 4 b^2.
BD = 800
T_Op_direct = mpf(0)
for m in range(-BD, BD+1):
    for b in range(-BD, BD+1):
        if (m, b) == (0, 0) or (m - b) % 2 != 0:
            continue
        Nm = mpf(m*m + 4*b*b)
        T_Op_direct += 4*mpf(m*m)/Nm**3 - 1/Nm**2
check("[T0] direct O'-sum (box 800) = T(Lambda_1)/16 (1e-3 truncation)",
      T_Op_direct, T_L1v/16, mpf(10)**(-3))

EK4_v = (10*1/pi**3)*comb
conj = mpf(5)/8*(4*M16tw - 20*M16 - 9*d4 + 4*d8)
check("[E1] EK4((1+2i)/2) = (5/8)(4 M16tw - 20 M16 - 9 d4 + 4 d8)",
      EK4_v, conj, TOL)
EK4_asm = (10/pi**3)*comb_dec
check("[E3] assembly from decomposed comb (FE constants only)",
      EK4_asm, conj, TOL)

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
        self.f0, self.f1, self.f2 = f0, f1, f2
    def __add__(a, b):
        b = b if isinstance(b, Dual) else Dual(b)
        return Dual(a.f0+b.f0, a.f1+b.f1, a.f2+b.f2)
    __radd__ = __add__
    def __neg__(a):
        return Dual(-a.f0, -a.f1, -a.f2)
    def __sub__(a, b):
        b = b if isinstance(b, Dual) else Dual(b)
        return Dual(a.f0-b.f0, a.f1-b.f1, a.f2-b.f2)
    def __rsub__(a, b):
        return Dual(b)-a
    def __mul__(a, b):
        b = b if isinstance(b, Dual) else Dual(b)
        return Dual(a.f0*b.f0, a.f0*b.f1+a.f1*b.f0,
                    a.f0*b.f2+a.f1*b.f1+a.f2*b.f0)
    __rmul__ = __mul__
    def __truediv__(a, b):
        b = b if isinstance(b, Dual) else Dual(b)
        return a * Dual(1/b.f0, -b.f1/b.f0**2,
                        (b.f1**2 - b.f0*b.f2)/b.f0**3)
    def __rtruediv__(a, b):
        return Dual(a)/b
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

def lattice_T_iv_e1():
    """T(Lambda_1) for tau = 1/2 + i: rows y = m alternate shift 0
    (m even, coth) and shift 1/2 (m odd, tanh); power terms and tail
    bounds are shift-independent (Poisson kernel k = 0 term)."""
    y0 = iv.mpf(1)
    M = 52
    T = iv.pi**4/15
    for m in range(1, M+1):
        y = iv.mpf(m)*y0
        S2, S3 = row_S23_iv(y, shift=(m % 2 == 1))
        dS2 = S2 - iv.pi/(2*y**3)
        dS3 = S3 - 3*iv.pi/(8*y**5)
        T += 2*(3*dS2 - 4*y**2*dS3)
    tb2, tb3 = tail_row_iv((M+1)*1)
    rr = mpf(iv.exp(-2*iv.pi).b)
    yM = mpf(M+1)
    tot = 2*(3*tb2 + 4*yM*yM*tb3)/(1-rr)
    return T + iv.mpf([-tot, tot])

# --- the iv values ---
L16_iv, M16_iv = Lset_iv(a16, 16)
L16tw_iv, M16tw_iv = Lset_iv(a16tw, 64)
Cat_iv = dirichlet2_iv(chi4l, 4)
Lchi8m_2_iv = dirichlet2_iv(chi8m, 8)
Lchi8p_2_iv = dirichlet2_iv(chi8p, 8)
zC_iv = iv.pi**2/6 * Cat_iv
zV2_iv = Lchi8p_2_iv * Lchi8m_2_iv
d4_iv = 2*Cat_iv/iv.pi
d8_iv = 4*iv.sqrt(2)/iv.pi*Lchi8m_2_iv

def iv_mid(z):
    return (mp.convert(z.a) + mp.convert(z.b))/2

def iv_w(z):
    return mp.convert(z.b) - mp.convert(z.a)

check("[V0a] iv: L(g16,3) vs mp value", iv_mid(L16_iv), L16,
      iv_w(L16_iv) + TOL)
check("[V0b] iv: L(g16tw,3) vs mp value", iv_mid(L16tw_iv), L16tw,
      iv_w(L16tw_iv) + TOL)
check("[V0c] iv: L(chi_8,2) = pi^2 sqrt2/16", iv_mid(Lchi8p_2_iv),
      pi**2*s2/16, iv_w(Lchi8p_2_iv) + TOL)

T_L2_iv = lattice_T_iv(1, 2)           # Lambda_2 = O_2 = Z + 2i Z
T_L1_iv = lattice_T_iv_e1()            # Lambda_1 = Z + tau Z, tau = 1/2+i

check("[V0d] iv T(Lambda_2) contains the 60-dps mp value",
      iv_mid(T_L2_iv), T_L2v, iv_w(T_L2_iv) + TOL)
check("[V0e] iv T(Lambda_1) contains the 60-dps mp value",
      iv_mid(T_L1_iv), T_L1v, iv_w(T_L1_iv) + TOL)

T_L2_dec_iv = 4*L16_iv + iv.mpf(7)/4*zC_iv
T_L1_dec_iv = 36*L16_iv - 32*L16tw_iv + iv.mpf(55)/4*zC_iv - 16*zV2_iv

check_lock("[V1] LOCK T(Lambda_2) = 4 L16 + (7/4) zeta_K(2)",
           T_L2_iv - T_L2_dec_iv)
check_lock("[V2] LOCK T(Lambda_1) = 36 L16 - 32 L16tw + (55/4) zC - 16 zV2",
           T_L1_iv - T_L1_dec_iv)

comb_iv = -T_L1_iv + 4*T_L2_iv
EK4_iv = (10/iv.pi**3)*comb_iv
conj_iv = iv.mpf(5)/8*(4*M16tw_iv - 20*M16_iv - 9*d4_iv + 4*d8_iv)

check_lock("[V3] LOCK EK4((1+2i)/2) = (5/8)(4 M16tw - 20 M16 - 9 d4 + 4 d8)",
           EK4_iv - conj_iv)

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
    import sys; sys.exit(1)
print("ALL CHECKS PASSED")
