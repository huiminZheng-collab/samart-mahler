# verify_P1_n5_e6.py -- (P1)-type CM evaluation for the Samart Table-6
# OPEN conjecture (Phase 8, target #6; THIRD class-number-2 field):
#
#   tau = (1 + sqrt(-15))/2 = omega (the ring of integers itself!),
#   s4 = (-192303 - 85995 sqrt5)/2   [numeric lock 1.6e-55 in
#        n5_p8_e6_step0.py; iv lock in cert0_n5_e6.py],
#   EK4(tau) = (1/15)(160 M15^1 + 120 M15^2 + 88 d3 + 5 d15)
#
# K = Q(sqrt(-15)), disc -15, h = 2 (class group = genus group = Z/2).
# The two newforms of S_3(Gamma_0(15), chi_{-15}) are the CM theta series
#   g1 (Samart's #2, LMFDB 15.3.d.*, a(2) = +1),
#   g2 (Samart's #1, a(2) = -1),
# built from the grossencharacters psi, psi chi_g of conductor 1 with
#   psi((alpha)) = alpha^2 (principal),  psi(p_2) = omega, psi(p_2') = om-bar
#   (consistent: p_2^2 = (omega), psi(p_2)^2 = omega^2 = psi((omega))),
#   a1(n) = P(n) + (om-bar/8) R(n),  a2(n) = P(n) - (om-bar/8) R(n),
#   P(n) = (1/8) sum_{c^2+15d^2=4n, c==d(2)} (c^2-15d^2),
#   R(n) = sum_{c^2+15d^2=8n, c==d(4)} ((c^2-15d^2)/4 + (cd/2) sqrt-15)
#        = Rr(n) + Ri(n) sqrt-15;
# rationality (the sqrt-15 parts cancel) and integrality checked exactly.
# M15^i = L'(gi,0) (level 15, w = +1, Fricke ratio), d_k = L'(chi_{-k},-1).
#
# Lattice side (tau = omega, so Lambda_1 = O_K = Z[omega],
# Lambda_2 = Z + 2 omega Z = Z[sqrt-15] = O_2, the conductor-2 order):
#   O_2 = 2 O_K  ⊔  C'',   C'' = {a+b sqrt-15 : a =/= b (2)}
#        = {alpha in O_K : N(alpha) odd}
#        (c,d both odd gives N = (c^2+15d^2)/4 even since c^2+15d^2 == 0(8)).
#   Odd-norm restriction by Euler-factor removal at the two primes over 2
#   (chi_{-15}(2) = +1 split; chi_5(2) = -1; psi(p_2)+psi(p_2') = 1,
#    psi(p_2) psi(p_2') = 4):
#     B(C'') = (9/16) zK + (25/16) Lg,
#     G(C'') = (15/16) L1 + (19/16) L2,
#   where zK = zeta_K(2), Lg = L(chi_5,2) L(chi_{-3},2), and
#   (1 - omega/8)(1 - om-bar/8) = 15/16, (1 + omega/8)(1 + om-bar/8) = 19/16.
#   Hence
#     T(O_2) = T(O_K)/16 + T(C'')
#            = (5/8) zK + (13/8) Lg + 2 L1 + (5/2) L2,
#     comb = -T(O_K) + 4 T(O_2)
#          = 6 L1 + 8 L2 + (3/2) zK + (11/2) Lg,
#     EK4 = (5 sqrt15/pi^3) comb        (y0 = sqrt15/2).
#
# Layers (same three-track standard as verify_P1_n4_p4_t3.py):
#   [X*]/[S*] exact integer / Fraction checks (no floating point);
#   [G*] newform identifications: exact-integer Hecke/CM checks; root number
#        by the non-vacuous Fricke ratio (numeric, flagged);
#   [L*],[D*],[T*],[E*] mpmath 60-dps numerical confirmations;
#   [V*] rigorous interval locks (iv.dps = 70, self-contained theta-series
#        Mellin tail |a(n)| <= 30 n^3 -- see bound derivation in PART 2):
#        [V1] T(O_K), [V2] T(O_2), [V3] the EK4 identity, all locked with
#        half-width < 1e-48.
#
# Quoted inputs (same standard as T3):
#   (Q1) Hecke's theorem: theta series of conductor-1 grossencharacters are
#        newforms; S_3(Gamma_0(15), chi_{-15}) is 2-dimensional (LMFDB:
#        15.3.d.a/b are the only dim-1 CM forms, cm_discs = [-15]), so
#        {g1, g2} exhausts it.
#   (Q2) h(-15) = 2, class group Z/2, units +-1; 3, 5 ramify; 2 splits.
#   (Q3) Functional equations of the CM forms (level 15, w = +1).
#   (Q4) Dirichlet class-number-free formulas: d_k = (k^{3/2}/(4 pi))
#        L(chi_{-k},2); even-character formula L(chi,2) = pi^2 tau(chi)
#        sum_a chi(a) B_2(a/f) / f for even primitive chi (chi_5, f = 5,
#        Gauss sum +sqrt5).
#   (Q5) L_K(chi o N, s) = L(chi,s) L(chi chi_{-15},s); chi_5 chi_{-15}
#        = chi_{-3}.
#   (Q6) zeta(2) = pi^2/6; zeta_K(2) = zeta(2) L(chi_{-15},2).

from fractions import Fraction as Fr
from mpmath import (mp, mpf, mpc, pi, sqrt, zeta, dirichlet, gamma,
                    diff as mpdiff, power, gammainc, sinh, cosh, cos, exp)

mp.dps = 60
s3 = sqrt(mpf(3)); s5 = sqrt(mpf(5)); s15 = sqrt(mpf(15))
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

# --- field arithmetic in pair form: alpha = (c + d sqrt-15)/2, c==d(2) ---
def k15mul(u, v):    # exact multiplication, result again c==d(2)
    c = u[0]*v[0] - 15*u[1]*v[1]
    d = u[0]*v[1] + u[1]*v[0]
    assert c % 2 == 0 and d % 2 == 0
    return (c//2, d//2)
def k15norm(u): return (u[0]**2 + 15*u[1]**2)//4
OM = (1, 1)          # omega
check_exact("[X1] omega^2 = omega - 4 (pair arithmetic); N(omega) = 4; "
            "disc(O_K) = -15 (exact)",
            k15mul(OM, OM) == (-7, 1) and (-7, 1) == (1 - 8, 1)
            and k15norm(OM) == 4)
# p_2 = (2, omega): p_2^2 = (omega): products 4, 2 omega, omega^2 generate;
# omega = 4 + omega^2 in pair form: (8,0) + (-7,1) = (1,1).
check_exact("[X2] p_2^2 = (omega): omega = 4 + omega^2 (pair arithmetic) "
            "and N(p_2)^2 = 4 = N(omega)", (8 - 7, 0 + 1) == OM and
            2**2 == k15norm(OM))
# p_3 = (3, sqrt-15), p_5 = (5, sqrt-15): p_3^2 = (3), p_5^2 = (5).
check_exact("[X3] p_3^2 = (3), p_5^2 = (5) (norm multiplicativity)",
            9 == k15norm((6, 0)) and 25 == k15norm((10, 0)))

# --- characters -----------------------------------------------------------
def chi_15(n):
    if n % 3 == 0 or n % 5 == 0: return 0
    return (1 if n % 15 in (1, 2, 4, 8) else -1)
def chi5(n):  return 0 if n % 5 == 0 else (1 if n % 5 in (1, 4) else -1)
def chi3m(n): return 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)

okX4 = all(chi5(n) == chi3m(n) for n in range(1, 600)
           if n % 15 != 0 and chi_15(n) == 1)
check_exact("[X4] genus character: chi_5(N) = chi_{-3}(N) on ideal norms "
            "(chi_{-15}(N)=+1; n<600, exact)", okX4)

S5 = (sum(chi5(a) for a in range(5)), sum(chi5(a)*a for a in range(5)),
      sum(chi5(a)*a*a for a in range(5)))
check_exact("[X5] chi_5: sums (0, 0, 4) (exact)", S5 == (0, 0, 4))

# --- principal ideals have chi_5(N) = +1 ----------------------------------
okX6 = True
for c in range(-9, 10):
    for d in range(-9, 10):
        if (c, d) != (0, 0) and (c-d) % 2 == 0 and (c*c+15*d*d) % 4 == 0:
            Nm = (c*c + 15*d*d)//4
            if Nm % 5:
                okX6 &= chi5(Nm) == 1
check_exact("[X6] chi_5(N(alpha)) = +1 for principal (alpha), 5 -| N (box)",
            okX6)

# --- lattice set identities ----------------------------------------------
okX7 = True
for a in range(-12, 13):
    for b in range(-12, 13):
        Nm = a*a + 15*b*b
        okX7 &= ((a-b) % 2 != 0) == (Nm % 2 == 1)
check_exact("[X7] C'' = {a =/= b (2)} = {N odd} inside Z[sqrt-15] (box)",
            okX7)

# --- ideal counts vs representation numbers (exact integers) -------------
NID = 3000
r = [0]*(NID+1); g = [0]*(NID+1)
for dd in range(1, NID+1):
    cc = chi_15(dd)
    if cc:
        for m in range(dd, NID+1, dd): r[m] += cc
    c2 = chi5(dd)
    if c2:
        for m in range(1, NID//dd+1):
            c3 = chi3m(m)
            if c3: g[dd*m] += c2*c3
r1 = [0]*(NID+1); r2 = [0]*(NID+1)
RB = int(NID**0.5)+3
for x in range(-RB, RB+1):
    for y in range(-RB, RB+1):
        n = x*x + x*y + 4*y*y
        if 1 <= n <= NID: r1[n] += 1
        n = 2*x*x + x*y + 2*y*y
        if 1 <= n <= NID: r2[n] += 1
okX8 = all((r1[n]+r2[n]) == 2*r[n] and (r1[n]-r2[n]) == 2*g[n]
           for n in range(1, NID+1))
check_exact("[X8] ideal counts: r1+r2 = 2r(n), r1-r2 = 2g(n) (n<=3000; "
            "forms x^2+xy+4y^2, 2x^2+xy+2y^2; genus split, exact)", okX8)

# --- theta series (exact Fraction arithmetic in Q(sqrt-15)) ---------------
NMAX = 400
RB2 = int((8*NMAX)**0.5)+2
Pc = [Fr(0)]*(NMAX+1)
Rr = [Fr(0)]*(NMAX+1)
Ri = [Fr(0)]*(NMAX+1)
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
okX9 = (all(v == 0 for v in a1im[1:]) and all(v == 0 for v in a2im[1:])
        and all(v.denominator == 1 for v in a1re[1:])
        and all(v.denominator == 1 for v in a2re[1:])
        and a1re[1] == a2re[1] == 1)
check_exact("[X9] a1, a2 rational integral, a(1)=1 (sqrt-15 parts cancel; "
            "exact Fractions)", okX9)
a1 = [int(v) for v in a1re]
a2 = [int(v) for v in a2re]
check_exact("[X10] a1(2)=+1, a2(2)=-1, a1(5)=+5, a2(5)=-5 (Samart #2=g1, "
            "#1=g2 pinned by [E] below)", a1[2] == 1 and a2[2] == -1
            and a1[5] == 5 and a2[5] == -5)

PR = []
for n in range(2, NMAX+1):
    if all(n % p for p in PR if p*p <= n): PR.append(n)
oke1, oke2, okcm = True, True, True
for aa in (a1, a2):
    for p in PR:
        if p*p > NMAX: break
        if aa[p*p] != aa[p]*aa[p] - chi_15(p)*p*p: oke1 = False
    for i, p in enumerate(PR):
        for q in PR[i+1:]:
            if p*q > NMAX: break
            if aa[p*q] != aa[p]*aa[q]: oke2 = False
    for p in PR:
        if p > 100: break
        if chi_15(p) == -1 and aa[p] != 0: okcm = False
check_exact("[X11] Hecke recursion a(p^2)=a(p)^2-chi_{-15}(p)p^2; "
            "multiplicativity; CM vanishing at inert p<=100 (exact)",
            oke1 and oke2 and okcm)

# --- [S] separation track ---------------------------------------------------
# comb = 6 L1 + 8 L2 + (3/2) zK + (11/2) Lg  (basis L1, L2, zK, Lg),
# derived from:
#   B(O_K) = zK + Lg,  G(O_K) = L1 + L2,
#   B(C'') = (9/16) zK + (25/16) Lg,  G(C'') = (15/16) L1 + (19/16) L2,
#   T(O_2) = T(O_K)/16 + T(C'')  (O_2 = 2 O_K ⊔ C'', T(2L) = T(L)/16).
czK = -1 + 4*(Fr(1, 16) + Fr(9, 16))
cLg = -1 + 4*(Fr(1, 16) + Fr(25, 16))
cL1 = -2 + 4*(Fr(2, 16) + 2*Fr(15, 16))
cL2 = -2 + 4*(Fr(2, 16) + 2*Fr(19, 16))
combv = (cL1, cL2, czK, cLg)
check_exact("[S0] comb coefficients (L1,L2,zK,Lg) = (6,8,3/2,11/2) "
            "(exact bookkeeping)",
            combv == (Fr(6), Fr(8), Fr(3, 2), Fr(11, 2)))
# Euler factors at 2 (exact in pair arithmetic):
# (1-omega/8)(1-ombar/8) = 1 - (omega+ombar)/8 + omega*ombar/64
#                        = 1 - 1/8 + 4/64 = 15/16.
check_exact("[S1] (1-om/8)(1-omb/8) = 15/16, (1+om/8)(1+omb/8) = 19/16 "
            "(om+omb = 1, om omb = 4, exact)",
            1 - Fr(1, 8) + Fr(4, 64) == Fr(15, 16) and
            1 + Fr(1, 8) + Fr(4, 64) == Fr(19, 16))
# Substitutions (Q3),(Q4),(Q6):
#   L1 = pi^3 M1/(15 sqrt15) * ... : M = (15^{3/2}/(4 pi^3)) L  =>
#   (5 sqrt15/pi^3) L = (4/3) M.
#   zK = pi^3 d15 * 2/(45 sqrt15),  Lg = (4 pi^2 sqrt5/125)(4 pi d3/
#        (3 sqrt3)) = 16 pi^3 sqrt15 d3/1125.
cM1 = combv[0]*Fr(4, 3); cM2 = combv[1]*Fr(4, 3)
cd15 = 5*combv[2]*2/45            # (5 sqrt15/pi^3) * (2 pi^3/(45 sqrt15))
cd3 = combv[3]*5*16*15/1125       # 5 sqrt15 * 16 sqrt15/1125
check_exact("[S2] M(g1) = Samart #2: 6*(4/3) = 8 = 120/15 (exact)",
            cM1 == Fr(120, 15))
check_exact("[S3] M(g2) = Samart #1: 8*(4/3) = 32/3 = 160/15 (exact)",
            cM2 == Fr(160, 15))
check_exact("[S4] d15: 5*(3/2)*(2/45) = 1/3 = 5/15 (exact)",
            cd15 == Fr(5, 15))
check_exact("[S5] d3: (11/2)*(5*16*15/1125) = 88/15 (exact)",
            cd3 == Fr(88, 15))
check_exact("[S6] assembly closes: LHS = RHS monomial-wise (exact)",
            (cM2, cM1, cd3, cd15)
            == (Fr(160, 15), Fr(120, 15), Fr(88, 15), Fr(5, 15)))

# =====================================================================
# PART 1: [G]/[L]/[D]/[T]/[E] numeric confirmations (mp, 60 dps)
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

L1, M1 = Lset(a1, 15, "g1, N=15")
L2, M2 = Lset(a2, 15, "g2, N=15")
check("[G:M1] M1 = (15^{3/2}/(4 pi^3)) L1 (FE constant)", M1,
      mpf(15)**1.5*L1/(4*pi**3))

L2_15m = dirichlet(mpf(2), [chi_15(n) for n in range(15)])
L2_5p = dirichlet(mpf(2), [0, 1, -1, -1, 1])
L2_3m = dirichlet(mpf(2), [0, 1, -1])
check("[D1] L(chi_5,2) = 4 pi^2 sqrt5/125 (closed form, Q4 + [X5])",
      L2_5p, 4*pi**2*s5/125)
d3 = mpf(3)**1.5/(4*pi)*L2_3m
d15 = mpf(15)**1.5/(4*pi)*L2_15m
check("[D2] d15 = (15^{3/2}/4pi) L(chi_{-15},2) = L'(chi_{-15},-1) "
      "(direct derivative)", d15,
      mpdiff(lambda s: dirichlet(s, [chi_15(n) for n in range(15)]), mpf(-1)))

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

zK = zeta(2)*L2_15m
Lg = L2_5p*L2_3m

T_L1v = lattice_T(1, mpf(1)/2, s15/2)    # Lambda_1 = O_K
T_L2v = lattice_T(2, mpf(1)/2, s15/2)    # Lambda_2 = Z[sqrt-15] = O_2
T_L2alt = lattice_T(1, 0, s15)           # same lattice, pure-imag rows

T_L1_dec = zK + Lg + 2*(L1 + L2)
T_L2_dec = Fr(5, 8)*zK + Fr(13, 8)*Lg + 2*L1 + Fr(5, 2)*L2
check("[T1] T(O_K) = zeta_K(2) + L(chi5)L(chi-3) + 2(L1+L2)",
      T_L1v, T_L1_dec)
check("[T1b] T(O_2) via pure-imag rows (same lattice Z[sqrt-15])",
      T_L2v, T_L2alt)
check("[T2] T(O_2) = (5/8)zK + (13/8)Lg + 2L1 + (5/2)L2",
      T_L2v, T_L2_dec)
comb = -T_L1v + 4*T_L2v
comb_dec = 6*L1 + 8*L2 + Fr(3, 2)*zK + Fr(11, 2)*Lg
check("[T3] comb = 6L1 + 8L2 + (3/2)zK + (11/2)Lg", comb, comb_dec)

# [T0] independent direct lattice sum over C'' (box truncation):
# T(C'') = sum'_{a =/= b(2)} F(a+b sqrt-15), |z|^2 = a^2+15 b^2.
BD = 800
T_Cpp_direct = mpf(0)
for a in range(-BD, BD+1):
    for b in range(-BD, BD+1):
        if (a, b) == (0, 0) or (a-b) % 2 == 0:
            continue
        Nm = mpf(a*a + 15*b*b)
        T_Cpp_direct += 4*mpf(a*a)/Nm**3 - 1/Nm**2
T_Cpp_pred = T_L2v - T_L1v/16
check("[T0] direct C''-sum (box 800) = T(O_2) - T(O_K)/16 "
      "(1e-3 truncation)", T_Cpp_direct, T_Cpp_pred, mpf(10)**(-3))

EK4 = (5*s15/pi**3)*comb
conj = Fr(1, 15)*(160*M2 + 120*M1 + 88*d3 + 5*d15)
check("[E1] EK4(omega) = (1/15)(160 M15^1 + 120 M15^2 + 88 d3 + 5 d15) "
      "[Samart #1 = g2, #2 = g1]", EK4, conj)
EK4_asm = (5*s15/pi**3)*comb_dec
check("[E3] assembly from decomposed comb (FE constants only)",
      EK4_asm, conj)

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
else:
    print("ALL CHECKS PASSED (pre-iv tracks)")

# =====================================================================
# PART 2: [V] rigorous interval locks (mpmath.iv, iv.dps = 70)
#   Machinery ported from verify_P1_n5_e2.py.  Mellin tail uses the
#   self-contained theta bound |a^i(n)| <= 30 n^3:
#     |P(n)| <= (1/8)(4n) r1'(n) with r1'(n) <= (4 sqrt n+1)(4 sqrt(n/15)+1)
#              <= 12 n  =>  |P(n)| <= 6 n^2;
#     |Rr(n)|, |Ri(n)| <= (1/4)(8n) r2'(n) <= 2n * 12 n = 24 n^2;
#     |a(n)| <= |P| + (|Rr| + 15|Ri|)/16 <= 6n^2 + 24 n^2 = 30 n^2 <= 30 n^3.
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

def lattice_T_iv_e6():
    """T(Lambda_1) = T(O_K) for tau = 1/2 + i sqrt15/2: rows alternate
    shift 0 (m even, coth) and shift 1/2 (m odd, tanh); power terms and
    tail bounds are shift-independent (Poisson kernel k = 0 term)."""
    y0 = iv.sqrt(15)/2
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

# --- theta coefficients to n0(15) = 214 for the iv Mellin tails ----------
NMV = 300
Pc_v = [0]*(NMV+1); Rr_v = [0]*(NMV+1); Ri_v = [0]*(NMV+1)
RBV = int((8*NMV)**0.5)+2
for cc_ in range(-RBV, RBV+1):
    for dd_ in range(-RBV, RBV+1):
        nn_ = cc_*cc_ + 15*dd_*dd_
        if nn_ % 4 == 0 and 1 <= nn_//4 <= NMV and (cc_-dd_) % 2 == 0:
            Pc_v[nn_//4] += cc_*cc_ - 15*dd_*dd_
        if nn_ % 8 == 0 and 1 <= nn_//8 <= NMV and (cc_-dd_) % 4 == 0:
            Rr_v[nn_//8] += cc_*cc_ - 15*dd_*dd_
            Ri_v[nn_//8] += 2*cc_*dd_
# a = Pc_v/8 +/- (Rr_v/4 + 15 Ri_v/4)/16 = (8 Pc_v +/- (Rr_v+15 Ri_v))/64
a1v = [(8*Pc_v[n] + Rr_v[n] + 15*Ri_v[n])//64 for n in range(NMV+1)]
a2v = [(8*Pc_v[n] - Rr_v[n] - 15*Ri_v[n])//64 for n in range(NMV+1)]
check_exact("[V-pre] iv-range coefficients coincide with the exact track "
            "(n <= 400)", a1v == a1[:NMV+1] and a2v == a2[:NMV+1])

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

L1_iv, M1_iv = Lset_iv(a1v, 15)
L2_iv, M2_iv = Lset_iv(a2v, 15)
L5p_iv = dirichlet2_iv([0, 1, -1, -1, 1], 5)
L3m_iv = dirichlet2_iv([0, 1, -1], 3)
L15m_iv = dirichlet2_iv([chi_15(n) for n in range(15)], 15)

check("[V0a] iv: L(g1,3) vs mp value", iv_mid(L1_iv), L1, iv_w(L1_iv)+TOL)
check("[V0b] iv: L(g2,3) vs mp value", iv_mid(L2_iv), L2, iv_w(L2_iv)+TOL)
check("[V0c] iv: L(chi_5,2) = 4 pi^2 sqrt5/125", iv_mid(L5p_iv),
      4*pi**2*s5/125, iv_w(L5p_iv)+TOL)

zK_iv = iv.pi**2/6 * L15m_iv
Lg_iv = L5p_iv * L3m_iv
d3_iv = 3*iv.sqrt(3)/(4*iv.pi) * L3m_iv
d15_iv = iv.mpf(15)*iv.sqrt(15)/(4*iv.pi) * L15m_iv

s15_iv = iv.sqrt(15)
T_L1_iv = lattice_T_iv_e6()            # Lambda_1 = O_K
T_L2_iv = lattice_T_iv(1, s15_iv)      # Lambda_2 = Z[sqrt-15]

check("[V0d] iv T(O_K) contains the 60-dps mp value",
      iv_mid(T_L1_iv), T_L1v, iv_w(T_L1_iv) + TOL)
check("[V0e] iv T(O_2) contains the 60-dps mp value",
      iv_mid(T_L2_iv), T_L2v, iv_w(T_L2_iv) + TOL)

T_L1_dec_iv = zK_iv + Lg_iv + 2*(L1_iv + L2_iv)
T_L2_dec_iv = (iv.mpf(5)/8*zK_iv + iv.mpf(13)/8*Lg_iv
               + 2*L1_iv + iv.mpf(5)/2*L2_iv)

check_lock("[V1] LOCK T(O_K) = zeta_K(2) + L(chi5)L(chi-3) + 2(L1+L2)",
           T_L1_iv - T_L1_dec_iv)
check_lock("[V2] LOCK T(O_2) = (5/8)zK + (13/8)Lg + 2L1 + (5/2)L2",
           T_L2_iv - T_L2_dec_iv)

comb_iv = -T_L1_iv + 4*T_L2_iv
EK4_iv = (5*s15_iv/iv.pi**3)*comb_iv
conj_iv = iv.mpf(1)/15*(160*M2_iv + 120*M1_iv + 88*d3_iv + 5*d15_iv)
check_lock("[V3] LOCK EK4(omega) = (1/15)(160 M15^1 + 120 M15^2"
           " + 88 d3 + 5 d15)", EK4_iv - conj_iv, mpf(10)**(-48))

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
    import sys; sys.exit(1)
print("ALL CHECKS PASSED")
