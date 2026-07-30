# verify_P1_n4_81.py
#
# (P1)-type exact CM evaluation for Samart's Table-6 target (n4 family)
#     EK4(tau2) = 40 M7,   tau2 = (7+sqrt(-7))/4,   M7 = L'(g7,0),
# g7 = eta(t)^3 eta(7t)^3.  Here EK4 is Samart's Prop 2.1(iii) expression
#     EK4(tau) = (10 Im tau/pi^3)(-T1 + 4 T2),
#     T_j(tau) = sum'_{m,n} [ 4 Re(j m tau + n)^2/|j m tau + n|^6 - |j m tau + n|^{-4} ].
#
# Lattice structure at tau2 (see P1_n4_81.md):
#   Lambda_2 = Z + Z 2 tau2 = Z + Z pi = O_K           (2 tau2 = 3 + pi)
#   Lambda_1 = Z + Z tau2  = Z + Z tau_w = (pibar)/2   (tau2 = tau_w + 2,
#             tau_w = (-1+sqrt(-7))/4 as in P1_n2pair.md)
# Both are O_K-ideals -- no imprimitive sums.  The zeta_K(2) terms cancel
# exactly in -T1 + 4 T2 = 28 L(g7,3), giving the single-term formula
#     EK4(tau2) = 70 sqrt7 L(g7,3)/pi^3 = 40 M7.
#
# Anchors for the EK4 implementation (Samart, arXiv:1205.4803, proved cases):
#   tau = i          (s4 = 648, Lemma 2.2):  EK4(i) = f4(648)
#                    = (160/pi^3) L(h,3) + (5/pi) L(chi_{-4},2)   [Thm 1.4 (1.9)]
#   tau = i sqrt3/2  (s2 = 256):             EK2 = f2(256)
#                    = (64 sqrt3/pi^3) L(g48,3) + (16/3 pi) L(chi_{-4},2)
#                                                                 [Thm 1.4 (1.6)]
#   tau = i sqrt2/2  (s4 = 256, boundary of the proved region):
#                    EK4(i/sqrt2) vs Rogers' 5F4 value of f4(256) [Prop 1.3(iii)]
#
# Every intermediate identity is verified numerically (60 dps working
# precision); each line prints PASS/FAIL.  Derivations: P1_n4_81.md.
#
# Checks:
#  [L1] g7: theta identity, root number +1, FE, M7 vs reference,
#       M7 = 7 sqrt7 L(g7,3)/(4 pi^3)
#  [S1] B(O_K)   = sum' |a|^-4          = 2 zeta_K(2)      (Poisson y-rows)
#  [S2] G(O_K)   = sum' abar^2 |a|^-6   = 2 L(g7,3)       (Poisson y-rows)
#  [S3] B((pibar)) = sum'_(pibar) |a|^-4 = zeta_K(2)/2
#  [S4] sum'_(pibar) abar^2 |a|^-6 = pi^2 L(g7,3)/4 (Re = -3L/8, Im = sqrt7 L/8)
#  [T1] T1(tau2) = -12 L(g7,3) + 8 zeta_K(2)   (direct Poisson-row T-sum)
#  [T2] T2(tau2) =   4 L(g7,3) + 2 zeta_K(2)
#  [T3] T1(tau2) = T1(tau_w)                    (tau2 = tau_w + 2)
#  [CB] -T1 + 4 T2 = 28 L(g7,3);  zeta_K(2) cancellation
#  [EK] EK4(tau2) = 70 sqrt7 L(g7,3)/pi^3 = 40 M7;  coefficient identities;
#       independent U-series form of EK4; periodicity EK4(tau2)=EK4(tau_w)
#  [A0] s4(tau2) = 81 (numeric), s4(i) = 648, s4(i/sqrt2) = 256,
#       s2(i sqrt3/2) = 256   [Samart Lemma 2.2 values]
#  [A1] h = eta(4t)^6 machinery: root number, FE, EK2(i/2) = f2(64)
#       = (128/pi^3) L(h,3)  [re-validation of samart_ek.py Step 1]
#  [A2] ANCHOR: EK4(i) (lattice & U-series) = f4(648) [Thm 1.4 (1.9)]
#       = Rogers 5F4 value = direct torus integration n4(648)
#  [A3] g48 machinery: a_n(g48) = chi_{-4}(n) a_n(g), a_n(g) from the form
#       m^2+3n^2 [Samart Lemma 2.7], exact eta-quotient series check;
#       Fricke root number at level 48;  EK2(i sqrt3/2) = f2(256)
#       [Thm 1.4 (1.6)] = Rogers 5F4(1/4)
#  [A4] EK4(i/sqrt2) (lattice & U) = Rogers f4(256) via 5F4 at z=1 (nsum);
#       comparison with the f2(256) closed form (coincidence test)

from mpmath import (mp, mpf, mpc, pi, sqrt, exp, sin, cos, sinh, cosh, zeta,
                    dirichlet, gamma, diff as mpdiff, hyper, log, gammainc,
                    power, quad, polyroots, catalan, nsum, inf, loggamma)

mp.dps = 60
s7 = sqrt(mpf(7))
s3 = sqrt(mpf(3))
FAILS = []

def check(name, got, want, tol):
    d = abs(got - want)
    ok = d < tol
    if not ok:
        FAILS.append(name)
    print("%-74s %s  (|diff| = %.2e)" % (name, "PASS" if ok else "FAIL", mpf(d)))

TOL = mpf(10) ** (-45)
TOLM = mpf(10) ** (-40)

# ================= g7 coefficients, L-values (as verify_P1_n2pair.py) ========
NMAX = 200
tri = []
j = 0
while j * (j + 1) // 2 < NMAX:
    tri.append(j * (j + 1) // 2)
    j += 1
a = [0] * (NMAX + 1)
for j, tj in enumerate(tri):
    for kk, tk in enumerate(tri):
        n = tj + 7 * tk + 1
        if n <= NMAX:
            a[n] += (-1) ** (j + kk) * (2 * j + 1) * (2 * kk + 1)

xN7 = s7

def I7(s):
    s = mpf(s)
    tot = mpf(0)
    for n in range(1, NMAX + 1):
        if a[n] == 0:
            continue
        tot += a[n] * power(2 * pi * n, -s) * gammainc(s, 2 * pi * n / xN7)
    return tot

def eta(tau, nterms=400):
    q = exp(2 * pi * 1j * tau)
    p = mpc(1)
    qn = q
    for n in range(1, nterms + 1):
        p *= (1 - qn)
        qn *= q
        if abs(qn) < mpf(10) ** (-65):
            break
    return exp(pi * 1j * tau / 12) * p

def g7f(tau):
    return eta(tau) ** 3 * eta(7 * tau) ** 3

y = mpf("0.6")
w7 = (g7f(1j / (7 * y)) / (y ** 3 * g7f(1j * y))) / xN7 ** 3
check("[L1] root number w(g7) = +1", w7, 1, TOLM)

def Lam7(s):
    return xN7 ** s * I7(s) + w7 * xN7 ** (3 - s) * I7(3 - s)

check("[L1] functional equation Lam(1) = Lam(2)", Lam7(1), Lam7(2), TOL)
Lam3 = Lam7(3)
M7 = w7 * Lam3
L3 = Lam3 * (2 * pi) ** 3 / (xN7 ** 3 * gamma(3))
M7_ref = mpf("0.10267160777890201121045659489829291399889482708922")
check("[L1] M7 = L'(g7,0) matches reference", M7, M7_ref, TOL)
check("[L1] FE identity M7 = 7 sqrt7 L(g7,3)/(4 pi^3)", M7, 7 * s7 * L3 / (4 * pi ** 3), TOL)

chi7 = [0, 1, 1, -1, 1, -1, -1]
Lchi7_2 = dirichlet(mpf(2), chi7)
zK2 = zeta(2) * Lchi7_2

# ================= [L1] theta identity =================
NT = 60
a_th = [0] * (NT + 1)
for N in range(1, NT + 1):
    rat, ompart = 0, 0
    B = int(sqrt(float(N / 2))) + 2
    Amax = int(sqrt(float(N))) + 2
    for n in range(-B, B + 1):
        for m in range(-Amax, Amax + 1):
            if m * m + m * n + 2 * n * n == N and (m, n) != (0, 0):
                rat += m * m - 2 * n * n
                ompart += 2 * m * n + n * n
    assert ompart == 0 and rat % 2 == 0
    a_th[N] = rat // 2
ok = all(a_th[n] == a[n] for n in range(1, NT + 1))
print("%-74s %s" % ("[L1] theta identity g7 = (1/2) sum' a^2 q^{N(a)} to q^60",
                    "PASS" if ok else "FAIL"))
if not ok:
    FAILS.append("theta")

# ================= Poisson-row machinery =================
def G_row(x, y):
    return (pi / y) * sinh(2 * pi * y) / (cosh(2 * pi * y) - cos(2 * pi * x))

def row_powers(x, y):
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y) / (2 * y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy) / (2 * yy), y) / (4 * y)
    return S2, S3

def lattice_T(d, tau):
    x0, y0 = tau.real, tau.imag
    B, M = mpf(0), mpf(0)
    for m in range(-300, 301):
        if m == 0:
            B += 2 * zeta(4)
            continue
        x = d * m * x0
        yy = abs(d * m) * y0
        S2, S3 = row_powers(x, yy)
        B += S2
        M += m * m * S3
        if abs(d * m) * y0 > 45 and m > 0:
            break
    return 3 * B - 4 * d * d * y0 * y0 * M

# ================= [S1]-[S4]: Hecke sums over O_K and (pibar) ================
# O_K = Z + Z pi, pi = (1+sqrt(-7))/2: elements x + y pi = (x + y/2) + i y sqrt7/2.
# B(O_K) = sum' |a|^-4,  G(O_K) = sum' abar^2 |a|^-6 = B - sum' 2 v^2/Q^3
# (v = y sqrt7/2, so 2 v^2 = 7 y^2/2);  G(O_K) is real (conjugation symmetry).
def sums_OK():
    # Poisson rows carry algebraic tails: S2(sh,v) = pi/(2 v^3) + exp-small,
    # S3(sh,v) = 3 pi/(8 v^5) + exp-small.  Subtract them per row and add the
    # analytic y-sums back (zeta(3) terms); otherwise the rows converge only
    # like 1/y^3 (in lattice_T these tails cancel in 3B - 4 d^2 y0^2 M).
    B = 2 * zeta(4)          # y = 0 row
    sub = mpf(0)
    for yy in range(1, 60):
        v = s7 * yy / 2
        sh = mpf(yy % 2) / 2
        S2, S3 = row_powers(sh, v)
        B += 2 * (S2 - pi / (2 * v ** 3))
        sub += 2 * (7 * yy * yy / 2) * (S3 - 3 * pi / (8 * v ** 5))
        if v > 45:
            break
    B += pi * (2 / s7) ** 3 * zeta(3)                 # 2 sum_y pi/(2 v^3)
    sub += (21 * pi / 8) * (2 / s7) ** 5 * zeta(3)    # sum_y 21 pi y^2/(8 v^5)
    return B, B - sub

B_OK, G_OK = sums_OK()
check("[S1] B(O_K) = sum'_{O_K} |a|^-4 = 2 zeta_K(2)", B_OK, 2 * zK2, TOL)
check("[S2] G(O_K) = sum' abar^2 |a|^-6 = 2 L(g7,3)", G_OK, 2 * L3, TOL)

# (pibar) = 2 Z + pibar Z: elements 2a + b pibar = (2a + b/2) - i b sqrt7/2.
# Row at fixed b: u = 2a + b/2, v = b sqrt7/2;  rescaling w = a + b/4 gives
# S_k(row) = 2^{-2k} S_k(b/4 mod 1, v/2).  b = 0 row: sum_{a!=0} (2a)^-4 = zeta(4)/8.
def sums_pibar():
    B = zeta(4) / 8
    sub_re, sub_im = mpf(0), mpf(0)
    for b in range(1, 100):
        v4 = s7 * b / 4          # v/2, v = b sqrt7/2
        sh = mpf(b % 4) / 4
        S2, S3 = row_powers(sh, v4)
        B += 2 * (S2 - pi / (2 * v4 ** 3)) / 16
        sub_re += 2 * (7 * b * b / 2) * (S3 - 3 * pi / (8 * v4 ** 5)) / 64
        # Im part: -2 sum u v/Q^3 with v = -b sqrt7/2 (row b>0), and
        # sum_a u/(u^2+v^2)^3 = (1/32) sum_w w/(w^2+(v/2)^2)^3
        #                     = (1/32) * (-(1/4) d/dx S2(x, v/2));
        # w/(w^2+c^2)^3 is odd and absolutely integrable, so this row sum
        # is exponentially convergent (no algebraic tail).
        dS2 = mpdiff(lambda x: row_powers(x, v4)[0], sh)
        sub_im += 2 * (s7 * b) * (mpf(1) / 32) * (-dS2 / 4)
        if v4 > 45:
            break
    B += 2 * (pi / 32) * (4 / s7) ** 3 * zeta(3)              # 2 sum_b pi/(32 v4^3)
    sub_re += 2 * (7 / 2) * (3 * pi / 512) * (4 / s7) ** 5 * zeta(3)
    return B, B - sub_re, sub_im

B_pi, G_pi_re, G_pi_im = sums_pibar()
pii = mpc(1, s7) / 2
check("[S3] B((pibar)) = sum'_(pibar) |a|^-4 = zeta_K(2)/2", B_pi, zK2 / 2, TOL)
check("[S4] Re sum'_(pibar) abar^2 |a|^-6 = Re(pi^2 L/4) = -3 L/8",
      G_pi_re, -3 * L3 / 8, TOL)
check("[S4] Im sum'_(pibar) abar^2 |a|^-6 = Im(pi^2 L/4) = sqrt7 L/8",
      G_pi_im, s7 * L3 / 8, TOL)
check("[S4] Re pi^2 = -3/2  (exact algebra)", (pii ** 2).real, mpf(-3) / 2, TOL)

# ================= [T1]-[CB]: T-sums at tau2 =================
tau2 = mpc(7, s7) / 4
tau_w = mpc(-1, s7) / 4
check("[T0] 2 tau2 = 3 + pi  (exact)", 2 * tau2, 3 + pii, TOL)
check("[T0] tau2 = tau_w + 2  (exact)", tau2, tau_w + 2, TOL)

T2v = lattice_T(2, tau2)
T1v = lattice_T(1, tau2)
check("[T2] T2(tau2) = 4 L(g7,3) + 2 zeta_K(2)   [Lambda_2 = O_K]",
      T2v, 4 * L3 + 2 * zK2, TOL)
check("[T1] T1(tau2) = -12 L(g7,3) + 8 zeta_K(2)   [Lambda_1 = (pibar)/2]",
      T1v, -12 * L3 + 8 * zK2, TOL)
T1w = lattice_T(1, tau_w)
check("[T3] T1(tau_w) = T1(tau2)  (same lattice)", T1w, T1v, TOL)
# T1 from the ideal decomposition directly: 16 [2 Re G((pibar)) + B((pibar))]
T1_ideal = 16 * (2 * G_pi_re + B_pi)
check("[T1] T1(tau2) = 16 [2 Re G((pibar)) + B((pibar))]  (ideal route)",
      T1_ideal, T1v, TOL)

comb = -T1v + 4 * T2v
check("[CB] -T1 + 4 T2 = 28 L(g7,3)   (zeta_K(2) cancels)", comb, 28 * L3, TOL)
check("[CB] zeta-coefficient identity -8 + 4*2 = 0  (exact integers)",
      mpf(-8 + 8), 0, TOL)

# ================= [EK]: final assembly =================
EK4_lattice = (10 * tau2.imag / pi ** 3) * comb
target = 40 * M7
check("[EK] EK4(tau2) = 40 M7", EK4_lattice, target, TOL)
check("[EK] EK4(tau2) = 70 sqrt7 L(g7,3)/pi^3", EK4_lattice, 70 * s7 * L3 / pi ** 3, TOL)
check("[EK] coefficient identity (5 sqrt7/2 pi^3)*28 = 40*(7 sqrt7/4 pi^3)",
      (5 * s7 / (2 * pi ** 3)) * 28, 40 * (7 * s7 / (4 * pi ** 3)), TOL)

# independent U-series form  EK4_U = Im[2 pi tau + (10/pi^3)(U1 - 2 U2)]
def U(j, tau):
    tot = mpc(0)
    for m in range(1, 3000):
        z = j * pi * m * tau
        term = cos(z) / sin(z) ** 3 / m
        tot += term
        if abs(term) < mpf(10) ** (-65):
            break
    return 2 * pi ** 3 * tot

def EK4_U(tau):
    return (2 * pi * tau + (10 / pi ** 3) * (U(1, tau) - 2 * U(2, tau))).imag

EK4u = EK4_U(tau2)
check("[EK] U-series EK4(tau2) = lattice EK4(tau2)", EK4u, EK4_lattice, TOL)
check("[EK] periodicity: EK4_U(tau_w) = EK4_U(tau2)", EK4_U(tau_w), EK4u, TOL)

# ================= [A0]: s4 / s2 values =================
def Delta(tau):
    return eta(tau) ** 24

def s4(tau):
    e1, e2, e4 = eta(tau), eta(2 * tau), eta(4 * tau)
    r = e1 * e4 ** 2 / e2 ** 3
    return (Delta(2 * tau) / Delta(tau)) * (16 * r ** 4 + r ** -4) ** 4

def s2(tau):
    return -Delta(tau + mpf(1) / 2) / Delta(2 * tau + 1)

check("[A0] s4(tau2) = 81  (numeric link to n4(81))", s4(tau2), 81, TOL)
check("[A0] s4(i) = 648  [Samart Lemma 2.2]", s4(mpc(0, 1)), 648, TOL)
check("[A0] s4(i/sqrt2) = 256  [Samart Lemma 2.2]", s4(mpc(0, 1) / sqrt(mpf(2))), 256, TOL)
check("[A0] s2(i sqrt3/2) = 256  [Samart Lemma 2.2]", s2(mpc(0, s3) / 2), 256, TOL)

# ================= [A1]: h = eta(4t)^6 machinery =================
NH = 300
triH = []
j = 0
while j * (j + 1) // 2 < NH:
    triH.append(j * (j + 1) // 2)
    j += 1
AH = [0] * (NH + 1)
for i_, ti in enumerate(triH):
    for jj, tj in enumerate(triH):
        n = ti + tj
        if n <= NH:
            AH[n] += (-1) ** (i_ + jj) * (2 * i_ + 1) * (2 * jj + 1)
ah = [0] * (4 * NH + 2)
for n in range(NH + 1):
    ah[4 * n + 1] = AH[n]

xN16 = mpf(4)

def I16(s):
    s = mpf(s)
    tot = mpf(0)
    for n in range(1, len(ah)):
        if ah[n] == 0:
            continue
        tot += ah[n] * power(2 * pi * n, -s) * gammainc(s, 2 * pi * n / xN16)
    return tot

def hf(tau):
    return eta(4 * tau) ** 6

wh1 = (hf(1j / (16 * mpf("0.35"))) / (mpf("0.35") ** 3 * hf(1j * mpf("0.35")))) / xN16 ** 3
wh2 = (hf(1j / (16 * mpf("0.53"))) / (mpf("0.53") ** 3 * hf(1j * mpf("0.53")))) / xN16 ** 3
check("[A1] root number w(h) = +1 (y = 0.35)", wh1, 1, TOLM)
check("[A1] root number w(h) = +1 (y = 0.53)", wh2, 1, TOLM)

def Lam16(s):
    return xN16 ** s * I16(s) + wh1 * xN16 ** (3 - s) * I16(3 - s)

check("[A1] FE check Lam_h(1) = Lam_h(2)", Lam16(1), Lam16(2), TOL)
Lh3 = Lam16(3) * (2 * pi) ** 3 / (xN16 ** 3 * gamma(3))

def EK2(tau):
    return (2 * tau.imag / pi ** 3) * (-lattice_T(1, tau) + 16 * lattice_T(4, tau))

check("[A1] EK2(i/2) = f2(64) = (128/pi^3) L(h,3)  [samart_ek.py Step 1]",
      EK2(mpc(0, mpf(1) / 2)), 128 / pi ** 3 * Lh3, TOL)

# ================= [A2]: ANCHOR at tau = i, f4(648) =================
Lchi4_2 = dirichlet(mpf(2), [0, 1, 0, -1])
check("[A2] L(chi_{-4},2) = Catalan  (sanity)", Lchi4_2, catalan, TOL)

f4_648 = 160 / pi ** 3 * Lh3 + 5 / pi * Lchi4_2
EK4i = (10 * mpc(0, 1).imag / pi ** 3) * (-lattice_T(1, mpc(0, 1)) + 4 * lattice_T(2, mpc(0, 1)))
check("[A2] EK4(i) [lattice] = f4(648) [Thm 1.4 (1.9)]", EK4i, f4_648, TOL)
check("[A2] EK4(i) [U-series] = f4(648)", EK4_U(mpc(0, 1)), f4_648, TOL)

f4_648_rogers = (log(648) - (mpf(24) / 648) *
                 hyper([mpf(5) / 4, mpf(3) / 2, mpf(7) / 4, 1, 1],
                       [2, 2, 2, 2], mpf(32) / 81)).real
check("[A2] Rogers 5F4(32/81) value of f4(648) agrees", f4_648_rogers, f4_648, TOL)

def mahler4(c, maxdegree):
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

c648 = mpf(648) ** mpf("0.25")
direct648 = 4 * mahler4(c648, 5)
check("[A2] direct torus integral n4(648) = f4(648)  (deg 5, loose tol)",
      direct648, f4_648, mpf(10) ** (-9))

# ================= [A3]: g48 machinery, f2(256) =================
NG = 400
ag = [0] * (NG + 1)
for N in range(1, NG + 1):
    s = 0
    Bn = int(sqrt(float(N / 3))) + 2
    Am = int(sqrt(float(N))) + 2
    for n in range(-Bn, Bn + 1):
        for m in range(-Am, Am + 1):
            if m * m + 3 * n * n == N:
                s += m * m - 3 * n * n
    assert s % 2 == 0
    ag[N] = s // 2

def chi4v(n):
    if n % 2 == 0:
        return 0
    return 1 if n % 4 == 1 else -1

a48 = [0] * (NG + 1)
for n in range(1, NG + 1):
    a48[n] = chi4v(n) * ag[n]
# spot values from Samart's paper: g48 = q + 3 q^3 - 2 q^7 + 9 q^9 - 22 q^13 - 26 q^19 - 6 q^21
spot = {1: 1, 3: 3, 7: -2, 9: 9, 13: -22, 19: -26, 21: -6}
okspot = all(a48[k] == v for k, v in spot.items())
print("%-74s %s" % ("[A3] g48 coefficients (form m^2+3n^2, twist chi_-4) match paper",
                    "PASS" if okspot else "FAIL"))
if not okspot:
    FAILS.append("g48spot")

# exact eta-quotient series check of g48 to O(q^100)
ORD = 100

def Pser(d, N):
    res = [0] * (N + 1)
    res[0] = 1
    k = 1
    while True:
        e1 = k * (3 * k - 1) // 2 * d
        e2 = k * (3 * k + 1) // 2 * d
        if e1 > N and e2 > N:
            break
        if e1 <= N:
            res[e1] += (-1) ** k
        if e2 <= N:
            res[e2] += (-1) ** k
        k += 1
    return res

def s_mul(A, B, N):
    C = [0] * (N + 1)
    for i in range(N + 1):
        if A[i] == 0:
            continue
        for j in range(N + 1 - i):
            if B[j]:
                C[i + j] += A[i] * B[j]
    return C

def s_pow(A, e, N):
    R = [0] * (N + 1)
    R[0] = 1
    for _ in range(e):
        R = s_mul(R, A, N)
    return R

def s_inv(A, N):
    assert A[0] == 1
    R = [0] * (N + 1)
    R[0] = 1
    for n in range(1, N + 1):
        R[n] = -sum(A[k] * R[n - k] for k in range(1, n + 1))
    return R

num = s_mul(s_pow(Pser(4, ORD), 9, ORD), s_pow(Pser(12, ORD), 9, ORD), ORD)
den = s_mul(s_mul(s_pow(Pser(2, ORD), 3, ORD), s_pow(Pser(6, ORD), 3, ORD), ORD),
            s_mul(s_pow(Pser(8, ORD), 3, ORD), s_pow(Pser(24, ORD), 3, ORD), ORD), ORD)
g48ser = s_mul(num, s_inv(den, ORD), ORD)   # = sum a48[n] q^{n-1}
oketa = all(g48ser[n - 1] == a48[n] for n in range(1, ORD))
print("%-74s %s" % ("[A3] g48 = eta(4)^9 eta(12)^9/(eta(2)^3 eta(6)^3 eta(8)^3 eta(24)^3) series check",
                    "PASS" if oketa else "FAIL"))
if not oketa:
    FAILS.append("g48eta")

xN48 = sqrt(mpf(48))

def I48(s):
    s = mpf(s)
    tot = mpf(0)
    for n in range(1, NG + 1):
        if a48[n] == 0:
            continue
        tot += a48[n] * power(2 * pi * n, -s) * gammainc(s, 2 * pi * n / xN48)
    return tot

def g48f(tau):
    return (eta(4 * tau) ** 9 * eta(12 * tau) ** 9 /
            (eta(2 * tau) ** 3 * eta(6 * tau) ** 3 * eta(8 * tau) ** 3 * eta(24 * tau) ** 3))

w48a = (g48f(1j / (48 * mpf("0.37"))) / (mpf("0.37") ** 3 * g48f(1j * mpf("0.37")))) / xN48 ** 3
w48b = (g48f(1j / (48 * mpf("0.51"))) / (mpf("0.51") ** 3 * g48f(1j * mpf("0.51")))) / xN48 ** 3
print("     [A3] g48 Fricke ratios:", w48a, w48b)
check("[A3] root number w(g48) = +1 (y = 0.37)", w48a, 1, TOLM)
check("[A3] root number w(g48) = +1 (y = 0.51)", w48b, 1, TOLM)

def Lam48(s):
    return xN48 ** s * I48(s) + w48a * xN48 ** (3 - s) * I48(3 - s)

L48_3 = Lam48(3) * (2 * pi) ** 3 / (xN48 ** 3 * gamma(3))

f2_256 = 64 * s3 / pi ** 3 * L48_3 + 16 / (3 * pi) * Lchi4_2
EK2_i32 = EK2(mpc(0, s3) / 2)
check("[A3] EK2(i sqrt3/2) = f2(256) [Thm 1.4 (1.6)]", EK2_i32, f2_256, TOL)
f2_256_rogers = (log(256) - (mpf(8) / 256) *
                 hyper([mpf(3) / 2, mpf(3) / 2, mpf(3) / 2, 1, 1],
                       [2, 2, 2, 2], mpf(1) / 4)).real
check("[A3] Rogers 5F4(1/4) value of f2(256) agrees", f2_256_rogers, f2_256, TOL)

# ================= [A4]: anchor at tau = i/sqrt2, f4(256) =================
taus2 = mpc(0, 1) / sqrt(mpf(2))
EK4_s2 = (10 * taus2.imag / pi ** 3) * (-lattice_T(1, taus2) + 4 * lattice_T(2, taus2))
check("[A4] EK4(i/sqrt2) [lattice] = [U-series]", EK4_U(taus2), EK4_s2, TOL)

# 5F4 at z = 1 converges algebraically (~ n^{-5/2}): sum the head by the term
# recurrence (exact small-number arithmetic) and extrapolate the tail with
# Euler--Maclaurin (Richardson/Shanks are poor for algebraic tails).
def hyp1term(n):
    n = mpf(n)
    return exp(loggamma(n + mpf(5) / 4) - loggamma(mpf(5) / 4)
               + loggamma(n + mpf(3) / 2) - loggamma(mpf(3) / 2)
               + loggamma(n + mpf(7) / 4) - loggamma(mpf(7) / 4)
               + loggamma(n + 1) - 4 * loggamma(n + 2))

N0 = 100000
t = mpf(1)
S_head = t
for n in range(0, N0):
    t *= (n + mpf(5) / 4) * (n + mpf(3) / 2) * (n + mpf(7) / 4) * (n + 1) / (n + 2) ** 4
    S_head += t
S1 = S_head + nsum(hyp1term, [N0 + 1, inf], method="euler-maclaurin")
f4_256_rogers = log(256) - (mpf(24) / 256) * S1
check("[A4] EK4(i/sqrt2) = Rogers f4(256) = log 256 - (24/256) 5F4(...;1)",
      EK4_s2, f4_256_rogers, TOL)
print("     [A4] EK4(i/sqrt2)      =", EK4_s2)
print("     [A4] f2(256) closed form=", f2_256, "  |diff| =", abs(EK4_s2 - f2_256))

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
else:
    print("ALL CHECKS PASSED")
