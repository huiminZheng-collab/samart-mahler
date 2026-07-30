# verify_P1_n2pair.py
#
# (P1)-type exact CM evaluation for Samart's Table-4 target
#     n2((47 +/- 45 sqrt(-7))/2) = (4/7)(54 M7 + d7),
# M7 = L'(g7,0), g7 = eta(t)^3 eta(7t)^3, d7 = L'(chi_{-7},-1).
#
# The evaluation is done at the Fricke partner tau_w = -1/(4 tau') = (-1+sqrt(-7))/4
# (Im tau_w = sqrt7/4 >= 1/2, inside Samart's domain), NOT at tau' = (1+sqrt(-7))/8,
# whose EK value is the wrong-sheet branch value (8/7)(44 M7 - d7) -- see P1_n2pair.md.
#
# Every intermediate identity is verified numerically (60 dps working precision);
# each line prints PASS/FAIL.  Derivations are in P1_n2pair.md.
#
# Checks:
#  [T1] theta identity  g7 = (1/2) sum'_{a in O_K} a^2 q^{N(a)}   (exact integers, q^60)
#  [L1] L'(g7,0) via incomplete-gamma series + root number; FE self-check Lam(1)=Lam(2)
#  [L2] d7 = L'(chi_{-7},-1): direct derivative vs functional equation (7 sqrt7/4pi) L(chi,2);
#       L(chi_{-7},-1) = 0;  zeta_K(2) = (2 pi^3/(21 sqrt7)) d7
#  [O1] B(O2) = sum'_{O2} |g|^{-4} = (5/4) zeta_K(2)
#       (a) direct Epstein zeta (closed form + Bessel tail)
#       (b) Euler-factor decomposition 2 zeta_K(2) [(1-2^-2)^2 + 2^-4]
#       (c) Glasser--Zucker 2(1-2^{1-s}+2^{1-2s}) zeta(s) L(chi_{-7},s) at s=2
#  [O2] G(O2) = sum'_{O2} g_bar^2 |g|^{-6} = 3 L(g7,3)
#       (a) direct lattice sum  (b) 2 L(g7,3) [(1-pi^2/8)(1-pibar^2/8) + 4^-2]
#  [O3] ray-class decomposition: sum'_{g = 1 mod 2} |g|^{-4} = (9/8) zeta_K(2),
#       sum'_{g = 1 mod 2} g^2 N(g)^{-3} = (23/8) L(g7,3)  (via shifted-coset sums
#       Sh_B = (3/4) zeta_K(2), Sh_G = -L(g7,3))
#  [T4] T4(tau_w) = 6 L(g7,3) + (5/4) zeta_K(2)   (Poisson-row lattice sum)
#  [T5] T1(tau_w) = -12 L(g7,3) + 8 zeta_K(2)
#  [T6] same at tau_w' = (1+sqrt(-7))/4 = pi/2
#  [EK] -T1 + 16 T4 = 108 L(g7,3) + 12 zeta_K(2);
#       EK(tau_w) = (2 Im tau_w/pi^3)(-T1+16T4) = (4/7)(54 M7 + d7);
#       independent EK via U-sum formula; f2(k) via Rogers' 5F4.
#  [FR] Fricke invariance s2(-1/(4t)) = s2(t) (numeric, several points);
#       s2(tau_w) = (47+45 sqrt(-7))/2, s2(tau_w') = conjugate (exact via cert0_s2_n2pair + FR).
#  [WS] wrong-sheet companion: EK(tau') = (8/7)(44 M7 - d7),
#       T4(tau') = 4 L3 + 2 zK2, T1(tau') = -288 L3 + 80 zK2.

from mpmath import (mp, mpf, mpc, pi, sqrt, exp, sin, cos, sinh, cosh, zeta,
                    dirichlet, gamma, diff as mpdiff, hyper, log)

mp.dps = 60
s7 = sqrt(mpf(7))
FAILS = []

def check(name, got, want, tol):
    d = abs(got - want)
    ok = d < tol
    if not ok:
        FAILS.append(name)
    print("%-72s %s  (|diff| = %.2e)" % (name, "PASS" if ok else "FAIL", mpf(d)))

TOL = mpf(10) ** (-45)
TOLM = mpf(10) ** (-40)

# ================= g7 coefficients, L-values =================
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

from mpmath import gammainc, power
xN = sqrt(mpf(7))

def I(s):
    s = mpf(s)
    tot = mpf(0)
    for n in range(1, NMAX + 1):
        if a[n] == 0:
            continue
        tot += a[n] * power(2 * pi * n, -s) * gammainc(s, 2 * pi * n / xN)
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

# root number via Fricke: f(-1/(7 tau)) = w 7^{3/2} (-i tau)^3 f(tau)
y = mpf("0.6")
w = (g7f(1j / (7 * y)) / (y ** 3 * g7f(1j * y))) / xN ** 3
check("[L1] root number w = +1", w, 1, TOLM)

def Lam(s):
    return xN ** s * I(s) + w * xN ** (3 - s) * I(3 - s)

check("[L1] functional equation Lam(1) = Lam(2)", Lam(1), Lam(2), TOL)
Lam3 = Lam(3)
M7 = w * Lam3                      # L'(g7,0) = w * Lambda(3)
L3 = Lam3 * (2 * pi) ** 3 / (xN ** 3 * gamma(3))   # L(g7,3) from Lambda(3)
M7_ref = mpf("0.10267160777890201121045659489829291399889482708922")
check("[L1] M7 = L'(g7,0) matches reference value", M7, M7_ref, TOL)
check("[L1] FE identity M7 = 7 sqrt7 L(g7,3) / (4 pi^3)", M7, 7 * s7 * L3 / (4 * pi ** 3), TOL)

chi = [0, 1, 1, -1, 1, -1, -1]
Lchi2 = dirichlet(mpf(2), chi)
zK2 = zeta(2) * Lchi2
check("[L2] L(chi_{-7},-1) = 0", dirichlet(mpf(-1), chi), 0, TOL)
d7 = mpdiff(lambda s: dirichlet(s, chi), mpf(-1))
d7_fe = (7 * s7 / (4 * pi)) * Lchi2
check("[L2] d7 = L'(chi_{-7},-1) = (7 sqrt7 / 4 pi) L(chi_{-7},2)", d7, d7_fe, TOL)
check("[L2] zeta_K(2) = (2 pi^3 / (21 sqrt7)) d7", zK2, (2 * pi ** 3 / (21 * s7)) * d7, TOL)

# ================= [T1] theta identity =================
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
print("%-72s %s" % ("[T1] theta identity g7 = (1/2) sum' a^2 q^{N(a)} to q^60", "PASS" if ok else "FAIL"))
if not ok:
    FAILS.append("T1")

# ================= O2 lattice sums (closed form + Bessel tail) =================
def K12(x):
    return sqrt(pi / (2 * x)) * exp(-x)

def K32(x):
    return K12(x) * (1 + 1 / x)

def K52(x):
    return K12(x) * (1 + 3 / x + 3 / x ** 2)

C7 = 2 * pi * s7   # 2 pi sqrt7

def epstein2():
    # sum'_{(x,y)} (x^2+7y^2)^-2
    D = mpf(0)
    for yy in range(1, 20):
        for n in range(1, 20):
            t = yy ** mpf("-1.5") * n ** mpf("1.5") * K32(C7 * n * yy)
            D += t
    return (2 * zeta(4) + 2 * sqrt(pi) * 7 ** mpf("-1.5") * gamma(mpf("1.5")) * zeta(3)
            + 8 * pi ** 2 * 7 ** mpf("-0.75") * D)   # factor 2: y and -y rows

def G_O2_direct():
    # sum' (x^2 - 7y^2)/(x^2+7y^2)^3 = epstein2 - 14 * C3
    D = mpf(0)
    for yy in range(1, 20):
        for n in range(1, 20):
            D += yy ** mpf("-0.5") * n ** mpf("2.5") * K52(C7 * n * yy)
    C3 = 3 * pi * zeta(3) / (4 * 7 ** mpf("2.5")) + 4 * pi ** 3 * 7 ** mpf("-1.25") * D  # factor 2: +/-y
    return epstein2() - 14 * C3

def shifted_B():
    # sum_{u,v in Z} ((u+1/2)^2 + 7(v+1/2)^2)^-2
    D = mpf(0)
    for v in range(0, 40):
        for n in range(1, 40):
            D += (v + mpf("0.5")) ** mpf("-1.5") * (-1) ** n * n ** mpf("1.5") * K32(C7 * n * (v + mpf("0.5")))
    return pi * zeta(3) / s7 + 8 * pi ** 2 * 7 ** mpf("-0.75") * D   # factor 2: +/-(v+1/2)

def shifted_G():
    # sum_{u,v} [((u+1/2)^2+7(v+1/2)^2)^-2 - 14 (v+1/2)^2 ((u+1/2)^2+7(v+1/2)^2)^-3]
    D = mpf(0)
    for v in range(0, 40):
        for n in range(1, 40):
            D += (v + mpf("0.5")) ** mpf("-0.5") * (-1) ** n * n ** mpf("2.5") * K52(C7 * n * (v + mpf("0.5")))
    C3sh = 3 * pi * zeta(3) / (4 * 7 ** mpf("1.5")) + 4 * pi ** 3 * 7 ** mpf("-1.25") * D  # factor 2
    return shifted_B() - 14 * C3sh

B_O2 = epstein2()
G_O2 = G_O2_direct()
check("[O1a] B(O2) = sum'_{O2} |g|^-4  (direct Epstein) = (5/4) zeta_K(2)", B_O2, mpf(5) / 4 * zK2, TOL)
check("[O1b] B(O2) = 2 zeta_K(2) [(1-2^-2)^2 + 2^-4]", B_O2, 2 * zK2 * ((1 - mpf(1)/4) ** 2 + mpf(1)/16), TOL)
check("[O1c] B(O2) = Glasser--Zucker 2(1-2^{1-s}+2^{1-2s}) zeta(s) L(chi,s) at s=2",
      B_O2, 2 * (1 - mpf(1)/2 + mpf(1)/8) * zK2, TOL)
check("[O2a] G(O2) = sum' g_bar^2 |g|^-6  (direct) = 3 L(g7,3)", G_O2, 3 * L3, TOL)

# exact algebraic factor: (1-pi^2/8)(1-pibar^2/8) = 23/16
pii = mpc(1, s7) / 2
pib = mpc(1, -s7) / 2
fac = (1 - pii ** 2 / 8) * (1 - pib ** 2 / 8)
check("[O2b] (1-pi^2/8)(1-pibar^2/8) = 23/16  (exact in K)", fac, mpf(23)/16, TOL)
check("[O2b] pi^2 + pibar^2 = -3, N(pi) = 2", pii ** 2 + pib ** 2, -3, TOL)
check("[O2b] G(O2) = 2 L(g7,3)[(1-pi^2/8)(1-pibar^2/8) + 4^-2]", G_O2, 2 * L3 * (fac + mpf(1)/16), TOL)

Sh_B = shifted_B()
Sh_G = shifted_G()
check("[O3] shifted coset Sh_B = (3/4) zeta_K(2)", Sh_B, mpf(3)/4 * zK2, TOL)
check("[O3] shifted coset Sh_G = -L(g7,3)", Sh_G, -L3, TOL)
ray_B = B_O2 - (B_O2 + Sh_B) / 16
ray_G = G_O2 - (G_O2 + Sh_G) / 16
check("[O3] ray class sum'_{g=1 mod 2} N(g)^-2 = (9/8) zeta_K(2) = 2 zeta_K(2)(1-2^-2)^2",
      ray_B, mpf(9)/8 * zK2, TOL)
check("[O3] ray class sum'_{g=1 mod 2} g^2 N(g)^-3 = (23/8) L(g7,3) = 2 L (1-pi^2/8)(1-pibar^2/8)",
      ray_G, mpf(23)/8 * L3, TOL)

# ================= T-sums via Poisson rows =================
def G_row(x, y):
    return (pi / y) * sinh(2 * pi * y) / (cosh(2 * pi * y) - cos(2 * pi * x))

def row_powers(x, y):
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y) / (2 * y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy) / (2 * yy), y) / (4 * y)
    return S2, S3

def lattice_T(tau, d):
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

tau_w = mpc(-1, s7) / 4
tau_wp = mpc(1, s7) / 4
tau_p = mpc(1, s7) / 8

T4w = lattice_T(tau_w, 4)
T1w = lattice_T(tau_w, 1)
check("[T4] T4(tau_w) = 6 L(g7,3) + (5/4) zeta_K(2)", T4w, 6 * L3 + mpf(5)/4 * zK2, TOL)
check("[T5] T1(tau_w) = -12 L(g7,3) + 8 zeta_K(2)", T1w, -12 * L3 + 8 * zK2, TOL)
combw = -T1w + 16 * T4w
check("[EK] -T1 + 16 T4 = 108 L(g7,3) + 12 zeta_K(2)", combw, 108 * L3 + 12 * zK2, TOL)

T4wp = lattice_T(tau_wp, 4)
T1wp = lattice_wp = lattice_T(tau_wp, 1)
check("[T6] T4(tau_w') = T4(tau_w)", T4wp, T4w, TOL)
check("[T6] T1(tau_w') = T1(tau_w)", T1wp, T1w, TOL)

# ================= EK: final assembly =================
EK_w = (2 * tau_w.imag / pi ** 3) * combw
target = mpf(4) / 7 * (54 * M7 + d7)
check("[EK] EK(tau_w) = (4/7)(54 M7 + d7)", EK_w, target, TOL)
# step-by-step FE conversions (exact algebra, numeric check):
check("[EK] (sqrt7/2pi^3)*108 = (4/7)*54 * 7 sqrt7/(4 pi^3)   [coefficient of L(g7,3)]",
      (s7 / (2 * pi ** 3)) * 108, mpf(4)/7 * 54 * (7 * s7 / (4 * pi ** 3)), TOL)
check("[EK] (sqrt7/2pi^3)*12 = (4/7) * 21 sqrt7/(2 pi^3)     [coefficient of zeta_K(2)]",
      (s7 / (2 * pi ** 3)) * 12, mpf(4)/7 * (21 * s7 / (2 * pi ** 3)), TOL)
check("[EK] zeta_K(2) * 21 sqrt7/(2 pi^3) = d7", zK2 * 21 * s7 / (2 * pi ** 3), d7, TOL)

# independent EK via U-sums (Samart's formula)
def U(j, tau):
    tot = mpc(0)
    for m in range(1, 2000):
        z = j * pi * m * tau
        term = cos(z) / sin(z) ** 3 / m
        tot += term
        if abs(term) < mpf(10) ** (-62):
            break
    return 2 * pi ** 3 * tot

def EK_direct(tau):
    return (2 * pi * tau + (2 / pi ** 3) * (U(1, tau) - 4 * U(4, tau))).imag

check("[EK] U-sum EK(tau_w) = T-sum EK(tau_w)", EK_direct(tau_w), EK_w, TOL)

# ================= Fricke invariance & s2 values =================
def s2(tau):
    return -eta(tau + mpf(1)/2) ** 24 / eta(2 * tau + 1) ** 24

for t in [mpc(1, s7) / 8, mpc("0.3", "0.9"), mpc(0, mpf(1)/3), mpc("0.7", "1.1")]:
    check("[FR] s2(-1/(4 t)) = s2(t) at t = %s" % t, s2(-1 / (4 * t)), s2(t), TOLM)

k = mpc(47, 45 * s7) / 2
kb = mpc(47, -45 * s7) / 2
check("[FR] s2(tau_w) = (47 + 45 sqrt(-7))/2", s2(tau_w), k, TOL)
check("[FR] s2(tau_w') = (47 - 45 sqrt(-7))/2", s2(tau_wp), kb, TOL)
check("[FR] tau_w = -1/(4 tau')", tau_w, -1 / (4 * tau_p), TOL)

# ================= Mahler-measure value (Rogers 5F4) =================
f2 = (log(k) - (8 / k) * hyper([mpf(3)/2, mpf(3)/2, mpf(3)/2, 1, 1],
                               [2, 2, 2, 2], 64 / k)).real
check("[EK] f2(k) = Re(log k - (8/k) 5F4(64/k)) = (4/7)(54 M7 + d7)", f2, target, TOLM)
check("[EK] EK(tau_w) = f2(s2(tau_w))  (Mahler interpretation)", EK_w, f2, TOLM)

# ================= wrong-sheet companion at tau' =================
T4p = lattice_T(tau_p, 4)
T1p = lattice_T(tau_p, 1)
check("[WS] T4(tau') = 4 L(g7,3) + 2 zeta_K(2)", T4p, 4 * L3 + 2 * zK2, TOL)
check("[WS] T1(tau') = -288 L(g7,3) + 80 zeta_K(2)", T1p, -288 * L3 + 80 * zK2, TOL)
combp = -T1p + 16 * T4p
check("[WS] -T1+16T4 at tau' = 352 L(g7,3) - 48 zeta_K(2)", combp, 352 * L3 - 48 * zK2, TOL)
EK_p = (2 * tau_p.imag / pi ** 3) * combp
ws_val = mpf(8) / 7 * (44 * M7 - d7)
check("[WS] EK(tau') = (8/7)(44 M7 - d7)  (wrong-sheet branch value)", EK_p, ws_val, TOL)
check("[WS] U-sum EK(tau') agrees", EK_direct(tau_p), EK_p, TOL)
print("     [WS] note: EK(tau') != f2(s2(tau')) = target;  |diff| =",
      abs(EK_p - target))

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
else:
    print("ALL CHECKS PASSED")
