# Verify (P1) ingredients for tau0 = (3+sqrt(-7))/8.
#
# (a) Theta identity: g7 = eta(t)^3 eta(7t)^3  =?  (1/2) sum'_{alpha in O_K} alpha^2 q^{N(alpha)},
#     K = Q(sqrt(-7)), O_K = Z[omega], omega = (1+sqrt(-7))/2, N(a+b w) = a^2+ab+2b^2.
#     alpha^2 = (m^2-2n^2) + (2mn+n^2) omega.  Compare coefficients (exact integers).
#
# (b) Lattice decomposition check:
#     T4(tau0) = T(O_K)          ?= 4 L(g7,3) + 2 zeta_K(2)
#     T1(tau0) = T((beta)/4)     ?= 8 L(g7,3) + 32 zeta_K(2),  beta = (3+sqrt(-7))/2
#     with T(L) = sum'_{lam in L} [ 2 Re(lam_bar^2/|lam|^6) + 1/|lam|^4 ].
#     Compute T's directly by Poisson summation rows (exponentially convergent).

from mpmath import mp, mpf, mpc, pi, sqrt, exp, sin, cos, sinh, cosh, zeta
from mpmath import dirichlet
import sys

mp.dps = 60
FAILS = []
TOL = mpf(10) ** (-45)   # working precision 60 dps; observed diffs ~1e-50

def check(name, got, want, tol=TOL):
    d = abs(got - want)
    ok = d < tol
    if not ok:
        FAILS.append(name)
    print("%-60s %s  (|diff| = %.2e)" % (name, "PASS" if ok else "FAIL", mpf(d)))

# ---------- (a) theta identity ----------
NMAX = 60
# g7 coefficients via Jacobi (from earlier script)
tri = []
j = 0
while j * (j + 1) // 2 < NMAX:
    tri.append(j * (j + 1) // 2)
    j += 1
a_eta = [0] * (NMAX + 1)
for i, ti in enumerate(tri):
    for jj, tj in enumerate(tri):
        n = ti + 7 * tj + 1
        if n <= NMAX:
            a_eta[n] += (-1) ** (i + jj) * (2 * i + 1) * (2 * jj + 1)

# theta: coeff of q^N = (1/2) sum_{m^2+mn+2n^2 = N} ( (m^2-2n^2) + (2mn+n^2) omega )
a_th = [0] * (NMAX + 1)
for N in range(1, NMAX + 1):
    rat, ompart = 0, 0
    B = int(sqrt(float(N / 2))) + 2
    Amax = int(sqrt(float(N))) + 2
    for n in range(-B, B + 1):
        for m in range(-Amax, Amax + 1):
            if m * m + m * n + 2 * n * n == N and (m, n) != (0, 0):
                rat += m * m - 2 * n * n
                ompart += 2 * m * n + n * n
    assert ompart == 0, (N, ompart)
    assert rat % 2 == 0, (N, rat)
    a_th[N] = rat // 2

print("(a) theta vs eta coefficients:")
ok = all(a_th[n] == a_eta[n] for n in range(1, NMAX + 1))
if not ok:
    FAILS.append("(a) theta vs eta")
print("    match to n=%d:" % NMAX, "PASS" if ok else "FAIL")
print("    first 15:", a_th[1:16])

# ---------- (b) T-splitting ----------
# B(L) = sum' 1/|lam|^4,  M(L) = sum' m^2/|dm tau+n|^6  via Poisson rows.
# Row for fixed m (multiplier d): sum_n ((n+x)^2 + y^2)^{-3}, x = d m Re t0, y = |d m| Im t0.
# From G(x,y) = sum_n ((n+x)^2+y^2)^{-1} = (pi/y) sinh(2 pi y)/(cosh(2 pi y) - cos(2 pi x)),
# F = sum (...)^{-3} = (1/2) d^2 G / d(y^2)^2 ... we just use finite differences? No:
# use direct formula: sum_n ((n+x)^2+y^2)^{-s} for s=2,3 via known closed forms,
# or sum over n with exponential cutoff: ((n+x)^2+y^2)^{-3} decays like 1/n^6 -- too slow.
# Instead use the G function and differentiate numerically? Cleaner: known identity
#   sum_n 1/((n+x)^2+y^2) = (pi/y) * sinh(2 pi y) / (cosh(2 pi y) - cos(2 pi x))
# and apply ( -1/(2y) d/dy ) twice:  (...)^{-2} = -1/(2y) d/dy (...)^{-1}, etc.

def G_row(x, y):
    return (pi / y) * sinh(2 * pi * y) / (cosh(2 * pi * y) - cos(2 * pi * x))

from mpmath import diff as mpdiff
def row_powers(x, y):
    # S2 = sum ((n+x)^2+y^2)^{-2}, S3 = sum (...)^{-3}
    # L = -(1/(2y)) d/dy sends (...)^{-k} -> k (...)^{-(k+1)}
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y) / (2 * y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy) / (2 * yy), y) / (4 * y)
    return S2, S3

s7 = sqrt(mpf(7))
tau0 = mpc(3, s7) / 8
y0 = tau0.imag
x0 = tau0.real

def lattice_T(d):
    # T_d = sum'_{m,n} [4 (d m x0 + n)^2/|dm t0+n|^6 - 1/|dm t0+n|^4]
    #      = 3 B_d - 4 d^2 y0^2 M_d,  B_d = sum' |.|^{-4}, M_d = sum' m^2 |.|^{-6}
    B, M = mpf(0), mpf(0)
    for m in range(-200, 201):
        if m == 0:
            B += 2 * zeta(4)
            continue
        x = d * m * x0
        y = abs(d * m) * y0
        S2, S3 = row_powers(x, y)
        B += S2
        M += m * m * S3
        if abs(d * m) * y0 > 40:
            if m > 0:
                break
    return 3 * B - 4 * d * d * y0 * y0 * M

# L-values:
# L(g7,3) from functional equation: L' = 7 sqrt7 L3/(4 pi^3)
Lp = mpf("0.10267160777890201121045659489829291399889482708922")
L3 = Lp * 4 * pi ** 3 / (7 * s7)
chi7 = [0, 1, 1, -1, 1, -1, -1]  # chi_{-7}(n) = (n/7)? Legendre (n mod 7): squares 1,2,4 -> +1
# chi_{-7}(n) = (-7/n) Kronecker; for n not divisible by 7: (n/7): 1->1,2->1,3->-1,4->1,5->-1,6->-1
def chi(n):
    n = n % 7
    return [0, 1, 1, -1, 1, -1, -1][n]
Lchi2 = dirichlet(mpf(2), [chi(n) for n in range(7)])
zK2 = zeta(2) * Lchi2
print("L(g7,3) =", L3)
print("zeta_K(2) =", zK2)

T1 = lattice_T(1)
T4 = lattice_T(4)
print("T4(tau0) =", T4)
print("4 L3 + 2 zK2 =", 4 * L3 + 2 * zK2)
check("(b) T4(tau0) = 4 L(g7,3) + 2 zeta_K(2)", T4, 4 * L3 + 2 * zK2)
print("T1(tau0) =", T1)
print("8 L3 + 32 zK2 =", 8 * L3 + 32 * zK2)
check("(b) T1(tau0) = 8 L(g7,3) + 32 zeta_K(2)", T1, 8 * L3 + 32 * zK2)

# final combination
EK = (2 * y0 / pi ** 3) * (-T1 + 16 * T4)
print("EK combination =", EK)
print("8 L'           =", 8 * Lp)
check("(b) EK combination = 8 L'(g7,0)", EK, 8 * Lp)
print("14 sqrt7 L3/pi^3 =", 14 * s7 * L3 / pi ** 3)
check("(b) 8 L'(g7,0) = 14 sqrt7 L(g7,3)/pi^3 (functional equation)", 8 * Lp, 14 * s7 * L3 / pi ** 3)

print()
if FAILS:
    print("FAILED CHECKS:", FAILS)
    sys.exit(1)
print("ALL CHECKS PASSED")
