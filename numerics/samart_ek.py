# Samart Prop 2.1(i) for f2(k) = 2 m( (x+1/x)(y+1/y)(z+1/z) + sqrt(k) ).
# Parametrization: s2(q) = -Delta(tau+1/2)/Delta(2 tau+1)  (= Weber f(2 tau)^24).
# For Im tau >= 1/2:
#   f2(s2(q(tau))) = Im[ 2 pi tau + (2/pi^3) ( U1(tau) - 4 U4(tau) ) ],
#   U_j(tau) = sum_{n in Z} sum_{m != 0} (1/m) (j m tau + n)^{-3}
#            = 2 pi^3 sum_{m>=1} (1/m) cos(j pi m tau) / sin^3(j pi m tau).
#
# Step 1: validate at tau = i/2  -> s2 = 64, expect f2(64) = 8 L'(h,0), h=eta(4t)^6,
#         and also compare with direct torus integration of 2 m(f + 8).
# Step 2: evaluate at disc -7 CM points where s2 = 1.

from mpmath import mp, mpf, mpc, pi, exp, log, sqrt, sin, cos, quad, binomial

mp.dps = 50

def eta(tau, nterms=300):
    q = exp(2 * pi * 1j * tau)
    p = mpc(1)
    qn = q
    for n in range(1, nterms + 1):
        p *= (1 - qn)
        qn *= q
        if abs(qn) < mpf(10) ** (-60):
            break
    return exp(pi * 1j * tau / 12) * p

def Delta(tau):
    return eta(tau) ** 24

def s2(tau):
    return -Delta(tau + mpf(1) / 2) / Delta(2 * tau + 1)

def U(j, tau):
    tot = mpc(0)
    for m in range(1, 2000):
        z = j * pi * m * tau
        term = cos(z) / sin(z) ** 3 / m
        tot += term
        if abs(term) < mpf(10) ** (-60):
            break
    return 2 * pi ** 3 * tot

def f2_EK(tau):
    val = 2 * pi * tau + (2 / pi ** 3) * (U(1, tau) - 4 * U(4, tau))
    return val.imag

# ---------- Step 1: tau = i/2 ----------
tau = mpc(0, mpf(1) / 2)
print("s2(i/2) =", s2(tau), " (expect 64)")
print("f2 via EK =", f2_EK(tau))

# direct: 2 m(f+8); |8 c1 c2 c3| <= 8 so Jensen branch log((8+sqrt(64-A^2))/2) everywhere
def inner8(t1):
    c = cos(t1)
    return quad(lambda t2: log((8 + sqrt(64 - 64 * c * c * cos(t2) ** 2)) / 2), [0, pi / 2])
m8 = 4 / pi ** 2 * quad(inner8, [0, pi / 2])
print("2 m(f+8) direct =", 2 * m8)

# L-value: h(tau) = eta(4 tau)^6, weight 3, level 16, char chi_{-4}
# a_n: h = q prod (1-q^{4n})^6 ; prod (1-q^n)^6 = (sum (-1)^j (2j+1) q^{T_j})^2
NMAX = 300
tri = []
j = 0
while j * (j + 1) // 2 < NMAX:
    tri.append(j * (j + 1) // 2)
    j += 1
A = [0] * (NMAX + 1)
for i, ti in enumerate(tri):
    for jj, tj in enumerate(tri):
        n = ti + tj
        if n <= NMAX:
            A[n] += (-1) ** (i + jj) * (2 * i + 1) * (2 * jj + 1)
a = [0] * (4 * NMAX + 2)
for n in range(NMAX + 1):
    a[4 * n + 1] = A[n]
print("h coeffs a1..a17:", [a[n] for n in range(1, 18)])

N16, k = 16, 3
xN = sqrt(mpf(N16))
def I(s):
    s = mpf(s)
    return mpf(sum(a[n] * (2 * pi * n) ** (-s) * float(0) for n in []))  # placeholder
from mpmath import gammainc, power
def Ireal(s):
    s = mpf(s)
    tot = mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0:
            continue
        x = 2 * pi * n / xN
        tot += a[n] * power(2 * pi * n, -s) * gammainc(s, x)
    return tot

# root number: Fricke f(-1/(16 tau)) = w 16^{3/2} (-i tau)^3 f(tau)
def hf(tau):
    return eta(4 * tau) ** 6
y = mpf("0.35")
C = hf(1j / (16 * y)) / (y ** 3 * hf(1j * y))
w = C / xN ** 3
print("root number w for h:", w)  # expect real +-1

def Lam(s):
    return xN ** s * Ireal(s) + w * xN ** (k - s) * Ireal(k - s)
print("FE check |Lam(1)-Lam(2)|:", abs(Lam(1) - Lam(2)))
Lam3 = Lam(3)
Lp0 = w * Lam3
print("L'(h,0) =", Lp0)
print("8 L'(h,0) =", 8 * Lp0)

# ---------- Step 2: disc -7 CM points ----------
s7 = sqrt(mpf(7))
for name, ta in [("(1+s-7)/4", mpc(1, s7) / 4), ("(3+s-7)/4", mpc(3, s7) / 4)]:
    sv = s2(ta)
    print(name, ": s2 =", sv)
    if abs(sv - 1) < mpf("1e-20"):
        print("   f2(1) via EK =", f2_EK(ta))
