# Compute L'(g7, 0) for g7(tau) = eta(tau)^3 eta(7 tau)^3
# g7 = q * prod (1-q^n)^3 (1-q^{7n})^3, weight k=3, level N=7, CM by Q(sqrt(-7)).
#
# Jacobi: prod (1-q^n)^3 = sum_{j>=0} (-1)^j (2j+1) q^{T_j}, T_j = j(j+1)/2.
# So a_n = sum_{T_j + 7 T_k = n-1} (-1)^{j+k} (2j+1)(2k+1).
#
# Mellin split at y = 1/sqrt(N), with Fricke f(-1/(N tau)) = w N^{k/2} tau^k f(tau):
#   Lambda(s) = N^{s/2} I(s) + w i^k N^{(k-s)/2} I(k-s),
#   I(s) = sum_n a_n (2 pi n)^{-s} Gamma(s, 2 pi n / sqrt(N)).
# Determine w from  Lambda(1) = w Lambda(2):
#   => w = I(1) / (sqrt(N) I(2))   (after cancelling (1+i) factors; see notes)
# Then Lambda(0) = w Lambda(3) and Gamma(s)~1/s forces L(0)=0, L'(0) = w Lambda(3).

from mpmath import mp, mpf, mpc, pi, sqrt, gammainc, power, fac
import itertools

mp.dps = 50
N = 7
k = 3

NMAX = 200  # plenty: incomplete gamma kills terms past n ~ 40 at 50 dps

# triangular numbers T_j < NMAX
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

# sanity: a_p = 0 for p inert mod 7? (CM check: a_2? 2 is inert in Q(sqrt-7)? (2)=split actually since -7 ≡ 1 mod 8... )
print("a_1..a_20:", a[1:21])
# multiplicativity spot check: a_2*a_3 == a_6 + chi(3)*3^{k-1}*a_2 ... skip; check a_4 vs a_2^2 - eps(p) p^2
print("check a4 - (a2^2 - chi(2)*4):", a[4] - (a[2] ** 2), "(should be chi(2)*4 if a2!=0; chi_7(2)=1 -> expect 4)")
print("check a9 - (a3^2 - chi(3)*9):", a[9] - (a[3] ** 2))

xN = sqrt(mpf(N))

def I(s):
    s = mpf(s)
    tot = mpf(0)
    for n in range(1, NMAX + 1):
        if a[n] == 0:
            continue
        x = 2 * pi * n / xN
        term = a[n] * power(2 * pi * n, -s) * gammainc(s, x)  # upper incomplete Gamma(s, x)
        tot += term
    return tot

I1, I2 = I(1), I(2)
print("I(1) =", I1)
print("I(2) =", I2)

# ---- determine root number w directly from the Fricke involution ----
# f(-1/(N*tau)) = w N^{k/2} tau^k f(tau).  Evaluate at tau = i*y with q-products.
def eta(tau, nterms=400):
    from mpmath import exp
    q = exp(2 * pi * 1j * tau)
    p = mpc(1)
    qn = q
    for n in range(1, nterms + 1):
        p *= (1 - qn)
        qn *= q
    return exp(pi * 1j * tau / 12) * p

def f(tau):
    return eta(tau) ** 3 * eta(N * tau) ** 3

y = mpf("0.6")
# correct Fricke convention: f(-1/(N*tau)) = w N^{k/2} (-i tau)^k f(tau)
C = f(1j / (N * y)) / ((y) ** k * f(1j * y))   # (-i)(i y) = y
w_unit = C / xN ** k
print("w =", w_unit, " (expect +1 or -1, real)")

# With this convention the Mellin split gives NO i^k factor:
#   Lambda(s) = N^{s/2} I(s) + w N^{(k-s)/2} I(k-s)
def Lam(s):
    return xN ** s * I(s) + w_unit * xN ** (k - s) * I(k - s)

Lam1, Lam2 = Lam(1), Lam(2)
print("Lambda(1) =", Lam1)
print("Lambda(2) =", Lam2)
print("FE check |Lam(1)-Lam(2)| =", abs(Lam1 - Lam2))

I0, I3 = I(0), I(3)
Lam3 = xN ** 3 * I3 + w_unit * I0
print("Lambda(3) =", Lam3)
Lp0 = w_unit * Lam3
print("L'(g7,0) =", Lp0)
print("8 L'(g7,0) =", 8 * Lp0)
