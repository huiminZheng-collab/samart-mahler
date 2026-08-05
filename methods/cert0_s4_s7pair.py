# Cert-0: rigorous proof that
#     s4(sqrt(-7))   = 8292456 + 3132675 sqrt7,
#     s4(sqrt(-7)/2) = 8292456 - 3132675 sqrt7,
# where s4(tau) = (eta(2t)/eta(t))^24 (16 W^4 + W^-4)^4,
# W = eta(t) eta(4t)^2 / eta(2t)^3  (Samart's signature-4 modular function).
#
# Chain:
#  (i)   s4 has integral q-expansion s4 = q^-2 + 104 q^-1 + ... in Z((q))
#        [formal series check to q^30 below, exact integer arithmetic].
#        For pure imaginary tau the value is real (real q-series).
#  (ii)  [QUOTED CM THEORY] s4 is an eta-quotient, hence a modular unit of
#        level 8 with integral q-expansions at all cusps; by Shimura's CM
#        theorem its values at CM points are algebraic integers, and
#        Shimura reciprocity identifies the Galois orbit.  In particular
#        s4(sqrt(-7)) and s4(sqrt(-7)/2) are conjugate algebraic integers
#        generating Q(sqrt7).  => their sum S and product P are rational
#        integers, and each is a root of X^2 - S X + P.
#  (iii) Interval arithmetic (mpmath.iv, iv.dps = 60, REAL intervals only:
#        q = e^{-2 pi Im tau} > 0 makes every eta factor real) with
#        rigorous product tails gives
#          s4(sqrt(-7)) + s4(sqrt(-7)/2) in 16584912 +- 1e-20,
#          s4(sqrt(-7)) * s4(sqrt(-7)/2) in 69257922561 +- 1e-10,
#        hence S = 16584912, P = 69257922561 exactly (integers).
#        Exact integer check: 8292456^2 - 7*3132675^2 = 69257922561,
#        so X^2 - S X + P = (X - r+)(X - r-) with
#        r+- = 8292456 +- 3132675 sqrt7, and r+ - r- > 1.6e7.
#  (iv)   Individual locks: |s4(sqrt(-7)) - r+| < 1e-20 and
#        |s4(sqrt(-7)/2) - r-| < 1e-20; since the two roots are 1.6e7
#        apart, each value equals its root.  QED (modulo the quoted
#        CM step (ii), which is the only non-self-contained input).

from fractions import Fraction as Fr
from math import comb
from mpmath import mp, mpf, sqrt, iv
import sys

mp.dps = 60
iv.dps = 60

FAILS = []
def report(name, cond):
    if not cond:
        FAILS.append(name)
    return 'PASS' if cond else 'FAIL'

def cv(z):
    return mp.convert(z)

# ---------------- (i) formal q-series of s4 to q^30 (exact integers) ------
N = 30
def Pser(d, N):
    """prod (1 - q^{d n}) as series."""
    res = [0]*(N+1)
    res[0] = 1
    k = 1
    while True:
        e1 = k*(3*k-1)//2*d
        e2 = k*(3*k+1)//2*d
        if e1 > N and e2 > N:
            break
        if e1 <= N: res[e1] += (-1)**k
        if e2 <= N: res[e2] += (-1)**k
        k += 1
    return res

def s_mul(A, B, N):
    C = [0]*(N+1)
    for i in range(N+1):
        if A[i] == 0: continue
        for j in range(N+1-i):
            if B[j]:
                C[i+j] += A[i]*B[j]
    return C

def s_pow(A, e, N):
    R = [0]*(N+1)
    R[0] = 1
    for _ in range(e):
        R = s_mul(R, A, N)
    return R

def s_inv(A, N):
    assert A[0] == 1
    R = [0]*(N+1)
    R[0] = 1
    for n in range(1, N+1):
        R[n] = -sum(A[k]*R[n-k] for k in range(1, n+1))
    return R

# W^4 = q prod (1-q^n)^4 (1-q^{4n})^8 / (1-q^{2n})^12
W4 = s_mul(s_mul(s_pow(Pser(1, N), 4, N), s_pow(Pser(4, N), 8, N), N),
           s_inv(s_pow(Pser(2, N), 12, N), N), N)          # series of W^4 / q
W4i = s_inv(W4, N)                                         # series of q / W^4
# 16 W^4 + W^-4 = q^-1 (16 q^2 W4 + W4i)
V1 = [16*W4[i-2] if i >= 2 else 0 for i in range(N+1)]
for i in range(N+1):
    V1[i] += W4i[i]
V = s_pow(V1, 4, N)
# s4 = q^2 prod((1-q^{2n})/(1-q^n))^24 * q^-4 V = q^-2 * s4ser
eta_ratio = s_mul(s_pow(Pser(2, N), 24, N), s_inv(s_pow(Pser(1, N), 24, N), N), N)
s4ser = s_mul(eta_ratio, V, N)
ok_i = all(isinstance(c, int) for c in s4ser) and s4ser[0] == 1
print('(i) s4 = q^-2 * (1 + %d q + %d q^2 + ...), integral coeffs to q^%d: %s'
      % (s4ser[1], s4ser[2], N, report('(i) integral q-series', ok_i)))

# ---------------- (ii)/(iii) exact minimal polynomial data -----------------
A0 = 8292456
B0 = 3132675
S_int = 2*A0
P_int = A0*A0 - 7*B0*B0
print('(iii) sum S =', S_int, ' product P =', P_int,
      ' (P = 3^12 * 19^4:', P_int == 3**12 * 19**4, ')')
ok_iii = (P_int == 69257922561)
print('     8292456^2 - 7*3132675^2 = 69257922561: %s'
      % report('(iii) exact norm', ok_iii))
# root separation: r+ - r- = 2*3132675*sqrt7 > 1.6e7
ok_sep = (2*B0)**2 * 7 > (16*10**6)**2
print('     r+ - r- = 2*3132675 sqrt7 > 1.6e7 (exact): %s'
      % report('(iii) root separation', ok_sep))

# ---------------- (iii)/(iv) interval evaluations --------------------------
def eta_iv(y, NT=80):
    """eta(i y) = e^{-pi y/12} prod (1 - q^n), q = e^{-2 pi y}, real iv with
    rigorous tail: |sum_{n>N} log(1-q^n)| <= sum_{n>N} q^n/(1-q^n)
    <= q^{N+1}/((1-q)(1-q^{N+1}))."""
    q = iv.exp(-2*iv.pi*y)
    r = q.b
    P = iv.mpf(1)
    qn = q
    for n in range(1, NT+1):
        P *= (1-qn)
        qn *= q
    E = (iv.mpf(r)**(NT+1)/((1-iv.mpf(r))*(1-iv.mpf(r)**(NT+1)))).b
    E = max(E, iv.mpf('1e-90'))
    return iv.exp(-iv.pi*y/12) * P * iv.exp(iv.mpf([-E, E]))

def s4_iv(y):
    e1 = eta_iv(y)
    e2 = eta_iv(2*y)
    e4 = eta_iv(4*y)
    W = e1*e4**2/e2**3
    return (e2/e1)**24 * (16*W**4 + W**(-4))**4

sA_iv = s4_iv(iv.sqrt(7))
sB_iv = s4_iv(iv.sqrt(7)/2)
def ivmid(z):
    return (cv(z.a) + cv(z.b))/2

print('(iii) s4(sqrt(-7))   ~', mp.nstr(ivmid(sA_iv), 25))
print('(iii) s4(sqrt(-7)/2) ~', mp.nstr(ivmid(sB_iv), 25))

sum_iv = sA_iv + sB_iv
prod_iv = sA_iv * sB_iv
dS = max(abs(cv(sum_iv.a) - S_int), abs(cv(sum_iv.b) - S_int))
dP = max(abs(cv(prod_iv.a) - P_int), abs(cv(prod_iv.b) - P_int))
print('     |sum - 16584912| <= %s ;  |product - 69257922561| <= %s'
      % (mp.nstr(dS, 3), mp.nstr(dP, 3)))
ok_lockSP = dS < mpf('1e-15') and dP < mpf('1e-5')
print('     => S = 16584912, P = 69257922561 in Z: %s'
      % report('(iii) sum/product integrality locks', ok_lockSP))

rp = A0 + B0*iv.sqrt(7)
rm = A0 - B0*iv.sqrt(7)
dA = cv(abs(sA_iv - rp).b)
dB = cv(abs(sB_iv - rm).b)
print('(iv)  |s4(sqrt(-7)) - r+|   <= %s' % mp.nstr(dA, 3))
print('(iv)  |s4(sqrt(-7)/2) - r-| <= %s' % mp.nstr(dB, 3))
ok_iv = dA < mpf('1e-20') and dB < mpf('1e-20')
print('      both << 1.6e7 root separation => s4 = the claimed roots: %s'
      % report('(iv) individual root locks', ok_iv))

print()
print('Chain summary: (i) exact integral q-series; (ii) quoted Shimura CM')
print('integrality + reciprocity (algebraic integers, conjugate in Q(sqrt7));')
print('(iii) interval locks force S, P in Z and equal 16584912, 69257922561;')
print('(iv) each value is pinned to its root (separation 1.6e7).  QED.')

print()
if FAILS:
    print('FAILED CHECKS:', FAILS)
    sys.exit(1)
print('ALL CHECKS PASSED')
