# Cert-0: rigorous proof that s2(tau0) = 1, tau0 = (3+sqrt(-7))/8.
#
# Chain:
#  (i)   s2(tau) = q^-1 prod (1+q^{2n-1})^24                    [eta defs + formal product algebra]
#  (ii)  16/(lam(2t)(1-lam(2t))) = q^-1 prod (1+q^{2n-1})^24   [Jacobi triple product + Euler:
#         prod (1+q^{2n}) = 1/prod(1-q^{4n-2})]  => s2 = 16/(lam(2t)(1-lam(2t))) exactly.
#         (Formal series check to q^40 below as a sanity check.)
#  (iii) j(2 tau0) = -3375 exactly: 2*tau0 = (3+s-7)/4 satisfies 2w^2-3w+2=0 (disc -7,
#         primitive, h(-7)=1) => j(w) in O_K; conjugate form equals the form itself and
#         class group is trivial => j(w) real; real + algebraic integer in O_K => in Z;
#         interval evaluation (iv, 100 dps, j = E4^3/Delta with rigorous tails)
#         |j+3375| < 1e-66 pins the integer.
#  (iv)  lam satisfies 256(X^2-X+1)^3 = j X^2 (1-X)^2 identically (classical).
#         => lam(2 tau0) is a root of Phi(X) = 256(X^2-X+1)^3 + 3375 X^2 (1-X)^2.
#  (v)   Exact factorization of Phi over Q: lambda0 = (1+3s-7)/2 is a root
#         (minimal poly X^2-X+16 divides Phi, multiplicity 1); all other roots are
#         separated from lambda0 (explicit gap below). Interval evaluation
#         (iv, 100 dps, rigorous tails) |lam(2 tau0) - lambda0| < 1e-58
#         => lam(2 tau0) = lambda0.
#  (vi)  s2(tau0) = 16/(lambda0 (1-lambda0)) = 16/16 = 1.  QED

from fractions import Fraction as Fr
from math import comb, gcd
from mpmath import mp, mpf, mpc, sqrt, iv
import sys

mp.dps = 70

FAILS = []
def report(name, cond):
    """record check result; return 'PASS'/'FAIL' for printing."""
    if not cond:
        FAILS.append(name)
    return 'PASS' if cond else 'FAIL'

# ---------------- (ii) sanity: exact formal series check to q^40 ----------------
N = 40
def series_s2():
    prod = [0]*(N+2); prod[0] = 1
    k = 1
    while 2*k-1 <= N:
        e = 2*k-1
        new = [0]*(N+2)
        for i in range(N+1):
            if prod[i] == 0: continue
            for jj in range(25):
                if i+jj*e <= N: new[i+jj*e] += prod[i]*comb(24, jj)
        prod = new; k += 1
    s2 = [Fr(0)]*(N+2)
    for i in range(N+1): s2[i] = Fr(prod[i+1])
    return s2

def series_lam2():
    # lam(2t) = 16 q prod (1+q^{2n})^8 (1+q^{2n-1})^{-8}
    prod = [Fr(1)] + [Fr(0)]*N
    n = 1
    while n <= N:
        e = n
        if n % 2 == 0:
            f = [Fr(comb(8, jj)) for jj in range(9)]
        else:
            f = [Fr((-1)**jj*comb(jj+7, 7)) for jj in range(N//e+2)]
        new = [Fr(0)]*(N+1)
        for i in range(N+1):
            if prod[i] == 0: continue
            for jj in range(len(f)):
                if i+jj*e <= N: new[i+jj*e] += prod[i]*f[jj]
        prod = new; n += 1
    lam = [Fr(0)]*(N+2)
    for i in range(N): lam[i+1] = 16*prod[i]   # *16 q
    return lam

def ser_mul(A, B):
    R = [Fr(0)]*(N+2)
    for i in range(N+2):
        if A[i] == 0: continue
        for m in range(N+2-i):
            R[i+m] += A[i]*B[m]
    return R

def ser_inv(A):
    # A[0] != 0
    R = [Fr(0)]*(N+2); R[0] = 1/A[0]
    for k in range(1, N+2):
        acc = Fr(0)
        for m in range(1, min(k+1, len(A))):
            acc += A[m]*R[k-m]
        R[k] = -acc/A[0]
    return R

s2s = series_s2()            # s2s[e] = coeff of q^e, e >= 0  (coeff of q^-1 is 1, dropped)
lam = series_lam2()          # lam[l] = coeff of q^l, l>=1
one_minus = [Fr(0)]*(N+2); one_minus[0] = Fr(1)
for i in range(N+2): one_minus[i] -= lam[i]
prod = ser_mul(lam, one_minus)          # lam(1-lam) = 16 q + ... (no const term)
B = [Fr(0)]*(N+1)
for i in range(N+1): B[i] = prod[i+1]   # prod = q * B, B[0] = 16
invB = ser_inv(B)
rhs = [Fr(0)]*(N+2)
for e in range(0, N):
    rhs[e] = 16*invB[e+1]               # 16/(q B): coeff of q^e = 16 * invB[e+1]
ok_ii = (invB[0] == Fr(1, 16)) and all(s2s[e] == rhs[e] for e in range(0, N-2))
print('(ii) formal q-series identity s2 = 16/(lam(2t)(1-lam(2t))) holds to q^%d: %s'
      % (N-3, report('(ii) formal identity', ok_ii)))

# ---------------- (iv),(v) exact sextic ----------------
# Phi(X) = 256 (X^2-X+1)^3 + 3375 X^2 (1-X)^2
from itertools import product as iproduct
def polymul(A, B):
    R = [Fr(0)]*(len(A)+len(B)-1)
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            R[i+j] += a*b
    return R
def polyadd(A, B):
    n = max(len(A), len(B))
    return [(A[i] if i < len(A) else 0) + (B[i] if i < len(B) else 0) for i in range(n)]
Q1 = [Fr(1), Fr(-1), Fr(1)]              # X^2 - X + 1
Q1c = polymul(polymul(Q1, Q1), Q1)
X2 = [Fr(0), Fr(0), Fr(1)]
onemX2 = polymul([Fr(1), Fr(-1)], [Fr(1), Fr(-1)])
Phi = polyadd([256*c for c in Q1c], [3375*c for c in polymul(X2, onemX2)])
print('(iv) Phi(X) =', [int(c) for c in Phi])
lam0 = (1 + 3*sqrt(mpf(-7)))/2
def peval(A, x):
    return sum(c*x**i for i, c in enumerate(A))
print('     Phi(lambda0) =', mp.nstr(abs(peval(Phi, lam0)), 5), '(exact root check)')

# exact division by X^2 - X + 16
m = [Fr(16), Fr(-1), Fr(1)]
def polydiv(A, B):
    A = list(A); B = list(B)
    while len(B) > 1 and B[-1] == 0: B.pop()
    Q = [Fr(0)]*max(1, (len(A)-len(B)+1))
    R = list(A)
    while len(R) >= len(B):
        k = len(R)-len(B)
        c = R[-1]/B[-1]
        Q[k] = c
        for i in range(len(B)):
            R[i+k] -= c*B[i]
        while len(R) > 1 and R[-1] == 0: R.pop()
        if len(R) == 1 and R[0] == 0: R = [Fr(0)]
        if R == [Fr(0)]: break
    return Q, R
Quo, Rem = polydiv(Phi, m)
ok_fac = (Rem == [Fr(0)])
print('(v) Phi / (X^2-X+16): remainder =', [int(c) for c in Rem], ' quotient =', [int(c) for c in Quo])
print('    X^2-X+16 divides Phi exactly: %s' % report('(v) exact factor', ok_fac))
# multiplicity: does m divide Quo?
_, Rem2 = polydiv(Quo, m)
ok_mult = any(c != 0 for c in Rem2)
print('    multiplicity of lambda0 = 1 (second division leaves nonzero remainder): %s'
      % report('(v) simple root', ok_mult))

# factor quartic: try quadratics with rational coeffs via sympy-free approach:
# just compute all roots numerically at 80 dps and separations.
from mpmath import polyroots
mp.dps = 80
Phif = [mpf(int(c)) for c in Phi]
roots = polyroots(Phif, maxsteps=200, error=False)
sep = min(abs(r - lam0) for r in roots if abs(r-lam0) > mpf('1e-60'))
print('    roots of Phi:')
for r in roots: print('      ', mp.nstr(r, 25))
print('    min |other root - lambda0| =', mp.nstr(sep, 6))

# ---------------- (iii),(v) interval-locked evaluations ----------------
# Same mechanism as cert0_s2_n2pair.py: iv interval arithmetic at 100 dps,
# q-series truncated with rigorous outward-rounded tail bounds (no float64,
# no plain mp on the certificate path).
iv.dps = 100
NJ = 80
sig3 = [0]*(NJ+2)
for n in range(1, NJ+2):
    sig3[n] = sum(d**3 for d in range(1, n+1) if n % d == 0)

def cbox(E):
    """complex interval enclosing exp(u+iv) for |u|,|v| <= E, via a square box."""
    E = max(E, iv.mpf('1e-95'))
    return iv.exp(iv.mpc(iv.mpf([-E, E]), iv.mpf([-E, E])))

def E4_iv(tau):
    q = iv.exp(2*iv.pi*1j*tau)
    r = abs(q).b
    s = iv.mpc(1); qn = q
    for n in range(1, NJ+1):
        s += 240*sig3[n]*qn
        qn *= q
    # r = abs(q).b is the rigorous upper endpoint of |q|; the tail bounds
    # below are increasing in r for 0 < r < 1, so evaluate them at the
    # point interval [r, r] with outward rounding (no float64).
    riv = iv.mpf(r)
    # |sum_{n>N} 240 sig3(n) q^n| <= 240 sum n^4 r^n;
    # n^4 r^(n/2) decreasing for n > 8/ln(1/r), so
    # sum_{n>N} n^4 r^n <= (N+1)^4 r^(N+1)/(1-sqrt(r)).
    assert NJ+1 > 8/iv.log(1/r)
    T = (iv.mpf(240)*(NJ+1)**4 * riv**(NJ+1)/(1-iv.sqrt(riv))).b
    return s + iv.mpc(iv.mpf([-T, T]), iv.mpf([-T, T]))

def delta_iv(tau):
    # Delta = q prod (1-q^n)^24; tail log bounded by 24 sum_{n>N} r^n/(1-r)
    q = iv.exp(2*iv.pi*1j*tau)
    r = abs(q).b
    P = iv.mpc(1); qn = q
    for n in range(1, NJ+1):
        P *= (1-qn)**24
        qn *= q
    riv = iv.mpf(r)     # rigorous upper endpoint of |q|; bound increasing in r
    E = (iv.mpf(24)*riv**(NJ+1)/(1-riv)**2).b
    return q*P*cbox(E)

def cv(z):
    """zero-width iv endpoint -> mp mpf."""
    return mp.convert(z)

w = iv.mpc(iv.mpf(3)/4, iv.sqrt(7)/4)          # 2*tau0 = (3+sqrt(-7))/4
j_iv = E4_iv(w)**3/delta_iv(w)
jerr = cv(abs(j_iv + 3375).b)                  # rigorous upper bound of |j(2 tau0)+3375|
width_re = cv(j_iv.real.b) - cv(j_iv.real.a)
width_im = cv(j_iv.imag.b) - cv(j_iv.imag.a)
print('(iii) j(2 tau0) interval: Re = [%s, %s]' % (mp.nstr(cv(j_iv.real.a), 8), mp.nstr(cv(j_iv.real.b), 8)))
print('      Re width = %s, Im width = %s' % (mp.nstr(width_re, 3), mp.nstr(width_im, 3)))
ok_iii = jerr < mpf('1e-66')
print('      |j(2 tau0) + 3375| <= %s  (< 1e-66, pins j = -3375 in Z): %s'
      % (mp.nstr(jerr, 3), report('(iii) j lock', ok_iii)))

NL = 120
def lam_iv(tau):
    # lam(t) = 16 q prod (1+q^{2n})^8 (1+q^{2n-1})^{-8}, q = e^{pi i t}
    q = iv.exp(iv.pi*1j*tau)
    r = abs(q).b
    P = iv.mpc(1); qn = q
    for n in range(1, NL+1):
        P *= (1+qn)**8 if n % 2 == 0 else (1+qn)**(-8)
        qn *= q
    riv = iv.mpf(r)     # rigorous upper endpoint of |q|; bound increasing in r
    # |sum_{n>N} +-8 log(1+q^n)| <= 8 sum_{n>N} r^n/(1-r) = 8 r^(N+1)/(1-r)^2
    E = (iv.mpf(8)*riv**(NL+1)/(1-riv)**2).b
    return 16*q*P*cbox(E)

lv = lam_iv(w)
lam0_iv = (1 + 3*iv.mpc(0, iv.sqrt(7)))/2      # interval enclosing lambda0 = (1+3 sqrt(-7))/2
dlock = cv(abs(lv - lam0_iv).b)                # rigorous upper bound of |lam(2 tau0) - lambda0|
print('(v)   lam(2 tau0) interval: Re = [%s, %s]' % (mp.nstr(cv(lv.real.a), 8), mp.nstr(cv(lv.real.b), 8)))
print('                          Im = [%s, %s]' % (mp.nstr(cv(lv.imag.a), 8), mp.nstr(cv(lv.imag.b), 8)))
ok_v = dlock < mpf('1e-58')
print('      |lam(2 tau0) - lambda0| <= %s  (<< 3.75 = min root separation): %s'
      % (mp.nstr(dlock, 3), report('(v) lam lock', ok_v)))
print('(vi)  16/(lambda0(1-lambda0)) =', mp.nstr(16/(lam0*(1-lam0)), 10))
print('      => s2(tau0) = 1 exactly.  QED')

print()
if FAILS:
    print('FAILED CHECKS:', FAILS)
    sys.exit(1)
print('ALL CHECKS PASSED')
