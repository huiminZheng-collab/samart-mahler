# Cert-0': rigorous proof that
#     s2(taup)  = (47 + 45 sqrt(-7))/2,   taup  = (1 + sqrt(-7))/8
#     s2(taupp) = (47 - 45 sqrt(-7))/2,   taupp = (-1 + sqrt(-7))/8 = -conj(taup)
# (the two arguments of Samart's Table 4 conjecture n2((47 +- 45 sqrt(-7))/2) = (4/7)(54 M7 + d7)).
#
# Chain (reuses the machinery of cert0_s2_eq_1.py):
#  (a)  s2(tau) = 16/(lam(2t)(1-lam(2t)))   [same formal product identity as cert0 (ii);
#       re-verified below as an integer q-series identity to q^37]
#  (b)  2*taup = pi/2 = (1+sqrt(-7))/4, root of the primitive form (2,-1,1) of disc -7
#       (NB: disc -7, the MAXIMAL order, h(-7) = 1 -- the same class group as cert0,
#       so j(pi/2) is again the unique singular modulus of disc -7; scaling the lattice
#       <1,pi/2> by 2 gives <2,pi>, whose endomorphism ring is O_K).
#       => j(pi/2) in Z. Interval evaluation (iv, 100 dps, rigorous tails) locks
#       j0 = -3375.
#       lam(pi/2) is then a root of Phi(X) = 256(X^2-X+1)^3 + 3375 X^2 (1-X)^2.
#       Phi splits COMPLETELY over Q(sqrt(-7)):
#           Phi = (X^2-X+16)(16X^2-X+1)(16X^2-31X+16)      [exact factorization below]
#       with six roots (1+-3s)/2, (1+-3s)/32, (31+-3s)/32, s = sqrt(-7).
#       Interval lock: lam(pi/2) lies within < 0.248 of (31+3s)/32, and every other
#       root is at distance >= 3 sqrt(7)/16 = 0.4960... from it, hence
#           lam(pi/2) = (31 + 3 sqrt(-7))/32   EXACTLY.
#       (Note: this is 1 - lam0 with lam0 = (1-3sqrt(-7))/32 the candidate root of
#       16X^2-X+1; lam and 1-lam give the same lam(1-lam), so the conclusion of the
#       original plan is unaffected. The direct numerical evaluation picks out
#       (31+3sqrt(-7))/32, NOT lam0 itself.)
#  (c)  exact algebra in Q(sqrt(-7)):  lam(1-lam) = (47 - 45 sqrt(-7))/512, hence
#       s2(taup) = 16/(lam(1-lam)) = (47 + 45 sqrt(-7))/2.
#  (d)  s2 has real q-coefficients, so s2(-x+iy) = conj(s2(x+iy)); taupp = -conj(taup)
#       gives s2(taupp) = (47 - 45 sqrt(-7))/2. Numerically verified to 40 digits.

from fractions import Fraction as Fr
from math import comb
import math
from mpmath import mp, mpf, mpc, pi, exp, sqrt, iv

def ok(cond):
    return 'PASS' if cond else 'FAIL'

# ================= exact arithmetic in Q(sqrt(-7)), s^2 = -7 =================
class QS:
    """a + b*sqrt(-7) with a, b Fractions."""
    def __init__(self, a, b=0):
        self.a, self.b = Fr(a), Fr(b)
    def __mul__(self, o):
        if not isinstance(o, QS): o = QS(o)
        return QS(self.a*o.a - 7*self.b*o.b, self.a*o.b + self.b*o.a)
    __rmul__ = __mul__
    def __add__(self, o):
        if not isinstance(o, QS): o = QS(o)
        return QS(self.a + o.a, self.b + o.b)
    __radd__ = __add__
    def __sub__(self, o):
        if not isinstance(o, QS): o = QS(o)
        return QS(self.a - o.a, self.b - o.b)
    def __rsub__(self, o): return QS(o) - self
    def __neg__(self):    return QS(-self.a, -self.b)
    def __truediv__(self, o):
        n = o.a*o.a + 7*o.b*o.b
        return QS((self.a*o.a + 7*self.b*o.b)/n, (self.b*o.a - self.a*o.b)/n)
    def __pow__(self, k):
        r = QS(1)
        for _ in range(k): r = r*self
        return r
    def __eq__(self, o): return self.a == o.a and self.b == o.b
    def __repr__(self):  return '(%s + %s*s)' % (self.a, self.b)

s7 = QS(0, 1)   # sqrt(-7)

def polymul(A, B):
    R = [A[0]*0]*(len(A)+len(B)-1)
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            R[i+j] = R[i+j] + a*b
    return R

# ================= (a) formal identity s2 = 16/(lam(2t)(1-lam(2t))) =================
# verbatim reuse of cert0_s2_eq_1.py section (ii): exact integer series check.
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
    for i in range(N): lam[i+1] = 16*prod[i]
    return lam

def ser_mul(A, B):
    R = [Fr(0)]*(N+2)
    for i in range(N+2):
        if A[i] == 0: continue
        for m in range(N+2-i):
            R[i+m] += A[i]*B[m]
    return R

def ser_inv(A):
    R = [Fr(0)]*(N+2); R[0] = 1/A[0]
    for k in range(1, N+2):
        acc = Fr(0)
        for m in range(1, min(k+1, len(A))):
            acc += A[m]*R[k-m]
        R[k] = -acc/A[0]
    return R

s2s = series_s2()
lam = series_lam2()
one_minus = [Fr(0)]*(N+2); one_minus[0] = Fr(1)
for i in range(N+2): one_minus[i] -= lam[i]
prod = ser_mul(lam, one_minus)
B = [Fr(0)]*(N+1)
for i in range(N+1): B[i] = prod[i+1]
invB = ser_inv(B)
rhs = [Fr(0)]*(N+2)
for e in range(0, N):
    rhs[e] = 16*invB[e+1]
ok_a = (invB[0] == Fr(1, 16)) and all(s2s[e] == rhs[e] for e in range(0, N-2))
print('(a) formal q-series identity s2 = 16/(lam(2t)(1-lam(2t))) to q^%d: %s' % (N-3, ok(ok_a)))

# ================= (b1) j(pi/2) locked to an integer by interval arithmetic =================
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

w = iv.mpc(iv.mpf(1)/4, iv.sqrt(7)/4)          # pi/2 = (1+sqrt(-7))/4
j_iv = E4_iv(w)**3/delta_iv(w)
def cv(z):
    """zero-width iv endpoint -> mp mpf."""
    return mp.convert(z)
jerr = cv(abs(j_iv + 3375).b)
width_re = cv(j_iv.real.b) - cv(j_iv.real.a)
width_im = cv(j_iv.imag.b) - cv(j_iv.imag.a)
digits_locked = -math.log10(float(max(width_re, width_im, mpf('1e-300'))))
ok_b1 = jerr < mpf('1e-60')
print('(b1) j(pi/2) interval: Re = [%s, %s]' % (mp.nstr(cv(j_iv.real.a), 8), mp.nstr(cv(j_iv.real.b), 8)))
print('     Re width = %s, Im width = %s' % (mp.nstr(width_re, 3), mp.nstr(width_im, 3)))
print('     |j + 3375| <= %s  (%d digits locked)' % (mp.nstr(jerr, 3), int(digits_locked)))
print('     disc -7, h = 1 => j in Z, locked j0 = -3375: %s' % ok(ok_b1))
j0 = -3375

# ================= (b2) exact sextic and its split over Q(sqrt(-7)) =================
# Phi(X) = 256 (X^2-X+1)^3 - j0 X^2 (1-X)^2  with j0 = -3375
Q1 = [QS(1), QS(-1), QS(1)]                    # X^2 - X + 1
Q1c = polymul(polymul(Q1, Q1), Q1)
X2omX2 = [QS(0), QS(0), QS(1), QS(-2), QS(1), QS(0), QS(0)]  # X^2 (1-X)^2, padded to deg 6
Phi = [256*c + 3375*d for c, d in zip(Q1c, X2omX2)]
print('(b2) Phi(X) =', [int(c.a) for c in Phi])

f1 = [QS(16), QS(-1), QS(1)]        # X^2 - X + 16
f2 = [QS(1), QS(-1), QS(16)]        # 16X^2 - X + 1
f3 = [QS(16), QS(-31), QS(16)]      # 16X^2 - 31X + 16
prod3 = polymul(polymul(f1, f2), f3)
ok_fac = all(p == q for p, q in zip(Phi, prod3))
print('     Phi = (X^2-X+16)(16X^2-X+1)(16X^2-31X+16) exactly: %s' % ok(ok_fac))

# each quadratic splits over Q(sqrt(-7)) (disc -63 = -7*3^2); verify linear factors
r1a, r1b = QS(Fr(1,2), Fr(3,2)),  QS(Fr(1,2), Fr(-3,2))    # (1+-3s)/2
r2a, r2b = QS(Fr(1,32), Fr(-3,32)), QS(Fr(1,32), Fr(3,32)) # (1-+3s)/32
r3a, r3b = QS(Fr(31,32), Fr(3,32)), QS(Fr(31,32), Fr(-3,32)) # (31+-3s)/32
def linfac_check(f, ra, rb, lead):
    # lead*(X-ra)(X-rb) == f ?
    g = polymul([-ra, QS(1)], [-rb, QS(1)])
    return all(lead*gi == fi for gi, fi in zip(g, f))
ok_l1 = linfac_check(f1, r1a, r1b, QS(1))
ok_l2 = linfac_check(f2, r2a, r2b, QS(16))
ok_l3 = linfac_check(f3, r3a, r3b, QS(16))
print('     linear factors over Q(sqrt(-7)): f1: %s, f2: %s, f3: %s'
      % (ok(ok_l1), ok(ok_l2), ok(ok_l3)))

roots = [r1a, r1b, r2a, r2b, r3a, r3b]
lam_star = r3a            # candidate value of lam(pi/2): (31+3s)/32 = 1 - lam0
def peval_QS(A, x):
    r = QS(0)
    for c in reversed(A): r = r*x + c
    return r
ok_roots = all(peval_QS(Phi, r) == QS(0) for r in roots)
print('     all six roots satisfy Phi exactly (Q(sqrt(-7)) arithmetic): %s' % ok(ok_roots))

# separations from lam_star (exact moduli: |a+bs| = sqrt(a^2+7b^2))
mp.dps = 60
def qsmod2(z): return float(z.a*z.a + 7*z.b*z.b)
seps = []
for r in roots:
    if r == lam_star: continue
    d = r - lam_star
    seps.append((sqrt(mpf(qsmod2(d))), d))
min_sep = min(s for s, _ in seps)
print('     separations from lam_star = (31+3s)/32:')
for sv, d in seps:
    print('       |lam_star - %s| = %s' % (d, mp.nstr(sv, 8)))
print('     min separation = 3 sqrt(7)/16 = %s' % mp.nstr(3*sqrt(mpf(7))/16, 10),
      ' check:', ok(abs(min_sep - 3*sqrt(mpf(7))/16) < mpf('1e-50')))

# ================= (b3) interval lock of lam(pi/2) =================
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

lam_star_mp = (31 + 3*sqrt(mpf(-7)))/32
lv = lam_iv(w)
dlock = cv(abs(lv - lam_star_mp).b)   # upper bound of |lam_iv - lam_star|
tol = 3*sqrt(mpf(7))/32               # half the min separation 3 sqrt(7)/16
ok_b3 = dlock < tol
print('(b3) lam(pi/2) interval: Re = [%s, %s]' % (mp.nstr(cv(lv.real.a), 8), mp.nstr(cv(lv.real.b), 8)))
print('                    Im = [%s, %s]' % (mp.nstr(cv(lv.imag.a), 8), mp.nstr(cv(lv.imag.b), 8)))
print('     |lam_iv - (31+3s)/32| <= %s  (< %s = half min separation): %s'
      % (mp.nstr(dlock, 3), mp.nstr(tol, 6), ok(ok_b3)))
print('     => lam(pi/2) = (31 + 3 sqrt(-7))/32 EXACTLY')
print('     NB: lam(pi/2) = 1 - lam0 with lam0 = (1-3s)/32 the root of 16X^2-X+1;')
print('         the numerical value selects (31+3s)/32, not lam0 (same lam(1-lam)).')

# ================= (c) exact s2(taup) =================
ll = lam_star
prod_exact = ll*(QS(1) - ll)
target_prod = QS(Fr(47,512), Fr(-45,512))
ok_c1 = (prod_exact == target_prod)
s2_exact = QS(16)/prod_exact
target_s2 = QS(Fr(47,2), Fr(45,2))
ok_c2 = (s2_exact == target_s2)
print('(c) lam(1-lam) = (47-45s)/512 exactly: %s' % ok(ok_c1))
print('    16/(lam(1-lam)) = (47+45s)/2 exactly: %s  => s2(taup) = (47+45 sqrt(-7))/2' % ok(ok_c2))

# numeric cross-check: direct q-product for s2(taup), 50 dps
def s2_mp(tau, Ntr=600):
    q = exp(2*pi*1j*tau); p = mpc(1); qn = q
    for n in range(1, Ntr+1):
        p *= (1+qn)**24; qn *= q*q
        if abs(qn) < mpf(10)**(-55): break
    return p/q
taup = mpc(1, sqrt(mpf(7)))/8
v = s2_mp(taup)
tgt = (47 + 45*sqrt(mpf(-7)))/2
d1 = abs(v - tgt)
ok_c3 = d1 < mpf('1e-40')
print('    numeric check |s2(taup) - (47+45s)/2| = %s (50 dps): %s' % (mp.nstr(d1, 3), ok(ok_c3)))

# ================= (d) conjugate point =================
# s2 = q^-1 prod (1+q^{2n-1})^24 has real coefficients; q(-x+iy) = conj(q(x+iy)),
# hence s2(-x+iy) = conj(s2(x+iy)). taupp = -1/8 + i sqrt(7)/8 = -conj(taup).
taupp = mpc(-1, sqrt(mpf(7)))/8
v2 = s2_mp(taupp)
tgt2 = (47 - 45*sqrt(mpf(-7)))/2
d2 = abs(v2 - tgt2)
d3 = abs(v2 - v.conjugate() if hasattr(v, 'conjugate') else v2 - mpc(v.real, -v.imag))
ok_d = (d2 < mpf('1e-40')) and (d3 < mpf('1e-40'))
print('(d) |s2(taupp) - (47-45s)/2| = %s,  |s2(taupp) - conj(s2(taup))| = %s (40+ digits): %s'
      % (mp.nstr(d2, 3), mp.nstr(d3, 3), ok(ok_d)))

print()
allok = all([ok_a, ok_b1, ok_fac, ok_l1, ok_l2, ok_l3, ok_roots, ok_b3, ok_c1, ok_c2, ok_c3, ok_d])
print('ALL CERT-0\' CHECKS PASSED' if allok else '*** SOME CHECK FAILED ***')
