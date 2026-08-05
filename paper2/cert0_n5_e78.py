# cert0_n5_e78.py -- Cert-0 (Phase 8, targets #7/#8): rigorous proof that
#     s4((3 + i sqrt21)/6)  = -893952 + 516096 sqrt3   =: r(+),
#     s4((1 + i sqrt21)/2)  = -893952 - 516096 sqrt3   =: r(-),
# the conjugate pair in Q(sqrt3), roots of X^2 + 1787904 X + c with
# c = 893952^2 - 3*516096^2 (computed exactly below).
#
# The two CM points (both with REAL negative q, so real values):
#   T7: tau = (3+i sqrt21)/6  (q = -e^{-pi sqrt21/3} < 0),
#   T8: tau = (1+i sqrt21)/2  (q = -e^{-pi sqrt21}   < 0).
#
# Chain (template cert0_n5_e56.py):
#  (i)   s4 = q^-1 (1 + ...) in Z((q)) [exact formal series to q^30];
#  (ii)  [QUOTED CM THEORY] s4 is an eta-quotient modular unit with
#        integral q-expansions at all cusps; by Shimura's CM theorem the
#        two values are algebraic integers, and they form ONE conjugate
#        pair => e1 = sum, e2 = product are rational integers.
#        This is the ONLY non-self-contained input of the proof.
#  (iii) Real interval arithmetic with rigorous product tails locks
#        e1 = -1787904 and e2 = c.
#  (iv)  Individual interval locks pin each value to its root
#        (root separation ~ 1.79e6 >> lock widths).
#  QED modulo the quoted CM step (ii), the only non-self-contained input.

from mpmath import mp, mpf, iv, sqrt
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
    res = [0]*(N+1)
    res[0] = 1
    k = 1
    while True:
        e1, e2 = k*(3*k-1)//2*d, k*(3*k+1)//2*d
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
            if B[j]: C[i+j] += A[i]*B[j]
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

Rq = s_mul(s_mul(s_pow(Pser(1, N), 4, N), s_pow(Pser(4, N), 8, N), N),
           s_inv(s_pow(Pser(2, N), 12, N), N), N)
Rqi = s_inv(Rq, N)
V1 = [16*Rq[i-1] if i >= 1 else 0 for i in range(N+1)]
for i in range(N+1):
    V1[i] += Rqi[i]
V = s_pow(V1, 4, N)
eta_ratio = s_mul(s_pow(Pser(2, N), 24, N), s_inv(s_pow(Pser(1, N), 24, N), N), N)
s4ser = s_mul(eta_ratio, V, N)   # s4 = q^-1 * s4ser
ok_i = all(isinstance(c, int) for c in s4ser) and s4ser[0] == 1
print('(i) s4 = q^-1 * (1 + %d q + %d q^2 + ...), integral coeffs to q^%d: %s'
      % (s4ser[1], s4ser[2], N, report('(i) integral q-series', ok_i)))

# ---------------- (ii-aux) exact symmetric functions of the roots ----------
# roots r(+-) = A +- B sqrt3 with A = -893952, B = 516096.
A, B = -893952, 516096
E1 = 2*A                        # sum = 2A = -1787904
E2 = A*A - 3*B*B                # product = c (exact integer arithmetic)
ok_ii = isinstance(E1, int) and isinstance(E2, int)
print('(ii-aux) minimal polynomial X^2 - (%d) X + (%d), discriminant '
      '%d = 3*%d^2 : %s' % (E1, E2, E1*E1 - 4*E2, 2*B,
      report('(ii-aux) symmetric functions integral',
             ok_ii and E1*E1 - 4*E2 == 3*(2*B)**2)))

# ---------------- (iii)/(iv) real-iv s4 evaluation ------------------------
NT = 60

def tailF(E_iv):
    E = max(E_iv.b, mpf(10)**-55)
    return iv.exp(iv.mpf([-E, E]))

def prod_pow(q, d_, e, NT):
    P = iv.mpf(1)
    qn = q**d_
    qd = q**d_
    for n in range(1, NT+1):
        P *= (1-qn)
        qn *= qd
    r = abs(q)
    E = abs(iv.mpf(e)) * r**(d_*(NT+1)) / ((1 - r**d_)*(1 - r))
    return (P * tailF(E))**e

def s4_iv_from_q(q):
    P24 = prod_pow(q, 2, 24, NT) / prod_pow(q, 1, 24, NT)
    R = (prod_pow(q, 1, 4, NT) * prod_pow(q, 4, 8, NT)) / prod_pow(q, 2, 12, NT)
    return P24 * (16*q*R + R**(-1))**4 / q

s21 = iv.sqrt(21)
v7 = s4_iv_from_q(-iv.exp(-iv.pi*s21/3))   # tau = (3 + i sqrt21)/6
v8 = s4_iv_from_q(-iv.exp(-iv.pi*s21))     # tau = (1 + i sqrt21)/2

s3iv = iv.sqrt(3)
rp = iv.mpf(A) + iv.mpf(B)*s3iv             # r(+) ~ 1.88
rm = iv.mpf(A) - iv.mpf(B)*s3iv             # r(-) ~ -1787905.9

# (iii) symmetric locks
sm = v7 + v8
pr = v7*v8
dev1 = max(abs(cv(sm.a) - E1), abs(cv(sm.b) - E1))
dev2 = max(abs(cv(pr.a) - E2), abs(cv(pr.b) - E2))
print('(iii) e1 = %d locked: |iv sum - target| <= %s : %s'
      % (E1, mp.nstr(dev1, 3),
         report('(iii) e1 integrality lock', dev1 < mpf('1e-20'))))
print('(iii) e2 = %d locked: |iv product - target| <= %s : %s'
      % (E2, mp.nstr(dev2, 3),
         report('(iii) e2 integrality lock', dev2 < mpf('1e-20'))))

# (iv) individual root locks
minsep = cv(abs(rp - rm).a)
okiv = True
for lb, v, r in (('T7 (3+i s21)/6 -> r(+)', v7, rp),
                 ('T8 (1+i s21)/2 -> r(-)', v8, rm)):
    dev = cv(abs(v - r).b)
    ok = dev < minsep/100
    okiv &= ok
    print('(iv) %-24s |s4 - root| <= %s (min separation %.4e) : %s'
          % (lb, mp.nstr(dev, 3), minsep, 'PASS' if ok else 'FAIL'))
report('(iv) individual root locks', okiv)
if not okiv:
    FAILS.append('(iv) individual root locks')

# ---------------- high-precision cross-check (not part of the proof) ------
def eta(tau, nterms=400):
    qq = mp.exp(2*mp.pi*1j*tau)
    p = mp.mpc(1)
    qn = qq
    for n in range(1, nterms+1):
        p *= (1-qn)
        qn *= qq
        if abs(qn) < mpf(10)**(-65):
            break
    return mp.exp(mp.pi*1j*tau/12)*p

def s4_hp(tau):
    e1, e2, e4 = eta(tau), eta(2*tau), eta(4*tau)
    W = e1*e4**2/e2**3
    return (e2/e1)**24 * (16*W**4 + W**(-4))**4

s3m = mp.sqrt(3); s21m = mp.sqrt(21)
d7 = s4_hp(mp.mpc(mp.mpf(1)/2, s21m/6)) - (mp.mpf(A) + mp.mpf(B)*s3m)
d8 = s4_hp(mp.mpc(mp.mpf(1)/2, s21m/2)) - (mp.mpf(A) - mp.mpf(B)*s3m)
print('(hp cross-check, 60 dps) |dev T7| = %s, |dev T8| = %s'
      % (mp.nstr(abs(d7), 3), mp.nstr(abs(d8), 3)))

print()
if FAILS:
    print('FAILED:', FAILS); sys.exit(1)
print('ALL CERT-0 CHECKS PASSED: s4((3+i s21)/6) = r(+), '
      's4((1+i s21)/2) = r(-) in Q(sqrt3)')
