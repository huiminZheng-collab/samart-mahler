# cert0_n4_p4_t3.py -- Cert-0 (Phase 4, target T3): rigorous proof that
#     s4(sqrt(-6))          = 1207368 + 853632 sqrt2 + 697680 sqrt3
#                             + 493272 sqrt6           =: r(+,+),
#     s4(i sqrt6/4)         = r(-,+),
#     s4((1+i sqrt6)/2)     = r(+,-),
#     s4(1/2 + i sqrt6/4)   = r(-,-),
# where r(e1,e2) = a + e1 b sqrt2 + e2 c sqrt3 + e1 e2 d sqrt6 with
# (a,b,c,d) = (1207368, 853632, 697680, 493272), generating Q(sqrt2, sqrt3).
#
# The four CM/level points (all with REAL q, so real values):
#   A: tau = i sqrt6      (disc -24,  q = e^{-2 pi sqrt6} > 0),
#   B: tau = i sqrt6/4    (disc -96,  q = e^{-pi sqrt6/2} > 0),
#   C: tau = (1+i sqrt6)/2 (disc -96, q = -e^{-pi sqrt6} < 0),
#   D: tau = 1/2 + i sqrt6/4 (disc -96, q = -e^{-pi sqrt6/2} < 0).
# [Point hunt 2026-08-05: i sqrt6/2, 2i sqrt6, 2i sqrt6/3 excluded
#  numerically (not conjugates); s4 is W_2-invariant, s4(-1/(2t)) = s4(t),
#  cf. scan notes in n4_phase4_report.md.]
#
# Chain (template cert0_n4_p3_t1t2.py):
#  (i)   s4 = q^-1 (1 + ...) in Z((q)) [exact formal series to q^30];
#        q real at all four points => real values.
#  (ii)  [QUOTED CM THEORY] s4 is an eta-quotient modular unit with
#        integral q-expansions at all cusps; by Shimura's CM theorem the
#        four values are algebraic integers; the four r(e1,e2) form ONE
#        full Galois orbit in the biquadratic (totally real) field
#        Q(sqrt2, sqrt3) => elementary symmetric functions e1..e4 in Z.
#  (iii) Real interval arithmetic with rigorous product tails locks
#        e1..e4 to the stated integers.
#  (iv)  Individual interval locks pin each value to its root
#        (minimum root separation ~ 305 >> lock widths).
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

# ---------------- exact symmetric functions of the four roots --------------
a, b, c, d = 1207368, 853632, 697680, 493272
E1 = 4*a
E2 = 6*a*a - 4*b*b - 6*c*c - 12*d*d
# E3, E4 computed symbolically via the u_i = e1 b s2 + e2 c s3 + e1 e2 d s6:
# use exact integer arithmetic in Z[s2,s3]: represent element as
# (w, x, y, z) = w + x s2 + y s3 + z s6.
def qmul(u, v):
    return (u[0]*v[0] + 2*u[1]*v[1] + 3*u[2]*v[2] + 6*u[3]*v[3],
            u[0]*v[1] + u[1]*v[0] + 3*u[2]*v[3] + 3*u[3]*v[2],
            u[0]*v[2] + 2*u[1]*v[3] + u[2]*v[0] + 2*u[3]*v[1],
            u[0]*v[3] + u[1]*v[2] + u[2]*v[1] + u[3]*v[0])
U = [(0, b, c, d), (0, -b, c, -d), (0, b, -c, -d), (0, -b, -c, d)]
R = [(a, b, c, d), (a, -b, c, -d), (a, b, -c, -d), (a, -b, -c, d)]
e3u = (0, 0, 0, 0)
for i in range(4):
    for j in range(i+1, 4):
        for k in range(j+1, 4):
            e3u = tuple(e3u[t] + qmul(qmul(U[i], U[j]), U[k])[t]
                        for t in range(4))
e4u = qmul(qmul(U[0], U[1]), qmul(U[2], U[3]))
e2u = -2*(2*b*b + 3*c*c + 6*d*d)
E3 = 4*a**3 + 2*a*e2u + e3u[0]
E4 = qmul(qmul(R[0], R[1]), qmul(R[2], R[3]))[0]
print('(ii-aux) minimal polynomial X^4 - (%d) X^3 + (%d) X^2 - (%d) X + (%d)'
      % (E1, E2, E3, E4))
print('         (exact integer symmetric functions via Z[s2,s3] arithmetic: %s)'
      % report('(ii-aux) symmetric functions integral',
               e3u[1:] == (0, 0, 0) and e4u[1:] == (0, 0, 0)))

# ---------------- (iii)/(iv) real-iv s4 evaluation -------------------------
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

s6 = iv.sqrt(6)
sA = s4_iv_from_q(iv.exp(-2*iv.pi*s6))          # tau = i sqrt6
sB = s4_iv_from_q(iv.exp(-iv.pi*s6/2))          # tau = i sqrt6/4
sC = s4_iv_from_q(-iv.exp(-iv.pi*s6))           # tau = (1+i sqrt6)/2
sD = s4_iv_from_q(-iv.exp(-iv.pi*s6/2))         # tau = 1/2 + i sqrt6/4

vals = [sA, sB, sC, sD]
labels = ['A (i s6) -> r(+,+)', 'B (i s6/4) -> r(-,+)',
          'C ((1+i s6)/2) -> r(+,-)', 'D (1/2+i s6/4) -> r(-,-)']

# (iii) symmetric-function locks
es = [E1, E2, E3, E4]
import itertools
got = [sum(vals), None, None, None]
s2iv = iv.sqrt(2); s3iv = iv.sqrt(3); s6iv = s6
roots = [a + s2iv*b + s3iv*c + s6iv*d, a - s2iv*b + s3iv*c - s6iv*d,
         a + s2iv*b - s3iv*c - s6iv*d, a - s2iv*b - s3iv*c + s6iv*d]
sym = [iv.mpf(0)]*4
for k in range(1, 5):
    tot = iv.mpf(0)
    for idxs in itertools.combinations(range(4), k):
        p = iv.mpf(1)
        for t in idxs:
            p *= vals[t]
        tot += p
    sym[k-1] = tot
    dev = max(abs(cv(tot.a) - es[k-1]), abs(cv(tot.b) - es[k-1]))
    scale = max(mpf(1), abs(mpf(es[k-1])))
    print('(iii) e%d = %d locked: |iv - target| <= %s : %s'
          % (k, es[k-1], mp.nstr(dev, 3),
             report('(iii) e%d integrality lock' % k,
                    dev < scale*mpf('1e-30') or dev < mpf('1e-25'))))

# (iv) individual root locks
seps = []
for i in range(4):
    for j in range(i+1, 4):
        seps.append(cv(abs(roots[i] - roots[j]).a))
minsep = min(seps)
okiv = True
for lb, v, r in zip(labels, vals, roots):
    dev = cv(abs(v - r).b)
    ok = dev < minsep/100
    okiv &= ok
    print('(iv) %-28s |s4 - root| <= %s (min separation %.4e) : %s'
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

s2m, s3m, s6m = mp.sqrt(2), mp.sqrt(3), mp.sqrt(6)
rt = lambda e1, e2: a + e1*b*s2m + e2*c*s3m + e1*e2*d*s6m
for tag, tau, tgt in (('A', mp.mpc(0, s6m), rt(1, 1)),
                      ('B', mp.mpc(0, s6m/4), rt(-1, 1)),
                      ('C', mp.mpc(mpf(1)/2, s6m/2), rt(1, -1)),
                      ('D', mp.mpc(mpf(1)/2, s6m/4), rt(-1, -1))):
    v = s4_hp(tau)
    print('(x) %s: |s4_hp - target| = %s' % (tag, mp.nstr(abs(v - tgt), 3)))
    if abs(v - tgt) > mpf(10)**(-38):
        FAILS.append('(x) ' + tag)

print()
print('Chain summary: (i) exact integral q-series; (ii) quoted Shimura CM')
print('integrality + one full Galois orbit in Q(sqrt2,sqrt3); (iii) real-iv')
print('locks force e1..e4 in Z; (iv) each value pinned to its root')
print('(min separation %.3e).  QED modulo (ii).' % minsep)
print()
if FAILS:
    print('FAILED CHECKS:', FAILS)
    sys.exit(1)
print('ALL CHECKS PASSED')
