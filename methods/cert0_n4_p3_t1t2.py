# cert0_n4_p3_t1t2.py -- Cert-0 (Phase 3, targets T1/T2): rigorous proof that
#     s4(sqrt(-2))      = 3656 + 2600 sqrt2,
#     s4((1+sqrt(-2))/2)= 3656 - 2600 sqrt2,
#     s4(2i)            = 143208 + 101574 sqrt2,
#     s4((1+2i)/2)      = 143208 - 101574 sqrt2,
# where s4(tau) = (eta(2t)/eta(t))^24 (16 W^4 + W^-4)^4,
# W = eta(t) eta(4t)^2 / eta(2t)^3  (Samart's signature-4 modular function).
#
# Chain (same template as cert0_s4_s7pair.py):
#  (i)   s4 has integral q-expansion s4 = q^-1 * (1 + ...) in Z((q))
#        [formal series check to q^30 below, exact integer arithmetic];
#        for the CM points below q = e^{2 pi i tau} is REAL (> 0 for the
#        pure imaginary point, < 0 for the Re = 1/2 point), so the values
#        are real.
#  (ii)  [QUOTED CM THEORY] s4 is an eta-quotient, hence a modular unit of
#        level 8 with integral q-expansions at all cusps; by Shimura's CM
#        theorem its values at CM points are algebraic integers, and Shimura
#        reciprocity identifies the Galois orbit.  The two points of each
#        pair have lattices with the same endomorphism ring
#        (T1: O_K, K = Q(sqrt(-2)), h(-8) = 1, both lattices homothetic to
#         O_K;  T2: O_2, the conductor-2 order of Q(i), h(-16) = 1), and
#        their s4-values form a conjugate pair generating Q(sqrt2) (real).
#        => sum S and product P are rational integers.
#  (iii) Interval arithmetic (mpmath.iv, REAL intervals only) with rigorous
#        product tails locks S and P to the stated integers; exact integer
#        check A0^2 - 2 B0^2 = P; root separation >> lock error.
#  (iv)  Individual locks pin each value to its root.  QED modulo the
#        quoted CM step (ii), the only non-self-contained input.

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

# W^4 = q^{1/2} R, R = prod (1-q^n)^4 (1-q^{4n})^8 / (1-q^{2n})^12;
# s4 = q^-1 prod((1-q^{2n})/(1-q^n))^24 (16 q R + R^-1)^4   [q = e^{2 pi i t}]
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

# ---------------- (iii)/(iv) real-iv s4 evaluation -------------------------
NT = 60

def tailF(E_iv):
    E = max(E_iv.b, mpf(10)**-55)
    return iv.exp(iv.mpf([-E, E]))

def prod_pow(q, d, e, NT):
    """prod_{n=1}^{NT} (1 - q^{d n})^e with rigorous tail (real q)."""
    P = iv.mpf(1)
    qn = q**d
    qd = q**d
    for n in range(1, NT+1):
        P *= (1-qn)
        qn *= qd
    r = abs(q)
    E = abs(iv.mpf(e)) * r**(d*(NT+1)) / ((1 - r**d)*(1 - r))
    return (P * tailF(E))**e

def s4_iv_from_q(q):
    """s4 = q^-1 P24 (16 q R + R^-1)^4, q real (either sign), |q| << 1."""
    P24 = prod_pow(q, 2, 24, NT) / prod_pow(q, 1, 24, NT)
    R = (prod_pow(q, 1, 4, NT) * prod_pow(q, 4, 8, NT)) / prod_pow(q, 2, 12, NT)
    return P24 * (16*q*R + R**(-1))**4 / q

def certify_pair(tag, yA, yB, A0, B0):
    """point A: tau = i yA (q > 0); point B: tau = 1/2 + i yB (q < 0).
    Target roots r+- = A0 +- B0 sqrt2."""
    qA = iv.exp(-2*iv.pi*yA)
    qB = -iv.exp(-iv.pi*2*yB)         # q = e^{2 pi i (1/2 + i yB)} = -e^{-2 pi yB}
    sA = s4_iv_from_q(qA)
    sB = s4_iv_from_q(qB)
    S_int = 2*A0
    P_int = A0*A0 - 2*B0*B0
    print('(%s) sum S = %d, product P = %d (= A0^2 - 2 B0^2: exact)' %
          (tag, S_int, P_int))
    sum_iv, prod_iv = sA + sB, sA * sB
    dS = max(abs(cv(sum_iv.a) - S_int), abs(cv(sum_iv.b) - S_int))
    dP = max(abs(cv(prod_iv.a) - P_int), abs(cv(prod_iv.b) - P_int))
    tolP = max(mpf('1e-5'), abs(P_int)*mpf('1e-20'))
    print('     |sum - S| <= %s ; |prod - P| <= %s' % (mp.nstr(dS, 3), mp.nstr(dP, 3)))
    print('     => S, P in Z locked: %s'
          % report('(%s) sum/product integrality locks' % tag,
                   dS < mpf('1e-15') and dP < tolP))
    rp = A0 + B0*iv.sqrt(2)
    rm = A0 - B0*iv.sqrt(2)
    sep = 2*B0*sqrt(2)
    dA = cv(abs(sA - rp).b)
    dB = cv(abs(sB - rm).b)
    print('     |s4(A) - r+| <= %s ; |s4(B) - r-| <= %s ; separation ~ %.3e'
          % (mp.nstr(dA, 3), mp.nstr(dB, 3), sep))
    print('     individual root locks: %s'
          % report('(%s) individual root locks' % tag,
                   dA < sep/100 and dB < sep/100))

# T1: tau_A = i sqrt2 (q = e^{-2 pi sqrt2}), tau_B = (1 + i sqrt2)/2
certify_pair('T1', iv.sqrt(2), iv.sqrt(2)/2, 3656, 2600)
# T2: tau_A = 2i (q = e^{-4 pi}), tau_B = (1 + 2i)/2
certify_pair('T2', iv.mpf(2), iv.mpf(1), 143208, 101574)

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

s2 = mp.sqrt(2)
for tag, tau, tgt in (('T1a', mp.mpc(0, s2), 3656 + 2600*s2),
                      ('T1b', mp.mpc(mpf(1)/2, s2/2), 3656 - 2600*s2),
                      ('T2a', mp.mpc(0, 2), 143208 + 101574*s2),
                      ('T2b', mp.mpc(mpf(1)/2, 1), 143208 - 101574*s2)):
    v = s4_hp(tau)
    print('(x) %s: |s4_hp - target| = %s' % (tag, mp.nstr(abs(v - tgt), 3)))
    if abs(v - tgt) > mpf(10)**(-40):
        FAILS.append('(x) ' + tag)

print()
print('Chain summary: (i) exact integral q-series; (ii) quoted Shimura CM')
print('integrality + reciprocity (conjugate algebraic integers in Q(sqrt2));')
print('(iii) real-iv locks force S, P in Z; (iv) each value pinned to its')
print('root (separations 7.4e3 and 2.9e5).  QED modulo (ii).')
print()
if FAILS:
    print('FAILED CHECKS:', FAILS)
    sys.exit(1)
print('ALL CHECKS PASSED')
