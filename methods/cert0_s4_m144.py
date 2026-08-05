# Cert-0 (n4, m144): rigorous proof that
#     s4(tau1) = -144,   tau1 = (1 + sqrt(-3))/2,
# where s4(tau) = (eta(2t)/eta(t))^24 (16 W^4 + W^-4)^4,
# W = eta(t) eta(4t)^2 / eta(2t)^3  (Samart's signature-4 modular function).
#
# Chain:
#  (i)   s4 has integral q-expansion s4 = q^-1 + ... in Z((q))
#        [formal series check to q^30 below, exact integer arithmetic];
#        in particular s4 is a function of q = e^{2 pi i tau}, hence
#        invariant under tau -> tau + 1.
#  (ii)  [QUOTED CM THEORY] s4 is an eta-quotient, hence a modular unit
#        with integral q-expansion; by Shimura's CM theorem its value at
#        the CM point tau1 (discriminant -3) is an algebraic integer lying
#        in the ring class field of the order Z[tau1].  Since h(-3) = 1 the
#        ring class field is K = Q(sqrt(-3)) itself:  s4(tau1) in K.
#        Reality: s4 has real q-coefficients, so
#            conj(s4(tau1)) = s4(-conj(tau1)) = s4(tau1 - 1) = s4(tau1),
#        using -conj(tau1) = -1/2 + i sqrt3/2 = tau1 - 1 and (i).
#        Hence s4(tau1) in K cap R = Q, and an algebraic integer in Q is a
#        rational integer:  s4(tau1) in Z.
#  (iii) q(tau1) = e^{2 pi i tau1} = -e^{-pi sqrt3} is REAL and negative,
#        so every eta product in the q-series form
#            s4 = q^-1 prod(1+q^n)^24 (16 q R + R^-1)^4,
#            R = prod (1-q^n)^4 (1-q^{4n})^8 / (1-q^{2n})^12
#        is a real product; mpmath.iv real interval arithmetic with
#        rigorous product tails locks s4(tau1) into (-144.5, -143.5).
#        Combined with (ii): s4(tau1) = -144 exactly.  QED (modulo the
#        quoted CM step (ii), the only non-self-contained input).
#
# As an independent numerical cross-check (not part of the proof):
# c(tau1) = s4(tau1)^{1/4} = 2 sqrt3 e^{-i pi/4}, i.e.
# c(tau1) = 2.449489743... - 2.449489743... i  (verified to 45 digits
# against an eta-quotient evaluation in high precision).

from mpmath import mp, mpf, iv

mp.dps = 60
iv.dps = 60

FAILS = []
def report(name, cond):
    if not cond:
        FAILS.append(name)
    return 'PASS' if cond else 'FAIL'

# ---------------- (i) formal q-series of s4 to q^30 (exact integers) ------
N = 30
def Pser(d, N):
    """prod (1 - q^{d n}) as series (Euler pentagonal)."""
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

# R(q) = prod (1-q^n)^4 (1-q^{4n})^8 / (1-q^{2n})^12  (series of W^4 / q^{1/2})
Rq = s_mul(s_mul(s_pow(Pser(1, N), 4, N), s_pow(Pser(4, N), 8, N), N),
           s_inv(s_pow(Pser(2, N), 12, N), N), N)
Rqi = s_inv(Rq, N)
# 16 q R + R^-1
V1 = [16*Rq[i-1] if i >= 1 else 0 for i in range(N+1)]
for i in range(N+1):
    V1[i] += Rqi[i]
V = s_pow(V1, 4, N)
# prod (1+q^n)^24 = prod (1-q^{2n})^24 / (1-q^n)^24
eta_ratio = s_mul(s_pow(Pser(2, N), 24, N), s_inv(s_pow(Pser(1, N), 24, N), N), N)
s4ser = s_mul(eta_ratio, V, N)   # s4 = q^-1 * s4ser
ok_i = all(isinstance(c, int) for c in s4ser) and s4ser[0] == 1
print('(i) s4 = q^-1 * (1 + %d q + %d q^2 + ...), integral coeffs to q^%d: %s'
      % (s4ser[1], s4ser[2], N, report('(i) integral q-series', ok_i)))

# ---------------- (iii) real interval evaluation at tau1 -------------------
# q = -e^{-pi sqrt3}: enclose it as a real interval.
q = -iv.exp(-iv.pi*iv.sqrt(3))
r_iv = abs(q)
assert r_iv.b < mpf('0.01'), '|q| too large for tail bounds'
NT = 40

def tailF(E_iv):
    """enclose exp(u) for |u| <= E (real)."""
    E = max(E_iv.b, mpf(10)**-55)
    return iv.exp(iv.mpf([-E, E]))

def prod_pow(q, d, e, NT):
    """prod_{n=1}^{NT} (1 - q^{d n})^e times rigorous tail factor.
    Tail: |sum_{n>NT} e log(1-q^{dn})| <= |e| r^{d(NT+1)}/((1-r^d)(1-r))."""
    P = iv.mpf(1)
    qn = q**d
    qd = q**d
    for n in range(1, NT+1):
        P *= (1-qn)
        qn *= qd
    r = abs(q)
    E = abs(iv.mpf(e)) * r**(d*(NT+1)) / ((1 - r**d)*(1 - r))
    return (P * tailF(E))**e if e >= 0 else (P * tailF(E))**e

# P24 = prod (1+q^n)^24 = prod (1-q^{2n})^24 / (1-q^n)^24
num = prod_pow(q, 2, 24, NT)
den = prod_pow(q, 1, 24, NT)
P24 = num/den
# R = prod (1-q^n)^4 (1-q^{4n})^8 / (1-q^{2n})^12
Rnum = prod_pow(q, 1, 4, NT) * prod_pow(q, 4, 8, NT)
Rden = prod_pow(q, 2, 12, NT)
R = Rnum/Rden
s4_iv = P24 * (16*q*R + R**(-1))**4 / q

def cv(z):
    return mp.convert(z)

mid = (cv(s4_iv.a) + cv(s4_iv.b))/2
w = cv(s4_iv.b) - cv(s4_iv.a)
print('(iii) s4(tau1) ~', mp.nstr(mid, 30), ' width', mp.nstr(w, 3))
d144 = max(abs(cv(s4_iv.a) + 144), abs(cv(s4_iv.b) + 144))
print('      |s4(tau1) + 144| <= %s' % mp.nstr(d144, 3))
ok_lock = d144 < mpf('0.5')
print('      enclosure inside (-144.5, -143.5): %s'
      % report('(iii) integer lock', ok_lock))

# ---------------- cross-check (not part of the proof) ----------------------
def eta(tau, nterms=300):
    qq = mp.exp(2*mp.pi*1j*tau)
    p = mp.mpc(1)
    qn = qq
    for n in range(1, nterms+1):
        p *= (1-qn)
        qn *= qq
        if abs(qn) < mpf(10)**(-65):
            break
    return mp.exp(mp.pi*1j*tau/12)*p

tau1 = mp.mpc(mpf(1)/2, mp.sqrt(3)/2)
A = eta(tau1)*eta(4*tau1)**2/eta(2*tau1)**3
c_hp = (eta(2*tau1)/eta(tau1))**6*(16*A**4 + A**-4)
s4_hp = c_hp**4
print('(x)  high-precision cross-check: s4(tau1) =', mp.nstr(s4_hp, 25))
print('     c(tau1) =', mp.nstr(c_hp, 25))
print('     |c_hp - 2 sqrt3 e^{-i pi/4}| =',
      mp.nstr(abs(c_hp - 2*mp.sqrt(3)*mp.exp(-mp.pi*1j/4)), 3))
ok_x = abs(s4_hp + 144) < mpf(10)**(-40) and cv(s4_iv.a) <= s4_hp.real <= cv(s4_iv.b)
print('     hp value inside iv enclosure: %s'
      % report('(x) hp containment', ok_x))

print()
print('Chain summary: (i) exact integral q-series (=> T-invariance);')
print('(ii) quoted Shimura CM: s4(tau1) in K cap R = Q, algebraic integer')
print('     => rational integer; (iii) real-iv lock |s4(tau1)+144| < 1/2.')
print('Hence s4((1+sqrt(-3))/2) = -144 exactly.  QED modulo (ii).')
print()
if FAILS:
    print('FAILED CHECKS:', FAILS)
    import sys; sys.exit(1)
print('ALL CHECKS PASSED')
