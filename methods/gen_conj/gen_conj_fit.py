# gen_conj_fit.py D -- per-disc fit machine for new n4 candidates.
# Usage: python gen_conj_fit.py D    (D in 7[control],11,19,27,43,67,163)
# Builds the CM newform of S_3(Gamma0(D), chi_{-D}) for the h=1 (order)
# disc -D via the Hecke-character theta sum over the principal form,
# checks Hecke/multiplicativity/CM-vanishing/Fricke=+1, computes
# M = L'(g,0) (via Lam3, w=+1 symmetry), d_D = D^{3/2}/(4 pi) L(chi_{-D},2),
# EK4(tau) at 60 dps, and pslq-fits EK4 against [M, dD] (+ spare d3, d4).
# Control run D=7 must reproduce Samart: EK4 = (10/7)(40 M7 + d7).
# New file; modifies nothing.
import sys
from mpmath import (mp, mpf, mpc, pi, zeta, sinh, cosh, cos, exp, sqrt,
                    diff as mpdiff, dirichlet, gamma, power, gammainc,
                    pslq, nstr, fabs)

D = int(sys.argv[1])
mp.dps = 60
sD = sqrt(mpf(D))
s3 = sqrt(mpf(3))

def eta_hp(tau, nterms=400):
    q = exp(2*pi*1j*tau)
    p = mpc(1); qn = q
    for n in range(1, nterms+1):
        p *= (1-qn); qn *= q
        if abs(qn) < mpf(10)**(-65): break
    return exp(pi*1j*tau/12)*p

def s4_hp(tau):
    e1, e2, e4 = eta_hp(tau), eta_hp(2*tau), eta_hp(4*tau)
    W = e1*e4**2/e2**3
    return (e2/e1)**24 * (16*W**4 + W**(-4))**4

def G_row(x, y):
    return (pi/y)*sinh(2*pi*y)/(cosh(2*pi*y) - cos(2*pi*x))

def row_powers(x, y):
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y)/(2*y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy)/(2*yy), y)/(4*y)
    return S2, S3

def lattice_T(d, x0, y0):
    B, M = mpf(0), mpf(0)
    for m in range(-300, 301):
        if m == 0:
            B += 2*zeta(4); continue
        S2, S3 = row_powers(d*m*x0, abs(d*m)*y0)
        B += S2; M += m*m*S3
        if abs(d*m)*y0 > 45 and m > 0: break
    return 3*B - 4*d*d*y0*y0*M

def jacobi(a, n):
    a %= n; t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): t = -t
            t = t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: t = -t
        a %= n
    return t if n == 1 else 0

if D == 27:
    chi = lambda n: 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)   # chi_{-3}
    per = 3
else:
    chi = lambda n: jacobi(n, D)   # (-D/n) = (n/D), D prime == 3 mod 4
    per = D

# --- theta coefficients of the Hecke character --------------------------------
NMAX = 2500
if D == 27:
    # order Z[3 omega] in Q(sqrt-3): norm a^2+3ab+9b^2, Re(alpha^2) =
    # ((2a+3b)^2 - 27 b^2)/4 ; a(n) = (1/8) sum over reps (+- pairing)
    cform = None
    BX = int(NMAX**0.5)+4; BY = int(2*(NMAX/27.0)**0.5)+4
    A = [0]*(NMAX+1)
    for b in range(-BY, BY+1):
        for a in range(-BX-3*abs(b), BX+3*abs(b)+1):
            n = a*a + 3*a*b + 9*b*b
            if 1 <= n <= NMAX:
                v = (2*a+3*b)**2 - 27*b*b
                assert v % 8 == 0 or True
                A[n] += v   # collected as /8 below
else:
    c = (D+1)//4
    BX = int(NMAX**0.5)+4; BY = int(2*(NMAX/D)**0.5)+4
    A = [0]*(NMAX+1)
    for y in range(-BY, BY+1):
        for x in range(-BX-abs(y), BX+abs(y)+1):
            n = x*x + x*y + c*y*y
            if 1 <= n <= NMAX:
                A[n] += (2*x+y)**2 - D*y*y
# a(n) = A[n]/8, must be integral
bad = [n for n in range(1, NMAX+1) if A[n] % 8 != 0]
print("integrality a(n)=A/8: bad count =", len(bad))
aa = [v//8 for v in A]
if D == 27:
    aa = [0 if n % 3 == 0 else v for n, v in enumerate(aa)]
print("a(1..12):", aa[1:13])

PR = []
for n in range(2, NMAX+1):
    if all(n % p for p in PR if p*p <= n): PR.append(n)

ok1 = ok2 = ok3 = True
for p in PR:
    if p*p > NMAX: break
    if aa[p*p] != aa[p]*aa[p] - chi(p)*p*p: ok1 = False
for i, p in enumerate(PR):
    for qq in PR[i+1:]:
        if p*qq > NMAX: break
        if aa[p*qq] != aa[p]*aa[qq]: ok2 = False
for p in PR:
    if p > 150: break
    if chi(p) == -1 and aa[p] != 0: ok3 = False
print("Hecke(p^2)=%s multiplicativity=%s CM-vanishing=%s" % (ok1, ok2, ok3))

def mellin_I(a, xN, s):
    s = mpf(s); tot = mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0: continue
        tot += a[n]*power(2*pi*n, -s)*gammainc(s, 2*pi*n/xN)
    return tot

def qser(a, tau):
    q = exp(2*pi*1j*tau)
    s = mpc(0); qn = q
    for n in range(1, len(a)):
        if a[n]: s += a[n]*qn
        qn *= q
        if abs(qn) < mpf(10)**(-65): break
    return s

def fricke(a, N, yv):
    xN = sqrt(mpf(N))
    return qser(a, 1j/(xN*yv))/(yv**3*qser(a, 1j*yv/xN))

f1 = fricke(aa, D, mpf("0.7")); f2 = fricke(aa, D, mpf("1.3"))
print("Fricke w: %s, %s" % (nstr(f1, 6), nstr(f2, 6)))

xN = sqrt(mpf(D))
I0, I3 = mellin_I(aa, xN, 0), mellin_I(aa, xN, 3)
Mv = xN**3*I3 + I0          # = Lambda(3) = Lambda(0) = L'(g,0)  (w=+1)
L3 = Mv*(2*pi)**3/(xN**3*gamma(3))
print("L(g,3) = %s" % nstr(L3, 15))
print("M = L'(g,0) = %s" % nstr(Mv, 18))

# --- Dirichlet values ----------------------------------------------------------
Lchi = dirichlet(mpf(2), [chi(n) for n in range(per)])
dD = mpf(D)**mpf("1.5")/(4*pi)*Lchi
L3m = dirichlet(mpf(2), [0, 1, -1])
Cat = dirichlet(mpf(2), [0, 1, 0, -1])
d3 = mpf(3)**mpf("1.5")/(4*pi)*L3m
d4 = 2*Cat/pi
print("d_D = %s   d3 = %s  d4 = %s" % (nstr(dD, 15), nstr(d3, 12), nstr(d4, 12)))

# --- s4 and EK4 ---------------------------------------------------------------
y0 = sD/2 if D != 27 else 3*s3/2
tau = mpc(mpf(1)/2, y0)
s4v = s4_hp(tau)
EK4 = (10*y0/pi**3)*(-lattice_T(1, mpf(1)/2, y0) + 4*lattice_T(2, mpf(1)/2, y0))
print("s4 = %s   |s4| = %s" % (nstr(s4v, 30), nstr(fabs(s4v), 8)))
print("EK4 = %s" % nstr(EK4, 25))

# --- fits ----------------------------------------------------------------------
def attempt(vec, names, tag):
    rel = pslq(vec, maxcoeff=10**7, maxsteps=10000)
    if rel is None:
        print("%s: NO FIT" % tag); return
    res = fabs(sum(r*v for r, v in zip(rel, vec)))
    print("%s: rel = %s  residual = %s" % (tag, rel, nstr(res, 3)))
    if res < mpf(10)**(-50):
        terms = ", ".join("%s*%s" % (r, n) for r, n in zip(rel, names) if r)
        print("   CANDIDATE (res<1e-50): 0 = %s  [%s]" % (terms, tag))

attempt([EK4, Mv, dD], ["EK4", "M", "dD"], "basis M,dD")
attempt([EK4, Mv, dD, d3, d4], ["EK4", "M", "dD", "d3", "d4"], "basis M,dD,d3,d4")

# control check for D=7: Samart says EK4 = (10/7)(40 M7 + d7)
if D == 7:
    r7 = EK4 - mpf(10)/7*(40*Mv + dD)
    print("CONTROL D=7: |EK4 - (10/7)(40M7+d7)| = %s" % nstr(fabs(r7), 3))
