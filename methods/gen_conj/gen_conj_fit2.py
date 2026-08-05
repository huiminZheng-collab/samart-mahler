# gen_conj_fit2.py D -- v2: adds twisted newforms g x chi4, g x chi8 and
# extended pslq bases.  Motivation: D == 3 mod 8 (2 inert) likely needs
# twists (cf. Q(i) case needing L16 and L16tw; Samart D=3 using M12).
# D=27: basis uses d3 (dD = 27 d3 trivially, which masked v1 fits).
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
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: t = -t
        a %= n
    return t if n == 1 else 0

if D == 27:
    chi = lambda n: 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)
    per = 3
else:
    chi = lambda n: jacobi(n, D)
    per = D
chi4 = lambda n: 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
chi8p = lambda n: 0 if n % 2 == 0 else (1 if n % 8 in (1, 7) else -1)
chi8m = lambda n: 0 if n % 2 == 0 else (1 if n % 8 in (1, 3) else -1)

# --- theta coefficients --------------------------------------------------------
NMAX = 2500
if D == 27:
    BX = int(NMAX**0.5)+4; BY = int(2*(NMAX/27.0)**0.5)+4
    A = [0]*(NMAX+1)
    for b in range(-BY, BY+1):
        for a in range(-BX-3*abs(b), BX+3*abs(b)+1):
            n = a*a + 3*a*b + 9*b*b
            if 1 <= n <= NMAX:
                A[n] += (2*a+3*b)**2 - 27*b*b
else:
    c = (D+1)//4
    BX = int(NMAX**0.5)+4; BY = int(2*(NMAX/D)**0.5)+4
    A = [0]*(NMAX+1)
    for y in range(-BY, BY+1):
        for x in range(-BX-abs(y), BX+abs(y)+1):
            n = x*x + x*y + c*y*y
            if 1 <= n <= NMAX:
                A[n] += (2*x+y)**2 - D*y*y
bad = [n for n in range(1, NMAX+1) if A[n] % 8 != 0]
aa = [v//8 for v in A]
if D == 27:
    aa = [0 if n % 3 == 0 else v for n, v in enumerate(aa)]
print("integrality bad =", len(bad), " a(1..12):", aa[1:13])

PR = []
for n in range(2, NMAX+1):
    if all(n % p for p in PR if p*p <= n): PR.append(n)

def hecke3(a):
    o1 = o2 = o3 = True
    for p in PR:
        if p*p > NMAX: break
        if a[p*p] != a[p]*a[p] - chi(p)*p*p: o1 = False
    for i, p in enumerate(PR):
        for qq in PR[i+1:]:
            if p*qq > NMAX: break
            if a[p*qq] != a[p]*a[qq]: o2 = False
    for p in PR:
        if p > 150: break
        if chi(p) == -1 and a[p] != 0: o3 = False
    return o1, o2, o3

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

def LamM(a, N):
    """M = Lambda(3) = w Lambda(0) = w L'(g,0); computed as xN^3 I3 + w I0."""
    xN = sqrt(mpf(N))
    w1 = fricke(a, N, mpf("0.7")); w2 = fricke(a, N, mpf("1.3"))
    w = 1 if abs(w1-1) < mpf("1e-40") else (-1 if abs(w1+1) < mpf("1e-40") else 0)
    I0, I3 = mellin_I(a, xN, 0), mellin_I(a, xN, 3)
    return xN**3*I3 + w*I0, w, (nstr(w1, 4), nstr(w2, 4))

Mv, w0, fw0 = LamM(aa, D)
print("g: Hecke=%s Fricke w=%d %s  M = %s" % (hecke3(aa), w0, fw0, nstr(Mv, 16)))

at4 = [v*chi4(n) for n, v in enumerate(aa)]
at8 = [v*chi8p(n) for n, v in enumerate(aa)]
Mt4, w4, fw4 = LamM(at4, 16*D)
Mt8, w8, fw8 = LamM(at8, 64*D)
print("gxchi4: Hecke=%s Fricke w=%d %s  Mt4 = %s" % (hecke3(at4), w4, fw4, nstr(Mt4, 16)))
print("gxchi8: Hecke=%s Fricke w=%d %s  Mt8 = %s" % (hecke3(at8), w8, fw8, nstr(Mt8, 16)))

# --- Dirichlet values ----------------------------------------------------------
Lchi = dirichlet(mpf(2), [chi(n) for n in range(per)])
dD = mpf(D)**mpf("1.5")/(4*pi)*Lchi
d3 = mpf(3)**mpf("1.5")/(4*pi)*dirichlet(mpf(2), [0, 1, -1])
d4 = 2*dirichlet(mpf(2), [0, 1, 0, -1])/pi
L8p = dirichlet(mpf(2), [chi8p(n) for n in range(8)])
L8m = dirichlet(mpf(2), [chi8m(n) for n in range(8)])
d8 = 64*L8p*L8m/pi**3        # d8 = 64 e4 (verify_P1_n5_e1 convention)
d8alt = 4*sqrt(mpf(2))/pi*L8m
print("dD = %s  d3 = %s  d4 = %s  d8 = %s (alt %s)"
      % (nstr(dD, 12), nstr(d3, 12), nstr(d4, 12), nstr(d8, 12), nstr(d8alt, 12)))

# --- EK4 -----------------------------------------------------------------------
y0 = sD/2 if D != 27 else 3*s3/2
tau = mpc(mpf(1)/2, y0)
EK4 = (10*y0/pi**3)*(-lattice_T(1, mpf(1)/2, y0) + 4*lattice_T(2, mpf(1)/2, y0))
print("EK4 = %s" % nstr(EK4, 25))

# --- fits ----------------------------------------------------------------------
dmain = d3 if D == 27 else dD
dname = "d3" if D == 27 else "dD"
def attempt(vec, names, tag):
    rel = pslq(vec, maxcoeff=10**7, maxsteps=10000)
    if rel is None:
        print("%s: NO FIT" % tag); return
    if rel[0] == 0:
        print("%s: trivial rel among basis %s (rejected)" % (tag, rel)); return
    res = fabs(sum(r*v for r, v in zip(rel, vec)))
    print("%s: rel = %s  residual = %s" % (tag, rel, nstr(res, 3)))
    if res < mpf(10)**(-50):
        terms = ", ".join("%s*%s" % (r, n) for r, n in zip(rel, names) if r)
        print("   CANDIDATE (res<1e-50): 0 = %s" % terms)

attempt([EK4, Mv, dmain], ["EK4", "M", dname], "B1 M,d")
attempt([EK4, Mv, Mt4, dmain, d4], ["EK4", "M", "Mt4", dname, "d4"], "B2 +tw4,d4")
attempt([EK4, Mv, Mt8, dmain, d8], ["EK4", "M", "Mt8", dname, "d8"], "B3 +tw8,d8")
attempt([EK4, Mv, Mt4, Mt8, dmain, d4, d8],
        ["EK4", "M", "Mt4", "Mt8", dname, "d4", "d8"], "B4 full")
attempt([EK4, Mt4, dmain, d4], ["EK4", "Mt4", dname, "d4"], "B5 tw4 only")
attempt([EK4, Mt8, dmain, d8], ["EK4", "Mt8", dname, "d8"], "B6 tw8 only")

if D == 7:
    print("CONTROL D=7: |EK4-(10/7)(40M+d7)| = %s"
          % nstr(fabs(EK4 - mpf(10)/7*(40*Mv + dmain)), 3))
