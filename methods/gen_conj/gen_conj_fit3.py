# gen_conj_fit3.py D -- v3: adds the h=3 conductor-2-order constants.
# For D == 3 mod 8 (2 inert), EK4 = (10 y0/pi^3)(-T(O_K) + 4 T(O_2)) and
# the G-side of T(O_2) involves L(Psi_1,3), L(Psi_2,3) (conjugate pair of
# Hecke characters of conductor (2), class group Z/3 of the order).
# Real combinations:
#   U = (L(Psi_1)+L(Psi_2))       = (2 S1 - Sw - Sw2)/2
#   V = (L(Psi_1)-L(Psi_2))/(i s3) = (Sw - Sw2)/2
# with S_r = sum_{alpha in r + 2 O_K} alpha^2 / N(alpha)^3, computed by the
# shifted G-row machine (absolutely convergent).  Self-check:
#   S1 + Sw + Sw2 = (15/8) L(g,3)   (2 inert, psi((2)) = 4).
# M-units: MU = (4D)^{3/2} U / (4 pi^3), MV likewise (Fricke factor of the
# level-4D pair, sign absorbed by the fit).
# D = 27: same with O_3 = Z[3 w] in Q(sqrt-3), twist level 108.
# New file; modifies nothing.
import sys
from mpmath import (mp, mpf, mpc, pi, zeta, sinh, cosh, cos, exp, sqrt,
                    diff as mpdiff, dirichlet, gamma, power, gammainc,
                    pslq, nstr, fabs)

D = int(sys.argv[1])
mp.dps = 60
sD = sqrt(mpf(D))
s3 = sqrt(mpf(3))

def G_row(x, y):
    return (pi/y)*sinh(2*pi*y)/(cosh(2*pi*y) - cos(2*pi*x))

def row_powers(x, y):
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y)/(2*y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy)/(2*yy), y)/(4*y)
    return S2, S3

def Gcoset(x0, y0, sx, st):
    """Re sum_{lambda in w+L} lambda^2/|lambda|^6, L = Z + tau Z,
    w = sx + st*tau, sx,st in {0,1/2}, st=0 => sx != 0.
    Row asymptotics: row_mm = -pi/(4 (|mm| y0)^3) + O(e^{-2 pi |mm| y0}),
    so an explicit zeta(3,.) tail is added (error ~ e^{-2 pi*45})."""
    tot = mpf(0)
    if st == 0:
        tot += zeta(4, sx) + zeta(4, 1-sx)
    M0 = 45/y0
    M = int(M0) + 1
    for k in range(-M, M+1):
        mm = k + st
        if mm == 0: continue
        S2, S3 = row_powers(mm*x0 + sx, abs(mm)*y0)
        tot += S2 - 2*(mm*y0)**2*S3
    # tail: rows |mm| > M (mm = k+st)
    kp = M + 1          # positive side: mm = kp+st, kp+1+st, ...
    tailp = zeta(3, mpf(kp) + st)
    tailn = zeta(3, mpf(kp) - st)         # negative side: |mm| = k-st
    tot -= (pi/(4*y0**3))*(tailp + tailn)
    return tot

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

# --- theta coefficients (as v2) -----------------------------------------------
NMAX = 2500
if D == 27:
    BX = int(NMAX**0.5)+4; BY = int(2*(NMAX/27.0)**0.5)+4
    A = [0]*(NMAX+1)
    for b in range(-BY, BY+1):
        for a in range(-BX-3*abs(b), BX+3*abs(b)+1):
            n = a*a + 3*a*b + 9*b*b
            if 1 <= n <= NMAX:
                A[n] += (2*a+3*b)**2 - 27*b*b
    x0L, y0L = mpf(1)/2, 3*s3/2
else:
    c = (D+1)//4
    BX = int(NMAX**0.5)+4; BY = int(2*(NMAX/D)**0.5)+4
    A = [0]*(NMAX+1)
    for y in range(-BY, BY+1):
        for x in range(-BX-abs(y), BX+abs(y)+1):
            n = x*x + x*y + c*y*y
            if 1 <= n <= NMAX:
                A[n] += (2*x+y)**2 - D*y*y
    x0L, y0L = mpf(1)/2, sD/2
bad = [n for n in range(1, NMAX+1) if A[n] % 8 != 0]
aa = [v//8 for v in A]
if D == 27:
    aa = [0 if n % 3 == 0 else v for n, v in enumerate(aa)]
print("integrality bad =", len(bad), " a(1..10):", aa[1:11])

def mellin_I(a, xN, s):
    s = mpf(s); tot = mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0: continue
        tot += a[n]*power(2*pi*n, -s)*gammainc(s, 2*pi*n/xN)
    return tot

xN = sqrt(mpf(D))
I0, I3 = mellin_I(aa, xN, 0), mellin_I(aa, xN, 3)
Mv = xN**3*I3 + I0
L3 = Mv*(2*pi)**3/(xN**3*gamma(3))
print("L(g,3) = %s   M = %s" % (nstr(L3, 16), nstr(Mv, 16)))

# --- coset G-sums ---------------------------------------------------------------
# D fundamental: cosets of 2O_K: 1 -> (1/2,0), w -> (0,1/2), 1+w -> (1/2,1/2)
# D=27 (O_3):    1 -> (1/2,0), 3w -> (1/2,1/2), 1+3w -> (0,1/2)
S1 = Gcoset(x0L, y0L, mpf(1)/2, mpf(0))/16
Sa = Gcoset(x0L, y0L, mpf(0), mpf(1)/2)/16
Sb = Gcoset(x0L, y0L, mpf(1)/2, mpf(1)/2)/16
Todd = S1 + Sa + Sb
print("S1 = %s" % nstr(S1, 16))
print("Sa = %s" % nstr(Sa, 16))
print("Sb = %s" % nstr(Sb, 16))
print("T_odd = %s   (15/8)L = %s   ratio = %s"
      % (nstr(Todd, 16), nstr(15*L3/8, 16), nstr(Todd/(15*L3/8), 8)))
U = (2*S1 - Sa - Sb)/2
V = (Sa - Sb)/2
lvl = 4*D if D != 27 else 108
fac = mpf(lvl)**mpf("1.5")/(4*pi**3)
MU = fac*U; MV = fac*V
print("U = %s   V = %s" % (nstr(U, 16), nstr(V, 16)))
print("MU = %s   MV = %s" % (nstr(MU, 16), nstr(MV, 16)))

# --- Dirichlet values -----------------------------------------------------------
Lchi = dirichlet(mpf(2), [chi(n) for n in range(per)])
dD = mpf(D)**mpf("1.5")/(4*pi)*Lchi
d3 = mpf(3)**mpf("1.5")/(4*pi)*dirichlet(mpf(2), [0, 1, -1])
d4 = 2*dirichlet(mpf(2), [0, 1, 0, -1])/pi
d8 = 64*dirichlet(mpf(2), [0,1,0,-1,0,-1,0,1])*dirichlet(mpf(2), [0,1,0,1,0,-1,0,-1])/pi**3
dmain = d3 if D == 27 else dD
dname = "d3" if D == 27 else "dD"

# --- EK4 ------------------------------------------------------------------------
def lattice_T(d, x0, y0):
    B, M = mpf(0), mpf(0)
    for m in range(-300, 301):
        if m == 0:
            B += 2*zeta(4); continue
        S2, S3 = row_powers(d*m*x0, abs(d*m)*y0)
        B += S2; M += m*m*S3
        if abs(d*m)*y0 > 45 and m > 0: break
    return 3*B - 4*d*d*y0*y0*M

EK4 = (10*y0L/pi**3)*(-lattice_T(1, mpf(1)/2, y0L) + 4*lattice_T(2, mpf(1)/2, y0L))
print("EK4 = %s" % nstr(EK4, 25))

# --- fits ------------------------------------------------------------------------
def attempt(vec, names, tag):
    rel = pslq(vec, maxcoeff=10**7, maxsteps=10000)
    if rel is None:
        print("%s: NO FIT" % tag); return
    if rel[0] == 0:
        print("%s: trivial basis rel %s (rejected)" % (tag, rel)); return
    res = fabs(sum(r*v for r, v in zip(rel, vec)))
    print("%s: rel = %s  residual = %s" % (tag, rel, nstr(res, 3)))
    if res < mpf(10)**(-50):
        terms = ", ".join("%s*%s" % (r, n) for r, n in zip(rel, names) if r)
        print("   CANDIDATE (res<1e-50): 0 = %s" % terms)

attempt([EK4, Mv, MU, dmain], ["EK4", "M", "MU", dname], "C1 M,MU,d")
attempt([EK4, Mv, MU, MV, dmain], ["EK4", "M", "MU", "MV", dname], "C2 +MV")
attempt([EK4, Mv, MU, MV, dmain, d4, d8],
        ["EK4", "M", "MU", "MV", dname, "d4", "d8"], "C3 +d4,d8")
attempt([EK4, Mv, MU, MV, dmain, d3, d4, d8],
        ["EK4", "M", "MU", "MV", dname, "d3", "d4", "d8"], "C4 full")
