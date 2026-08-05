# gen_conj_fit4.py D -- v4: exact B/G anchor decomposition for
# tau = (1+sqrt(-D))/2, D == 3 mod 8, h(-D) = 1.
#   EK4 = (10 y0/pi^3)(-T(O_K) + 4 T(O_2)),  O_2 = 2 O_K ⊔ (1 + 2 O_K)
#   comb = zeta_K(2) + (8/3) BU + 2 L(g,3) + (16/3) U
# with U  = (2 E1 - Ea - Eb)/2   (G-side ray twists, conductor-(2) Hecke
#      characters, class group Z/3 of O_2),
#      BU = (2 B1 - Ba - Bb)/2   (B-side partial ray zetas).
# Predicted identity (M = L'(g,0), MU = (4D)^{3/2} U/(4 pi^3),
# DU = BU * 3 D^{3/2}/(2 pi^3)):
#   EK4 = (40/D) M + (40/(3D)) MU + (10/(3D)) d_D + (80/(9D)) DU.
# Self-checks: E-odd = 15 L/8, B-odd = 15 zeta_K/8, G(O_K) = 4L,
# B(O_K) = 2 zeta_K.  New file; modifies nothing.
import sys
from mpmath import (mp, mpf, mpc, pi, zeta, sinh, cosh, cos, exp, sqrt,
                    diff as mpdiff, dirichlet, gamma, power, gammainc,
                    pslq, nstr, fabs)

D = int(sys.argv[1])
mp.dps = 60
sD = sqrt(mpf(D))

def G_row(x, y):
    return (pi/y)*sinh(2*pi*y)/(cosh(2*pi*y) - cos(2*pi*x))

def row_powers(x, y):
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y)/(2*y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy)/(2*yy), y)/(4*y)
    return S2, S3

def coset_BG(x0, y0, sx, st):
    """(B, E) of the coset w + L, L = Z + tau Z, w = sx + st tau:
    B = sum 1/|lam|^4,  E = sum Re(lam^2)/|lam|^6.
    Row asymptotics S2 ~ pi/(2 y^3), S3 ~ 3 pi/(8 y^5): explicit zeta tails."""
    B = mpf(0); E = mpf(0)
    if st == 0:
        B += zeta(4, sx) + zeta(4, 1-sx)
        E += zeta(4, sx) + zeta(4, 1-sx)
    M0 = 45/y0; M = int(M0) + 1
    for k in range(-M, M+1):
        mm = k + st
        if mm == 0: continue
        S2, S3 = row_powers(mm*x0 + sx, abs(mm)*y0)
        B += S2; E += S2 - 2*(mm*y0)**2*S3
    kp = M + 1
    tl = zeta(3, mpf(kp) + st) + zeta(3, mpf(kp) - st)
    B += (pi/(2*y0**3))*tl
    E -= (pi/(4*y0**3))*tl
    return B, E

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

chi = lambda n: jacobi(n, D)

# --- newform a(n) ---------------------------------------------------------------
NMAX = 2500
c = (D+1)//4
BX = int(NMAX**0.5)+4; BY = int(2*(NMAX/D)**0.5)+4
A = [0]*(NMAX+1)
for y in range(-BY, BY+1):
    for x in range(-BX-abs(y), BX+abs(y)+1):
        n = x*x + x*y + c*y*y
        if 1 <= n <= NMAX:
            A[n] += (2*x+y)**2 - D*y*y
aa = [v//8 for v in A]

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

# --- Dirichlet/zeta values -------------------------------------------------------
Lchi = dirichlet(mpf(2), [chi(n) for n in range(D)])
dD = mpf(D)**mpf("1.5")/(4*pi)*Lchi
zK = zeta(2)*Lchi

# --- coset sums over 2 O_K -------------------------------------------------------
y0 = sD/2; x0 = mpf(1)/2
B1, E1 = coset_BG(x0, y0, mpf(1)/2, mpf(0))
Ba, Ea = coset_BG(x0, y0, mpf(0), mpf(1)/2)
Bb, Eb = coset_BG(x0, y0, mpf(1)/2, mpf(1)/2)
B1, E1 = B1/16, E1/16
Ba, Ea = Ba/16, Ea/16
Bb, Eb = Bb/16, Eb/16
print("self-check E-odd = 15L/8:",
      nstr(fabs(E1+Ea+Eb - 15*L3/8), 3))
print("self-check B-odd = 15 zK/8:",
      nstr(fabs(B1+Ba+Bb - 15*zK/8), 3))
U = (2*E1 - Ea - Eb)/2
BU = (2*B1 - Ba - Bb)/2
MU = mpf(4*D)**mpf("1.5")/(4*pi**3)*U
DU = BU*3*mpf(D)**mpf("1.5")/(2*pi**3)
print("L = %s  M = %s" % (nstr(L3, 14), nstr(Mv, 14)))
print("U = %s  MU = %s" % (nstr(U, 14), nstr(MU, 14)))
print("BU = %s  DU = %s" % (nstr(BU, 14), nstr(DU, 14)))
print("zK = %s  dD = %s" % (nstr(zK, 14), nstr(dD, 14)))

# --- EK4 and the predicted identity ---------------------------------------------
EK4 = (10*y0/pi**3)*(-lattice_T(1, mpf(1)/2, y0) + 4*lattice_T(2, mpf(1)/2, y0))
pred = (mpf(40)/D)*Mv + (mpf(40)/(3*D))*MU + (mpf(10)/(3*D))*dD + (mpf(80)/(9*D))*DU
print("EK4 = %s" % nstr(EK4, 25))
print("pred = %s" % nstr(pred, 25))
print("|EK4 - pred| = %s" % nstr(fabs(EK4-pred), 3))

# pslq cross-check on the same basis
rel = pslq([EK4, Mv, MU, dD, DU], maxcoeff=10**7, maxsteps=10000)
if rel and rel[0] != 0:
    res = fabs(sum(r*v for r, v in zip(rel, [EK4, Mv, MU, dD, DU])))
    print("pslq [EK4,M,MU,dD,DU] = %s  residual %s" % (rel, nstr(res, 3)))
