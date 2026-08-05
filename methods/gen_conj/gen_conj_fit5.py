# gen_conj_fit5.py -- D=27 branch: tau = (1+3 i s3)/2, O_3 = Z[3 w] in
# K = Q(sqrt-3), h(-27) = 1.  Lattice structure:
#   Lambda_1 = O_3,  Lambda_2 = O_6 (conductor 6),  O_6 = 2 O_3 ⊔ (1+2 O_3)
#   EK4 = (10 y0/pi^3)(-T(O_3)+4 T(O_6))
#   comb = (1/2) B(O_3) + (8/3) BU + 2 L27 + (16/3) U,
#   B(O_3) = 6 zeta_K(2) - 2 X,  X = sum_{gamma in w+O_3} 1/N^2.
# Predicted (M27 = L'(g27,0), MU = 108^{3/2} U/(4 pi^3),
# DU = BU*9 s3/(2 pi^3), DX = X*9 s3/(2 pi^3), d3 = L'(chi_{-3},-1)):
#   EK4 = (40/27) M27 + (40/81) MU + 10 d3 - (10/3) DX + (80/9) DU.
# (DX turns out to equal (56/27) d3, so this simplifies to
#  EK4 = (40/27) M27 + (40/81) MU + (250/81) d3 + (80/9) DU.)
# New file; modifies nothing.
from mpmath import (mp, mpf, mpc, pi, zeta, sinh, cosh, cos, exp, sqrt,
                    diff as mpdiff, dirichlet, gamma, power, gammainc,
                    pslq, nstr, fabs)
mp.dps = 60
s3 = sqrt(mpf(3))

def G_row(x, y):
    return (pi/y)*sinh(2*pi*y)/(cosh(2*pi*y) - cos(2*pi*x))

def row_powers(x, y):
    Gv = lambda yy: G_row(x, yy)
    S2 = -mpdiff(Gv, y)/(2*y)
    S3 = -mpdiff(lambda yy: -mpdiff(Gv, yy)/(2*yy), y)/(4*y)
    return S2, S3

def coset_BG(x0, y0, sx, st):
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

chi3 = lambda n: 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)

# newform of level 27 (order disc -27)
NMAX = 2500
BX = int(NMAX**0.5)+4; BY = int(2*(NMAX/27.0)**0.5)+4
A = [0]*(NMAX+1)
for b in range(-BY, BY+1):
    for a in range(-BX-3*abs(b), BX+3*abs(b)+1):
        n = a*a + 3*a*b + 9*b*b
        if 1 <= n <= NMAX:
            A[n] += (2*a+3*b)**2 - 27*b*b
aa = [0 if n % 3 == 0 else v//8 for n, v in enumerate(A)]

def mellin_I(a, xN, s):
    s = mpf(s); tot = mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0: continue
        tot += a[n]*power(2*pi*n, -s)*gammainc(s, 2*pi*n/xN)
    return tot

xN = sqrt(mpf(27))
I0, I3 = mellin_I(aa, xN, 0), mellin_I(aa, xN, 3)
M27 = xN**3*I3 + I0
L27 = M27*(2*pi)**3/(xN**3*gamma(3))

d3 = mpf(3)**mpf("1.5")/(4*pi)*dirichlet(mpf(2), [0, 1, -1])
zK = zeta(2)*dirichlet(mpf(2), [0, 1, -1])

x0, y0 = mpf(1)/2, 3*s3/2
# mod-2 cosets of O_3: 1 -> (1/2,0); 3w -> (1/2,1/2); 1+3w -> (0,1/2)
B1, E1 = coset_BG(x0, y0, mpf(1)/2, mpf(0))
Ba, Ea = coset_BG(x0, y0, mpf(1)/2, mpf(1)/2)
Bb, Eb = coset_BG(x0, y0, mpf(0), mpf(1)/2)
B1, E1 = B1/16, E1/16; Ba, Ea = Ba/16, Ea/16; Bb, Eb = Bb/16, Eb/16
print("self-check E-odd' = 15 L27/8:", nstr(fabs(E1+Ea+Eb - 15*L27/8), 3))
# mod-3 cosets: w + O_3 -> (1/3,1/3); 2w + O_3 -> (2/3,2/3)
Xa, _ = coset_BG(x0, y0, mpf(1)/3, mpf(1)/3)
Xb, _ = coset_BG(x0, y0, mpf(2)/3, mpf(2)/3)
print("self-check X(w) = X(2w):", nstr(fabs(Xa-Xb), 3), " (diff)")
X = Xa
B_O3 = 6*zK - 2*X
print("self-check B-odd' = 15 B(O_3)/16:",
      nstr(fabs(B1+Ba+Bb - 15*B_O3/16), 3))

U = (2*E1 - Ea - Eb)/2
BU = (2*B1 - Ba - Bb)/2
MU = mpf(108)**mpf("1.5")/(4*pi**3)*U
DU = BU*9*s3/(2*pi**3)
DX = X*9*s3/(2*pi**3)
print("L27 = %s  M27 = %s" % (nstr(L27, 14), nstr(M27, 14)))
print("U = %s  MU = %s" % (nstr(U, 14), nstr(MU, 14)))
print("BU = %s  DU = %s  DX = %s" % (nstr(BU, 14), nstr(DU, 14), nstr(DX, 14)))

EK4 = (10*y0/pi**3)*(-lattice_T(1, mpf(1)/2, y0) + 4*lattice_T(2, mpf(1)/2, y0))
pred = mpf(40)/27*M27 + mpf(40)/81*MU + 10*d3 - mpf(10)/3*DX + mpf(80)/9*DU
print("EK4 = %s" % nstr(EK4, 25))
print("pred = %s" % nstr(pred, 25))
print("|EK4 - pred| = %s" % nstr(fabs(EK4-pred), 3))
rel = pslq([EK4, M27, MU, d3, DX, DU], maxcoeff=10**7, maxsteps=10000)
if rel and rel[0] != 0:
    res = fabs(sum(r*v for r, v in zip(rel, [EK4, M27, MU, d3, DX, DU])))
    print("pslq = %s  residual %s" % (rel, nstr(res, 3)))
