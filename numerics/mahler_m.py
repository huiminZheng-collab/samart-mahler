# Mahler measure m = m( (x+1/x)(y+1/y)(z+1/z) + 1 )
# On the torus x=e^{i t1}, y=e^{i t2}, z=e^{i t3}:  P = 8 cos t1 cos t2 cos t3 + 1.
# m = E[ log|1 + 8 c1 c2 c3| ].
#
# Integrate over t3 first (Jensen):
#   J(A) = (1/2pi) int_0^{2pi} log|1 + A cos t| dt
#        = log((1+sqrt(1-A^2))/2)  if |A|<=1
#        = log(|A|/2)              if |A|>=1
# Symmetry reduces to [0, pi/2]^2:
#   m = (4/pi^2) int_0^{pi/2} int_0^{pi/2} J(8 cos t1 cos t2) dt2 dt1
# Integrand has a sqrt-type cusp along 8 cos t1 cos t2 = 1: split there.

from mpmath import mp, mpf, pi, log, sqrt, cos, acos, quad

mp.dps = 40

def J(A):
    a2 = A * A
    if a2 <= 1:
        return log((1 + sqrt(1 - a2)) / 2)
    else:
        return log(abs(A) / 2)

def inner(t1):
    c = cos(t1)
    cusp = 1 / (8 * c)
    if cusp >= 1:            # no cusp in t2 range
        return quad(lambda t2: J(8 * c * cos(t2)), [0, pi / 2])
    t2star = acos(cusp)
    return quad(lambda t2: J(8 * c * cos(t2)), [0, t2star, pi / 2])

# outer integral; for t1 > acos(1/8) there is no cusp at all
t1split = acos(mpf(1) / 8)
outer = quad(inner, [0, t1split, pi / 2], maxdegree=8)
m = 4 * outer / pi ** 2
print("m =", m)

# comparison targets from the L-value computation:
Lp = mpf("0.10267160777890201121045659489829291399889482708922")
print("4 L'(g7,0) =", 4 * Lp, " diff:", m - 4 * Lp)
print("8 L'(g7,0) =", 8 * Lp, " diff:", m - 8 * Lp)
