# gen_conj_minpoly.py -- minimal-polynomial hunt for s4((1+sqrt-D)/2),
# D = 11, 19, 27, 43, 67, 163 (h=1 discs NOT in Samart Table 6).
# s4 values are real, non-integer, near j(D)-640 (numerically).  Try
# degree <= 4 pslq relations at 110 dps.  New file; modifies nothing.
from mpmath import mp, mpf, mpc, pi, exp, sqrt, pslq, nstr, fabs

mp.dps = 110

def eta_hp(tau):
    q = exp(2*pi*1j*tau)
    p = mpc(1); qn = q
    for n in range(1, 600):
        p *= (1-qn); qn *= q
        if abs(qn) < mpf(10)**(-115): break
    return exp(pi*1j*tau/12)*p

def s4_hp(tau):
    e1, e2, e4 = eta_hp(tau), eta_hp(2*tau), eta_hp(4*tau)
    W = e1*e4**2/e2**3
    return (e2/e1)**24 * (16*W**4 + W**(-4))**4

s3 = sqrt(mpf(3))
for D in (11, 19, 27, 43, 67, 163):
    y0 = sqrt(mpf(D))/2 if D != 27 else 3*s3/2
    s = s4_hp(mpc(mpf(1)/2, y0)).real
    print("D = %d" % D)
    print("  s4 = %s" % nstr(s, 80))
    found = False
    for deg, maxc in ((3, 10**8), (4, 10**6)):
        vec = [s**k for k in range(deg, 0, -1)] + [mpf(1)]
        rel = pslq(vec, maxcoeff=maxc, maxsteps=20000)
        if rel is not None and rel[0] != 0:
            res = fabs(sum(r*v for r, v in zip(rel, vec)))
            print("  deg-%d rel = %s  residual = %s" % (deg, rel, nstr(res, 3)))
            found = True
            break
    if not found:
        print("  no deg<=4 relation found")
