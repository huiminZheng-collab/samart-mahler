# Cert-2 (revised cut): rigorous interval-arithmetic certification that the
# extended path
#     i --(axis)--> 9i/16 --(quarter arc)--> 1/16 + i/2 --(horiz)-->
#     3/8 + i/2 --(vert)--> tau0 = (3 + sqrt(-7))/8
# stays in V = { tau : s2(tau) notin [0,64] } except for the endpoint tau0,
# where s2(tau0) = 1 (Cert-0) and the cap argument gives Im s2 < 0.
#
# Why [0,64]: the correct statement of Samart's theorem is
#     f2(s2(tau)) = 2 m(f + sqrt(s2(tau))) = EK(tau)   (shift sqrt(k), not k),
# so the branch cut to avoid in k-space is [0,64] -- the preimage of the
# c-plane cut [-8,8] under k = c^2. (At k=1, sqrt(1)=1, so the old [-8,8]
# formulation was accidentally equivalent and the error was invisible.
# Numerical check: EK(i/2) = 3.9803828365... = 2m(f+8), while 2m(f-64) ~
# 8.316.) Since s2(i/2) = 64 lies ON the new cut, the path starts at i
# (|s2| > 261 > 64, elementary) and detours around i/2 on a quarter circle
# of radius 1/16.
#
# Segments:
#   (i)   axis:     x = 0, y in [9/16, 1]
#   (ii)  arc:      tau = i/2 + (1/16) exp(i theta), theta in [0, pi/2]
#   (iii) horizon:  y = 1/2, x in [1/16, 3/8]
#   (iv)  vertical: x = 3/8, y in [y0+1e-3, 1/2],   y0 = sqrt(7)/8
#   (v)   sub-cap:  x = 3/8, y in [y0+eps0, y0+1e-3], eps0 = 1/64000
#   (vi)  cap:      y in (y0, y0+eps0], Taylor argument as in cert2_path.py
#                   (Im s2(tau0 + i eps) <= -13.549...*eps < 0, hence
#                   s2 misses the real interval [0,64] on the cap).
#
# Machinery identical to cert2_path.py (strict iv end to end; no float(...),
# no math.*, constants via outward iv enclosures). Cut predicate:
#     Z certified to miss [0,64]  iff  0 notin Im Z  or  Re Z cap [0,64] = {}.
#
# Arc evaluation rigor (mpmath 1.3.0, libmpi.py): exp(i theta) for a real
# interval theta is mpci_exp on the box 0 + i*theta; it computes
# r = mpi_exp([0,0]) (outward rounding) and (c,s) = mpi_cos_sin(theta)
# (directed-rounding cos/sin enclosure, valid for intervals of any width),
# then re = r*c, im = r*s with outward interval multiplication. Hence
# iv.exp(iv.mpc(0, theta_iv)) is a guaranteed rectangular enclosure of
# {exp(i t) : t in theta_iv}, and tau_iv = i/2 + (1/16)*that is a guaranteed
# enclosure of the arc block. Verified empirically in SELF-TEST below.

from mpmath import iv, mp, mpf

iv.dps = 40
mp.dps = 50   # must precede constants: mpf(iv-endpoint) rounding at default
              # 15 dps would cross the 40-dps iv endpoints (coverage gaps)
NTR = 120

# ---------------------------------------------------------------- constants
y0_iv = iv.sqrt(iv.mpf(7))/8       # interval containing sqrt(7)/8 exactly
Y0W_iv = y0_iv - y0_iv.a           # outward enclosure of the rounding width
Y0W = mpf(Y0W_iv.b)                # mpf upper bound of |y0_iv - sqrt(7)/8|
y0 = mpf(y0_iv.b)                  # mpf >= sqrt(7)/8 (safe side for the cap)

EPS0_iv = 1/iv.mpf(64000)          # interval containing 1/64000 exactly
eps0 = mpf(EPS0_iv.b)              # mpf >= 1/64000
EPSUB_iv = iv.mpf(eps0) + iv.mpf([0, Y0W])
eps_ub = mpf(EPSUB_iv.b)

RHO2_iv = (iv.mpf(2)/100)**2       # interval containing 0.02^2 exactly

PI2_iv = iv.pi/2                   # interval containing pi/2 exactly
pi2 = mpf(PI2_iv.b)                # mpf >= pi/2 (arc theta range covers pi/2)

ZERO = iv.mpf(0)
HALF = iv.mpf(1)/2                 # exact
R16 = iv.mpf(1)/16                 # exact
CENTER = iv.mpc(ZERO, HALF)        # i/2, exact

def cbox(E_iv):
    """encloses exp(u+iv) for |u|,|v| <= max(E_iv), all in iv arithmetic."""
    Eb = max(mpf(E_iv.b), mpf(10)**-38)
    box = iv.mpc(iv.mpf([-Eb, Eb]), iv.mpf([-Eb, Eb]))
    return iv.exp(box)

def s2_iv_tau(tau):
    """tau: complex interval. Returns complex interval Z with s2(t) in Z for
    every t in tau. Pure iv arithmetic throughout."""
    q = iv.exp(2*iv.pi*1j*tau)
    R_iv = abs(q)                  # iv enclosing |q| for all tau in the block
    assert R_iv.b < mpf('0.9'), '|q| too large for tail bound'
    P = iv.mpc(1)
    qodd = q
    for n in range(1, NTR+1):
        P *= (1 + qodd)**24        # rigorous complex integer power
        qodd *= q*q
    E_iv = 24*R_iv**(2*NTR+1)/((1 - R_iv)*(1 - R_iv**2))
    return P*cbox(E_iv)/q

def s2_iv(x, y):
    """x, y: real intervals (or exact mpf)."""
    return s2_iv_tau(iv.mpc(x, y))

def arc_tau(th_iv):
    """guaranteed enclosure of { i/2 + (1/16) exp(i t) : t in th_iv }."""
    return CENTER + R16*iv.exp(iv.mpc(ZERO, th_iv))

# ------------------------------------------------- avoidance predicate
KLO, KHI = 0, 64                   # the cut [0,64]

def margin(z):
    """Certified distance of z from [0,64], using ONLY iv endpoint
    comparisons; None if avoidance is not certified.
    min(|Im z|) if 0 notin Im z; else distance of Re z from [0,64]."""
    if 0 not in z.imag:
        return min(abs(mpf(z.imag.a)), abs(mpf(z.imag.b)))
    if z.real.b < KLO:
        return KLO - mpf(z.real.b)
    if z.real.a > KHI:
        return mpf(z.real.a) - KHI
    return None

def width(z):
    return max(mpf(z.real.b) - mpf(z.real.a), mpf(z.imag.b) - mpf(z.imag.a))

def certify(name, Zfn, a, b, maxdepth=40):
    """Certify s2 misses [0,64] on the blocks Zfn(t), t a bisection of [a,b].
    Returns (pieces, min_margin, max_depth, records)."""
    stack = [(iv.mpf([a, b]), 0)]
    count = 0
    min_marg = None
    deepest = 0
    recs = []
    while stack:
        t, depth = stack.pop()
        z = Zfn(t)
        m = margin(z)
        if m is not None:
            count += 1
            deepest = max(deepest, depth)
            min_marg = m if min_marg is None else min(min_marg, m)
            recs.append((m, width(z), t, z, depth))
            continue
        assert depth < maxdepth, 'subdivision exhausted on %s -- possible hit!' % name
        mid = (t.a + t.b)/2
        stack.append((iv.mpf([t.a, mid]), depth+1))
        stack.append((iv.mpf([mid, t.b]), depth+1))
    print('%s: certified in %d pieces, max depth %d, min margin %s'
          % (name, count, deepest, mp.nstr(min_marg, 8)))
    return count, min_marg, deepest, recs

print('mpmath version:', __import__('mpmath').__version__, '; iv.dps =', iv.dps, '; NTR =', NTR)
print('cut to avoid: [0,64] (k = c^2 preimage of [-8,8]); path starts at i')
print('eps0 =', mp.nstr(eps0, 20), ' (>= 1/64000); eps_ub =', mp.nstr(eps_ub, 8))
print()

print('--- axis start (elementary): s2(i) >= e^{2 pi} ---')
mp.dps = 50
print('e^{2 pi} =', mp.nstr(mp.e**(2*mp.pi), 15), '> 64:', mp.e**(2*mp.pi) > 64)
print()

print('--- (i) axis segment: x = 0, y in [9/16, 1] ---')
c1, m1, d1, r1 = certify('axis', lambda t: s2_iv(ZERO, t), mpf(9)/16, mpf(1))
print()
print('--- (ii) quarter arc: tau = i/2 + (1/16) e^{i theta}, theta in [0, pi/2] ---')
c2, m2, d2, r2 = certify('arc', lambda th: s2_iv_tau(arc_tau(th)), mpf(0), pi2)
print()
print('--- (iii) horizontal: y = 1/2, x in [1/16, 3/8] ---')
c3, m3, d3, r3 = certify('horizontal', lambda t: s2_iv(t, HALF), mpf(1)/16, mpf(3)/8)
print()
print('--- (iv) vertical: x = 3/8, y in [y0+1e-3, 1/2] ---')
c4, m4, d4, r4 = certify('vertical', lambda t: s2_iv(iv.mpf(mpf(3)/8), t),
                         y0 + mpf('1e-3'), mpf(1)/2)
print()
print('--- (v) sub-cap: x = 3/8, y in [y0+eps0, y0+1e-3] ---')
c5, m5, d5, r5 = certify('sub-cap', lambda t: s2_iv(iv.mpf(mpf(3)/8), t),
                         y0 + eps0, y0 + mpf('1e-3'))
print()

# ------------------------------------------------- certificate quality report
segnames = ['(i) axis', '(ii) arc', '(iii) horiz', '(iv) vert', '(v) sub-cap']
allrecs = [(nm, r) for nm, rr in zip(segnames, [r1, r2, r3, r4, r5]) for r in rr]
print('--- certificate quality (per segment, closest enclosure) ---')
for nm, rr in zip(segnames, [r1, r2, r3, r4, r5]):
    m, w, t, z, depth = min(rr, key=lambda rec: rec[0])
    print('%s: closest to [0,64] at t in [%s, %s] (depth %d)'
          % (nm, mp.nstr(mpf(t.a), 12), mp.nstr(mpf(t.b), 12), depth))
    print('    Z = Re [%s, %s]' % (mp.nstr(mpf(z.real.a), 10), mp.nstr(mpf(z.real.b), 10)))
    print('        Im [%s, %s]' % (mp.nstr(mpf(z.imag.a), 10), mp.nstr(mpf(z.imag.b), 10)))
    print('    width(Z) = %s ; dist(Z, [0,64]) = %s ; margin/width = %s'
          % (mp.nstr(w, 6), mp.nstr(m, 6), mp.nstr(m/w, 6)))
worst = min(allrecs, key=lambda w: w[1][0])
widest = max(allrecs, key=lambda w: w[1][1])
print('global closest: segment %s, dist = %s, width = %s, margin/width = %s'
      % (worst[0], mp.nstr(worst[1][0], 8), mp.nstr(worst[1][1], 6),
         mp.nstr(worst[1][0]/worst[1][1], 6)))
print('global widest : segment %s, width = %s, dist = %s'
      % (widest[0], mp.nstr(widest[1][1], 6), mp.nstr(widest[1][0], 8)))
gmin = worst[1][0]
print('GLOBAL minimum margin over all %d pieces: %s' % (len(allrecs), mp.nstr(gmin, 8)))
print()

print('--- (vi) cap (y0, y0+eps0]: Taylor argument (same as cert2_path.py) ---')
bx = iv.mpf(mpf(3)/8) + iv.mpf([mpf('-0.03'), mpf('0.03')])
by = y0_iv + iv.mpf([mpf('-0.03'), mpf('0.03')])
M0_iv = abs(s2_iv(bx, by))
M0 = mpf(M0_iv.b)
print('M0 =', mp.nstr(M0, 20))
M2_iv = 2*M0_iv/RHO2_iv
M2 = mpf(M2_iv.b)
print('M2 =', mp.nstr(M2, 20))
Y_iv = y0_iv + iv.mpf(eps0)
ze = s2_iv(iv.mpf(mpf(3)/8), Y_iv)
imdiff = (ze - 1).imag
EPSPOS_iv = iv.mpf(eps0) + (y0_iv - y0_iv.b)   # iv containing every eps' > 0
ERR0 = (M2_iv*EPSUB_iv*EPSUB_iv/2).b           # bounds M2 eps'^2/2
S2P_iv = (imdiff + iv.mpf([-ERR0, ERR0]))/EPSPOS_iv
lo, hi = mpf(S2P_iv.a), mpf(S2P_iv.b)
print('Re s2\'(tau0) in [%s, %s]' % (mp.nstr(lo, 15), mp.nstr(hi, 15)))
ERRUB_iv = M2_iv*EPSUB_iv/2
C = mpf((S2P_iv.b + ERRUB_iv).b)
print('coefficient bound hi + M2*eps_ub/2 <=', mp.nstr(C, 15))
print('=> Im s2(tau0 + i eps) <= eps * (%s) < 0 for all eps in (0, %s]'
      % (mp.nstr(C, 12), mp.nstr(eps_ub, 6)))
assert C < 0, 'cap argument FAILED'
print('cap OK: Im s2 < 0 strictly, so s2 misses the real interval [0,64] on the cap')
print()

# ------------------------------------------------- SELF-TEST
mp.dps = 80
def s2_hp(x, y):
    """High-precision (non-interval) s2; tail error < 1e-200 << 1e-80."""
    tau = mp.mpc(x, y)
    q = mp.exp(2*mp.pi*1j*tau)
    P = mp.mpc(1)
    qodd = q
    for n in range(1, NTR+1):
        P *= (1 + qodd)**24
        qodd *= q*q
    return P/q

pts = []
for j in range(4):                 # (i) axis
    pts.append((mpf(0), mpf(9)/16 + (mpf(1) - mpf(9)/16)*j/3))
for j in range(6):                 # (ii) arc, 6 sample points
    th = mp.pi*j/10
    pts.append((mp.cos(th)/16, mpf(1)/2 + mp.sin(th)/16))
for j in range(4):                 # (iii) horizontal
    pts.append((mpf(1)/16 + (mpf(3)/8 - mpf(1)/16)*j/3, mpf(1)/2))
for j in range(4):                 # (iv) vertical
    pts.append((mpf(3)/8, y0 + mpf('1e-3') + (mpf(1)/2 - y0 - mpf('1e-3'))*j/3))
pts += [(mpf(3)/8, y0 + eps0), (mpf(3)/8, y0 + (eps0 + mpf('1e-3'))/2)]
assert len(pts) == 20
ok = 0
for x, y in pts:
    z = s2_iv(iv.mpf(x), iv.mpf(y))
    tv = s2_hp(x, y)
    if z.real.a <= tv.real <= z.real.b and z.imag.a <= tv.imag <= z.imag.b:
        ok += 1
    else:
        print('  CONTAINMENT FAILURE at', mp.nstr(x, 6), mp.nstr(y, 6))
print('SELF-TEST 1 (iv enclosure contains 80-dps truth): %d/20 passed' % ok)
assert ok == 20

# Enclosure monotonicity on the arc: point enclosure inside block enclosure.
ok2 = 0
for j in range(5):
    a = mpf(j)/10
    b = a + mpf(1)/20
    zb = s2_iv_tau(arc_tau(iv.mpf([a, b])))
    thp = a + mpf(1)/40
    zp = s2_iv_tau(arc_tau(iv.mpf(thp)))
    if zb.real.a <= zp.real.a and zp.real.b <= zb.real.b \
       and zb.imag.a <= zp.imag.a and zp.imag.b <= zb.imag.b:
        ok2 += 1
print('SELF-TEST 2 (arc point enclosure inside block enclosure): %d/5 passed' % ok2)
assert ok2 == 5
print()

# ------------------------------------------------- final verdict
print('--- conclusion ---')
print('Every accepted piece carries a strict enclosure s2(I) subset Z_I with')
print('Z_I cap [0,64] = EMPTY, decided by iv endpoint comparisons only:')
for (nm, c, m) in zip(segnames, [c1, c2, c3, c4, c5], [m1, m2, m3, m4, m5]):
    print('  %-11s: %3d pieces -- PASS (min margin %s)' % (nm, c, mp.nstr(m, 6)))
print('  (vi) cap   : Taylor coefficient %s < 0, Im s2 < 0 -- PASS' % mp.nstr(C, 8))
print()
print('ALL CERT-2 (cut [0,64]) CHECKS PASSED (fully rigorous interval certificate)')
