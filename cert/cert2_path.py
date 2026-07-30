# Cert-2: rigorous interval-arithmetic certification that the path
#     gamma:  i/2 --(horiz)--> 3/8 + i/2 --(vert)--> tau0 = (3 + sqrt(-7))/8
# stays in V = { tau : s2(tau) notin [-8,8] } except for the endpoint tau0,
# where s2(tau0) = 1 (Cert-0).
#
# Also: elementary bound s2(it) >= e^{2 pi t} > 8 for t >= 1/2 (imaginary axis
# lies in V) -- needed for Step 1 of (P2).
#
# Method: s2(tau) = q^-1 prod_{n>=1} (1 + q^{2n-1})^24,  q = e^{2 pi i tau}.
# Truncate at N; the tail factor prod_{n>N}(1+q^{2n-1})^24 = exp(S) with
#     |S| <= E := 24 r^{2N+1} / ((1-r)(1-r^2)),  r = upper bound of |q|,
# since |log(1+z)| <= -log(1-|z|) <= |z|/(1-r) for |z| <= r < 1.
# The tail is enclosed as exp of the COMPLEX box |u|,|v| <= E (bounding the
# argument as well as the modulus), backported from cert2_path_n2pair.py.
#
# STRICTNESS POLICY (revision): the certificate path uses mpmath.iv interval
# arithmetic END TO END. No float(...) conversions, no math.* functions, no
# ordinary double arithmetic anywhere on the path from the enclosure of q to
# the final inequalities. Constants that are not exactly representable
# (1/64000, 0.02, sqrt(7)/8) enter through outward-rounded iv constructions,
# so the rounding of constants is accounted for inside the enclosures.
#
# mpmath.iv rigor (verified for mpmath 1.3.0, see site-packages):
#   * real interval arithmetic rounds endpoints outward (MPFR-style directed
#     rounding, round_floor / round_ceiling) -- ctx_iv.py, libmpi.py;
#   * monotonic functions (exp, log, sqrt, atan) are evaluated at the two
#     endpoints with floor resp. ceiling rounding (libmpi.py mpi_exp etc.);
#   * complex intervals are rectangles; arithmetic/multiplication rounds
#     outward; complex exp uses exp(re)*cis(im) with guard digits and outward
#     rounding (libmpi.py mpci_exp); integer powers use binary exponentiation
#     at prec+20 guard bits with a final outward round (mpci_pow_int);
#   * iv.pi is an interval enclosing pi (ivmpf_constant, round_floor/ceiling).
# The empirical self-test at the end (SELF-TEST section) verifies containment
# against 80-dps high-precision values on 20 sample points, plus enclosure
# monotonicity (point enclosure inside block enclosure), 20/20 and 5/5.

from mpmath import iv, mp, mpf

iv.dps = 40
mp.dps = 50   # must precede constants: mpf(iv-endpoint) rounding at default
              # 15 dps would cross the 40-dps iv endpoints (coverage gaps)
NTR = 120

# ---------------------------------------------------------------- constants
# Every non-representable constant enters through an outward iv enclosure.

y0_iv = iv.sqrt(iv.mpf(7))/8       # interval containing sqrt(7)/8 exactly
Y0W_iv = y0_iv - y0_iv.a           # outward enclosure of the rounding width
Y0W = mpf(Y0W_iv.b)                # mpf upper bound of |y0_iv - sqrt(7)/8|
y0 = mpf(y0_iv.b)                  # mpf >= sqrt(7)/8 (safe side for the cap)

EPS0_iv = 1/iv.mpf(64000)          # interval containing 1/64000 exactly
eps0 = mpf(EPS0_iv.b)              # mpf >= 1/64000 (delta/64, delta = 1e-3)
# eps_ub bounds every eps' = y - sqrt(7)/8 with y <= y0_iv.b + eps0:
EPSUB_iv = iv.mpf(eps0) + iv.mpf([0, Y0W])
eps_ub = mpf(EPSUB_iv.b)

RHO2_iv = (iv.mpf(2)/100)**2       # interval containing 0.02^2 exactly

def cbox(E_iv):
    """encloses exp(u+iv) for |u|,|v| <= max(E_iv), all in iv arithmetic."""
    Eb = E_iv.b                    # mpf upper endpoint of the iv bound
    Eb = max(Eb, mpf(10)**-38)
    box = iv.mpc(iv.mpf([-Eb, Eb]), iv.mpf([-Eb, Eb]))
    return iv.exp(box)

def s2_iv(x, y):
    """x, y: real intervals (or exact mpf). Returns a complex interval Z with
    s2(tau) in Z for every tau in x + i*y. Pure iv arithmetic throughout."""
    tau = iv.mpc(x, y)
    q = iv.exp(2*iv.pi*1j*tau)
    R_iv = abs(q)                  # iv enclosing |q| for all tau in the block
    assert R_iv.b < mpf('0.9'), '|q| too large for tail bound'
    P = iv.mpc(1)
    qodd = q
    for n in range(1, NTR+1):
        P *= (1 + qodd)**24        # rigorous complex integer power
        qodd *= q*q
    # Tail bound E = 24 r^{2N+1}/((1-r)(1-r^2)), fully in iv arithmetic;
    # f(r) is increasing on [0,1), so E_iv.b bounds E at r = max |q|.
    E_iv = 24*R_iv**(2*NTR+1)/((1 - R_iv)*(1 - R_iv**2))
    return P*cbox(E_iv)/q

# ------------------------------------------------- avoidance predicate
def margin(z):
    """Certified distance of z from [-8,8], using ONLY iv endpoint
    comparisons; None if avoidance is not certified.
    min(|Im z|) if 0 notin Im z; else distance of Re z from [-8,8]."""
    if 0 not in z.imag:
        return min(abs(mpf(z.imag.a)), abs(mpf(z.imag.b)))
    if z.real.b < -8:
        return -8 - mpf(z.real.b)
    if z.real.a > 8:
        return mpf(z.real.a) - 8
    return None

def misses_K(z):
    return margin(z) is not None

def width(z):
    return max(mpf(z.real.b) - mpf(z.real.a), mpf(z.imag.b) - mpf(z.imag.a))

def certify_segment(name, axis, a, b, fixed, maxdepth=40):
    """Certify s2 misses [-8,8] on { axis = t in [a,b], other coord = fixed }.
    Adaptive bisection. Returns (pieces, min_margin, max_depth, records)."""
    stack = [(iv.mpf([a, b]), 0)]
    count = 0
    min_marg = None
    deepest = 0
    recs = []
    while stack:
        t, depth = stack.pop()
        z = s2_iv(t, fixed) if axis == 'x' else s2_iv(fixed, t)
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
print('eps0 =', mp.nstr(eps0, 20), ' (>= 1/64000); eps_ub =', mp.nstr(eps_ub, 8))
print()

print('--- axis bound (Step 1 of P2) ---')
mp.dps = 50
epi = mp.e**mp.pi
print('s2(it) >= e^{2 pi t} >= e^pi =', mp.nstr(epi, 15))
print('e^pi > 8:', epi > 8)
print()

print('--- horizontal segment (A): y = 1/2, x in [0, 3/8] ---')
c1, m1, d1, r1 = certify_segment('horizontal', 'x', mpf(0), mpf(3)/8, mpf(1)/2)
print()
print('--- vertical segment (B): x = 3/8, y in [y0+1e-3, 1/2] ---')
c2, m2, d2, r2 = certify_segment('vertical', 'y', y0 + mpf('1e-3'), mpf(1)/2, mpf(3)/8)
print()
print('--- sub-cap (B2): x = 3/8, y in [y0+eps0, y0+1e-3] ---')
c3, m3, d3, r3 = certify_segment('sub-cap', 'y', y0 + eps0, y0 + mpf('1e-3'), mpf(3)/8)
print()

# ------------------------------------------------- certificate quality report
allrecs = [('A', r) for r in r1] + [('B', r) for r in r2] + [('B2', r) for r in r3]
worst = min(allrecs, key=lambda w: w[1][0])          # closest to [-8,8]
widest = max(allrecs, key=lambda w: w[1][1])         # widest enclosure
print('--- certificate quality ---')
for tag, rec in [('closest to [-8,8]', worst), ('widest enclosure', widest)]:
    seg, (m, w, t, z, depth) = rec
    print('%s: segment %s, t in [%s, %s] (depth %d)'
          % (tag, seg, mp.nstr(mpf(t.a), 12), mp.nstr(mpf(t.b), 12), depth))
    print('    Z = Re [%s, %s]' % (mp.nstr(mpf(z.real.a), 10), mp.nstr(mpf(z.real.b), 10)))
    print('        Im [%s, %s]' % (mp.nstr(mpf(z.imag.a), 10), mp.nstr(mpf(z.imag.b), 10)))
    print('    width(Z) = %s ; dist(Z, [-8,8]) = %s ; margin/width = %s'
          % (mp.nstr(w, 6), mp.nstr(m, 6), mp.nstr(m/w, 6)))
gmin = worst[1][0]
print('GLOBAL minimum margin over all %d pieces: %s' % (len(allrecs), mp.nstr(gmin, 8)))
print()

print('--- cap (y0, y0+eps0]: Taylor argument ---')
# M0 = upper bound of |s2| on the box |Re tau - 3/8| <= 0.03,
# |Im tau - sqrt(7)/8| <= 0.03 (y0 rounding absorbed by the iv box).
bx = iv.mpf(mpf(3)/8) + iv.mpf([mpf('-0.03'), mpf('0.03')])
by = y0_iv + iv.mpf([mpf('-0.03'), mpf('0.03')])
zbox = s2_iv(bx, by)
M0_iv = abs(zbox)
M0 = mpf(M0_iv.b)                  # strict mpf upper bound, straight from iv
print('M0 =', mp.nstr(M0, 20))
# Cauchy: |s2''(w)| <= 2 M0 / rho^2 for |w - tau0| <= 0.01, rho = 0.02;
# the disc of radius 0.02 about any such w lies inside the 0.03-box.
M2_iv = 2*M0_iv/RHO2_iv
M2 = mpf(M2_iv.b)
print('M2 =', mp.nstr(M2, 20))
# Re s2'(tau0):  s2(tau0 + i eps') - 1 = i eps' s2'(tau0) + R,
# |R| <= M2 eps'^2/2, where eps' = Im w - sqrt(7)/8 in
# Eps = eps0 + (y0_iv - sqrt(7)/8)  (contains the true eps').
Y_iv = y0_iv + iv.mpf(eps0)        # evaluation ordinate block (width ~1e-40)
ze = s2_iv(iv.mpf(mpf(3)/8), Y_iv)
imdiff = (ze - 1).imag
EPSPOS_iv = iv.mpf(eps0) + (y0_iv - y0_iv.b)   # iv containing every eps' > 0
ERR0 = (M2_iv*EPSUB_iv*EPSUB_iv/2).b           # bounds M2 eps'^2/2
S2P_iv = (imdiff + iv.mpf([-ERR0, ERR0]))/EPSPOS_iv
lo, hi = mpf(S2P_iv.a), mpf(S2P_iv.b)
print('Re s2\'(tau0) in [%s, %s]' % (mp.nstr(lo, 15), mp.nstr(hi, 15)))
# For eps in (0, eps_ub]:  Im s2(tau0 + i eps) = eps Re s2'(tau0) + Im R,
# |Im R| <= M2 eps_ub/2 * eps, hence Im s2 <= eps*(hi + ERR_UB).
ERRUB_iv = M2_iv*EPSUB_iv/2
C_iv = S2P_iv.b + ERRUB_iv         # iv enclosing the coefficient
C = mpf(C_iv.b)
print('coefficient bound hi + M2*eps_ub/2 <=', mp.nstr(C, 15))
print('=> Im s2(tau0 + i eps) <= eps * (%s) for all eps in (0, %s]'
      % (mp.nstr(C, 12), mp.nstr(eps_ub, 6)))
print('   upper bound of eps*(hi + M2*eps0/2): eps_ub*C =',
      mp.nstr(eps_ub*C, 12), '(scale -13.5 * eps)')
assert C < 0, 'cap argument FAILED'
print('cap OK: coefficient strictly negative, so s2 misses [-8,8] on the cap')
print()

# ------------------------------------------------- SELF-TEST
# Empirical verification that iv enclosures contain 80-dps "truth" values.
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
for j in range(8):                 # horizontal segment
    pts.append((mpf(3)*j/56, mpf(1)/2))
for j in range(7):                 # vertical segment
    pts.append((mpf(3)/8, y0 + mpf('1e-3') + (mpf(1)/2 - y0 - mpf('1e-3'))*j/6))
pts += [(mpf(3)/8, y0 + eps0), (mpf(3)/8, y0 + (eps0 + mpf('1e-3'))/2),
        (mpf(3)/8, y0 + mpf('1e-3')),
        (mpf(3)/8, y0 + eps0/4), (mpf(3)/8, y0 + 3*eps0/4)]
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

# Enclosure monotonicity: point enclosure inside enclosing block enclosure.
ok2 = 0
for j in range(5):
    a = mpf(j)/20
    b = a + mpf(1)/40
    zb = s2_iv(iv.mpf([a, b]), iv.mpf(mpf(1)/2))
    zp = s2_iv(iv.mpf(a + mpf(1)/80), iv.mpf(mpf(1)/2))
    if zb.real.a <= zp.real.a and zp.real.b <= zb.real.b \
       and zb.imag.a <= zp.imag.a and zp.imag.b <= zb.imag.b:
        ok2 += 1
print('SELF-TEST 2 (point enclosure inside block enclosure): %d/5 passed' % ok2)
assert ok2 == 5
print()

# ------------------------------------------------- final verdict
print('--- conclusion ---')
print('Every accepted piece carries a strict enclosure s2(I) subset Z_I with')
print('Z_I cap [-8,8] = EMPTY, decided by iv endpoint comparisons only:')
print('  (A) horizontal : %3d pieces -- PASS (min margin %s)' % (c1, mp.nstr(m1, 6)))
print('  (B) vertical   : %3d pieces -- PASS (min margin %s)' % (c2, mp.nstr(m2, 6)))
print('  (B2) sub-cap   : %3d pieces -- PASS (min margin %s)' % (c3, mp.nstr(m3, 6)))
print('  (C) cap        : Taylor bound coefficient %s < 0 -- PASS' % mp.nstr(C, 8))
print()
print('ALL CERT-2 CHECKS PASSED (fully rigorous interval certificate)')
