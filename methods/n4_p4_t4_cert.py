# Cert-2 (n4, T4): rigorous interval-arithmetic certification that the path
#     gamma:  i --(horiz, y=1)--> 1/2 + i --(vert, x=1/2)--> tau4 = (1+sqrt(-2))/2
# stays in the good region
#     V4 = { tau in H : c(tau) NOTIN Ast } ,
# where c(tau) = (eta(2t)/eta(t))^6 (16 A^4 + A^-4), A = eta(t)eta(4t)^2/eta(2t)^3,
# and Ast = { c in C : |Re c|^{2/3} + |Im c|^{2/3} <= 4^{2/3} } is the closed
# ASTROID DISK = critical image P(T^3), P = (x^4+y^4+z^4+1)/(xyz), of the n4
# family.  Avoidance of Ast implies avoidance of every bad set of the n4
# machine (the slit cross [-4,4] u i[-4,4] lies inside Ast).
#
# POSITION OF tau4 (60-dps pre-scan, this file's SELF-TEST, and step-0):
#   c(tau4) = 1.5128932208...(1+i);  astroid functional = 2.6357... >
#   4^{2/3} = 2.5198..., margin ~ 0.116 > 0.  Since c is holomorphic (open
#   mapping) and Ast is closed, tau4 is an INTERIOR point of V4.
#   Along x = 1/2 the functional decreases monotonically from 4.38 (y=1);
#   it crosses the threshold near y ~ 0.693 (margin +0.057 at y = 0.702,
#   -0.011 at y = 0.690).  Hence the vertical leg certified down to
#   y0 = 0.702 < 1/sqrt2 = 0.70711 stays strictly inside V4, and tau4 is
#   an interior point of the certified leg (0.702 < 1/sqrt2 < 1).
#
# CONSEQUENCE (no boundary-continuity step needed): although tau4 lies on
# the boundary {Im tau = 1/sqrt2} of Samart's Prop 2.1(iii) statement and
# the Prop-2.1 proof does NOT cover tau4 (q = -e^{-pi sqrt2} < 0 real and
# |s4| = 20.96 < 256, both hypotheses fail), the certified path lands tau4
# in the SAME connected component W of V4 as the anchor disk D around i
# (where the identity n4(s4(tau)) = EK4(tau) holds by Prop 2.1(iii):
# q = e^{-2 pi} in (0,1), |s4| > 256, Im = 1 >= 1/sqrt2).  The propagation
# lemma (formal.tex Section 3, lem:propagation + certified-path criterion,
# subsec:criterion) then gives the identity at tau4 itself.
#
# Method: identical to cert2_path_n4_m144.py (ported tail-bound machinery);
# the only changes are the leg-B lower endpoint (0.702 instead of
# sqrt3/2 ~ 0.866) and the |q| cap: |q| <= e^{-2 pi * 0.702} < 0.0122, so
# the tail bound assert is relaxed from 0.01 to 0.02 (r^{NTR+1} <
# 0.0122^41 < 1e-78, still far below the iv.dps = 40 working precision).
#
# STRICTNESS POLICY: mpmath.iv interval arithmetic end to end (same rigor
# notes as cert2_path.py: outward rounding, complex rectangles, exp of boxes,
# integer powers with guard bits; verified for mpmath 1.3.0).  No float()
# conversions on the certification path.  SELF-TEST at the end checks
# containment against 80-dps eta-quotient values (20 sample points) and
# enclosure monotonicity (5/5).

from mpmath import iv, mp, mpf

iv.dps = 40
mp.dps = 50   # must precede constants (mpf rounding of iv endpoints)
NTR = 40      # |q| <= e^{-2 pi * 0.702} < 0.0122 on the whole path;
              # r^{41} < 1e-78

# ---------------------------------------------------------------- constants
THR_iv = iv.mpf(4)**(iv.mpf(2)/3)      # interval enclosing 4^{2/3}
TWO3_iv = iv.mpf(2)/3

def cbox(E_iv):
    """encloses exp(u+iv) for |u|,|v| <= max(E_iv), all in iv arithmetic."""
    Eb = E_iv.b
    Eb = max(Eb, mpf(10)**-38)
    box = iv.mpc(iv.mpf([-Eb, Eb]), iv.mpf([-Eb, Eb]))
    return iv.exp(box)

def prod_pow_iv(q, d, e, R_iv):
    """prod_{n=1}^{NTR} (1 - q^{d n})^e with rigorous complex tail.
    Tail log-bound E = |e| r^{d(NTR+1)}/((1-r^d)(1-r)), f increasing in r."""
    P = iv.mpc(1)
    qd = q**d
    qn = qd
    for n in range(1, NTR+1):
        P *= (1 - qn)
        qn *= qd
    E_iv = abs(e) * R_iv**(d*(NTR+1)) / ((1 - R_iv**d)*(1 - R_iv))
    return (P * cbox(E_iv))**e

def c_iv(x, y):
    """x, y: real intervals. Returns a complex interval Z with c(tau) in Z
    for every tau in x + i*y. Pure iv arithmetic throughout."""
    tau = iv.mpc(x, y)
    q = iv.exp(2*iv.pi*1j*tau)
    R_iv = abs(q)
    assert R_iv.b < mpf('0.02'), '|q| too large for tail bound'
    # P6 = prod (1+q^n)^6 = prod (1-q^{2n})^6 / (1-q^n)^6
    P6 = prod_pow_iv(q, 2, 6, R_iv) / prod_pow_iv(q, 1, 6, R_iv)
    # R = prod (1-q^n)^4 (1-q^{4n})^8 / (1-q^{2n})^{12}
    Rv = (prod_pow_iv(q, 1, 4, R_iv) * prod_pow_iv(q, 4, 8, R_iv)
          / prod_pow_iv(q, 2, 12, R_iv))
    qm4 = iv.exp(-iv.pi*1j*tau/2)          # q^{-1/4} = e^{-pi i tau/2}
    return qm4 * P6 * (16*q*Rv + Rv**(-1))

# ------------------------------------------------- avoidance predicate
def ast_margin(z):
    """Certified lower bound of dist to Ast in the astroid functional:
    min_{Z} (|x|^{2/3}+|y|^{2/3}) - 4^{2/3} if positive; else None."""
    rx, ry = z.real, z.imag
    mx = mpf(0) if 0 in rx else min(abs(mpf(rx.a)), abs(mpf(rx.b)))
    my = mpf(0) if 0 in ry else min(abs(mpf(ry.a)), abs(mpf(ry.b)))
    lhs = iv.mpf(mx)**TWO3_iv + iv.mpf(my)**TWO3_iv
    if lhs.a > THR_iv.b:
        return mpf(lhs.a) - mpf(THR_iv.b)
    return None

def width(z):
    return max(mpf(z.real.b) - mpf(z.real.a), mpf(z.imag.b) - mpf(z.imag.a))

def certify_segment(name, axis, a, b, fx, fy, maxdepth=40):
    """Certify c misses Ast on { axis = t in [a,b] } with the other coordinate
    given by the mpf functions fx(t), fy(t). Adaptive bisection."""
    stack = [(iv.mpf([a, b]), 0)]
    count = 0
    min_marg = None
    deepest = 0
    recs = []
    while stack:
        t, depth = stack.pop()
        # affine coordinate boxes: co(t) = c0 + d*t over the block t (d in {0,1})
        x = iv.mpf(fx.c0) + fx.d*iv.mpf(t)
        y = iv.mpf(fy.c0) + fy.d*iv.mpf(t)
        z = c_iv(x, y)
        m = ast_margin(z)
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

class Lin:
    """linear coordinate function co(t) = c0 + d*(t - t0), exact in iv."""
    def __init__(self, c0, d):
        self.c0, self.d = c0, d
    def __call__(self, tm):
        return self.c0

print('mpmath version:', __import__('mpmath').__version__, '; iv.dps =', iv.dps, '; NTR =', NTR)
print('4^{2/3} in [%s, %s]' % (mp.nstr(mpf(THR_iv.a), 15), mp.nstr(mpf(THR_iv.b), 15)))
print()

# --- anchor disk D around i: certify c(D) outside Ast and |c| > 4 on D ----
print('--- anchor disk D: |Re tau| <= 1/32, |Im tau - 1| <= 1/32 ---')
bx = iv.mpf([mpf(-1)/32, mpf(1)/32])
by = iv.mpf([mpf(31)/32, mpf(33)/32])
zD = c_iv(bx, by)
mD = ast_margin(zD)
print('c(D) subset Z = Re [%s, %s]' % (mp.nstr(mpf(zD.real.a), 8), mp.nstr(mpf(zD.real.b), 8)))
print('              Im [%s, %s]' % (mp.nstr(mpf(zD.imag.a), 8), mp.nstr(mpf(zD.imag.b), 8)))
print('astroid margin on D:', mp.nstr(mD, 8) if mD is not None else 'NOT CERTIFIED')
assert mD is not None, 'anchor disk not certified outside Ast'
absD = abs(zD)
print('|c| on D >= %s (>4: %s); |s4| = |c|^4 >= %s (>256: %s)'
      % (mp.nstr(mpf(absD.a), 8), absD.a > 4,
         mp.nstr(mpf((absD**4).a), 8), (absD**4).a > 256))
assert absD.a > 4 and (absD**4).a > 256
print('D certified: D subset V4 and |s4| > 256 on D (Rogers convergence).')
print()

# --- leg (A): horizontal, y = 1, x in [0, 1/2] --------------------------
print('--- leg (A): y = 1, x in [0, 1/2] ---')
fxA, fyA = Lin(None, 1), Lin(1, 0)
fxA.c0 = 0
cA, mA, dA, rA = certify_segment('leg A (horiz)', 'x', mpf(0), mpf(1)/2, fxA, fyA)
print()

# --- leg (B): vertical, x = 1/2, y in [0.702, 1] -------------------------
# Lower endpoint 0.702 < 1/sqrt2 = 0.70710678..., so tau4 = 1/2 + i/sqrt2 is
# an INTERIOR point of the certified leg.  Pre-scan (60 dps): astroid
# functional at y = 0.702 is ~ 2.577, margin ~ +0.057 > 0; the crossing
# into Ast happens near y ~ 0.693, strictly below the leg.
print('--- leg (B): x = 1/2, y in [0.702, 1] ---')
fxB, fyB = Lin(mpf(1)/2, 0), Lin(0, 1)
cB, mB, dB, rB = certify_segment('leg B (vert)', 'y', mpf('0.702'), mpf(1),
                                 fxB, fyB)
print()

# ------------------------------------------------- certificate quality report
allrecs = [('A', r) for r in rA] + [('B', r) for r in rB]
worst = min(allrecs, key=lambda w: w[1][0])
widest = max(allrecs, key=lambda w: w[1][1])
print('--- certificate quality ---')
for tag, rec in [('closest to Ast', worst), ('widest enclosure', widest)]:
    seg, (m, w, t, z, depth) = rec
    print('%s: leg %s, t in [%s, %s] (depth %d)'
          % (tag, seg, mp.nstr(mpf(t.a), 12), mp.nstr(mpf(t.b), 12), depth))
    print('    Z = Re [%s, %s]' % (mp.nstr(mpf(z.real.a), 10), mp.nstr(mpf(z.real.b), 10)))
    print('        Im [%s, %s]' % (mp.nstr(mpf(z.imag.a), 10), mp.nstr(mpf(z.imag.b), 10)))
    print('    width(Z) = %s ; astroid margin = %s ; margin/width = %s'
          % (mp.nstr(w, 6), mp.nstr(m, 6), mp.nstr(m/w, 6)))
gmin = worst[1][0]
print('GLOBAL minimum margin over all %d pieces: %s' % (len(allrecs), mp.nstr(gmin, 8)))
print()

# ------------------------------------------------- SELF-TEST
mp.dps = 80
def eta_hp(tau, nterms=300):
    q = mp.exp(2*mp.pi*1j*tau)
    p = mp.mpc(1)
    qn = q
    for n in range(1, nterms+1):
        p *= (1-qn)
        qn *= q
        if abs(qn) < mpf(10)**(-85):
            break
    return mp.exp(mp.pi*1j*tau/12)*p

def c_hp(x, y):
    tau = mp.mpc(x, y)
    A = eta_hp(tau)*eta_hp(4*tau)**2/eta_hp(2*tau)**3
    return (eta_hp(2*tau)/eta_hp(tau))**6*(16*A**4 + A**-4)

pts = []
for j in range(10):                # leg A
    pts.append((mpf(j)/20, mpf(1)))
for j in range(10):                # leg B
    pts.append((mpf(1)/2, 1 - (1 - mpf('0.702'))*mpf(j)/9))
pts.append((mpf(1)/2, 1/mp.sqrt(2)))   # tau4 itself
assert len(pts) == 21
ok = 0
for x, y in pts:
    z = c_iv(iv.mpf(x), iv.mpf(y))
    tv = c_hp(x, y)
    if z.real.a <= tv.real <= z.real.b and z.imag.a <= tv.imag <= z.imag.b:
        ok += 1
    else:
        print('  CONTAINMENT FAILURE at', mp.nstr(x, 6), mp.nstr(y, 6))
print('SELF-TEST 1 (iv enclosure contains 80-dps eta value): %d/21 passed' % ok)
assert ok == 21

# tau4 strictly inside leg B's y-interval and strictly outside Ast:
z4 = c_iv(iv.mpf(mpf(1)/2), iv.mpf(1/iv.sqrt(2)))
m4 = ast_margin(z4)
print('tau4 enclosure: astroid margin = %s (must be > 0)'
      % (mp.nstr(m4, 8) if m4 is not None else 'NOT CERTIFIED'))
assert m4 is not None and m4 > 0
assert mpf('0.702') < 1/mp.sqrt(2) < 1

ok2 = 0
for j in range(5):
    a = mpf(j)/20
    b = a + mpf(1)/40
    zb = c_iv(iv.mpf([a, b]), iv.mpf(1))
    zp = c_iv(iv.mpf(a + mpf(1)/80), iv.mpf(1))
    if zb.real.a <= zp.real.a and zp.real.b <= zb.real.b \
       and zb.imag.a <= zp.imag.a and zp.imag.b <= zb.imag.b:
        ok2 += 1
print('SELF-TEST 2 (point enclosure inside block enclosure): %d/5 passed' % ok2)
assert ok2 == 5
print()

# ------------------------------------------------- final verdict
print('--- conclusion ---')
print('Every accepted block carries a strict enclosure c(I) subset Z_I with')
print('Z_I cap Ast = EMPTY, decided by iv endpoint comparisons only:')
print('  (D) anchor disk : 1 box -- PASS (margin %s, |s4|>256)' % mp.nstr(mD, 6))
print('  (A) horizontal  : %3d pieces -- PASS (min margin %s)' % (cA, mp.nstr(mA, 6)))
print('  (B) vertical    : %3d pieces -- PASS (min margin %s)' % (cB, mp.nstr(mB, 6)))
print('Endpoint tau4 = (1+sqrt(-2))/2 lies INSIDE V4 (c(tau4) strictly outside')
print('Ast, certified margin above), so the path lands tau4 in W itself;')
print('no boundary-continuity step is needed.')
print()
print('ALL CERT-2 (n4, T4) CHECKS PASSED (fully rigorous interval certificate)')
