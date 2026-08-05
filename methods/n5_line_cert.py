# Cert-2 (n4, Phase 8 line certificate): rigorous interval-arithmetic
# certification that the path
#     gamma:  i --(horiz, y=1)--> 1/2 + i --(vert, x=1/2)--> 1/2 + i sqrt(21)/2
# stays in the good region
#     V4 = { tau in H : c(tau) NOTIN Ast } ,
# c(tau) = (eta(2t)/eta(t))^6 (16 A^4 + A^-4), A = eta(t)eta(4t)^2/eta(2t)^3,
# Ast = { c : |Re c|^{2/3} + |Im c|^{2/3} <= 4^{2/3} } (critical image).
#
# ONE certificate for THREE open entries of Samart's Table 6, whose CM
# points all lie on the line x = 1/2:
#   #2: tau = (1+i sqrt6)/2  = (1/2, sqrt6/2  = 1.2247...)
#   #6: tau = (1+i sqrt15)/2 = (1/2, sqrt15/2 = 1.9365...)
#   #8: tau = (1+i sqrt21)/2 = (1/2, sqrt21/2 = 2.2913...)
# (Entry #1 at (1/2, 1) is the top corner of the legs; entry #7 at
# (1/2, sqrt21/6 = 0.7638...) lies on the lower part of leg B; both were
# already covered by n4_p4_t4_cert.py, whose leg B is the y in [0.702, 1]
# part of this one.  Entry #8's endpoint is certified interior as well.)
#
# CONSEQUENCE: the whole path lands in the same connected component W of
# V4 as the anchor disk D around i (where n4(s4(tau)) = EK4(tau) holds by
# Samart Prop 2.1(iii)), so the propagation theorem gives the identity at
# every point of the path, in particular at the three targets.  Each
# target is certified to be a strict INTERIOR point of V4 (positive
# astroid margin on its point enclosure), so no boundary-continuity step
# is needed for any of them.
#
# Method: identical to n4_p4_t4_cert.py (same tail-bound machinery,
# NTR = 40; the binding |q| cap is at the bottom y = 0.702:
# |q| <= e^{-2 pi * 0.702} < 0.0122, r^41 < 1e-78).  The extension upward
# is strictly safer: larger y means smaller |q| and larger astroid margin.
#
# STRICTNESS POLICY: mpmath.iv end to end, no float() conversions on the
# certification path.  SELF-TESTs: 23/23 point containments against
# 80-dps eta values; 5/5 enclosure monotonicity; interior-point verdicts
# for the three targets.

from mpmath import iv, mp, mpf

iv.dps = 40
mp.dps = 50   # must precede constants (mpf rounding of iv endpoints)
NTR = 40      # |q| <= e^{-2 pi * 0.702} < 0.0122 on the whole path;
              # r^{41} < 1e-78

# ---------------------------------------------------------------- constants
THR_iv = iv.mpf(4)**(iv.mpf(2)/3)      # interval enclosing 4^{2/3}
TWO3_iv = iv.mpf(2)/3
SQ21 = mpf(mp.sqrt(21))/2              # top of leg B

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
    P6 = prod_pow_iv(q, 2, 6, R_iv) / prod_pow_iv(q, 1, 6, R_iv)
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
    """Certify c misses Ast on { axis = t in [a,b] } with the other
    coordinate given by the linear functions fx(t), fy(t)."""
    stack = [(iv.mpf([a, b]), 0)]
    count = 0
    min_marg = None
    deepest = 0
    recs = []
    while stack:
        t, depth = stack.pop()
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
        assert depth < maxdepth, \
            'subdivision exhausted on %s -- possible hit!' % name
        mid = (t.a + t.b)/2
        stack.append((iv.mpf([t.a, mid]), depth+1))
        stack.append((iv.mpf([mid, t.b]), depth+1))
    print('%s: certified in %d pieces, max depth %d, min margin %s'
          % (name, count, deepest, mp.nstr(min_marg, 8)))
    return count, min_marg, deepest, recs

class Lin:
    """linear coordinate function co(t) = c0 + d*t, exact in iv."""
    def __init__(self, c0, d):
        self.c0, self.d = c0, d

print('mpmath version:', __import__('mpmath').__version__, '; iv.dps =',
      iv.dps, '; NTR =', NTR)
print('4^{2/3} in [%s, %s]' % (mp.nstr(mpf(THR_iv.a), 15),
                               mp.nstr(mpf(THR_iv.b), 15)))
print()

# --- anchor disk D around i ------------------------------------------
print('--- anchor disk D: |Re tau| <= 1/32, |Im tau - 1| <= 1/32 ---')
bx = iv.mpf([mpf(-1)/32, mpf(1)/32])
by = iv.mpf([mpf(31)/32, mpf(33)/32])
zD = c_iv(bx, by)
mD = ast_margin(zD)
print('astroid margin on D:', mp.nstr(mD, 8) if mD is not None
      else 'NOT CERTIFIED')
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
fxA, fyA = Lin(0, 1), Lin(1, 0)
cA, mA, dA, rA = certify_segment('leg A (horiz)', 'x', mpf(0), mpf(1)/2,
                                 fxA, fyA)
print()

# --- leg (B): vertical, x = 1/2, y in [0.702, sqrt(21)/2] ---------------
print('--- leg (B): x = 1/2, y in [0.702, sqrt(21)/2 = %s] ---'
      % mp.nstr(SQ21, 10))
fxB, fyB = Lin(mpf(1)/2, 0), Lin(0, 1)
cB, mB, dB, rB = certify_segment('leg B (vert)', 'y', mpf('0.702'), SQ21,
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
    print('    width(Z) = %s ; astroid margin = %s ; margin/width = %s'
          % (mp.nstr(w, 6), mp.nstr(m, 6), mp.nstr(m/w, 6)))
gmin = worst[1][0]
print('GLOBAL minimum margin over all %d pieces: %s'
      % (len(allrecs), mp.nstr(gmin, 8)))
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
    pts.append((mpf(1)/2, mpf('0.702') + (SQ21 - mpf('0.702'))*mpf(j)/9))
targets = [('#2', mpf(mp.sqrt(6))/2), ('#6', mpf(mp.sqrt(15))/2),
           ('#8', SQ21)]
for _, yy in targets:
    pts.append((mpf(1)/2, yy))
assert len(pts) == 23
ok = 0
for x, y in pts:
    z = c_iv(iv.mpf(x), iv.mpf(y))
    tv = c_hp(x, y)
    if z.real.a <= tv.real <= z.real.b and z.imag.a <= tv.imag <= z.imag.b:
        ok += 1
    else:
        print('  CONTAINMENT FAILURE at', mp.nstr(x, 6), mp.nstr(y, 6))
print('SELF-TEST 1 (iv enclosure contains 80-dps eta value): %d/23 passed'
      % ok)
assert ok == 23

# each target strictly outside Ast (interior-point verdict):
for tag, yy in targets:
    zt = c_iv(iv.mpf(mpf(1)/2), iv.mpf(yy))
    mt = ast_margin(zt)
    print('target %s (y = %s): astroid margin = %s (must be > 0)'
          % (tag, mp.nstr(yy, 10),
             mp.nstr(mt, 8) if mt is not None else 'NOT CERTIFIED'))
    assert mt is not None and mt > 0

ok2 = 0
for j in range(5):
    a = mpf(j)/20
    b = a + mpf(1)/40
    zb = c_iv(iv.mpf([a, b]), iv.mpf(1))
    zp = c_iv(iv.mpf(a + mpf(1)/80), iv.mpf(1))
    if zb.real.a <= zp.real.a and zp.real.b <= zb.real.b \
       and zb.imag.a <= zp.imag.a and zp.imag.b <= zb.imag.b:
        ok2 += 1
print('SELF-TEST 2 (point enclosure inside block enclosure): %d/5 passed'
      % ok2)
assert ok2 == 5
print()

# ------------------------------------------------- final verdict
print('--- conclusion ---')
print('Every accepted block carries a strict enclosure c(I) subset Z_I with')
print('Z_I cap Ast = EMPTY, decided by iv endpoint comparisons only:')
print('  (D) anchor disk : 1 box -- PASS (margin %s, |s4|>256)'
      % mp.nstr(mD, 6))
print('  (A) horizontal  : %3d pieces -- PASS (min margin %s)'
      % (cA, mp.nstr(mA, 6)))
print('  (B) vertical    : %3d pieces -- PASS (min margin %s)'
      % (cB, mp.nstr(mB, 6)))
print('Targets #2, #6, #8 lie strictly INSIDE V4 on the certified leg B,')
print('hence in W; the propagation theorem gives n4(s4(tau)) = EK4(tau)')
print('at all three CM points.  No boundary-continuity step needed.')
print()
print('ALL CERT-2 (n4, Phase 8 line) CHECKS PASSED'
      ' (fully rigorous interval certificate)')
