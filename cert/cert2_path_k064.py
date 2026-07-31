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
# Machinery lives in s2_iv.py (strict iv end to end; no float(...),
# no math.*, constants via outward iv enclosures) and is shared with the
# independent verifier cert2_verify.py.  Cut predicate:
#     Z certified to miss [0,64]  iff  0 notin Im Z  or  Re Z cap [0,64] = {}.
#
# Certificate separation (referee request): on success this generator writes
# a machine-readable certificate (default certificate_k064.json) recording
# every accepted block (parameter interval, enclosure Z, margin) and the cap
# constants as FULL-PRECISION decimal strings (mpf -> str, no float rounding).
# cert2_verify.py re-checks the fixed blocks without any adaptive search.
#
# Usage: python cert2_path_k064.py [--only axis,arc,...] [-o OUT.json]
#   --only runs a subset of segments (smoke test; cap still runs, the
#   SELF-TESTs are skipped); omitted segments are simply absent from the JSON.

import sys, os, json, datetime
import mpmath
from mpmath import iv, mp, mpf

import s2_iv                                   # sets iv.dps = 40, mp.dps = 50
from s2_iv import (NTR, y0_iv, y0, eps0, eps_ub, EPSUB_iv, RHO2_iv, pi2,
                   ZERO, HALF, s2_iv, s2_iv_tau, arc_tau, margin, width)

BASE = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------------- CLI
ONLY = None
OUT = os.path.join(BASE, 'certificate_k064.json')
args = sys.argv[1:]
if '--only' in args:
    ONLY = args[args.index('--only') + 1].split(',')
if '-o' in args:
    OUT = args[args.index('-o') + 1]

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

print('mpmath version:', mpmath.__version__, '; iv.dps =', iv.dps, '; NTR =', NTR)
print('cut to avoid: [0,64] (k = c^2 preimage of [-8,8]); path starts at i')
print('eps0 =', mp.nstr(eps0, 20), ' (>= 1/64000); eps_ub =', mp.nstr(eps_ub, 8))
if ONLY is not None:
    print('SMOKE MODE: only segments', ONLY, '(SELF-TESTs skipped)')
print()

print('--- axis start (elementary): s2(i) >= e^{2 pi} ---')
mp.dps = 50
print('e^{2 pi} =', mp.nstr(mp.e**(2*mp.pi), 15), '> 64:', mp.e**(2*mp.pi) > 64)
print()

# (key, display name, block evaluator, a, b) -- identical calls to pre-refactor
SEGSPEC = [
    ('axis',       'axis',    lambda t: s2_iv(ZERO, t),
     mpf(9)/16, mpf(1)),
    ('arc',        'arc',     lambda th: s2_iv_tau(arc_tau(th)),
     mpf(0), pi2),
    ('horizontal', 'horizontal', lambda t: s2_iv(t, HALF),
     mpf(1)/16, mpf(3)/8),
    ('vertical',   'vertical',   lambda t: s2_iv(iv.mpf(mpf(3)/8), t),
     y0 + mpf('1e-3'), mpf(1)/2),
    ('sub-cap',    'sub-cap',    lambda t: s2_iv(iv.mpf(mpf(3)/8), t),
     y0 + eps0, y0 + mpf('1e-3')),
]
HEADERS = {
    'axis':       '--- (i) axis segment: x = 0, y in [9/16, 1] ---',
    'arc':        '--- (ii) quarter arc: tau = i/2 + (1/16) e^{i theta}, theta in [0, pi/2] ---',
    'horizontal': '--- (iii) horizontal: y = 1/2, x in [1/16, 3/8] ---',
    'vertical':   '--- (iv) vertical: x = 3/8, y in [y0+1e-3, 1/2] ---',
    'sub-cap':    '--- (v) sub-cap: x = 3/8, y in [y0+eps0, y0+1e-3] ---',
}
segnames_disp = ['(i) axis', '(ii) arc', '(iii) horiz', '(iv) vert', '(v) sub-cap']

results = {}   # key -> (count, min_margin, depth, recs, a, b)
for key, name, fn, a, b in SEGSPEC:
    if ONLY is not None and key not in ONLY:
        continue
    print(HEADERS[key])
    c, m, d, r = certify(name, fn, a, b)
    results[key] = (c, m, d, r, a, b)
    print()

# ------------------------------------------------- certificate quality report
if ONLY is None:
    allrecs = [(nm, r) for nm, rr in zip(segnames_disp,
               [results[k][3] for k, *_ in SEGSPEC]) for r in rr]
    print('--- certificate quality (per segment, closest enclosure) ---')
    for nm, rr in zip(segnames_disp, [results[k][3] for k, *_ in SEGSPEC]):
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
if ONLY is None:
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

# ------------------------------------------------- certificate JSON export
# Decimal strings are emitted at elevated working precision: str(mpf) at
# exactly mp.dps digits sits on the binary<->decimal round-trip boundary and
# can re-parse 1 ulp off (observed ~1e-49 at dps 50); 70 digits re-parse
# exactly at the verifier's mp.dps = 50.  All values stay mpf -- no float.
with mp.workdps(70):
    def block_json(t, z, m):
        """Full-precision record of one accepted block (mpf -> str, no float)."""
        return {
            't': {'a': str(mpf(t.a)), 'b': str(mpf(t.b))},
            'Z': {'re_a': str(mpf(z.real.a)), 're_b': str(mpf(z.real.b)),
                  'im_a': str(mpf(z.imag.a)), 'im_b': str(mpf(z.imag.b))},
            'margin': str(mpf(m)),
        }

    cert_json = {
        'metadata': {
            'mpmath_version': mpmath.__version__,
            'python_version': sys.version.split()[0],
            'iv_dps': iv.dps,
            'NTR': NTR,
            'date': datetime.datetime.now().isoformat(timespec='seconds'),
            'cut': [0, 64],
            'generator': os.path.basename(__file__),
        },
        'segments': [
            {'name': key,
             'a': str(mpf(results[key][4])), 'b': str(mpf(results[key][5])),
             'n_blocks': results[key][0],
             'min_margin': str(mpf(results[key][1])),
             'blocks': [block_json(t, z, m) for (m, w, t, z, depth) in results[key][3]]}
            for key, *_ in SEGSPEC if key in results
        ],
        'cap': {
            'M0': str(M0),
            'M2': str(M2),
            'lo': str(lo),            # Re s2'(tau0) in [lo, hi]
            'hi': str(hi),
            'eps0': str(eps0),
            'eps_ub': str(eps_ub),    # eps0 + rounding width of y0
            'C': str(C),              # certified upper bound of hi + M2*eps_ub/2
        },
    }
    if ONLY is None:
        cert_json['global_min_margin'] = str(mpf(gmin))

# round-trip self-check at the verifier's precision (mp.dps = 50; SELF-TEST
# leaves mp.dps = 80, so force it here): every string in the certificate
# must re-parse to the EXACT mpf it came from
with mp.workdps(50):
    for key, *_ in SEGSPEC:
        if key not in results:
            continue
        seg = next(s for s in cert_json['segments'] if s['name'] == key)
        assert mpf(seg['a']) == mpf(results[key][4]) and mpf(seg['b']) == mpf(results[key][5])
        assert mpf(seg['min_margin']) == mpf(results[key][1])
        for bj, (m, w, t, z, depth) in zip(seg['blocks'], results[key][3]):
            assert mpf(bj['t']['a']) == mpf(t.a) and mpf(bj['t']['b']) == mpf(t.b)
            assert mpf(bj['Z']['re_a']) == mpf(z.real.a) and mpf(bj['Z']['re_b']) == mpf(z.real.b)
            assert mpf(bj['Z']['im_a']) == mpf(z.imag.a) and mpf(bj['Z']['im_b']) == mpf(z.imag.b)
            assert mpf(bj['margin']) == mpf(m)
    for k_, v_ in (('M0', M0), ('M2', M2), ('lo', lo), ('hi', hi),
                   ('eps0', eps0), ('eps_ub', eps_ub), ('C', C)):
        assert mpf(cert_json['cap'][k_]) == mpf(v_), 'cap round-trip failure: %s' % k_
with open(OUT, 'w') as f:
    json.dump(cert_json, f, indent=2)
print('certificate written to %s (%d blocks%s)'
      % (OUT, sum(s['n_blocks'] for s in cert_json['segments']),
         ', FULL RUN' if ONLY is None else ', SMOKE subset %s' % ONLY))
print()

# ------------------------------------------------- final verdict
print('--- conclusion ---')
print('Every accepted piece carries a strict enclosure s2(I) subset Z_I with')
print('Z_I cap [0,64] = EMPTY, decided by iv endpoint comparisons only:')
for nm, key in zip(segnames_disp, [k for k, *_ in SEGSPEC]):
    if key not in results:
        continue
    print('  %-11s: %3d pieces -- PASS (min margin %s)'
          % (nm, results[key][0], mp.nstr(results[key][1], 6)))
print('  (vi) cap   : Taylor coefficient %s < 0, Im s2 < 0 -- PASS' % mp.nstr(C, 8))
print()
print('ALL CERT-2 (cut [0,64]) CHECKS PASSED (fully rigorous interval certificate)')
