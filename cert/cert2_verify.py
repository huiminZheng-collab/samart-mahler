# Independent verifier for the Cert-2 (cut [0,64]) interval certificate
# produced by cert2_path_k064.py (default: certificate_k064.json).
#
# Referee-requested separation of concerns: this script performs NO adaptive
# subdivision and NO search.  It reads the fixed certificate data and
# re-verifies, using the SAME shared evaluation core (s2_iv.py):
#   1. for every recorded block, recompute the strict iv enclosure Z of s2 on
#      the recorded parameter interval and assert the recorded Z encloses the
#      recomputed one (in practice they coincide exactly);
#   2. assert Z cap [0,64] = EMPTY with a positive margin, and that the
#      recorded margin matches the recomputation;
#   3. assert the blocks of each segment chain end-to-end (no gaps) and cover
#      the whole certified parameter range [a, b] of that segment;
#   4. re-check the cap inequality: recompute hi + M2*eps_ub/2 from the JSON
#      constants (outward iv) and assert it is strictly negative and no
#      larger than the recorded coefficient bound C;
#   5. metadata consistency (mpmath version, iv.dps, NTR).
# Prints per-segment block counts and margins; exits 0 with
# 'CERTIFICATE VERIFIED' on success, nonzero (AssertionError) on any failure.
#
# Usage: python cert2_verify.py [certificate.json]

import sys, os, json
import mpmath
from mpmath import iv, mp, mpf

import s2_iv                                   # sets iv.dps = 40, mp.dps = 50
from s2_iv import NTR, ZERO, HALF, s2_iv, s2_iv_tau, arc_tau, margin, KLO, KHI

BASE = os.path.dirname(os.path.abspath(__file__))

# block evaluators, identical to the generator's SEGSPEC functions
BUILDERS = {
    'axis':       lambda t: s2_iv(ZERO, t),
    'arc':        lambda th: s2_iv_tau(arc_tau(th)),
    'horizontal': lambda t: s2_iv(t, HALF),
    'vertical':   lambda t: s2_iv(iv.mpf(mpf(3)/8), t),
    'sub-cap':    lambda t: s2_iv(iv.mpf(mpf(3)/8), t),
}

def main(path):
    with open(path) as f:
        data = json.load(f)

    meta = data['metadata']
    print('verifying certificate:', path)
    print('metadata: mpmath %s, python %s, iv.dps %s, NTR %s, date %s'
          % (meta['mpmath_version'], meta['python_version'], meta['iv_dps'],
             meta['NTR'], meta['date']))
    assert meta['mpmath_version'] == mpmath.__version__, 'mpmath version mismatch'
    assert int(meta['iv_dps']) == iv.dps, 'iv.dps mismatch: rebuild with iv.dps %s' % meta['iv_dps']
    assert int(meta['NTR']) == NTR, 'NTR mismatch'
    print()

    gmin = None
    total = 0
    for seg in data['segments']:
        name = seg['name']
        fn = BUILDERS[name]
        a, b = mpf(seg['a']), mpf(seg['b'])
        blocks = sorted(seg['blocks'], key=lambda blk: mpf(blk['t']['a']))
        assert len(blocks) == int(seg['n_blocks']), '%s: block count mismatch' % name

        # coverage: first block starts at/below a, last ends at/above b,
        # consecutive blocks touch or overlap (bisection tiling -> no gaps)
        assert mpf(blocks[0]['t']['a']) <= a, '%s: blocks do not cover left endpoint' % name
        assert mpf(blocks[-1]['t']['b']) >= b, '%s: blocks do not cover right endpoint' % name

        mmin = None
        prev_b = None
        for blk in blocks:
            ta, tb = mpf(blk['t']['a']), mpf(blk['t']['b'])
            assert ta < tb, '%s: degenerate block' % name
            if prev_b is not None:
                assert ta <= prev_b, '%s: GAP between consecutive blocks at t=%s' % (name, mp.nstr(ta, 12))
            prev_b = tb

            # recompute the strict enclosure on the recorded parameter interval
            t = iv.mpf([ta, tb])
            z = fn(t)
            rz = blk['Z']
            # recorded Z must enclose the recomputed enclosure
            assert mpf(rz['re_a']) <= z.real.a and z.real.b <= mpf(rz['re_b']), \
                '%s: recorded Re Z is not an outer enclosure' % name
            assert mpf(rz['im_a']) <= z.imag.a and z.imag.b <= mpf(rz['im_b']), \
                '%s: recorded Im Z is not an outer enclosure' % name

            # avoidance of the cut [0,64], from iv endpoint comparisons only
            assert (0 not in z.imag) or z.real.b < KLO or z.real.a > KHI, \
                '%s: enclosure INTERSECTS [0,64]' % name
            m = margin(z)
            assert m is not None and m > 0, '%s: non-positive margin' % name
            assert mpf(blk['margin']) == m, '%s: recorded margin mismatch' % name
            mmin = m if mmin is None else min(mmin, m)

        total += len(blocks)
        gmin = mmin if gmin is None else min(gmin, mmin)
        assert mpf(seg['min_margin']) == mmin, '%s: recorded min_margin mismatch' % name
        print('  %-11s: %3d blocks re-verified -- PASS (min margin %s)'
              % (name, len(blocks), mp.nstr(mmin, 6)))

    print()
    print('  total %d blocks, GLOBAL minimum margin %s' % (total, mp.nstr(gmin, 8)))
    if 'global_min_margin' in data:
        assert mpf(data['global_min_margin']) == gmin, 'recorded global min margin mismatch'

    # ------------------------------------------------ cap (Taylor argument)
    cap = data.get('cap')
    if cap is not None:
        M0 = mpf(cap['M0'])
        M2 = mpf(cap['M2'])
        lo, hi = mpf(cap['lo']), mpf(cap['hi'])
        eps0 = mpf(cap['eps0'])
        eps_ub = mpf(cap['eps_ub'])
        C = mpf(cap['C'])
        assert M0 > 0 and M2 > 0, 'cap: non-positive M0/M2'
        assert lo <= hi, 'cap: empty Re s2\'(tau0) enclosure'
        assert eps_ub >= eps0 > 0, 'cap: bad eps bounds'
        # recompute the coefficient bound from the JSON constants (outward iv)
        C_check = mpf((iv.mpf(hi) + iv.mpf(M2)*iv.mpf(eps_ub)/2).b)
        print('  cap: recomputed hi + M2*eps_ub/2 = %s (recorded bound %s)'
              % (mp.nstr(C_check, 12), mp.nstr(C, 12)))
        assert C_check <= C, 'cap: recorded C is not an upper bound'
        assert C_check < 0, 'cap: coefficient bound NOT strictly negative'
        print('  cap        : Im s2(tau0 + i eps) <= eps * (%s) < 0 on (0, %s] -- PASS'
              % (mp.nstr(C_check, 8), mp.nstr(eps_ub, 6)))

    print()
    print('CERTIFICATE VERIFIED')
    return 0

if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, 'certificate_k064.json')
    sys.exit(main(path))
