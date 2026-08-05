# gen_conj_s4.py -- numeric reconnaissance: s4 at class-number-1 CM points
# NOT in Samart Table 6.  Candidates: odd fundamental discs D = 11, 19, 43,
# 67, 163 at tau = (1+sqrt(-D))/2, and the order disc -27 at tau=(1+3i s3)/2.
# (All pure-imaginary h=1 order discs -4,-8,-12,-16,-28 are already in
# Table 6, hence excluded by the dedup rule.)
# For each point: 60-dps s4, integrality/minpoly guess, |s4|>256 check, and a
# numerical scan of |s4(1/2+iy)| above the point (evidence that tau lies in
# the i-infty component of {|s4|>256}, Lemma 2.13 anchor).
# New file; modifies nothing.
from mpmath import mp, mpf, mpc, pi, exp, sqrt, pslq, nstr, fabs

mp.dps = 60

def eta_hp(tau, nterms=400):
    q = exp(2*pi*1j*tau)
    p = mpc(1); qn = q
    for n in range(1, nterms+1):
        p *= (1-qn); qn *= q
        if abs(qn) < mpf(10)**(-65): break
    return exp(pi*1j*tau/12)*p

def s4_hp(tau):
    e1, e2, e4 = eta_hp(tau), eta_hp(2*tau), eta_hp(4*tau)
    W = e1*e4**2/e2**3
    return (e2/e1)**24 * (16*W**4 + W**(-4))**4

# Table 6 s-values (for dedup reporting only)
T6 = ["256", "3656+/-2600s2", "26856+/-15300s3", "-144", "648",
      "143208+/-101574s2", "-1024", "2304", "1207368+... (4 vals)",
      "8292456+/-3132675s7", "81", "-3969", "-12288", "20736", "-82944",
      "-192303+/-85995s5 /2 (2 vals)", "614656", "-893952+/-516096s3",
      "347648256+/-141926400s6"]

CANDS = [(11, "fund disc -11, h=1"), (19, "fund disc -19, h=1"),
         (43, "fund disc -43, h=1"), (67, "fund disc -67, h=1"),
         (163, "fund disc -163, h=1")]
s3 = sqrt(mpf(3))

print("=== s4 values at (1+sqrt(-D))/2, D class-number-1 ===")
for D, note in CANDS:
    tau = mpc(mpf(1)/2, sqrt(mpf(D))/2)
    s = s4_hp(tau).real
    # integer guess, then quadratic guess
    r = pslq([s, 1])
    rel = pslq([s, sqrt(mpf(D))]) if r is None else None
    qrel = pslq([s*s, s, 1], maxcoeff=10**6) if r is None else None
    print("D=%-4d %s" % (D, note))
    print("  s4 = %s" % nstr(s, 55))
    print("  |s4| = %s  (>256: %s)" % (nstr(fabs(s), 8), fabs(s) > 256))
    if r is not None:
        print("  pslq [s4, 1]: %s  => s4 ~ %s" % (r, nstr(-mpf(r[1])/r[0], 12)))
    if rel is not None:
        print("  pslq [s4, sqrt D]: %s" % (rel,))
    if qrel is not None:
        print("  pslq [s4^2, s4, 1]: %s" % (qrel,))

print("=== order disc -27: tau = (1+3 i sqrt3)/2 ===")
tau = mpc(mpf(1)/2, 3*s3/2)
s = s4_hp(tau).real
r = pslq([s, 1])
qrel = pslq([s*s, s, 1], maxcoeff=10**6) if r is None else None
print("  s4 = %s" % nstr(s, 55))
print("  |s4| = %s (>256: %s)  pslq[s,1]: %s  pslq[s^2,s,1]: %s"
      % (nstr(fabs(s), 8), fabs(s) > 256, r, qrel))

# --- anchor component evidence: |s4(1/2+iy)| along the vertical ray --------
mp.dps = 40
print("=== scan |s4(1/2+iy)|, y from y_D upward (Lemma 2.13 component) ===")
for D, _ in CANDS + [(27, "")]:
    y0 = sqrt(mpf(D))/2 if D != 27 else 3*sqrt(mpf(3))/2
    mn = None
    for k in range(31):
        y = y0 + k*(mpf(8)/30)
        v = fabs(s4_hp(mpc(mpf(1)/2, y)))
        mn = v if mn is None else min(mn, v)
    print("D=%-4d min|s4| on [y0, y0+8] = %s  (>256: %s)"
          % (D, nstr(mn, 8), mn > 256))
