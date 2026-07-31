# Shared rigorous interval-arithmetic core for s2(tau), extracted verbatim
# from cert2_path_k064.py so that the certificate generator
# (cert2_path_k064.py) and the independent verifier (cert2_verify.py)
# evaluate s2 through ONE shared implementation.  Numerical behaviour is
# identical to the original script: strict iv end to end; no float(...),
# no math.*, constants via outward iv enclosures.
#
# s2(tau) = q^{-1} prod_{n>=1} (1 + q^{2n-1})^24,  q = e^{2 pi i tau},
# evaluated by truncating the eta-type product at n = NTR and enclosing the
# tail factor prod_{n>NTR} (1 + q^{2n-1})^24 in exp(box) with
# box = [-E,E] + i[-E,E],
#     E = 24 R^{2 NTR + 1} / ((1 - R)(1 - R^2)),  R = |q|  (tail lemma).
#
# Arc evaluation rigor (mpmath 1.3.0, libmpi.py): exp(i theta) for a real
# interval theta is mpci_exp on the box 0 + i*theta; it computes
# r = mpi_exp([0,0]) (outward rounding) and (c,s) = mpi_cos_sin(theta)
# (directed-rounding cos/sin enclosure, valid for intervals of any width),
# then re = r*c, im = r*s with outward interval multiplication. Hence
# iv.exp(iv.mpc(0, theta_iv)) is a guaranteed rectangular enclosure of
# {exp(i t) : t in theta_iv}, and tau_iv = i/2 + (1/16)*that is a guaranteed
# enclosure of the arc block.

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
