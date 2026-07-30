# Find tau0 with s2(q(tau0)) = 1 via the signature-2 hypergeometric relation:
#   s2(q2(alpha)) = 16/(alpha(1-alpha)),  q2(alpha) = exp(-pi F(1-alpha)/F(alpha)),
#   F = 2F1(1/2,1/2;1; . ),  q = e^{2 pi i tau}  =>  tau = (i/2) F(1-alpha)/F(alpha).
# s2 = 1 => alpha(1-alpha) = 16 => alpha = (1 +- 3 sqrt(-7))/2.
# Then evaluate f2 via the EK formula at tau0 and compare with 8 L'(g7,0).

from mpmath import mp, mpf, mpc, pi, exp, sqrt, sin, cos, hyper

mp.dps = 60

def F(z):
    return hyper([mpf(1) / 2, mpf(1) / 2], [mpf(1)], z)

alpha = (1 + 3 * sqrt(mpc(-7))) / 2
print("alpha =", alpha)
tau0 = (1j / 2) * F(1 - alpha) / F(alpha)
print("tau0 =", tau0)
print("Im tau0 =", tau0.imag)

def eta(tau, nterms=400):
    q = exp(2 * pi * 1j * tau)
    p = mpc(1)
    qn = q
    for n in range(1, nterms + 1):
        p *= (1 - qn)
        qn *= q
        if abs(qn) < mpf(10) ** (-70):
            break
    return exp(pi * 1j * tau / 12) * p

def s2(tau):
    return -eta(tau + mpf(1) / 2) ** 24 / eta(2 * tau + 1) ** 24

print("s2(tau0) =", s2(tau0), " (expect 1)")

def U(j, tau):
    tot = mpc(0)
    for m in range(1, 3000):
        z = j * pi * m * tau
        term = cos(z) / sin(z) ** 3 / m
        tot += term
        if abs(term) < mpf(10) ** (-65):
            break
    return 2 * pi ** 3 * tot

val = 2 * pi * tau0 + (2 / pi ** 3) * (U(1, tau0) - 4 * U(4, tau0))
print("holomorphic mtilde-value (complex) =", val)
print("Im part (= f2 via EK) =", val.imag)
target = mpf("0.82137286223121608968365275918634331199115861671374")
print("8 L'(g7,0)     =", target)
print("diff           =", val.imag - target)
print("2*true m       =", target)
