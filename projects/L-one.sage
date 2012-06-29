# some of this code is probably wrong

from psage.number_fields.sqrt5 import *
from psage.ellcurve.lseries.aplist_sqrt5 import aplist
import mpmath

sigma2, sigma = F.embeddings(RR)
phi = sigma(a)
phibar = sigma2(a)

def make_tot_pos(n):
    """
    Given an element n of Q(sqrt5), return a totally positive unit multiple
    of n; i.e. return m such that m/n is a unit and m is taken to a positive
    number under both embeddings of K into R.
    """
    if sigma(n) > 0:
        if sigma2(n) < 0:
            return n * a
        else:
            return n
    else:
        if sigma2(n) < 0:
            return n * a * (1 - a)
        else:
            return n * (1 - a)

def apn(a, p, n):
    """
    This is a little silly, and we should do this better.
    """
    if n == 0:
        return 0
    if n == 1:
        return a
    if n == 2:
        return a^2 - p
    if n == 3:
        return a^3 - 2*a*p
    if n == 4:
        return a^4 - 3*a^2*p + p^2
    if n == 5:
        return a^5 - 4*a^3*p + 3*a*p^2
    if n == 6:
        return a^6 - 5*a^4*p + 6*a^2*p^2 - p^3
    if n == 7:
        return a^7 - 6*a^5*p + 10*a^3*p^2 - 4*a*p^3
    if n == 8:
        return a^8 - 7*a^6*p + 15*a^4*p^2 - 10*a^2*p^3 + p^4
    if n == 9:
        return a^9 - 8*a^7*p + 21*a^5*p^2 - 20*a^3*p^3 + 5*a*p^4
    if n == 10:
        return a^10 - 9*a^8*p + 28*a^6*p^2 - 35*a^4*p^3 + 15*a^2*p^4 - p^5


def anlist_(E, bound=2000):
    level = E.conductor()
    
    primes = [make_tot_pos(p.sage_ideal().gens_reduced()[0]) for p in primes_of_bounded_norm(bound)]
    aps = aplist(E, bound)

    bad_primes = [ [p.norm(), sigma(p), sigma2(p), ap] for p, ap in zip(primes, aps) if F.ideal(p).divides(level)]
    primes = [ [p.norm(), sigma(p), sigma2(p), ap] for p, ap in zip(primes, aps)]
    prime_powers = [ [p,] for p in primes]
    for P in prime_powers:
        norm = P[0][0]
        ap = P[0][3]
        if P[0] in bad_primes:
            while P[0][0] * P[-1][0] < bound:
                P.append( [x*y for (x,y) in zip(P[0], P[-1])] )
        else:
            while P[0][0] * P[-1][0] < bound:
                P.append( [x*y for (x,y) in zip(P[0], P[-1])] )
            for n in range(len(P)):
                P[n][3] = apn(ap, norm, n+1)

    L = [[1, 1.0, 1.0, 1]]
    for P in prime_powers:
        L2 = [l for l in L]
        for p in P:
            for n in L:
                if n[0] * p[0] < bound:
                    L2.append([x*y for (x,y) in zip(p,n)])
                else:
                    break
        L2.sort()
        L = L2
    return L

def Lambdaprime_one(E, bound=300):
    N = make_tot_pos(E.conductor().gens_reduced()[0])
    anlist = anlist_(E, bound)
    rootN = RR(sqrt(sigma(N)))
    rootNbar = RR(sqrt(sigma2(N)))
    rootD = RR(sqrt(5))
    PI = RR(pi)

    S = RR(0)
    #E1 = exponential_integral_1
    E1 = lambda t : RR(mpmath.expint(1, t))

    for (norm, n, nbar, an) in anlist:
        X = 0
        for k in srange(-10, 10):
            A = RR(2 * PI * n * phi^(2*k + 1)/(rootN * rootD))
            B = RR(2 * PI * nbar * phibar^(2*k + 1)/(rootNbar * rootD))
            
            x = E1(A) * ( (B/phi).exp() - (B*phi).exp() )
            y = (-A).exp() * ( log(phi) * ((B*phi).exp() + (B/phi).exp()) - E1(-B/phi) + E1(-B*phi) )
            
            X = X + (-x + y)/(A*B)

        S = S + an*X
        print norm, S
        sys.stdout.flush()

    return 2 * S

def Lprime_one(E, bound=300):
    # compute L'(E, 1) assuming that L(E, 1) = 0 and
    # L(E,s) has odd sign.

    Q = RR(5 * sqrt(E.conductor().norm())/(4 * pi^2))
    return Lambdaprime_one(E, bound)/Q
