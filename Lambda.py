from sage.misc.all import cputime
from sage.rings.all import Integer, ZZ
from psage.number_fields.sqrt5 import *
from psage.ellcurve.lseries.aplist_sqrt5 import aplist
from psage.modform.hilbert.sqrt5.hmf import *
import mpmath

R=RealField(170)
sigma2, sigma = F.embeddings(RR)
phi = R(sigma(a))
phibar = R(sigma2(a))

#V is an EllipticCurveFactor

def Lambda(V,r, bound=300):
    """
    INPUT:
    
    - V is an EllipticCurveFactor
    - r is a whole number
    - bound?
    
    OUTPUT:
    
    - value of the completed L-function at 1
    """
    
    if r==1:
        return Lambda1(V,bound)
    elif r==2:
        return Lambda2(V,bound)
    
def Lambda1(V, bound):    
    N = make_tot_pos(V.conductor().gens_reduced()[0])
    anlist = anlist_(V, bound)
    rootN = R(sqrt(sigma(N)))
    rootNbar = R(sqrt(sigma2(N)))
    rootD = R(sqrt(5)) #specific to our F
    D = 5 #specific to our F
    PI = R(pi)

    S = R(0)
    
    for (norm, n, nbar, an) in anlist:
        X = 0
        for k in srange(-10, 10):
            A = R(2 * PI * n * phi^(2*k + 1)/(rootN * rootD))
            B = R(-2 * PI * nbar * phibar^(2*k + 1)/(rootNbar * rootD))          
            
            X = X + (-log(phi)*G(0,A)*G(0,B/phi)+G(0,A)*G(1,B/phi)+G(1,A)*G(0,B/phi)-log(phi)*G(0,A)*G(0,B*phi)-G(0,A)*G(1,B*phi)-G(1,A)*G(0,B*phi))/(A*B)

        S = S + an*X
        #print norm, S
        #sys.stdout.flush()

    return 2 * S
    
def Lambda2(V, bound):    
    N = make_tot_pos(V.conductor().gens_reduced()[0])
    anlist = anlist_(V, bound)
    rootN = R(sqrt(sigma(N)))
    rootNbar = R(sqrt(sigma2(N)))
    rootD = R(sqrt(5)) #specific to our F
    D = 5 #specific to our F
    PI = R(pi)

    S = R(0)
    
    for (norm, n, nbar, an) in anlist:
        X = 0
        for k in srange(-15, 15):
            A = R(2 * PI * n * phi^(2*k + 1)/(rootN * rootD))
            B = R(-2 * PI * nbar * phibar^(2*k + 1)/(rootNbar * rootD))            
            
            X = X + (-log(phi)*G(0,A)*G(1,B/phi)-log(phi)*G(1,A)*G(0,B/phi)+G(0,A)*G(2,B/phi)+G(1,A)*G(1,B/phi)+G(2,A)*G(0,B/phi)-log(phi)*G(0,A)*G(1,B*phi)-log(phi)*G(1,A)*G(0,B*phi)-G(0,A)*G(2,B*phi)-G(1,A)*G(1,B*phi)-G(2,A)*G(0,B/phi))/(A*B)

        S = S + an*X
        #print norm, S
        #sys.stdout.flush()

    return 4 * S
        

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


def anlist_(V, bound=2000):
    level = V.conductor()
    
    primes = [make_tot_pos(p.sage_ideal().gens_reduced()[0]) for p in primes_of_bounded_norm(bound)]
    aps = V.aplist(bound)

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

#I found some code in Gr_series.py

def Qr_poly(r, R=RDF):
    """
    Q_r(t) polynomial from AMEC p.44
    """
    Rt,t = PolynomialRing(R, 't').objgen()
    if r==0:
        return Rt(1)
    if r==1:
        return t
    tt = PowerSeriesRing(R, 'tt').gen()
    ps = sum((-1)**n * zeta(R(n))/n * tt**n for n in range(2,r+1)).exp(r+1)
    return sum( ps[r-i] * t**i/factorial(i) for i in range(0,r+1) )

def Pr_poly(r, R=RDF):
    """
    P_r(t) polynomial from AMEC p.44
    """
    Qr = Qr_poly(r, R)
    t = parent(Qr).gen()
    gamma = R(euler_gamma)
    return Qr(t-gamma)

#This was in Gr_series.py, but Cohen says not to use this for x>10, instead I implemented what Cremona implemented in C++
#def Gr_series(r, R=RDF, prec=None, MAX=20):
#    """
#    power series part of G_r(x), defined as
#    .. math::
#    
#        CG_r(x) = \sum_{n=1}^{\infty} \frac{{-1}^{n-r}}{n^r\cdot n!}\,x^n
#    
#    see AMEC p.44
#    """
#    t = cputime()
#    MAX = R(MAX)
#    x = PowerSeriesRing(R, 'x').gen()
#    if prec is None:
#        ser = 0
#        ser_max = R(0)
#        for n in range(1, 1000):
#            ser += (-x)**n / (-n)**r / factorial(n)
#            next = ser(MAX)
#            # test if the last term makes a difference for evaluating at MAX
#            if next == ser_max:
#                verbose("Computed Gr_series(%d) to precision %d" % (r,n), t)
#                return ser + O(x**n)
#           ser_max = next
#       else:
#            raise RuntimeError, "Reached end of loop in Gr_series!"
#    else:
#        return sum( (-x)**n / (-n)**r / factorial(n)   for n in range(1,prec)) + O(x**prec)


def is_approx_zero(x):
   return abs(x) < 10^(-20)

def CG(r,x): #Cohen's G_r(x), used for small x
    emx=exp(-x)
    ans=x
    term=x
    Av=[1]
    for j in [1..r]:
        Av.append(1)    
    n=1
    while (not is_approx_zero(emx*term*Av[r])):
        n += 1
        for j in [1..r]:
            Av[j] += (Av[j-1]/n) #update Av vector
        term *= (x/n)
        ans += (Av[r]*term)
    return emx*ans


# Slow...
def Gsmall(r, x):
    """
    G_r function, defined as G_r(x) = P_r(-log(x)) + CG_r(x)
    """ 
    R = parent(x)
    return Pr_poly(r, R) (x) + GC(r,x)

@cached_function
def Gsmall_fast_float(r, MAX=20):
    from sage.ext.fast_eval import fast_float_arg, fast_float_constant, fast_float_func
    x = fast_float_arg(0)
    return Pr_poly(r, RDF) (-log(x)) + CG(r,x)

# precache G_0, G_1, G_2, G_3, G_4
Gsmall_fast_float(0)
Gsmall_fast_float(1)
Gsmall_fast_float(2)
Gsmall_fast_float(3)
Gsmall_fast_float(4)

def G(r,x):
    """
    INPUT:
    
    - r is a whole number, G(r,x) is used to compute the r^th derivative
    - x is a real number
    
    OUTPUT:
    
    - the value of Cohen's Gamma_r(1,x), choosing how to compute it depending on how large x is and using an implemented function when r=0, 1
    """
    if r==0:
        return (-x).exp()
    elif r==1:
        E1 = lambda t : RR(mpmath.expint(1, t)) #E1 = exponential_integral_1
        return E1(x)
    elif 0 < x < 50:
        return Gsmall(r,x)
    elif x >= 50:
        return Glarge(r,x)
    else:
        raise ValueError, "Everything should be totally positive, but there is a negative value of x coming up"
        
#The choice of x < 50 should be appropriate to get Gamma_r to 18 significant digits, see discussion on page 582.
     
def Glarge(r,x):  #Cohen's Gamma_r(x) for large x
    emx=exp(-x)
    ans=0
    term=-1/x
    Cv=[1]
    for j in [1..r]:
        Cv.append(0)
    n=0
    while (not is_approx_zero(emx*term)):
        n = n + 1
        for j in [1..r]:
            Cv[-j] = Cv[-j] + (Cv[-j-1]/n) #update Cv vector
        term = term * (-n/x)
        ans = ans + (Cv[r]*term)
    return 2*emx*ans
    









