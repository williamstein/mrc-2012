#######################################
#
# Fast computation of G_r(x) functions
#

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

def Gr_series(r, R=RDF, prec=None, MAX=20):
    """
    power series part of G_r(x), defined as
    .. math::
    
        CG_r(x) = \sum_{n=1}^{\infty} \frac{{-1}^{n-r}}{n^r\cdot n!}\,x^n
    
    see AMEC p.44
    """
    t = cputime()
    MAX = R(MAX)
    x = PowerSeriesRing(R, 'x').gen()
    if prec is None:
        ser = 0
        ser_max = R(0)
        for n in range(1, 1000):
            ser += (-x)**n / (-n)**r / factorial(n)
            next = ser(MAX)
            # test if the last term makes a difference for evaluating at MAX
            if next == ser_max:
                verbose("Computed Gr_series(%d) to precision %d" % (r,n), t)
                return ser + O(x**n)
            ser_max = next
        else:
            raise RuntimeError, "Reached end of loop in Gr_series!"
    else:
        return sum( (-x)**n / (-n)**r / factorial(n)   for n in range(1,prec)) + O(x**prec)

# Slow...
def Gr(r, x):
    """
    G_r function, defined as G_r(x) = P_r(-log(x)) + CG_r(x)
    """ 
    R = parent(x)
    return Pr_poly(r, R) (x) + Gr_series(r, R, MAX=x) (x)

@cached_function
def Gr_fast_float(r, MAX=20):
    from sage.ext.fast_eval import fast_float_arg, fast_float_constant, fast_float_func
    x = fast_float_arg(0)
    return Pr_poly(r, RDF) (-log(x)) + Gr_series(r, RDF, MAX=MAX) (x)

# precache G_0, G_1, G_2, G_3, G_4
Gr_fast_float(0)
Gr_fast_float(1)
Gr_fast_float(2)
Gr_fast_float(3)
Gr_fast_float(4)
