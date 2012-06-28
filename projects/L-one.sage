# some of this code is probably wrong

from psage.number_fields.sqrt5 import *
from psage.ellcurve.lseries.aplist_sqrt5 import aplist

sigma2, sigma = F.embeddings(RR)
phi = sigma(a)
phibar = sigma2(a)

def make_tot_pos(n):
    """
    Given an element n of Q(sqrt5), return a totally positive unit multiple
    of n; i.e. return m such that m/n is a unit and m is taken to a positive
    number under both embeddings of K into R.

    EXAMPLES::
        
        sage: from psage.modform.hilbert.sqrt5.ellcurve import K, a, phi, phi2, make_tot_pos
        sage: N = 15 + 30 * a
        sage: M = make_tot_pos(N); M
        45*a + 30
        sage: M/N
        a
        sage: phi(M) > 0
        True
        sage: phi2(M) > 0
        True
        sage: K.ideal(N) == K.ideal(M)
        True
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

def build_an_list(E, bound = 2000, verbose=0):
    import sys
    level = E.conductor()

    an_list = {}
    ap_dict = {}

    #print "finding ap values and creating prime ideals"
    sys.stdout.flush()

    ideals = F.ideals_of_bdd_norm(bound)
    prime_ideals = primes_of_bounded_norm(bound)
    ap_list = aplist(E, bound)
    for prime, ap in zip(prime_ideals, ap_list):
        ap_dict[prime.sage_ideal()] = ap
    #for n in range(1, bound):
    #    if verbose:
    #        print "looking for primes of norm", n
    #    for p in ideals[n]:
    #        p = K.ideal(p.gens_reduced())
    #        #if is_prime(p):
    #        #    ap_dict[p] = psage.ellcurve.lseries.aplist.ap(E, p)
    #        if is_prime(p):
    #            if p.divides(E.discriminant()):
    #                #print "computing ap for p =", p, "n =", n
    #                ap_dict[p] = psage.ellcurve.lseries.aplist.ap(E, p)
    #            else:
    #                #print "computing ap for p =", p, "n =", n
    #                ap_dict[p] = E.change_ring(p.residue_field()).trace_of_frobenius()
                
    #print "finding an values"
    sys.stdout.flush()

    if 0: # bug somewhere in this code...
        for n in range(1, bound):
            for I in ideals[n]:
                an_list[I] = 1


        for p, ap in zip(prime_ideals, ap_list):
            p = p.sage_ideal()
            pk = 1
            for k in range(1, floor(log(bound)/log(p.norm()))):
                pk = pk * p
                if p.divides(level):
                    apk = ap**k
                else:
                    if k == 1:
                        apk = ap
                    else:
                        apk = ap * an_list[pk/p] - p.norm() * an_list[pk/(p * p)]
                for n in range(1, floor(bound/pk.norm())):
                    for I in ideals[n]:
                        if not p.divides(I):
                            print (pk*I).norm()
                            an_list[pk * I] *= apk

        # finally, since we want our an list to be indexed by totally positive generators,
        # and in slightly different form, we put it in that form

        an_dict = {}
        for n, an in an_list.iteritems():
            N = make_tot_pos(n.gens_reduced()[0])
            an_dict[N] = (sigma(N), sigma2(N), an)

        return an_dict
    else:
        for n in range(1, bound):
            if verbose:
                print "computing an for norm", n
            for I in ideals[n]:
                #print "computing an for norm", n
                an = 1
                f = I.factor()
                for (p, e) in f:
                    p = F.ideal(p.gens_reduced())
                    ap = ap_dict[p]
                    if p.divides(level):
                        an *= ap**e
                    else:
                        if e == 1:
                            an *= ap
                        elif e == 2:
                            an *= (ap**2 - p.norm())
                        elif e == 3:
                            an *= (ap**3 - 2 * p.norm() * ap)

                        elif e == 4:
                            an *= (ap**4 - 4 * ap**2 * p.norm() + p.norm()**2)
                        elif e == 5:
                            an *= ap**5 - 4 * ap**3 * p.norm() + 3 * ap * p.norm()**2
                        elif e == 6:
                            an *= ap**6 - 5*ap**4*p.norm() + 6*ap**2*p.norm()**2 - p.norm()**3
                        elif e == 7:
                            an *= ap**7 - 6*ap**5*p.norm() + 10*ap**3*p.norm()**2 - 4*ap*p.norm()**3
                        else:
                            raise ValueError("not ready to handle things that are divisible by higher than cubic powers")
                N = make_tot_pos(I.gens_reduced()[0])
                an_list[N] = (sigma(N), sigma2(N), an)
        return an_list

def Lambdaprime_one(E, bound=300):
    # something here is broken
    N = make_tot_pos(E.conductor().gens_reduced()[0])
    anlist = build_an_list(E, bound)
    rootN = RR(sqrt(sigma(N)))
    rootNbar = RR(sqrt(sigma2(N)))
    rootD = RR(sqrt(5))
    PI = RR(pi)

    S = RR(0)
    #E1 = exponential_integral_1
    E1 = lambda t : RR(mpmath.expint(1, t))

    ideals = anlist.keys()
    ideals.sort(cmp = lambda a,b: -1r if a.norm() < b.norm() else 0r if a.norm() == b.norm else 1r)
    for n in ideals:
        an = anlist[n][2]
        X = 0
        for k in srange(-10, 10):
            A = RR(2 * PI * sigma(n) * phi^(2*k + 1)/(rootN * rootD))
            B = RR(2 * PI * sigma2(n) * phibar^(2*k + 1)/(rootNbar * rootD))
            
            x = E1(A) * ( (B/phi).exp() - (B*phi).exp() )
            y = (-A).exp() * ( log(phi) * ((B*phi).exp() + (B/phi).exp()) - E1(-B/phi) + E1(-B*phi) )
            
            X = X + (x + y)/(A*B)

        S = S + an*X
        print S
        sys.stdout.flush()

    return S
