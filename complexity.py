from psage.modform.hilbert.sqrt5.hmf import *

def dim(N):
    """
    This seems to be an upper bound on the dimension for large N
    """
    return N.norm()/30.0

def kernel_compute_time(N,p):
    hasse = 2*Integer(p.norm()).isqrt()+1
    return 2.0*hasse*(dim(N)**2)/(10**6)

def coeff_compute_time(Y):
    """
    Time it takes to compute the a_p up to norm Y
    """
    return Y**2 #need to put in correct time!

def compute_time(X):
    t1 = 0
    #Ideals = F.ideals_of_bdd_norm(X)
    #for i in range(1,X):
    #    if Ideals[i] != 0:
    #        for N in Ideals[i]:
    #            t1 = t1 + kernel_compute_time(N,F.prime_above(2))
    Primes = primes_of_bounded_norm(X)
    for p in Primes:
        t1 = t1 + kernel_compute_time(p,F.prime_above(2))


    Y = 1000 #this should depend on X            
    t2 = 1.0*X*coeff_compute_time(Y)
    return t1+t2, t1, t2
